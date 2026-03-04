from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from esgpull.models import Query

from climateset.download.esgpull_downloader import (
    EsgpullDownloader,
    _apply_facet_fallback,
    _apply_version_fallback,
)


@pytest.fixture
def mock_esg_context():
    with patch("climateset.download.esgpull_downloader.isolated_esgpull_context") as mock_isolated:
        mock_esg = MagicMock()
        mock_isolated.return_value.__enter__.return_value = mock_esg
        yield mock_esg


def test_apply_facet_fallback_preferred():
    mock_esg = MagicMock()
    query = Query()
    mock_esg.context.hints.return_value = [{"grid_label": {"gn": 10, "gr": 5}}]
    logger = MagicMock()

    _apply_facet_fallback(mock_esg, query, "grid_label", "gn", logger)
    assert query.selection["grid_label"] == ["gn"]


def test_apply_facet_fallback_not_preferred():
    mock_esg = MagicMock()
    query = Query()
    mock_esg.context.hints.return_value = [{"grid_label": {"gr": 5}}]
    logger = MagicMock()

    _apply_facet_fallback(mock_esg, query, "grid_label", "gn", logger)
    assert query.selection["grid_label"] == ["gr"]


def test_apply_version_fallback_latest():
    mock_esg = MagicMock()
    query = Query()
    logger = MagicMock()

    _apply_version_fallback(mock_esg, query, "latest", logger)
    assert query.options.latest.name == "true"


@patch("climateset.download.esgpull_downloader._download_and_move_files")
def test_search_and_download_esgf_raw_single_var(mock_download_and_move, mock_esg_context, tmp_path):
    downloader = EsgpullDownloader()
    mock_esg_context.context.hints.side_effect = [
        [{"grid_label": {"gn": 10}}],
        [{"nominal_resolution": {"100 km": 5}}],
        [{"frequency": {"mon": 10}}],
        [{"target_mip": {"CMIP": 10}}],
        [{"version": {"v2020": 10}}],
    ]
    mock_esg_context.context.search.return_value = ["file1", "file2"]
    mock_download_and_move.return_value = ["path/to/file1.nc", "path/to/file2.nc"]

    files = downloader.search_and_download_esgf_raw_single_var(
        variable="tas",
        institution_id="INST",
        project="input4MIPs",
        default_grid_label="gn",
        default_frequency="mon",
        preferred_version="latest",
        data_dir=tmp_path,
    )

    assert files == ["path/to/file1.nc", "path/to/file2.nc"]
    mock_esg_context.context.search.assert_called_once()
    mock_download_and_move.assert_called_once_with(
        mock_esg_context,
        ["file1", "file2"],
        tmp_path / "input4MIPs" / "raw_input_vars" / "INST" / "tas",
        downloader.logger,
    )

    # Assert query passed to search
    called_query = mock_esg_context.context.search.call_args[0][0]
    assert called_query.selection["project"] == ["input4MIPs"]
    assert called_query.selection["variable"] == ["tas"]
    assert called_query.selection["institution_id"] == ["INST"]
    assert called_query.selection["grid_label"] == ["gn"]
    assert called_query.selection["frequency"] == ["mon"]
    assert called_query.selection["target_mip"] == ["CMIP"]
    assert called_query.options.latest.name == "true"


@patch("climateset.download.esgpull_downloader._download_and_move_files")
def test_search_and_download_esgf_model_single_var(mock_download_and_move, mock_esg_context, tmp_path):
    downloader = EsgpullDownloader()
    mock_esg_context.context.hints.side_effect = [
        [{"frequency": {"mon": 10}}],
        [{"grid_label": {"gn": 10}}],
        [{"variant_label": {"r1i1p1f1": 10, "r2i1p1f1": 10}}],
        [{"version": {"v2020": 10}}],
    ]
    mock_esg_context.context.search.return_value = ["file1"]
    mock_download_and_move.return_value = ["path/to/file1.nc"]

    files = downloader.search_and_download_esgf_model_single_var(
        model="Model-1",
        variable="tas",
        experiment="historical",
        project="CMIP6",
        default_grid_label="gn",
        default_frequency="mon",
        preferred_version="v2020",
        max_ensemble_members=1,
        ensemble_members=[],
        data_dir=tmp_path,
    )

    assert files == ["path/to/file1.nc"]
    mock_download_and_move.assert_called_once_with(
        mock_esg_context, ["file1"], tmp_path / "CMIP6" / "Model-1" / "tas", downloader.logger
    )
    called_query = mock_esg_context.context.search.call_args[0][0]
    assert called_query.selection["source_id"] == ["Model-1"]
    assert called_query.selection["experiment_id"] == ["historical"]
    assert len(called_query.selection["variant_label"]) == 1
    assert called_query.selection["version"] == ["v2020"]


@pytest.mark.integration
def test_esgpull_downloader_integration_search(tmp_path):
    """
    Integration test that performs a real search against ESGF using esgpull.

    Requires network access.
    """
    downloader = EsgpullDownloader(distrib=True)

    # Do a very specific search to limit results and ensure we get something predictable
    files = downloader.search_and_download_esgf_model_single_var(
        model="CanESM5",
        variable="tas",
        experiment="historical",
        project="CMIP6",
        default_grid_label="gn",
        default_frequency="mon",
        preferred_version="latest",
        max_ensemble_members=1,
        ensemble_members=["r1i1p1f1"],
        data_dir=tmp_path,
    )

    # Should at least find something
    assert files is not None
    assert len(files) > 0

    # Check that returned objects are Paths to the downloaded chunk files
    assert isinstance(files[0], Path)

    # Verify the file name matches our query
    filename = files[0].name
    assert filename.endswith(".nc")
    assert "CanESM5" in filename
    assert "tas" in filename
    assert "historical" in filename
