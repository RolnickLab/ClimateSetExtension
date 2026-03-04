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


def test_search_and_download_esgf_raw_single_var(mock_esg_context, tmp_path):
    downloader = EsgpullDownloader()
    mock_esg_context.context.hints.side_effect = [
        [{"grid_label": {"gn": 10}}],
        [{"nominal_resolution": {"100 km": 5}}],
        [{"frequency": {"mon": 10}}],
        [{"target_mip": {"CMIP": 10}}],
        [{"version": {"v2020": 10}}],
    ]
    mock_esg_context.context.search.return_value = ["file1", "file2"]

    files = downloader.search_and_download_esgf_raw_single_var(
        variable="tas",
        institution_id="INST",
        project="input4MIPs",
        default_grid_label="gn",
        default_frequency="mon",
        preferred_version="latest",
        data_dir=tmp_path,
    )

    assert files == ["file1", "file2"]
    mock_esg_context.context.search.assert_called_once()

    # Assert query passed to search
    called_query = mock_esg_context.context.search.call_args[0][0]
    assert called_query.selection["project"] == ["input4MIPs"]
    assert called_query.selection["variable"] == ["tas"]
    assert called_query.selection["institution_id"] == ["INST"]
    assert called_query.selection["grid_label"] == ["gn"]
    assert called_query.selection["frequency"] == ["mon"]
    assert called_query.selection["target_mip"] == ["CMIP"]
    assert called_query.options.latest.name == "true"


def test_search_and_download_esgf_model_single_var(mock_esg_context, tmp_path):
    downloader = EsgpullDownloader()
    mock_esg_context.context.hints.side_effect = [
        [{"frequency": {"mon": 10}}],
        [{"grid_label": {"gn": 10}}],
        [{"variant_label": {"r1i1p1f1": 10, "r2i1p1f1": 10}}],
        [{"version": {"v2020": 10}}],
    ]
    mock_esg_context.context.search.return_value = ["file1"]

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

    assert files == ["file1"]
    called_query = mock_esg_context.context.search.call_args[0][0]
    assert called_query.selection["source_id"] == ["Model-1"]
    assert called_query.selection["experiment_id"] == ["historical"]
    assert len(called_query.selection["variant_label"]) == 1
    assert called_query.selection["version"] == ["v2020"]
