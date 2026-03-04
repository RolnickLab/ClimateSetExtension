from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from esgpull.models import Query

from climateset.download.esgpull_downloader import (
    EsgpullDownloader,
    _apply_facet_fallback,
    _apply_version_fallback,
)
from climateset.download.utils import isolated_esgpull_context


@pytest.fixture
def real_esg_context(tmp_path):
    with isolated_esgpull_context(tmp_path) as esg:
        yield esg


@pytest.fixture(autouse=True)
def mock_esgpull_download():
    """Intercept Esgpull.download globally for these tests to prevent actual data transfer."""
    with patch("esgpull.esgpull.Esgpull.download", new_callable=AsyncMock) as mock_download:
        mock_download.return_value = ([], [])
        yield mock_download


def test_apply_facet_fallback_preferred_real(real_esg_context):
    """Verifies that facet fallback works with real esgpull context (network hit for hints)."""
    query = Query()
    query.selection["project"] = ["CMIP6"]
    query.selection["variable"] = ["tas"]
    logger = MagicMock()

    # Hit real ESGF for hints
    _apply_facet_fallback(real_esg_context, query, "grid_label", "gn", logger)

    # We expect 'gn' to be available for tas in CMIP6
    assert query.selection.grid_label == ["gn"]


def test_apply_version_fallback_latest():
    mock_esg = MagicMock()
    query = Query()
    logger = MagicMock()

    _apply_version_fallback(mock_esg, query, "latest", logger)
    assert query.options.latest.name == "true"


def test_search_and_download_esgf_model_single_var_real_search(tmp_path, mock_esgpull_download):
    """Verifies that search_and_download_esgf_model_single_var performs a real search but uses the mocked download
    function."""
    downloader = EsgpullDownloader(distrib=True)

    # This will hit real ESGF nodes for the search phase
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

    # Verification:
    # 1. Search should have found something
    # (Since we mocked download to return empty, and move_files relies on cached files,
    # it might return [] if we don't mock more, but we want to check if search was called)

    # In search_and_download_esgf_model_single_var, it returns _download_and_move_files(...)
    # which returns [] because nothing was actually downloaded.
    assert files == []

    # 2. Esgpull.download should have been called
    assert mock_esgpull_download.called


def test_search_and_download_esgf_raw_single_var_real_search(tmp_path, mock_esgpull_download):
    downloader = EsgpullDownloader(distrib=True)

    files = downloader.search_and_download_esgf_raw_single_var(
        variable="tas",
        institution_id="MRI",  # Known variable for MRI in input4MIPs
        project="input4MIPs",
        default_grid_label="gn",
        default_frequency="3hrPt",
        preferred_version="latest",
        data_dir=tmp_path,
    )

    assert files == []
    assert mock_esgpull_download.called


@pytest.mark.integration
def test_esgpull_downloader_integration_search_real(tmp_path):
    """This is effectively redundant now that all tests do real searches, but we keep it to verify the whole flow."""
    downloader = EsgpullDownloader(distrib=True)

    with patch("esgpull.esgpull.Esgpull.download", new_callable=AsyncMock) as mock_download:
        mock_download.return_value = ([], [])

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

        assert files is not None
        assert mock_download.called
