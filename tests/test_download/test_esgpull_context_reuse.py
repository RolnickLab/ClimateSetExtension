import shutil
from unittest.mock import AsyncMock, patch

import pytest
import sqlalchemy as sa
from esgpull.models import File

from climateset.download.esgpull_utils import (
    esgpull_search_and_download_esgf_model_single_var,
    isolated_esgpull_context,
)


@pytest.fixture
def clean_tmp_path(tmp_path):
    yield tmp_path
    if tmp_path.exists():
        shutil.rmtree(tmp_path, ignore_errors=True)


class TestEsgpullContextReuse:
    @pytest.mark.integration
    @patch("climateset.download.esgpull_utils.Esgpull.download", new_callable=AsyncMock)
    def test_context_reuse_prevents_facet_bleed(self, mock_download, clean_tmp_path):
        """
        Verify that reusing a single esg context across sequential calls does not cause search facets (like variable or
        experiment) to bleed between queries.

        We mock the actual download to save time/bandwidth, but let the search hit the network.
        """
        mock_download.return_value = ([], [])

        with isolated_esgpull_context(clean_tmp_path) as esg:
            # First call: Search for 'tas' in 'historical'
            _ = esgpull_search_and_download_esgf_model_single_var(
                esg=esg,
                model="NorESM2-LM",
                variable="tas",
                experiment="historical",
                project="CMIP6",
                default_grid_label="gn",
                default_frequency="mon",
                preferred_version="latest",
                max_ensemble_members=1,
                ensemble_members=[],
                data_dir=clean_tmp_path,
                distrib=True,
            )

            # Retrieve files added to the DB in the first call
            files_1 = list(esg.db.session.scalars(sa.select(File)))
            assert len(files_1) > 0, "First search should return results"

            # Verify constraints were respected
            for file in files_1:
                assert "tas" in file.file_id
                assert "historical" in file.file_id

            # Clear DB to cleanly assert on the second query
            for file in files_1:
                esg.db.delete(file)

            # Second call: Search for 'pr' in 'ssp126'
            _ = esgpull_search_and_download_esgf_model_single_var(
                esg=esg,
                model="NorESM2-LM",
                variable="pr",
                experiment="ssp126",
                project="CMIP6",
                default_grid_label="gn",
                default_frequency="mon",
                preferred_version="latest",
                max_ensemble_members=1,
                ensemble_members=[],
                data_dir=clean_tmp_path,
                distrib=True,
            )

            # Retrieve files added to the DB in the second call
            files_2 = list(esg.db.session.scalars(sa.select(File)))
            assert len(files_2) > 0, "Second search should return results"

            # Crucial assertion: ensure no bleed from the first query ('tas', 'historical')
            for file in files_2:
                assert "pr" in file.file_id
                assert "ssp126" in file.file_id
                assert "tas" not in file.file_id, "Facet bleed detected: 'tas' found in 'pr' query results"
                assert (
                    "historical" not in file.file_id
                ), "Facet bleed detected: 'historical' found in 'ssp126' query results"

    def test_isolated_context_lifecycle(self, clean_tmp_path):
        """Verify that isolated_esgpull_context creates exactly one UUID directory and cleans it up afterwards."""
        jobs_dir = clean_tmp_path / ".esgpull_jobs"

        assert not jobs_dir.exists() or len(list(jobs_dir.iterdir())) == 0

        with isolated_esgpull_context(clean_tmp_path) as esg:
            assert jobs_dir.exists(), ".esgpull_jobs directory should be created"

            uuid_dirs = list(jobs_dir.iterdir())
            assert len(uuid_dirs) == 1, "Exactly one UUID directory should be created"

            uuid_dir = uuid_dirs[0]
            assert uuid_dir.is_dir(), "The item should be a directory"
            # In the latest esgpull the root path might not be exposed as `root` directly,
            # but rather data/db paths are inside it.
            assert uuid_dir.name in str(esg.config.paths.data), "Esgpull data path should be inside the UUID dir"

        # After the context manager exits, the UUID directory should be removed
        assert not uuid_dir.exists(), "UUID directory should be cleaned up on exit"
