from pathlib import Path

from climateset.download.abstract_downloader import AbstractDownloader
from climateset.download.esgpull_utils import (
    esgpull_search_and_download_esgf_biomass_single_var,
    esgpull_search_and_download_esgf_model_single_var,
    esgpull_search_and_download_esgf_raw_single_var,
)
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class EsgpullDownloader(AbstractDownloader):
    def __init__(self, config=None, distrib: bool = False):
        self.config = config
        self.distrib = distrib
        self.logger = LOGGER

    def download(self):
        # Dispatch based on config type (or could be an abstract base)
        # Note: EsgpullDownloader executes searches via esgpull.
        # The actual download logic via esgpull is in task 4,
        # so for now `download` can just invoke search to satisfy the interface.
        pass

    def esgpull_search_and_download_esgf_raw_single_var(
        self,
        variable: str,
        institution_id: str,
        project: str,
        default_grid_label: str,
        default_frequency: str,
        preferred_version: str,
        data_dir: Path | str,
    ):
        return esgpull_search_and_download_esgf_raw_single_var(
            variable=variable,
            institution_id=institution_id,
            project=project,
            default_grid_label=default_grid_label,
            default_frequency=default_frequency,
            preferred_version=preferred_version,
            data_dir=data_dir,
            distrib=self.distrib,
            logger=self.logger,
        )

    def esgpull_search_and_download_esgf_biomass_single_var(
        self,
        variable: str,
        variable_id: str,
        institution_id: str,
        project: str,
        default_grid_label: str,
        default_frequency: str,
        preferred_version: str,
        data_dir: Path | str,
    ):
        return esgpull_search_and_download_esgf_biomass_single_var(
            variable=variable,
            variable_id=variable_id,
            institution_id=institution_id,
            project=project,
            default_grid_label=default_grid_label,
            default_frequency=default_frequency,
            preferred_version=preferred_version,
            data_dir=data_dir,
            distrib=self.distrib,
            logger=self.logger,
        )

    def esgpull_search_and_download_esgf_model_single_var(
        self,
        model: str,
        variable: str,
        experiment: str,
        project: str,
        default_grid_label: str,
        default_frequency: str,
        preferred_version: str,
        max_ensemble_members: int,
        ensemble_members: list[str],
        data_dir: Path | str,
    ):
        return esgpull_search_and_download_esgf_model_single_var(
            model=model,
            variable=variable,
            experiment=experiment,
            project=project,
            default_grid_label=default_grid_label,
            default_frequency=default_frequency,
            preferred_version=preferred_version,
            max_ensemble_members=max_ensemble_members,
            ensemble_members=ensemble_members,
            data_dir=data_dir,
            distrib=self.distrib,
            logger=self.logger,
        )
