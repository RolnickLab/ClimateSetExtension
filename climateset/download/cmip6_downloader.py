from pathlib import Path

from esgpull import Esgpull

from climateset.download.abstract_downloader import AbstractDownloader
from climateset.download.constants.esgf import CMIP6
from climateset.download.downloader_config import (
    CMIP6DownloaderConfig,
    create_cmip6_downloader_config_from_file,
)
from climateset.download.esgpull_utils import (
    esgpull_search_and_download_esgf_model_single_var,
    isolated_esgpull_context,
)
from climateset.download.utils import search_and_download_esgf_model_single_var
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class CMIP6Downloader(AbstractDownloader):
    def __init__(self, config: CMIP6DownloaderConfig):
        self.logger = LOGGER
        self.config = config

    def download(self):
        """
        Function handling the download of all variables that are associated with a model's output.

        Searches for all files associated with the respected variables and experiment that the downloader
        was initialized with.

        A search connection is established and the search is iteratively constraint to meet all specifications.
        Data is downloaded and stored in a separate file for each year. The default format is netCDF4.

        Resulting hierarchy:

        `CMIPx/model_id/ensemble_member/experiment/variable/nominal_resolution/frequency/year.nc`

        If the constraints cannot be met, per default behaviour for the downloader to select first other
        available value
        """
        for model in self.config.models:
            self.logger.info(f"Downloading data for model: [{model}]")
            for variable in self.config.variables:
                self.logger.info(f"Downloading data for variable: [{variable}]")
                for experiment in self.config.experiments:
                    self.logger.info(f"Downloading data for experiment: [{experiment}]")
                    self.download_from_model_single_var(
                        model=model, project=self.config.project, variable=variable, experiment=experiment
                    )

    def download_from_model_single_var(
        self,
        model: str,
        variable: str,
        experiment: str,
        project: str = CMIP6,
        default_frequency: str = "mon",
        preferred_version: str = "latest",
        default_grid_label: str = "gn",
    ):
        """
        Function handling the download of a single variable-experiment pair that is associated with a model's output
        (CMIP data).

        Args:
            model (str): The model ID
            variable: variable ID
            experiment: experiment ID
            project: umbrella project id e.g. CMIPx
            default_frequency: default frequency to download
            preferred_version: data upload version, if 'latest', the newest version will get selected always
            default_grid_label: default gridding method in which the data is provided
        """
        results_list = search_and_download_esgf_model_single_var(
            model=model,
            variable=variable,
            experiment=experiment,
            project=project,
            default_frequency=default_frequency,
            default_grid_label=default_grid_label,
            preferred_version=preferred_version,
            ensemble_members=self.config.ensemble_members,
            max_ensemble_members=self.config.max_ensemble_members,
            base_path=self.config.data_dir,
        )
        self.logger.info(f"Download results: {results_list}")


class CMIP6DownloaderV2(AbstractDownloader):
    def __init__(self, config: CMIP6DownloaderConfig, distrib: bool = True):
        self.logger = LOGGER
        self.distrib = distrib
        self.config = config

    def download(self):
        """
        Function handling the download of all variables that are associated with a model's output.

        Searches for all files associated with the respected variables and experiment that the downloader
        was initialized with.

        A search connection is established and the search is iteratively constraint to meet all specifications.
        Data is downloaded and stored in a separate file for each year. The default format is netCDF4.

        Resulting hierarchy:

        `CMIPx/model_id/ensemble_member/experiment/variable/nominal_resolution/frequency/year.nc`

        If the constraints cannot be met, per default behaviour for the downloader to select first other
        available value
        """
        with isolated_esgpull_context(self.config.data_dir) as esg:
            for model in self.config.models:
                self.logger.info(f"Downloading data for model: [{model}]")
                for variable in self.config.variables:
                    self.logger.info(f"Downloading data for variable: [{variable}]")
                    for experiment in self.config.experiments:
                        self.logger.info(f"Downloading data for experiment: [{experiment}]")
                        self.download_from_model_single_var(
                            esg=esg, model=model, project=self.config.project, variable=variable, experiment=experiment
                        )

    def download_from_model_single_var(
        self,
        esg: Esgpull,
        model: str,
        variable: str,
        experiment: str,
        project: str = CMIP6,
        default_frequency: str = "mon",
        preferred_version: str = "latest",
        default_grid_label: str = "gn",
    ):
        """
        Function handling the download of a single variable-experiment pair that is associated with a model's output
        (CMIP data).

        Args:
            model (str): The model ID
            variable: variable ID
            experiment: experiment ID
            project: umbrella project id e.g. CMIPx
            default_frequency: default frequency to download
            preferred_version: data upload version, if 'latest', the newest version will get selected always
            default_grid_label: default gridding method in which the data is provided
        """
        results_list = esgpull_search_and_download_esgf_model_single_var(
            esg=esg,
            model=model,
            variable=variable,
            experiment=experiment,
            project=project,
            default_frequency=default_frequency,
            default_grid_label=default_grid_label,
            preferred_version=preferred_version,
            ensemble_members=self.config.ensemble_members,
            max_ensemble_members=self.config.max_ensemble_members,
            data_dir=Path(self.config.data_dir),
            distrib=self.distrib,
        )
        self.logger.info(f"Download results: {results_list}")


def cmip6_download_from_config(config):
    config_object = create_cmip6_downloader_config_from_file(config)
    downloader = CMIP6DownloaderV2(config=config_object)
    downloader.download()
