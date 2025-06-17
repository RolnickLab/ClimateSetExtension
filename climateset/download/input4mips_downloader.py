from climateset.download.abstract_downloader import AbstractDownloader
from climateset.download.constants.esgf import INPUT4MIPS
from climateset.download.downloader_config import (
    Input4mipsDownloaderConfig,
    create_input4mips_downloader_config_from_file,
)
from climateset.download.utils import (
    search_and_download_esgf_biomass_single_var,
    search_and_download_esgf_raw_single_var,
)
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class Input4MipsDownloader(AbstractDownloader):
    def __init__(self, config: Input4mipsDownloaderConfig):
        self.config: Input4mipsDownloaderConfig = config
        self.logger = LOGGER

    def download(self):
        for variable in self.config.variables:
            if variable.endswith("openburning"):
                institution_id = "IAMC"
            else:
                institution_id = "PNNL-JGCRI"
            self.logger.info(f"Downloading data for variable: {variable}")
            self.download_raw_input_single_var(variable=variable, institution_id=institution_id)

        if self.config.download_biomass_burning & ("historical" in self.config.experiments):
            for variable in self.config.biomass_vars:
                self.logger.info(f"Downloading biomassburing data for variable: {variable}")
                self.download_raw_input_single_var(variable=variable, institution_id="VUA")

        if self.config.download_metafiles:
            for variable in self.config.meta_vars_percentage:
                # percentage are historic and have no scenarios
                self.logger.info(f"Downloading meta percentage data for variable: {variable}")
                self.download_meta_historic_biomassburning_single_var(variable=variable, institution_id="VUA")
            for variable in self.config.meta_vars_share:
                self.logger.info(f"Downloading meta openburning share data for variable: {variable}")
                self.download_raw_input_single_var(variable=variable, institution_id="IAMC")

    def download_raw_input_single_var(
        self,
        variable: str,
        project: str = INPUT4MIPS,
        institution_id: str = "PNNL-JGCRI",
        default_frequency: str = "mon",
        preferred_version: str = "latest",
        default_grid_label: str = "gn",
    ):
        """
        Function handling the download of all input4mips data associated with a single variable.

        Args:
            variable: variable ID
            project: umbrella project, here "input4MIPs"
            institution_id: id of the institution that provides the data
            default_frequency: default frequency to download
            preferred_version: data upload version, if 'latest', the newest version will get selected always
            default_grid_label: default gridding method in which the data is provided
        """
        self.logger.info("Using download_raw_input_single_var() function")

        # Search context is sensitive to order and sequence, which is why
        # it's done in different steps instead of putting everything in `new_context`
        results_list = search_and_download_esgf_raw_single_var(
            variable=variable,
            project=project,
            institution_id=institution_id,
            default_grid_label=default_grid_label,
            default_frequency=default_frequency,
            preferred_version=preferred_version,
            data_dir=self.config.data_dir,
        )
        self.logger.info(f"Download results: {results_list}")

    def download_meta_historic_biomassburning_single_var(
        self,
        variable: str,
        institution_id: str,
        project: str = INPUT4MIPS,
        default_grid_label: str = "gn",
        default_frequency: str = "mon",
        preferred_version: str = "latest",
    ):
        """
        Function handling the download of all metadata associated with a single input4mips variable.

        Args:
            variable: variable ID
            project: umbrella project
            institution_id: id of the institution that provides the data
            default_grid_label: default gridding method in which the data is provided
            default_frequency: default frequency to download
            preferred_version: data upload version, if 'latest', the newest version will get selected always
        """
        variable_id = variable.replace("_", "-")
        variable_search = f"percentage_{variable_id.replace('-', '_').split('_')[-1]}"
        self.logger.info(variable, variable_id, institution_id)

        # Search context is sensitive to order and sequence, which is why
        # it's done in different steps instead of putting everything in `new_context`
        results = search_and_download_esgf_biomass_single_var(
            variable=variable_search,
            variable_id=variable_id,
            project=project,
            institution_id=institution_id,
            default_grid_label=default_grid_label,
            default_frequency=default_frequency,
            preferred_version=preferred_version,
            base_path=self.config.data_dir,
        )
        self.logger.info(f"Download results: {results}")


def input4mips_download_from_config(config):
    config_object = create_input4mips_downloader_config_from_file(config)
    downloader = Input4MipsDownloader(config=config_object)
    downloader.download()
