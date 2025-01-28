import logging
from abc import ABC
from pathlib import Path
from typing import Union

import yaml

from climateset import CONFIGS, RAW_DATA
from climateset.download.constants.esgf import ESGF_PROJECTS, ESGF_PROJECTS_CONSTANTS
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class AbstractDownloaderConfig(ABC):
    def __init__(
        self,
        project: str,
        data_dir: Union[str, Path] = RAW_DATA,
        experiments: list[str] = None,
        variables: list[str] = None,
        overwrite: bool = False,
        logger: logging.Logger = LOGGER,
    ):
        self.logger = logger

        self.project = ""
        uppercase_project = project.upper()
        for p in ESGF_PROJECTS:
            if p.upper() == uppercase_project:
                self.project = p

        if self.project not in ESGF_PROJECTS:
            self.logger.error(f"Project {self.project} has not been implemented in the Downloader yet.")
            raise ValueError(
                f"Project {self.project} is not recognized. Consider adding a constant class in download/constants and "
                f"the esgf.py file."
            )

        if isinstance(data_dir, str):
            data_dir = Path(data_dir)
        self.data_dir = data_dir

        self.experiments = experiments
        self.variables = variables
        self.overwrite = overwrite

        # init shared constants
        self.proj_constants = ESGF_PROJECTS_CONSTANTS[self.project]
        self.node_link = self.proj_constants.NODE_LINK
        self.avail_variables = self.proj_constants.VAR_SOURCE_LOOKUP
        self.avail_experiments = self.proj_constants.SUPPORTED_EXPERIMENTS

    def generate_config_file(self, config_file_name: str, config_path: Union[str, Path] = CONFIGS) -> None:
        if isinstance(config_path, str):
            config_path = Path(config_path)
        if not config_file_name.endswith(".yaml"):
            config_file_name = f"{config_file_name}.yaml"

        config_full_path = config_path / config_file_name
        data = {self.project: {}}
        for key, value in self.__dict__.items():
            if key not in ["project", "logger"] and not callable(value):
                data[self.project][key] = value
        with open(config_full_path, "w") as config_file:
            yaml.dump(data, config_file, indent=2)


class Input4mipsDownloaderConfig(AbstractDownloaderConfig):
    def __init__(
        self,
        project: str,
        data_dir: str = RAW_DATA,
        experiments: list[str] = None,
        variables: list[str] = None,
        download_biomassburning: bool = True,  # get biomassburning data for input4mips
        download_metafiles: bool = True,  # get input4mips meta files
        use_plain_emission_vars: bool = True,
        overwrite: bool = False,
        logger: logging.Logger = LOGGER,
    ):
        super().__init__(project, data_dir, experiments, variables, overwrite, logger)

        self.download_metafiles: bool = download_metafiles  # TODO infer automatically from vars
        self.download_biomass_burning: bool = download_biomassburning  # TODO infer automatically from vars
        self.use_plain_emission_vars: bool = use_plain_emission_vars

        self.emissions_endings = self.proj_constants.EMISSIONS_ENDINGS
        self.meta_endings_prc = self.proj_constants.META_ENDINGS_PRC
        self.meta_endings_share = self.proj_constants.META_ENDINGS_SHAR
        self.mip_area = self.proj_constants.MIP_ERA
        self.target_mip = self.proj_constants.TARGET_MIP

        # Attributes that are going to be retrieved / set within this class for
        ## (all)
        self.vars: list[str] = variables
        ## (climate model inputs)
        self.biomass_vars: list[str] = []
        self.meta_vars_percentage: list[str] = []
        self.meta_vars_share: list[str] = []

        self._handle_emission_variables(
            variables=variables,
        )

    def _handle_emission_variables(self, variables: list[str]):
        self.vars = []
        self._generate_raw_emission_vars(variables=variables)
        self._generate_plain_emission_vars()
        self.logger.info(f"Emission variables to download: {self.vars}")
        if self.download_biomass_burning:
            self.logger.info(f"Biomass burning vars to download: {self.biomass_vars}")
        if self.download_metafiles:
            self.logger.info(
                f"Meta emission vars to download:\n\t{self.meta_vars_percentage}\n\t{self.meta_vars_share}"
            )

    def _generate_raw_emission_vars(self, variables: list[str]):
        if variables is None:
            # variables = ["tas", "pr", "SO2_em_anthro", "BC_em_anthro"]
            raise ValueError("No variables have been given to the downloader. Variables must be given for downloader.")
        variables = [v.replace(" ", "_").replace("-", "_") for v in variables]
        self.logger.info(f"Cleaned variables : {variables}")
        for v in variables:
            self.vars.append(v)

    def _generate_plain_emission_vars(self):
        if self.use_plain_emission_vars:
            # plain vars are biomass vars
            self.biomass_vars = self.vars
            self.meta_vars_percentage = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in self.meta_endings_prc
            ]
            self.meta_vars_share = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in self.meta_endings_share
            ]

            self.vars = [
                variable + emission_ending for variable in self.vars for emission_ending in self.emissions_endings
            ]
            # be careful with CO2
            if "CO2_em_openburning" in self.vars:
                self.vars.remove("CO2_em_openburning")
        else:
            # get plain input4mips vars = biomass vars for historical
            self.biomass_vars = list({v.split("_")[0] for v in self.vars})
            # remove biomass vars from normal vars list
            for b in self.biomass_vars:
                try:
                    self.vars.remove(b)
                except Exception as error:
                    self.logger.warning(f"Caught the following exception but continuing : {error}")

            self.meta_vars_percentage = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in self.meta_endings_prc
            ]
            self.meta_vars_share = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in self.meta_endings_share
            ]


class CMIP6DownloaderConfig(AbstractDownloaderConfig):
    def __init__(
        self,
        project: str,
        data_dir: str = RAW_DATA,
        experiments: list[str] = None,
        variables: list[str] = None,
        overwrite: bool = False,
        logger: logging.Logger = LOGGER,
    ):
        super().__init__(project, data_dir, experiments, variables, overwrite, logger)

        self.avail_models = self.proj_constants.MODEL_SOURCES
