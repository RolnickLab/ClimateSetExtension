import copy
import inspect
import logging
from abc import ABC
from pathlib import Path

import yaml

from climateset import CONFIGS, RAW_DATA
from climateset.download.constants.esgf import (
    CMIP6,
    ESGF_PROJECTS,
    ESGF_PROJECTS_CONSTANTS,
    INPUT4MIPS,
)
from climateset.download.utils import match_key_in_list
from climateset.utils import create_logger, get_yaml_config

LOGGER = create_logger(__name__)

AVAILABLE_CONFIGS = frozenset([CMIP6, INPUT4MIPS])


class AbstractDownloaderConfig(ABC):
    def __init__(
        self,
        project: str,
        data_dir: str | Path = RAW_DATA,
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
        self.config_is_valid = True

        self._validate_item_list(
            item_list=self.variables, available_items=self.avail_variables, name_of_item="variable"
        )
        self._validate_item_list(
            item_list=self.experiments, available_items=self.avail_experiments, name_of_item="experiment"
        )

    def _validate_item_list(self, item_list: list[str], available_items: list[str], name_of_item: str) -> None:
        """
        This small function checks that the given items (variables, models, experiments, etc.) are valid for their given
        project (Input4MIPs, CMIP6, etc.).

        Also remove unvalid items from the list of items as to not.

        Args:
            item_list: List of items to check (like self.variables, self.experiments, etc.)
            available_items: List of available items against which to check (like self.avail_variables, etc.)
            name_of_item: Name of item to check. Write lowercase and singular: ie. variable, experiment, etc.

        Returns:
            None
        """
        error_in_item_list = False
        for e in item_list:
            if e not in available_items:
                self.logger.error(f"{name_of_item.capitalize()} [{e}] not supported.")
                item_list.remove(e)
                error_in_item_list = True
        if error_in_item_list:
            self.logger.error(f"Some, or all submitted {name_of_item}s were not found found - Please verify")
            self.logger.error(f"Available {name_of_item}s: {available_items}")
            self.logger.warning(f"List of valid submitted {name_of_item}s: {available_items}")
            self.config_is_valid = False

    @staticmethod
    def _handle_yaml_config_path(config_file_name, config_path):
        if isinstance(config_path, str):
            config_path = Path(config_path)
        if not config_file_name.endswith(".yaml"):
            config_file_name = f"{config_file_name}.yaml"
        config_full_path = config_path / config_file_name
        return config_full_path

    def generate_config_dict(self):
        init_params = inspect.signature(self.__init__).parameters
        init_args = set(init_params.keys()) - {"self"}
        config_dict = {self.project: {}}
        for key, value in self.__dict__.items():
            if key in init_args and key not in ["project", "logger"] and not callable(value):
                config_dict[self.project][key] = value
        return config_dict

    def generate_config_file(self, config_file_name: str, config_path: str | Path = CONFIGS) -> None:
        config_full_path = self._handle_yaml_config_path(config_file_name, config_path)
        data = self.generate_config_dict()
        with open(config_full_path, "w", encoding="utf-8") as config_file:
            yaml.dump(data, config_file, indent=2)

    def add_to_config_file(self, config_file_name: str, config_path: str | Path = CONFIGS) -> None:
        config_full_path = self._handle_yaml_config_path(config_file_name, config_path)
        existing_config = {}
        if config_full_path.exists():
            existing_config = get_yaml_config(config_full_path)
            existing_config.update(existing_config)
        new_config = self.generate_config_dict()
        existing_config.update(new_config)
        with open(config_full_path, "w", encoding="utf-8") as config_file:
            yaml.dump(existing_config, config_file, indent=2)


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
        # (all)
        # (climate model inputs)
        self.biomass_vars: list[str] = []
        self.meta_vars_percentage: list[str] = []
        self.meta_vars_share: list[str] = []

        self._handle_emission_variables()

    def _handle_emission_variables(self):
        self._generate_raw_emission_vars()
        self._generate_plain_emission_vars()
        self.logger.info(f"Emission variables to download: {self.variables}")
        if self.download_biomass_burning:
            self.logger.info(f"Biomass burning vars to download: {self.biomass_vars}")
        if self.download_metafiles:
            self.logger.info(
                f"Meta emission vars to download:\n\t{self.meta_vars_percentage}\n\t{self.meta_vars_share}"
            )

    def _generate_raw_emission_vars(self):
        variables = copy.deepcopy(self.variables)
        if variables is None:
            raise ValueError("No variables have been given to the downloader. Variables must be given for downloader.")
        self.variables = [v.replace(" ", "_").replace("-", "_") for v in variables]
        self.logger.info(f"Cleaned variables : {self.variables}")

    def _generate_plain_emission_vars(self):
        if self.use_plain_emission_vars:
            # plain vars are biomass vars
            self.biomass_vars = self.variables
            self.variables = [
                variable + emission_ending for variable in self.variables for emission_ending in self.emissions_endings
            ]
            # be careful with CO2
            if "CO2_em_openburning" in self.variables:
                self.variables.remove("CO2_em_openburning")
        else:
            # get plain input4mips vars = biomass vars for historical
            self.biomass_vars = list({v.split("_")[0] for v in self.variables})
            # remove biomass vars from normal vars list
            for b in self.biomass_vars:
                try:
                    self.variables.remove(b)
                except Exception as error:  # pylint: disable=W0718
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
        models: list[str] = None,
        experiments: list[str] = None,
        ensemble_members: list[str] = None,  # preferred ensemble members used, if None not considered
        max_ensemble_members: int = 10,  # if -1 take all
        variables: list[str] = None,
        overwrite: bool = False,
        logger: logging.Logger = LOGGER,
    ):
        super().__init__(project, data_dir, experiments, variables, overwrite, logger)

        if not models:
            models = ["NorESM2-LM"]
        if isinstance(models, str):
            models = [models]
        self.models: list[str] = models
        self.avail_models = self.proj_constants.MODEL_SOURCES
        self.ensemble_members: list[str] = ensemble_members
        self.max_ensemble_members: int = max_ensemble_members

        self._validate_item_list(item_list=self.models, available_items=self.avail_models, name_of_item="model")


def _get_config_from_file(config_file, config_id, config_class, logger=LOGGER):
    configs = get_yaml_config(config_file)
    config_key = config_id
    if config_key not in configs:
        config_key = match_key_in_list(config_key, list(configs.keys()))
    if not config_key:
        logger.error(f"Config key [{config_id}] not found in config file [{config_file}]")
    class_configs = configs[config_key]
    config_object = config_class(project=config_id, **class_configs)
    return config_object


def create_input4mips_downloader_config_from_file(config_file) -> Input4mipsDownloaderConfig:
    config_object = _get_config_from_file(config_file, INPUT4MIPS, Input4mipsDownloaderConfig)
    return config_object


def create_cmip6_downloader_config_from_file(config_file) -> CMIP6DownloaderConfig:
    config_object = _get_config_from_file(config_file, CMIP6, CMIP6DownloaderConfig)
    return config_object
