import logging
import pathlib
from typing import Union

from pyesgf.search import SearchConnection

from climateset import RAW_DATA
from climateset.download.cmip6_downloader import CMIP6Downloader
from climateset.download.constants.esgf import (
    CMIP6,
    ESGF_MODEL_OUTPUT_LIST,
    ESGF_PROJECTS_CONSTANTS,
    INPUT4MIPS,
)
from climateset.download.downloader_config import (
    AVAILABLE_CONFIGS,
    create_cmip6_downloader_config_from_file,
    create_input4mips_downloader_config_from_file,
    match_project_key,
)
from climateset.download.input4mips_downloader import Input4MipsDownloader
from climateset.download.utils import (
    download_metadata_variable,
    download_model_variable,
    download_raw_input_variable,
    get_upload_version,
    handle_base_search_constraints,
)
from climateset.utils import create_logger, get_yaml_config

LOGGER = create_logger(__name__)


class Downloader:
    """
    Class handling the downloading of the data.

    It communicates with the esgf nodes to search and download the specified data.
    """

    def __init__(  # noqa: C901
        self,
        project: str = "CMIP6",  # default as in ClimateBench
        model: Union[str, None] = "NorESM2-LM",  # default as in ClimateBench
        experiments: list[str] = None,  # sub-selection of ClimateBench default
        variables: list[str] = None,
        data_dir: str = RAW_DATA,
        max_ensemble_members: int = 10,  # if -1 take all
        ensemble_members: list[str] = None,  # preferred ensemble members used, if None not considered
        overwrite: bool = False,  # flag if files should be overwritten
        download_biomassburning: bool = True,  # get biomassburning data for input4mips
        download_metafiles: bool = True,  # get input4mips meta files
        use_plain_emission_vars: bool = True,  # specifies if plain variables are given and rest is inferred
        logger: logging.Logger = LOGGER,
    ):
        """
        Init method for the Downloader.

        Args:
            project (str): Which categorie the data belongs to. Can be: CMIP6, CMIP6Plus, E3SM, input4mips, obs4mips, and more.
                To date, only CMIP6, and input4mips are supported.
            model: Model ID from which output should be downloaded. A list of all supported model ids can
                be found in parameters.constants.MODEL_SOURCES. Model data only.
            experiments:  List of simulations from which data should be downloaded. Model data only.
            experiments: List of variables for which data should be downloaded. Both model and raw data.
            data_dir: Relative or absolute path to the directory where data should be stored. Will be created
                if not yet existent.
            meta_dir: Relative or absolute path to the directory where the metadata should be sored. Will be
                created if not yet existent.
            overwrite: Flag if files should be overwritten, if they already exist.
            download_biomassburning: Flag if biomassburning data for input4mips variables should be downloaded.
            download_metafiles: Flag if metafiles for input4mips variables should be downloaded.
        """
        # Args init for
        ## (all)
        self.logger = logger
        self.project: str = project
        self.data_dir: Union[str, pathlib.Path] = data_dir
        self.overwrite: bool = overwrite
        ## (climate model output) (e.g. cmip6)
        self.model: str = model
        self.experiments: list[str] = experiments
        self.ensemble_members: list[str] = ensemble_members
        self.max_ensemble_members: int = max_ensemble_members
        ## (climate model input) (e.g. input4mips)
        self.download_metafiles: bool = download_metafiles  # TODO infer automatically from vars
        self.download_biomass_burning: bool = download_biomassburning  # TODO infer automatically from vars
        self.use_plain_emission_vars: bool = use_plain_emission_vars  # TODO infer automatically from vars
        # ESGF project constants for
        ## (all)
        self.node_link: str = None
        self.avail_variables: list[str] = None
        self.avail_experiments: list[str] = None
        ## (climate model output)
        self.avail_models: list[str] = None
        ## (climate model input)
        self.emissions_endings: list[str] = None
        self.meta_endings_prc: list[str] = None
        self.meta_endings_share: list[str] = None
        self.mip_area: str = None
        self.target_mip: str = None
        # Attributes that are going to be retrieved / set within this class for
        ## (all)
        self.vars: list[str] = variables
        ## (climate model inputs)
        self.biomass_vars: list[str] = []
        self.meta_vars_percentage: list[str] = []
        self.meta_vars_share: list[str] = []

        self._init_project_constants()

        # if max ensemble member number is too large --> we are relying on the server to complain?

        # adapt variables in case of input4mips
        if self.project == "input4MIPs":
            self._handle_emission_variables(
                variables=variables,
            )

        self._check_desired_params()

    # @Francis (JK) I am still not happy about this function. Please let me know if I can improve this bit
    def _init_project_constants(self):
        """Assign/init attributed depending on the project."""
        if self.project not in ESGF_PROJECTS_CONSTANTS:
            self.logger.info(f"Project {self.project} has not been implemented in the Downloader yet.")
            raise ValueError(
                f"Project {self.project} is not recognized. Consider adding a constant class in download/constants and the esgf.py file."
            )
        proj_constants = ESGF_PROJECTS_CONSTANTS[self.project]

        # init shared constants
        self.node_link = proj_constants.NODE_LINK
        self.avail_variables = proj_constants.VAR_SOURCE_LOOKUP
        self.avail_experiments = proj_constants.SUPPORTED_EXPERIMENTS

        # init climate model output constants
        if self.project in ESGF_MODEL_OUTPUT_LIST:
            self.avail_models = proj_constants.MODEL_SOURCES

        # init input4mips constants
        if self.project == "input4MIPs":
            self.emissions_endings = proj_constants.EMISSIONS_ENDINGS
            self.meta_endings_prc = proj_constants.META_ENDINGS_PRC
            self.meta_endings_share = proj_constants.META_ENDINGS_SHAR
            self.mip_area = proj_constants.MIP_ERA
            self.target_mip = proj_constants.TARGET_MIP

    def _check_desired_params(self):
        """Check if the desired params exist."""
        # check model
        if (self.model is not None) and (self.model not in self.avail_models):
            self.logger.warning(f"Model {self.model} unknown.")
            raise ValueError(
                f"Model {self.model} is not in the list of supported models. Check for typos and consider adding it manually."
            )

        # check experiments
        for exp in self.experiments:
            if exp not in self.avail_experiments:
                self.logger.warning(f"Experiment {exp} unknown.")
                raise ValueError(
                    f"Experiment {exp} is not in the list of supported experiments. Check for typos and consider adding it manually."
                )

        # check variables
        for var in self.vars:
            if var not in self.avail_variables:
                self.logger.warning(f"Variable {var} unknown.")
                raise ValueError(
                    f"Variable {var} is not in the list of supported variables. Check for typos and consider adding it manually."
                )

    # check variables
    def _handle_emission_variables(self, variables: list[str]):
        if self.project == "input4MIPs":
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

    # def _check_models(self):
    #     # check if model, variable, and experiment exists
    #     try:
    #         self.node_link = self.avail_models[self.model]["node_link"]
    #     except KeyError:
    #         if self.model is not None:
    #             self.logger.info(f"WARNING: Model {self.model} unknown.")
    #             raise ValueError(
    #                 "Model {} is not in the list of supported models. Consider adding manually to esgf_server.py".format(
    #                     self.model
    #                 )
    #             )
    #         self.node_link = self.avail_models[self.model]["node_link"]

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

    def _generate_raw_emission_vars(self, variables: list[str]):
        if variables is None:
            # variables = ["tas", "pr", "SO2_em_anthro", "BC_em_anthro"]
            raise ValueError("No variables have been given to the downloader. Variables must be given for downloader.")
        variables = [v.replace(" ", "_").replace("-", "_") for v in variables]
        self.logger.info(f"Cleaned variables : {variables}")
        for v in variables:
            self.vars.append(v)

    #
    # Class functions
    #
    def download_from_model_single_var(  # noqa: C901
        self,
        variable: str,
        experiment: str,
        project: str = "CMIP6",
        default_frequency: str = "mon",
        preferred_version: str = "latest",
        default_grid_label: str = "gn",
    ):
        """
        Function handling the download of a single variable-experiment pair that is associated with a model's output
        (CMIP data).

        Args:
            variable: variable ID
            experiment: experiment ID
            project: umbrella project id e.g. CMIPx
            default_frequency: default frequency to download
            preferred_version: data upload version, if 'latest', the newest version will get selected always
            default_grid_label: default gridding method in which the data is provided
        """

        facets = (
            "project,experiment_id,source_id,variable,frequency,variant_label,variable, nominal_resolution, "
            "version, grid_label, experiment_id"
        )

        self.logger.info("Using download_from_model_single_var() function")

        # Search context is sensitive to order and sequence, which is why
        # it's done in different steps instead of putting everything in `new_context`
        conn = SearchConnection(url=self.node_link, distrib=False)
        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            source_id=self.model,
            variable=variable,
            facets=facets,
        )

        ctx = handle_base_search_constraints(ctx, default_frequency, default_grid_label)

        variants = list(ctx.facet_counts["variant_label"])

        if len(variants) < 1:
            self.logger.info(
                "No items were found for this request. Please check on the esgf server if the combination of your model/scenarios/variables exists."
            )
            raise ValueError(
                "Downloader did not find any items on esgf for your request with: Project {project}, Experiment {experiment}, Model {self.model}, Variable {variable}."
            )

        self.logger.info(f"Available variants : {variants}\n")
        self.logger.info(f"Length : {len(variants)}")

        # TODO refactor logic of if/else
        if not self.ensemble_members:
            if self.max_ensemble_members > len(variants):
                self.logger.info("Less ensemble members available than maximum number desired. Including all variants.")
                ensemble_member_final_list = variants
            else:
                self.logger.info(
                    f"{len(variants)} ensemble members available than desired (max {self.max_ensemble_members}. "
                    f"Choosing only the first {self.max_ensemble_members}.)."
                )
                ensemble_member_final_list = variants[: self.max_ensemble_members]
        else:
            self.logger.info(f"Desired list of ensemble members given: {self.ensemble_members}")
            ensemble_member_final_list = list(set(variants) & set(self.ensemble_members))
            if len(ensemble_member_final_list) == 0:
                self.logger.info("WARNING: no overlap between available and desired ensemble members!")
                self.logger.info("Skipping.")
                return None

        for ensemble_member in ensemble_member_final_list:
            self.logger.info(f"Ensembles member: {ensemble_member}")
            ctx_ensemble = ctx.constrain(variant_label=ensemble_member)

            version = get_upload_version(context=ctx, preferred_version=preferred_version)
            if version:
                ctx_ensemble = ctx_ensemble.constrain(version=version)

            results = ctx_ensemble.search()

            self.logger.info(f"Result len {len(results)}")

            download_model_variable(
                model_id=self.model, search_results=results, variable=variable, base_path=self.data_dir
            )

    def download_raw_input_single_var(  # noqa: C901
        self,
        variable: str,
        project: str = "input4mips",
        institution_id: str = "PNNL-JGCRI",
        default_frequency: str = "mon",
        preferred_version: str = "latest",
        default_grid_label: str = "gn",
    ):
        """
        Function handling the download of all input4mips data associated with a single variable.

        Args:
            variable: variable ID
            project: umbrella project, here "input4mips"
            institution_id: id of the institution that provides the data
            default_frequency: default frequency to download
            preferred_version: data upload version, if 'latest', the newest version will get selected always
            default_grid_label: default gridding method in which the data is provided
        """
        self.logger.info("Using download_raw_input_single_var() function")

        facets = "project,frequency,variable,nominal_resolution,version,target_mip,grid_label"

        # Search context is sensitive to order and sequence, which is why
        # it's done in different steps instead of putting everything in `new_context`
        conn = SearchConnection(url=self.node_link, distrib=False)
        ctx = conn.new_context(
            project=project,
            variable=variable,
            institution_id=institution_id,
            facets=facets,
        )
        ctx = handle_base_search_constraints(ctx, default_frequency, default_grid_label)

        mips_targets = list(ctx.facet_counts["target_mip"])
        self.logger.info(f"Available target mips: {mips_targets}")

        for target in mips_targets:
            ctx_target = ctx.constrain(target_mip=target)
            version = get_upload_version(context=ctx_target, preferred_version=preferred_version)
            if version:
                ctx_target = ctx_target.constrain(version=version)

            results = ctx_target.search()
            self.logger.info(f"Result len  {len(results)}")
            if len(results) > 0:
                download_raw_input_variable(
                    institution_id=institution_id, search_results=results, variable=variable, base_path=self.data_dir
                )

    def download_meta_historic_biomassburning_single_var(
        self,
        variable: str,
        institution_id: str,
        project: str = "input4mips",
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
        facets = "nominal_resolution,version"

        # Search context is sensitive to order and sequence, which is why
        # it's done in different steps instead of putting everything in `new_context`
        conn = SearchConnection(url=self.node_link, distrib=False)
        ctx = conn.new_context(
            project=project,
            variable=variable_search,
            variable_id=variable_id,
            institution_id=institution_id,
            target_mip="CMIP",
            facets=facets,
        )
        ctx = handle_base_search_constraints(ctx, default_frequency, default_grid_label)

        version = get_upload_version(context=ctx, preferred_version=preferred_version)
        if version:
            ctx = ctx.constrain(version=version)

        results = ctx.search()
        self.logger.info(f"Result len  {len(results)}")

        result_list = [r.file_context().search() for r in results]
        self.logger.info(f"List of results :\n{result_list}")

        download_metadata_variable(
            institution_id=institution_id, search_results=results, variable=variable, base_path=self.data_dir
        )

    def download_from_model(self):
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
        for variable in self.vars:
            self.logger.info(f"Downloading data for variable: {variable}")
            for experiment in self.experiments:
                if experiment not in self.avail_experiments:
                    self.logger.info(
                        f"Chosen experiment {experiment} not supported. All supported experiments: "
                        f"{self.avail_experiments}. Skipping."
                    )
                    continue
                self.logger.info(f"Downloading data for experiment: {experiment}")
                self.download_from_model_single_var(project=self.project, variable=variable, experiment=experiment)

    def download_raw_input(self):
        """
        Function handling the download of all variables that are associated with a model's input (input4mips).

        Searches for all files associated with the respected variables that the downloader was initialized with.
        A search connection is established and the search is iteratively constraint to meet all specifications.
        Data is downloaded and stored in a separate file for each year. The default format is netCDF4.

        Resulting hierarchy:

        `input4mips/experiment/variable/nominal_resolution/frequency/year.nc`

        If the constraints cannot be met, the default behaviour for the downloader is to select first other
        available value.
        """
        for variable in self.vars:
            institution_id = "PNNL-JGCRI"
            if variable.endswith("openburning"):
                institution_id = "IAMC"
            self.logger.info(f"Downloading data for variable: {variable}")
            self.download_raw_input_single_var(variable=variable, institution_id=institution_id)

        if self.download_biomass_burning & ("historical" in self.experiments):
            for variable in self.biomass_vars:
                self.logger.info(f"Downloading biomassburing data for variable: {variable}")
                self.download_raw_input_single_var(variable=variable, institution_id="VUA")

        if self.download_metafiles:
            for variable in self.meta_vars_percentage:
                # percentage are historic and have no scenarios
                self.logger.info(f"Downloading meta percentage data for variable: {variable}")
                self.download_meta_historic_biomassburning_single_var(variable=variable, institution_id="VUA")
            for variable in self.meta_vars_share:
                self.logger.info(f"Downloading meta openburning share data for variable: {variable}")
                self.download_raw_input_single_var(variable=variable, institution_id="IAMC")


def download_from_config_file(config_file: Union[str, pathlib.Path], logger: logging.Logger = LOGGER):
    """
    This function downloads variables automatically from input config file
    Args:
        config_file: Path to a configuration yaml file
        logger: Logging instance
    """
    if isinstance(config_file, str):
        config_file = pathlib.Path(config_file)
    config_dict = get_yaml_config(config_file)

    downloader_factory = {
        INPUT4MIPS: {"configs": create_input4mips_downloader_config_from_file, "downloader": Input4MipsDownloader},
        CMIP6: {"configs": create_cmip6_downloader_config_from_file, "downloader": CMIP6Downloader},
    }

    verified_config_keys = []
    for config_key in config_dict:
        verified_key = match_project_key(input_key=config_key, key_list=AVAILABLE_CONFIGS)
        if verified_key:
            verified_config_keys.append(verified_key)

    for config_key in verified_config_keys:
        configs = downloader_factory[config_key]["configs"](config_file=config_file)
        downloader = downloader_factory[config_key]["downloader"](config=configs)
        downloader.download()
