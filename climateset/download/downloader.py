import logging
import pathlib
from typing import Union

from pyesgf.search import SearchConnection

import climateset.download.constants.cmip6_constants as cmip6_constants
import climateset.download.constants.cmip6plus_constants as cmip6plus_constants
import climateset.download.constants.input4mips_constants as input4mips_constants
from climateset import RAW_DATA
from climateset.download.constants.data_constants import (
    EMISSIONS_ENDINGS,
    META_ENDINGS_PRC,
    META_ENDINGS_SHAR,
)
from climateset.download.utils import (
    _handle_base_search_constraints,
    download_metadata_variable,
    download_model_variable,
    download_raw_input_variable,
    get_upload_version,
)
from climateset.utils import create_logger, get_keys_from_value, get_yaml_config

LOGGER = create_logger(__name__)


class Downloader:
    """
    Class handling the downloading of the data.

    It communicates with the esgf nodes to search and download the specified data.
    """

    # TODO Fix complexity issue
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
        # Args init
        self.logger = logger
        # init global variables depending on project type
        self._init_globs(project)
        self.project: str = project
        self.model: str = model
        self.model_node_link: str = ""
        #  to reduce computation cost
        self.experiments: list[str] = experiments
        self.raw_vars: list[str] = []
        self.model_vars: list[str] = []
        self.biomass_vars: list[str] = []
        self.meta_vars_percentage: list[str] = []
        self.meta_vars_share: list[str] = []
        self.data_dir: Union[str, pathlib.Path] = data_dir
        self.ensemble_members: list[str] = ensemble_members
        self.max_ensemble_members: int = max_ensemble_members
        self.overwrite: bool = overwrite
        self.download_metafiles: bool = download_metafiles
        self.download_biomass_burning: bool = download_biomassburning
        self.use_plain_emission_vars: bool = use_plain_emission_vars
        self.model_node_link = self.NODE_LINK

        # if max ensemble member number is too large --> we are relying on the server to complain?

        self._check_desired_params()

        # Args processing
        self._handle_variables(
            variables=variables,
        )
        # self._handle_model_params()

    # TODO we need to make the downloader an abstract parent class
    # each project needs its own constant file + downloader function, the rest stays the same
    # this function should not be done this way, this is the first naive approach
    def _init_globs(self, project: str):
        """Load globs depending on project."""
        if project == "CMIP6":
            self.MODEL_SOURCES = cmip6_constants.MODEL_SOURCES
            self.SUPPORTED_EXPERIMENTS = cmip6_constants.SUPPORTED_EXPERIMENTS
            self.VAR_SOURCE_LOOKUP = cmip6_constants.VAR_SOURCE_LOOKUP
            self.NODE_LINK = cmip6_constants.NODE_LINK
        elif project == "input4mips":
            self.NODE_LINK = input4mips_constants.NODE_LINK
        elif project == "CMIP6Plus":
            self.MODEL_SOURCES = cmip6plus_constants.MODEL_SOURCES
            self.SUPPORTED_EXPERIMENTS = cmip6plus_constants.SUPPORTED_EXPERIMENTS
            self.VAR_SOURCE_LOOKUP = cmip6plus_constants.VAR_SOURCE_LOOKUP
            self.NODE_LINK = cmip6plus_constants.NODE_LINK
        else:
            self.logger.info(f"Project {project} has not been implemented in the Downloader yet.")
            raise NotImplementedError(f"Project {project} has not been implemented in the downloader.")

    def _check_desired_params(self):
        """Check if the desired params exist."""
        # check model
        if self.model not in self.MODEL_SOURCES:
            self.logger.info(f"WARNING: Model {self.model} unknown.")
            raise ValueError(
                f"Model {self.model} is not in the list of supported models. Consider adding manually to esgf_server.py"
            )

        # check experiments
        # loop over experiments and check for each experiment in the list

        # check variables

    def _handle_variables(self, variables: list[str]):
        self._generate_variables(variables=variables)
        self._generate_plain_emission_vars()
        self.logger.info(f"Raw variables to download: {self.raw_vars}")
        self.logger.info(f"Model predicted vars to download: {self.model_vars}")
        if self.download_biomass_burning:
            self.logger.info(f"Download biomass burning vars: {self.biomass_vars}")
        if self.download_metafiles:
            self.logger.info(f"Downloading meta vars:\n\t{self.meta_vars_percentage}\n\t{self.meta_vars_share}")

    def _handle_model_params(self):
        # check if model, variable, and experiment exists
        try:
            self.model_node_link = self.MODEL_SOURCES[self.model]["node_link"]
        except KeyError:
            if self.model is not None:
                self.logger.info(f"WARNING: Model {self.model} unknown.")
                raise ValueError(
                    "Model {} is not in the list of supported models. Consider adding manually to esgf_server.py".format(
                        self.model
                    )
                )
            self.model_node_link = self.MODEL_SOURCES[self.model]["node_link"]

    def _generate_plain_emission_vars(self):
        if self.use_plain_emission_vars:
            # plain vars are biomass vars
            self.biomass_vars = self.raw_vars
            self.meta_vars_percentage = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in META_ENDINGS_PRC
            ]
            self.meta_vars_share = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in META_ENDINGS_SHAR
            ]

            self.raw_vars = [
                variable + emission_ending for variable in self.raw_vars for emission_ending in EMISSIONS_ENDINGS
            ]
            # be careful with CO2
            if "CO2_em_openburning" in self.raw_vars:
                self.raw_vars.remove("CO2_em_openburning")
        else:
            # get plain input4mips vars = biomass vars for historical
            self.biomass_vars = list({v.split("_")[0] for v in self.raw_vars})
            # remove biomass vars from normal raw vars list
            for b in self.biomass_vars:
                try:
                    self.raw_vars.remove(b)
                except Exception as error:
                    self.logger.warning(f"Caught the following exception but continuing : {error}")

            self.meta_vars_percentage = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in META_ENDINGS_PRC
            ]
            self.meta_vars_share = [
                biomass_var + ending
                for biomass_var in self.biomass_vars
                if biomass_var != "CO2"
                for ending in META_ENDINGS_SHAR
            ]

    def _generate_variables(self, variables: list[str]):
        if variables is None:
            variables = ["tas", "pr", "SO2_em_anthro", "BC_em_anthro"]
        variables = [v.replace(" ", "_").replace("-", "_") for v in variables]
        self.logger.info(f"Cleaned variables : {variables}")
        for v in variables:
            t = get_keys_from_value(d=self.VAR_SOURCE_LOOKUP, val=v, logger=self.logger)
            if t == "model":
                self.model_vars.append(v)
            elif t == "raw":
                self.raw_vars.append(v)

            else:
                self.logger.info(f"WARNING: unknown source type for var {v}. Not supported. Skipping.")

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
        conn = SearchConnection(url=self.model_node_link, distrib=False)

        facets = (
            "project,experiment_id,source_id,variable,frequency,variant_label,variable, nominal_resolution, "
            "version, grid_label, experiment_id"
        )

        self.logger.info("Using download_from_model_single_var() function")

        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            source_id=self.model,
            variable=variable,
            facets=facets,
        )

        ctx = _handle_base_search_constraints(ctx, default_frequency, default_grid_label)

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
        conn = SearchConnection(url=self.model_node_link, distrib=False)

        ctx = conn.new_context(
            project=project,
            variable=variable,
            institution_id=institution_id,
            facets=facets,
        )

        ctx = _handle_base_search_constraints(ctx, default_frequency, default_grid_label)

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
        conn = SearchConnection(url=self.model_node_link, distrib=False)
        facets = "nominal_resolution,version"
        ctx = conn.new_context(
            project=project,
            variable=variable_search,
            variable_id=variable_id,
            institution_id=institution_id,
            target_mip="CMIP",
            facets=facets,
        )

        ctx = _handle_base_search_constraints(ctx, default_frequency, default_grid_label)

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

        for variable in self.model_vars:
            self.logger.info(f"Downloading data for variable: {variable}")
            for experiment in self.experiments:
                if experiment in self.SUPPORTED_EXPERIMENTS:
                    self.logger.info(f"Downloading data for experiment: {experiment}")
                    self.download_from_model_single_var(project=self.project, variable=variable, experiment=experiment)
                else:
                    self.logger.info(
                        f"Chosen experiment {experiment} not supported. All supported experiments: "
                        f"{self.SUPPORTED_EXPERIMENTS}. Skipping."
                    )

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
        for variable in self.raw_vars:
            if variable.endswith("openburning"):
                institution_id = "IAMC"
            else:
                institution_id = "PNNL-JGCRI"
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


def download_from_config_file(config: str, logger: logging.Logger = LOGGER):
    """
    This function downloads variables automatically from input config file
    Args:
        config: Can be a dictionary containing configurations or a path to a configuration yaml file
        logger: Logging instance
    """
    if not isinstance(config, dict):
        if isinstance(config, str):
            config = pathlib.Path(config)
        config = get_yaml_config(config)
    try:
        project = config["project"]
    except KeyError as e:
        logger.warning(
            f"No project specified. Assuming CMIP6 data should be downloaded. Caught the following exception: {e}"
        )
        project = "CMIP6"
    try:
        models = config["models"]
    except KeyError as e:
        logger.warning(f"Caught the following exception but continuing : {e}")
        logger.info("No climate models specified. Assuming only input4mips data should be downloaded.")
        models = [None]
    downloader_kwargs = config["downloader_kwargs"]
    logger.info(f"Downloader kwargs : {downloader_kwargs}")

    # TODO @Francis I think we need to implement an abstract Downloader.
    # Each project should get its own Downloader: CMIPXDownloader, input4mipsDownloader, etc.
    # These classes only need to implement the different downloading functions needed for their specific datasets.
    # Here, I am just doing the naive way with the stuff we have right now:
    if project == "input4mips":
        downloader = Downloader(project=project, model=models, **downloader_kwargs, logger=logger)
        downloader.download_raw_input()
    elif (project == "CMIP6") or (project == "CMIP6Plus"):
        for m in models:
            downloader = Downloader(project=project, model=m, **downloader_kwargs, logger=logger)
            downloader.download_from_model()
    else:
        logger.info(
            f"Project {project} is not supported. Consider implementing your own downloader childclass for this."
        )
        raise ValueError(
            f"Project {project} is not supported. Currently supported projects are: CMIP6, CMIP6Plus, input4mips."
        )
