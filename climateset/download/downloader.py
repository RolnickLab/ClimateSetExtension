import logging
import pathlib
from typing import Union

import pandas as pd
from pyesgf.search import SearchConnection

from climateset import RAW_DATA
from climateset.download.constants.data_constants import (
    EMISSIONS_ENDINGS,
    META_ENDINGS_PRC,
    META_ENDINGS_SHAR,
)
from climateset.download.constants.esgf_server import (
    MODEL_SOURCES,
    SUPPORTED_EXPERIMENTS,
    VAR_SOURCE_LOOKUP,
)
from climateset.download.utils import (
    _handle_base_search_constraints,
    download_metadata_variable,
    download_model_variable,
    download_raw_input_variable,
    get_max_ensemble_member_number,
    get_select_model_scenarios,
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
        self.model: str = model
        self.model_node_link: str = ""
        self.model_source_center: str = ""
        if experiments is None:
            experiments = [
                "historical",
                "ssp370",
                "hist-GHG",
                "piControl",
                "ssp434",
                "ssp126",
            ]
        # TODO: have a list of supported experiments before trying to look for them on the node
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

        # Args processing
        selected_scenarios = get_select_model_scenarios()
        self._hande_max_possible_member_number(
            df_model_source=selected_scenarios, max_ensemble_members=max_ensemble_members
        )
        self._handle_variables(
            variables=variables,
        )
        self._handle_model_params()

    #
    # Internal helper functions for class init
    #
    def _hande_max_possible_member_number(self, df_model_source: pd.DataFrame, max_ensemble_members: int):
        max_possible_member_number = get_max_ensemble_member_number(
            df_model_source=df_model_source, experiments=self.experiments, model=self.model
        )
        if max_ensemble_members == -1:
            self.logger.info("Trying to take all ensemble members available.")
            self.max_ensemble_members = max_possible_member_number
            # verify that we have enough members for wanted experiments
            # else choose the smallest available for all
        if max_ensemble_members > max_possible_member_number:
            self.logger.info("Not enough members available. Choosing smallest maximum.")
            self.max_ensemble_members = max_possible_member_number
        self.logger.info(f"Downloading data for {self.max_ensemble_members} members.")

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
        try:
            self.model_node_link = MODEL_SOURCES[self.model]["node_link"]
            self.model_source_center = MODEL_SOURCES[self.model]["center"]
        except KeyError:
            self.model = next(iter(MODEL_SOURCES))
            if self.model is not None:
                self.logger.info(f"WARNING: Model {self.model} unknown. Using default instead.")
                self.logger.info(f"Using : {self.model}")
            self.model_node_link = MODEL_SOURCES[self.model]["node_link"]
            self.model_source_center = MODEL_SOURCES[self.model]["center"]

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
            t = get_keys_from_value(d=VAR_SOURCE_LOOKUP, val=v, logger=self.logger)
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

        variants = list(ctx.facet_counts["variant_label"].keys())

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

        mips_targets = list(ctx.facet_counts["target_mip"].keys())
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

            # files_list = temp_download_path.glob("*.nc")
            #
            # for f in files_list:
            #     experiment = extract_target_mip_exp_name(str(f), target)
            #     self.logger.info(f"Experiment : {experiment}")
            #
            #     # make sure to only download data for wanted scenarios
            #     if experiment in self.experiments:
            #         self.logger.info(f"Saving data for experiment : {experiment}")
            #     else:
            #         self.logger.info(
            #             f"Experiment {experiment} not in wanted experiments ({self.experiments}). Skipping"
            #         )
            #         continue
            #
            #     try:
            #         self.logger.info(f"Opening dataset [{f}]")
            #         with xr.open_dataset(f) as ds:
            #             dataset = ds
            #     except OSError as os_error:
            #         self.logger.error(f"Having problems opening the dateset [{f}]. Original file will not be")
            #         self.logger.error(os_error)
            #         continue
            #
            #     if nominal_resolution == "none":
            #         nominal_resolution = infer_nominal_resolution(dataset, nominal_resolution)
            #
            #     years = np.unique(dataset.time.dt.year.to_numpy())
            #     self.logger.info(f"Data covering years: {years[0]} to {years[-1]}")
            #     year_tag = f"{years[0]}_{years[-1]}"
            #
            #     if variable in self.biomass_vars:
            #         variable = f"{variable}_em_biomassburning"
            #     nominal_resolution = nominal_resolution.strip()
            #     nominal_resolution = nominal_resolution.replace(" ", "_")
            #     # Check whether the specified path exists or not
            #     base_filename = f"{experiment}_{variable}_{nominal_resolution}_{frequency}_{grid_label}_{year_tag}.nc"
            #     if save_to_meta:
            #         # if meta, we have future openburning stuff
            #
            #         out_dir = (
            #             f"future-openburning/{experiment}/{variable.split('_')[0]}/{nominal_resolution}/{frequency}/"
            #         )
            #         out_name = f"future_openburning_{base_filename}"
            #         path = os.path.join(self.meta_dir_parent, out_dir)
            #     else:
            #         out_dir = f"{project}/{experiment}/{variable}/{nominal_resolution}/{frequency}/"
            #         out_name = f"{project}_{base_filename}"
            #         path = os.path.join(self.data_dir_parent, out_dir)
            #
            #     os.makedirs(path, exist_ok=True)
            #     outfile = path + out_name
            #
            #     if (not self.overwrite) and os.path.isfile(outfile):
            #         self.logger.info(f"File {outfile} already exists, skipping.")
            #     else:
            #         self.logger.info("Writing file")
            #         self.logger.info(outfile)
            #         chunk_size = RES_TO_CHUNKSIZE[frequency]
            #         dataset = dataset.chunk({"time": chunk_size})
            #         dataset.to_netcdf(outfile, engine="h5netcdf")

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
        #
        # files_list = temp_download_path.glob("*.nc")
        # self.logger.info(f"List of files downloaded : \n{files_list}")
        #
        # for f in files_list:
        #     # find out chunking dependent on resolution
        #     chunk_size = RES_TO_CHUNKSIZE[frequency]
        #     self.logger.info(f"Chunksize : {chunk_size}")
        #
        #     # replacing spaces for file naming
        #     nominal_resolution = nominal_resolution.replace(" ", "_")
        #
        #     try:
        #         dataset = xr.open_dataset(f, chunks={"time": chunk_size})
        #     except OSError:
        #         self.logger.info("Having problems downloading the dataset. The server might be down. Skipping")
        #         continue
        #
        #     years = np.unique(dataset.time.dt.year.to_numpy())
        #     self.logger.info(f"Data covering years: {years[0]} to {years[-1]}")
        #     year_tag = f"{years[0]}_{years[-1]}"
        #
        #     out_dir = f"historic-biomassburning/{variable_save}/{nominal_resolution}/{frequency}/"
        #
        #     # Check whether the specified path exists or not
        #     path = os.path.join(self.meta_dir_parent, out_dir)
        #     os.makedirs(path, exist_ok=True)
        #
        #     base_file_name = f"{variable}_{nominal_resolution}_{frequency}_{grid_label}_{year_tag}.nc"
        #     outfile = path + base_file_name
        #
        #     if (not self.overwrite) and os.path.isfile(outfile):
        #         self.logger.info(f"File {outfile} already exists, skipping.")
        #     else:
        #         self.logger.info("Writing file")
        #         self.logger.info(outfile)
        #         dataset = dataset.chunk({"time": chunk_size})
        #         dataset.to_netcdf(outfile, engine="h5netcdf")

        # for i, files in enumerate(files_list):
        #     file_names = [files[i].opendap_url for i in range(len(files))]
        #     self.logger.info(f"File {i} names : {file_names}")
        #
        #     # find out chunking dependent on resolution
        #     chunksize = RES_TO_CHUNKSIZE[frequency]
        #     self.logger.info(f"Chunksize : {chunksize}")
        #
        #     # replacing spaces for file naming
        #     nominal_resolution = nominal_resolution.replace(" ", "_")
        #
        #     for f in file_names:
        #         try:
        #             ds = xr.open_dataset(f, chunks={"time": chunksize})
        #         except OSError:
        #             self.logger.info("Having problems downloading the dataset. The server might be down. Skipping")
        #             continue
        #
        #         years = np.unique(ds.time.dt.year.to_numpy())
        #         self.logger.info(f"Data covering years: {years[0]} to {years[-1]}")
        #
        #         for y in years:
        #             y = str(y)
        #             out_dir = f"historic-biomassburning/{variable_save}/{nominal_resolution}/{frequency}/{y}/"
        #
        #             # Check whether the specified path exists or not
        #             path = os.path.join(self.meta_dir_parent, out_dir)
        #             os.makedirs(path, exist_ok=True)
        #
        #             out_name = f"{variable}_{nominal_resolution}_{frequency}_{grid_label}_{y}.nc"
        #             outfile = path + out_name
        #
        #             if (not self.overwrite) and os.path.isfile(outfile):
        #                 self.logger.info(f"File {outfile} already exists, skipping.")
        #             else:
        #                 self.logger.info(f"Selecting specific year : {y}")
        #                 ds_y = ds.sel(time=y)
        #                 self.logger.info(ds_y)
        #
        #                 self.logger.info("Writing file")
        #                 self.logger.info(outfile)
        #                 ds_y.to_netcdf(outfile)

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
                if experiment in SUPPORTED_EXPERIMENTS:
                    self.logger.info(f"Downloading data for experiment: {experiment}")
                    self.download_from_model_single_var(variable=variable, experiment=experiment)
                else:
                    self.logger.info(
                        f"Chosen experiment {experiment} not supported. All supported experiments: "
                        f"{SUPPORTED_EXPERIMENTS}. Skipping."
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

    Args:
        config: Can be a dictionary containing configurations or a path to a configuration yaml file
        logger: Logging instance
    """
    if not isinstance(config, dict):
        if isinstance(config, str):
            config = pathlib.Path(config)
        config = get_yaml_config(config)
    try:
        models = config["models"]
    except KeyError as e:
        logger.warning(f"Caught the following exception but continuing : {e}")
        logger.info("No climate models specified. Assuming only input4mips data should be downloaded.")
        models = [None]
    downloader_kwargs = config["downloader_kwargs"]
    logger.info(f"Downloader kwargs : {downloader_kwargs}")
    for m in models:
        downloader = Downloader(model=m, **downloader_kwargs, logger=logger)
        downloader.download_raw_input()
        if m is not None:
            downloader.download_from_model()
