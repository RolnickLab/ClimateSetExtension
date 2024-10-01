import argparse
from typing import Union

from pyesgf.search import SearchConnection

from climateset import META_DATA, RAW_DATA
from climateset.download.constants.data_constants import (
    DATA_CSV,
    EMISSIONS_ENDINGS,
    META_ENDINGS_PRC,
    META_ENDINGS_SHAR,
)
from climateset.download.constants.esgf_server import (  # RES_TO_CHUNKSIZE,
    MODEL_SOURCES,
    SUPPORTED_EXPERIMENTS,
    VAR_SOURCE_LOOKUP,
)
from climateset.download.utils import (  # extract_target_mip_exp_name,; filter_download_script,; infer_nominal_resolution,
    download_metadata_variable,
    download_model_variable,
    download_raw_input_variable,
    get_frequency,
    get_grid_label,
    get_max_ensemble_member_number,
    get_nominal_resolution,
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
        model: Union[str, None] = "NorESM2-LM",  # defaul as in ClimateBench
        experiments: list[str] = None,  # sub-selection of ClimateBench defaul
        variables: list[str] = None,
        data_dir: str = RAW_DATA,
        meta_dir: str = META_DATA,
        max_ensemble_members: int = 10,  # if -1 take all
        ensemble_members: list[str] = None,  # preferred ensemble members used, if None not considered
        overwrite=False,  # flag if files should be overwritten
        download_biomassburning=True,  # get biomassburning data for input4mips
        download_metafiles=True,  # get input4mips meta files
        plain_emission_vars=True,  # specifies if plain variabsle for emissions data are given and rest is inferred
        # or if variables are specified
        start_year=None,
        end_year=None,
        logger=LOGGER,
    ):
        """
        Init method for the Downloader.

        Args:
            model (str): Id of the model from which output should be downloaded. A list of all supported model ids can
                be found in parameters.constants.MODEL_SOURCES. Model data only.
            experiments ([str]):  List of simulations from which data should be downloaded. Model data only.
            experiments ([str]): List of variables for which data should be downloaded. Both model and raw data.
            data_dir: (str): Relative or absolute path to the directory where data should be stored. Will be created
                if not yet existent.
            meta_dir: (str) Relative or absolute path to the directory where the meta data should be sored. Will be
                created if not yet existent.
            overwrite: (bool) Flag if files should be overwritten, if they already exist.
            download_biomassburning: (bool) Flag if biomassburning data for input4mips variables should be downloaded.
            download_metafiles: (bool) Flag if metafiles for input4mips variables should be downloaded.
        """
        self.logger = logger
        self.model = model
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
        self.experiments = experiments
        # assign vars to either target or raw source
        self.raw_vars = []
        self.model_vars = []

        self.ensemble_members = ensemble_members
        self.max_ensemble_members = max_ensemble_members
        self.year_max = 2100

        self.overwrite = overwrite

        # csv with supported models and sources
        df_model_source = DATA_CSV
        self.start_year = start_year
        self.end_year = end_year

        max_possible_member_number = get_max_ensemble_member_number(df_model_source, experiments, model)
        # if taking all ensemble members

        if max_ensemble_members == -1:
            self.logger.info("Trying to take all ensemble members available.")
            self.max_ensemble_members = max_possible_member_number
            # verify that we have enough members for wanted experiments
            # else choose the smallest available for all
        if max_ensemble_members > max_possible_member_number:
            self.logger.info("Not enough members available. Choosing smallest maximum.")
            self.max_ensemble_members = max_possible_member_number
        self.logger.info(f"Downloading data for {self.max_ensemble_members} members.")

        # determine if we are on a slurm cluster cluster = "none" if "SLURM_TMPDIR" in
        # os.environ: cluster = "slurm"
        #
        # if cluster == "slurm":
        #     data_dir=f"{os.environ['SLURM_TMPDIR']}/causalpaca/data/"
        # else:
        #     data_dir = str(ROOT_DIR) + "/data"

        self._generate_variables(variables)

        self._generate_plain_emission_vars(download_biomassburning, plain_emission_vars)

        self.download_metafiles = download_metafiles
        self.download_biomass_burning = download_biomassburning

        self.logger.info(f"Raw variables to download: {self.raw_vars}")
        self.logger.info(f"Model predicted vars to download: {self.model_vars}")

        if self.download_biomass_burning:
            self.logger.info(f"Download biomass burning vars: {self.biomass_vars}")
        if self.download_metafiles:
            self.logger.info("Downloading meta vars:")
            self.logger.info(self.meta_vars_percentage)
            self.logger.info(self.meta_vars_share)

        try:
            self.model_node_link = MODEL_SOURCES[self.model]["node_link"]
            self.model_source_center = MODEL_SOURCES[self.model]["center"]
        except KeyError:
            self.model = next(iter(MODEL_SOURCES))
            if model is not None:
                self.logger.info(f"WARNING: Model {self.model} unknown. Using default instead.")
                self.logger.info(f"Using : {self.model}")
            # else None but we still need the links
            self.model_node_link = MODEL_SOURCES[self.model]["node_link"]
            self.model_source_center = MODEL_SOURCES[self.model]["center"]

        self.data_dir_parent = data_dir
        self.meta_dir_parent = meta_dir

    def _generate_plain_emission_vars(self, download_biomassburning, plain_emission_vars):
        if plain_emission_vars:
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

            # create var names from emissions
            self.raw_vars = [v + e for v in self.raw_vars for e in EMISSIONS_ENDINGS]
            # be careful with CO2
            if "CO2_em_openburning" in self.raw_vars:
                self.raw_vars.remove("CO2_em_openburning")
        else:
            # get plain input4mips vars = biomass vars for historical
            self.biomass_vars = list(set([v.split("_")[0] for v in self.raw_vars]))
            # remove biomass vars from normal raw vars list
            for b in self.biomass_vars:
                try:
                    self.raw_vars.remove(b)
                except Exception as error:
                    self.logger.warning(f"Caught the following exception but continuing : {error}")

            self.download_biomass_burning = download_biomassburning

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

    def _generate_variables(self, variables):
        if variables is None:
            variables = ["tas", "pr", "SO2_em_anthro", "BC_em_anthro"]
        # take care of var mistype (node takes no spaces or '-' only '_')
        variables = [v.replace(" ", "_").replace("-", "_") for v in variables]
        self.logger.info(f"Cleaned variables : {variables}")
        for v in variables:
            t = get_keys_from_value(VAR_SOURCE_LOOKUP, v, self.logger)
            if t == "model":
                self.model_vars.append(v)
            elif t == "raw":
                self.raw_vars.append(v)

            else:
                self.logger.info(f"WARNING: unknown source type for var {v}. Not supported. Skipping.")

    # TODO Fix complexity issue
    # TODO Refactor download to use wget scripts instead of Opendap, and setup test case
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
            variable (str): variable Id
            experiment (str): experiment Id
            project (str): umbrella project id e.g. CMIPx
            default_frequency (str): default frequency to download
            preferred_version (str): data upload version, if 'latest', the newest version will get selected always
            default_grid_label (str): default gridding method in which the data is provided
        """
        conn = SearchConnection(self.model_node_link, distrib=False)

        facets = (
            "project,experiment_id,source_id,variable,frequency,variant_label,variable, nominal_resolution, "
            "version, grid_label, experiment_id"
        )
        """"

        # extracting available facets
        ctx = conn.new_context(project=project, source_id=self.model)
        available_facets=ctx.facet_counts
        for k in available_facets.keys():
            self.logger.info(f"\n facet {k}")
            vs=[str(v) for v in available_facets[k].keys()]
            self.logger.info(vs)
        raise RuntimeError

        """
        self.logger.info("Using download_from_model_single_var() function")

        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            source_id=self.model,
            variable=variable,
            facets=facets,
        )

        # dealing with grid labels
        grid_label = get_grid_label(ctx, default_grid_label)
        if grid_label:
            ctx = ctx.constrain(grid_label=grid_label)

        try:
            nominal_resolutions = list(ctx.facet_counts["nominal_resolution"].keys())
            self.logger.info(f"Available nominal resolution : {nominal_resolutions}")
            #  deal with multipl nom resolutions availabe
            if len(nominal_resolutions) > 1:
                self.logger.info(
                    "Multiple nominal resolutions exist, choosing smallest_nominal resolution (trying), "
                    "please do a check up"
                )

            nominal_resolution = nominal_resolutions[-1]
            self.logger.info(f"Choosing nominal resolution : {nominal_resolution}")
        except IndexError:
            self.logger.info("No nominal resolution")

        # dealing with frequencies
        self.logger.info(f"Available frequencies : {ctx.facet_counts['frequency'].keys()}")
        frequency = "mon"  # list(ctx.facet_counts['frequency'].keys())[-1]
        self.logger.info(f"choosing frequency : {frequency}")

        ctx_origin = ctx.constrain(frequency=frequency, nominal_resolution=nominal_resolution)

        variants = list(ctx.facet_counts["variant_label"].keys())

        self.logger.info(f"Available variants : {variants}\n")
        self.logger.info(f"Length : {len(variants)}")

        if self.ensemble_members is None:
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

        for i, ensemble_member in enumerate(ensemble_member_final_list):
            self.logger.info(f"Ensembles member: {ensemble_member}")
            ctx = ctx_origin.constrain(variant_label=ensemble_member)
            version = get_upload_version(ctx, preferred_version)

            if version:
                ctx = ctx.constrain(version=version)

            results = ctx.search()

            self.logger.info(f"Result len {len(results)}")

            download_model_variable(self.model, results, variable)

    # TODO: test, improve and cleanup download part
    def download_meta_historic_biomassburning_single_var(
        self,
        variable: str,
        institution_id: str,
        project="input4mips",
        frequency="mon",
        preferred_version="latest",
        grid_label="gn",
    ):
        """
        Function handling the download of all meta data associated with a single input4mips variable.

        Args:
            variable (str): variable Id
            project (str): umbrella project
            institution_id (str): id of the institution that provides the data
            frequency (str): default frequency to download
            preferred_version (str): data upload version, if 'latest', the newest version will get selected always
            grid_label (str): default gridding method in which the data is provided
        """
        variable_id = variable.replace("_", "-")
        # variable_save = variable.split("_")[0]
        variable_search = f"percentage_{variable_id.replace('-', '_').split('_')[-1]}"
        self.logger.info(variable, variable_id, institution_id)
        conn = SearchConnection(self.model_node_link, distrib=False)
        facets = "nominal_resolution,version"
        ctx = conn.new_context(
            project=project,
            variable=variable_search,
            variable_id=variable_id,
            institution_id=institution_id,
            target_mip="CMIP",
            grid_label=grid_label,
            frequency=frequency,
            facets=facets,
        )

        # choose nominal resolution if existent
        try:
            nominal_resolution = get_nominal_resolution(ctx)
        except IndexError:
            self.logger.info("No nominal resolution")
            nominal_resolution = "none"
        ctx = ctx.constrain(nominal_resolution=nominal_resolution)

        ctx_origin_v = ctx

        version = get_upload_version(ctx, preferred_version)

        if preferred_version:
            ctx = ctx_origin_v.constrain(version=version)

        results = ctx.search()

        self.logger.info(f"Result len  {len(results)}")

        result_list = [r.file_context().search() for r in results]
        self.logger.info(f"List of results :\n{result_list}")

        download_metadata_variable(institution_id, results, variable)
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

    # TODO Fix complexity issue
    def download_raw_input_single_var(  # noqa: C901
        self,
        variable,
        project="input4mips",
        institution_id="PNNL-JGCRI",  # make sure that we have the correct data
        default_frequency="mon",
        preferred_version="latest",
        default_grid_label="gn",
        save_to_meta=False,
    ):
        """
        Function handling the download of a all input4mips data associated with a single variable.

        Args:
            variable (str): variable Id
            project (str): umbrella project, here "input4mips"
            institution_id (str): id of the institution that provides the data
            default_frequency (str): default frequency to download
            preferred_version (str): data upload version, if 'latest', the newest version will get selected always
            default_grid_label (str): default gridding method in which the data is provided
            save_to_meta (bool): if data should be saved to the meta folder instead of the input4mips folder
        """
        self.logger.info("Using download_raw_input_single_var() function")
        conn = SearchConnection(self.model_node_link, distrib=False)

        facets = "project,frequency,variable,nominal_resolution,version,target_mip,grid_label"

        # basic constraining (project, var, institution)

        ctx = conn.new_context(
            project=project,
            variable=variable,
            institution_id=institution_id,
            facets=facets,
        )
        grid_label = get_grid_label(ctx, default_grid_label)
        if grid_label:
            ctx = ctx.constrain(grid_label=grid_label)

        nominal_resolution = get_nominal_resolution(ctx)
        ctx = ctx.constrain(nominal_resolution=nominal_resolution)

        frequency = get_frequency(ctx, default_frequency)
        if frequency:
            ctx = ctx.constrain(frequency=frequency)

        mips_targets = list(ctx.facet_counts["target_mip"].keys())
        self.logger.info(f"Available target mips: {mips_targets}")

        for target in mips_targets:
            ctx_target = ctx.constrain(target_mip=target)
            version = get_upload_version(ctx_target, preferred_version)
            if version:
                ctx_target = ctx_target.constrain(version=version)

            results = ctx_target.search()
            self.logger.info(f"Result len  {len(results)}")
            if len(results) > 0:
                download_raw_input_variable(institution_id, results, variable)

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
            #     base_file_name = f"{experiment}_{variable}_{nominal_resolution}_{frequency}_{grid_label}_{year_tag}.nc"
            #     if save_to_meta:
            #         # if meta, we have future openburning stuff
            #
            #         out_dir = (
            #             f"future-openburning/{experiment}/{variable.split('_')[0]}/{nominal_resolution}/{frequency}/"
            #         )
            #         out_name = f"future_openburning_{base_file_name}"
            #         path = os.path.join(self.meta_dir_parent, out_dir)
            #     else:
            #         out_dir = f"{project}/{experiment}/{variable}/{nominal_resolution}/{frequency}/"
            #         out_name = f"{project}_{base_file_name}"
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

        # iterate over respective vars
        for v in self.model_vars:
            self.logger.info(f"Downloading data for variable: {v}")
            # iterate over experiments
            for e in self.experiments:
                # check if experiment is availabe
                if e in SUPPORTED_EXPERIMENTS:
                    self.logger.info(f"Downloading data for experiment: {e}")
                    self.download_from_model_single_var(v, e)
                else:
                    self.logger.info(
                        f"Chosen experiment {e} not supported. All supported experiments: {SUPPORTED_EXPERIMENTS}. "
                        "Skipping."
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
        for v in self.raw_vars:
            if v.endswith("openburning"):
                institution_id = "IAMC"
            else:
                institution_id = "PNNL-JGCRI"
            self.logger.info(f"Downloading data for variable: {v}")
            self.download_raw_input_single_var(v, institution_id=institution_id)

        # if download historical + openburning
        if self.download_biomass_burning & ("historical" in self.experiments):
            for v in self.biomass_vars:
                self.logger.info(f"Downloading biomassburing data for variable: {v}")
                self.download_raw_input_single_var(v, institution_id="VUA")

        if self.download_metafiles:
            for v in self.meta_vars_percentage:
                # percentage are historic and have no scenarios
                self.logger.info(f"Downloading meta percentage data for variable: {v}")
                self.download_meta_historic_biomassburning_single_var(variable=v, institution_id="VUA")
            for v in self.meta_vars_share:
                self.logger.info(f"Downloading meta openburning share data for variable: {v}")
                self.download_raw_input_single_var(v, institution_id="IAMC", save_to_meta=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", help="Path to config file.")
    args = parser.parse_args()
    cli_logger = create_logger("Runtime")

    cfg = get_yaml_config(args.cfg)

    try:
        models = cfg["models"]
    except Exception as error:
        cli_logger.warning(f"Caught the following exception but continuing : {error}")
        cli_logger.info("No climate models specified. Assuming only input4mips data should be downloaded.")
        models = [None]
    downloader_kwargs = cfg["downloader_kwargs"]
    cli_logger.info(f"Downloader kwargs : {downloader_kwargs}")

    # one downloader per climate model
    for m in models:
        downloader = Downloader(model=m, **downloader_kwargs, logger=cli_logger)
        downloader.download_raw_input()
        if m is not None:
            downloader.download_from_model()
