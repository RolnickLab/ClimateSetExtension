import os
import os.path

import numpy as np
import xarray as xr
from pyesgf.search import SearchConnection

from climateset import DATA_DIR
from climateset.download.constants.downloader import SCENARIOS, VARS
from climateset.download.constants.esgf_server import (
    MODEL_SOURCES,
    RES_TO_CHUNKSIZE,
    SUPPORTED_EXPERIMENTS,
    VAR_SOURCE_LOOKUP,
)
from climateset.utils import create_logger, get_keys_from_value

# Question @Charlie (by Julia): Why are we not doing the downloading process
# as it is done here? https://github.com/aditya-grover/climate-learn/blob/main/src/climate_learn/data/download.py
# we could use subprocess.check_call() with the parameter "--no-check-certificate" to circumvent the cookie situation

LOGGER = create_logger(__name__)


class Downloader:
    """
    Class handling the downloading of the data. It communicates with the esgf nodes to search and download the specified data.
    """

    def __init__(
        self,
        model: str = "NorESM2-LM",  # default as in ClimateBench
        experiments: list[str] = None,
        variables: list[str] = None,
        data_dir: str = DATA_DIR,
        max_ensemble_members: int = 10,  # max ensemble members
        ensemble_members: list[str] = None,  # preferred ensemble members used, if None not considered
        overwrite: bool = False,
        logger=LOGGER,
    ):
        """Init method for the Downloader

        Parameters
        ----------
        model :
            Id of the model from which output should be downloaded. A list of all supported model ids can
            be find in parameters.constants.MODEL_SOURCES. Model data only.
        experiments :
            List of simulations from which data should be downloaded. Model data only.
        variables :
            List of variables for which data should be downloaded. Both model and raw data.
        data_dir :
            Relative or absolute path to the directory where data should be stored.
            Will be created if not yet existent.
        """
        self.logger = logger
        self.model = model
        # TODO: Have a list of supported experiments before trying to look for them on the node to reduce computation cost
        self.experiments = experiments
        if self.experiments is None:
            self.experiments = [
                "historical",
                "ssp370",
                "hist-GHG",
                "piControl",
                "ssp434",
                "ssp126",
            ]
        # assign vars to either target or raw source
        self.raw_vars = []
        self.model_vars = []
        self.max_ensemble_members = max_ensemble_members
        self.ensemble_members = ensemble_members

        # take care of var mistype (node takes no spaces or '-' only '_')
        if variables is None:
            variables = ["tas", "pr", "SO2", "BC"]
        cleaned_vars = [v.replace(" ", "_").replace("-", "_") for v in variables]
        self.logger.info("Cleaned vars", cleaned_vars)
        for v in cleaned_vars:
            t = get_keys_from_value(VAR_SOURCE_LOOKUP, v)
            if t == "model":
                self.model_vars.append(v)
            elif t == "raw":
                self.raw_vars.append(v)

            else:
                self.logger.warning(f"WARNING: unknown source type for var {v}. Not supported. Skipping.")

        self.logger.info(f"Raw variables to download: {self.raw_vars}")
        self.logger.info(f"Model predicted vars to download: {self.model_vars}")

        try:
            self.model_node_link = MODEL_SOURCES[self.model]["node_link"]
            self.model_source_center = MODEL_SOURCES[self.model]["center"]
        except KeyError:
            self.logger.warning(f"WARNING: Model {self.model} unknown. Using default instead.")
            self.model = next(iter(MODEL_SOURCES))
            self.model_node_link = MODEL_SOURCES[self.model]["node_link"]
            self.model_source_center = MODEL_SOURCES[self.model]["center"]
            self.logger.warning("Using:", self.model)
        self.logger.info("model node link:", self.model_node_link)

        # log on Manager
        # self.lm = LogonManager()
        # self.lm.logoff()
        # self.lm.logon_with_openid(openid=OPENID, password=PASSWORD, bootstrap=True)
        # self.logger.info("Log In to Node:", self.lm.is_logged_on())

        self.data_dir_parent = data_dir
        self.overwrite = overwrite

        # TODO: create folder hierachy / check if existent make new if not
        # TODO: more checkups?

    def download_from_model_single_var(
        self,
        variable: str,
        experiment: str,
        project: str = "CMIP6",
        default_frequency: str = "mon",
        default_version: str = "latest",
        default_grid_label: str = "gn",
        overwrite: bool = None,
    ):
        """Function handling the download of a single variable-experiment pair that is associated wtih a model's output (CMIP data).

        Parameters
        ----------
            variable :
                variable Id
            experiment (str):
                experiment Id
            project :
                umbrella project id e.g. CMIPx
            default_frequency :
                default frequency to download
            default_version :
                data upload version, if 'latest', the newest version will get selected always
            defaul_grid_label :
                default gridding method in which the data is provided

        """
        conn = SearchConnection(self.model_node_link, distrib=False)

        facets = "project,experiment_id,source_id,variable,frequency,variant_label,variable, nominal_resolution, version, grid_label, experiment_id"

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

        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            source_id=self.model,
            variable=variable,
            facets=facets,
        )

        if overwrite is None:
            overwrite = self.overwrite

        # dealing with grid labels
        grid_labels = list(ctx.facet_counts["grid_label"].keys())
        self.logger.info("Available grid labels:", grid_labels)
        if default_grid_label in grid_labels:
            self.logger.info("Choosing grid:", default_grid_label)
            grid_label = default_grid_label
        else:
            self.logger.info("Default grid label not available.")
            grid_label = grid_labels[0]
            self.logger.info(f"Choosing grid {grid_label} instead.")
        ctx = ctx.constrain(grid_label=grid_label)

        try:
            nominal_resolutions = list(ctx.facet_counts["nominal_resolution"].keys())
            self.logger.info("Available nominal resolution:", nominal_resolutions)
            #  deal with multipl nom resolutions availabe
            if len(nominal_resolutions) > 1:
                self.logger.info(
                    "Multiple nominal resolutions exist, choosing smallest_nominal resolution (trying), please do a check up"
                )

            nominal_resolution = nominal_resolutions[-1]
            self.logger.info("Choosing nominal resolution", nominal_resolution)
        except IndexError:
            self.logger.warning("No nominal resolution")

        # dealing with frequencies
        self.logger.info("Available frequencies: ", ctx.facet_counts["frequency"].keys())
        frequency = "mon"  # list(ctx.facet_counts['frequency'].keys())[-1]
        self.logger.info("choosing frequency: ", frequency)

        ctx_origin = ctx.constrain(frequency=frequency, nominal_resolution=nominal_resolution)

        variants = list(ctx.facet_counts["variant_label"].keys())

        # download_files = {}

        self.logger.info("Available variants:", variants, "\n")

        if self.ensemble_members is None:
            if self.max_ensemble_members > len(variants):
                self.logger.warning(
                    "Less ensemble members available than maximum number desired. Including all variants."
                )
                ensemble_member_final_list = variants
            else:
                self.logger.warning(
                    f"{len(variants)} ensemble members available than desired (max {self.max_ensemble_members}. Choosing only the first {self.max_ensemble_members}.)."
                )
                ensemble_member_final_list = variants[: self.max_ensemble_members]
        else:
            self.logger.info(f"Desired list of ensemble members given: {self.ensemble_members}")
            ensemble_member_final_list = list(set(variants) & set(self.ensemble_members))
            if len(ensemble_member_final_list) == 0:
                self.logger.warning("WARNING: no overlap between available and desired ensemble members!")
                self.logger.warning("Skipping.")
                return None

        for i, ensemble_member in enumerate(ensemble_member_final_list):
            self.logger.info(f"Ensembles member: {ensemble_member}")
            ctx = ctx_origin.constrain(variant_label=ensemble_member)

            # pick a version
            versions = list(ctx.facet_counts["version"].keys())
            self.logger.info("Available versions:", versions)

            if default_version == "latest":
                version = versions[0]
                self.logger.info("Chooosing latetst version:", version)
            else:
                try:
                    version = versions[default_version]
                except KeyError:
                    self.logger.warning(f"Preferred version {default_version} does not exist.")
                    version = versions[0]
                    self.logger.warning(f"Resuming with latest version:{version}")

            ctx = ctx.constrain(version=version)

            result = ctx.search()

            self.logger.info(f"Result len: {len(result)}")

            files_list = [r.file_context().search() for r in result]

            for i, files in enumerate(files_list):
                file_names = [files[i].opendap_url for i in range(len(files))]
                self.logger.info(f"File {i} names: ", file_names)

                # num_files = len(file_names)

                chunksize = RES_TO_CHUNKSIZE[frequency]
                self.logger.info("Chunksize", chunksize)

                nominal_resolution = nominal_resolution.replace(" ", "_")

                for f in file_names:
                    # try to opend datset
                    try:
                        ds = xr.open_dataset(f, chunks={"time": chunksize}, engine="netcdf4")

                    except OSError:
                        self.logger.error("Having problems downloading the dataset. The server might be down. Skipping")
                        continue

                    years = np.unique(ds.time.dt.year.to_numpy())
                    self.logger.info(f"Data covering years: {years[0]} to {years[-1]}")

                    for y in years:
                        y = str(y)
                        out_dir = f"{project}/{self.model}/{ensemble_member}/{experiment}/{variable}/{nominal_resolution}/{frequency}/{y}/"

                        # check if path is existent
                        path = self.data_dir_parent + out_dir
                        isExist = os.path.exists(path)

                        if not isExist:
                            # Create a new directory because it does not exist
                            os.makedirs(path)
                            self.logger.info("The new directory is created!")

                        out_name = f"{project}_{self.model}_{ensemble_member}_{experiment}_{variable}_{nominal_resolution}_{frequency}_{grid_label}_{y}.nc"
                        outfile = path + out_name

                        if (not overwrite) and os.path.isfile(outfile):
                            self.logger.info(f"File {outfile} already exists, skipping.")
                        else:
                            self.logger.info("Selecting specific year", y)
                            ds_y = ds.sel(time=y)
                            self.logger.info(ds_y)
                            self.logger.info("writing file")
                            self.logger.info(outfile)
                            ds_y.to_netcdf(outfile)

    def download_raw_input_single_var(
        self,
        variable,
        project="input4mips",
        institution_id="PNNL-JGCRI",  # make sure that we have the correct data
        default_frequency="mon",
        default_version="latest",
        default_grid_label="gn",
        overwrite: bool = None,
    ):
        """Function handling the download of a all input4mips data associated with a single variable. A

        Parameters
        ----------
            variable :
                variable Id
            project :
                umbrella project, here "input4mips"
            institution_id :
                id of the institution that provides the data
            default_frequency :
                default frequency to download
            default_version :
                data upload version, if 'latest', the newest version will get selected always
            defaul_grid_label :
                default gridding method in which the data is provided

        """
        if overwrite is None:
            overwrite = self.overwrite

        conn = SearchConnection(self.model_node_link, distrib=False)

        facets = "project,frequency,variable,nominal_resolution,version,target_mip,grid_label"

        # basic constraining (project, var, institution)
        ctx = conn.new_context(
            project=project,
            variable=variable,
            institution_id=institution_id,
            facets=facets,
        )

        # dealing with grid labels
        grid_labels = list(ctx.facet_counts["grid_label"].keys())
        self.logger.info("Available grid labels:", grid_labels)

        if default_grid_label in grid_labels:
            self.logger.info("Choosing grid:", default_grid_label)
            grid_label = default_grid_label
        else:
            self.logger.info("Default grid label not available.")
            grid_label = grid_labels[0]
            # grid_label = "gn"
            self.logger.info(f"Choosing grid {grid_label} instead.")
        ctx = ctx.constrain(grid_label=grid_label)

        # choose nominal resolution if existent
        try:
            nominal_resolutions = list(ctx.facet_counts["nominal_resolution"].keys())
            self.logger.info("Available nominal resolution:", nominal_resolutions)

            # deal with mulitple nominal resoulitions, taking smalles one as default
            if len(nominal_resolutions) > 1:
                self.logger.info(
                    "Multiple nominal resolutions exist, choosing smallest_nominal resolution (trying), please do a check up"
                )
            nominal_resolution = nominal_resolutions[0]
            self.logger.info("Choosing nominal resolution", nominal_resolution)
            ctx = ctx.constrain(nominal_resolution=nominal_resolution)

        except IndexError:
            self.logger.warning("No nominal resolution")
            nominal_resolution = "none"

        # choose default frequency if wanted
        frequencies = list(ctx.facet_counts["frequency"].keys())
        self.logger.info("Available frequencies: ", frequencies)

        if default_frequency in frequencies:
            frequency = default_frequency
            self.logger.info("Choosing default frequency", frequency)
        else:
            frequency = frequencies[0]
            self.logger.warning(
                "Default frequency not available, choosing first available one instead: ",
                frequency,
            )
        ctx = ctx.constrain(frequency=frequency)

        # target mip group
        target_mips = list(ctx.facet_counts["target_mip"].keys())
        self.logger.info(f"Available target mips: {target_mips}")
        ctx_origin = ctx

        self.logger.info("\n")
        for t in target_mips:
            self.logger.info(f"Target mip: {t}")
            ctx = ctx_origin.constrain(target_mip=t)

            versions = list(ctx.facet_counts["version"].keys())
            self.logger.info("Available versions", versions)
            ctx_origin_v = ctx

            # deal with different versions
            if default_version == "latest":
                version = versions[0]
                self.logger.info("Chooosing latetst version:", version)
            else:
                try:
                    version = versions[default_version]
                except KeyError:
                    self.logger.warning(f"Preferred version {default_version} does not exist.")
                    version = versions[0]
                    self.logger.warning(f"Resuming with latest version: {version}")

            ctx = ctx_origin_v.constrain(version=version)

            result = ctx.search()

            self.logger.info(f"Result len: {len(result)}")

            files_list = [r.file_context().search() for r in result]

            for i, files in enumerate(files_list):
                file_names = [files[i].opendap_url for i in range(len(files))]
                self.logger.info(f"File {i} names: ", file_names)
                # num_files = len(file_names)

                # find out chunking dependent on resolution
                chunksize = RES_TO_CHUNKSIZE[frequency]
                self.logger.info("Chunksize", chunksize)

                # replacing spaces for file naming
                nominal_resolution = nominal_resolution.replace(" ", "_")

                for f in file_names:
                    experiment = self.extract_target_mip_exp_name(f, t)

                    # make sure to only download data for wanted scenarios
                    if experiment in self.experiments:
                        self.logger.info("Downloading data for experiment:", experiment)
                    else:
                        self.logger.warning(
                            f"Experiment {experiment} not in wanted experiments ({self.experiments}). Skipping"
                        )
                        continue

                    try:
                        ds = xr.open_dataset(f, chunks={"time": chunksize})
                    except OSError:
                        self.logger.error("Having problems downloading the dataset. The server might be down. Skipping")
                        continue

                    years = np.unique(ds.time.dt.year.to_numpy())
                    self.logger.info(f"Data covering years: {years[0]} to {years[-1]}")

                    for y in years:
                        y = str(y)
                        out_dir = f"{project}/{experiment}/{variable}/{nominal_resolution}/{frequency}/{y}/"

                        # Check whether the specified path exists or not
                        path = self.data_dir_parent + out_dir
                        isExist = os.path.exists(path)

                        if not isExist:
                            # Create a new directory because it does not exist
                            os.makedirs(path)
                            self.logger.info("The new directory is created!")

                        out_name = (
                            f"{project}_{experiment}_{variable}_{nominal_resolution}_{frequency}_{grid_label}_{y}.nc"
                        )
                        outfile = path + out_name

                        if (not overwrite) and os.path.isfile(outfile):
                            self.logger.info(f"File {outfile} already exists, skipping.")
                        else:
                            self.logger.info("Selecting specific year ", y)
                            ds_y = ds.sel(time=y)
                            self.logger.info(ds_y)

                            self.logger.info("Writing file")
                            self.logger.info(outfile)
                            ds_y.to_netcdf(outfile)

    def extract_target_mip_exp_name(self, filename: str, target_mip: str):
        """Helper function extracting the target experiment name from a given file name and the target's umbrella MIP.
        supported target mips: "CMIP" "ScenarioMIP", "DAMIP", "AerChemMIP"

        params:
            filename (str): name of the download url to extract the information from
            target_mip (str): name of the umbreall MIP

        """
        year_end = filename.split("_")[-1].split("-")[1].split(".")[0][:4]
        # self.logger.info(f'years from {year_from} to {year_end}')

        if (target_mip == "ScenarioMIP") or (target_mip == "DAMIP"):
            # extract ssp experiment from file name
            experiment = f"ssp{filename.split('ssp')[-1][:3]}"
            if "covid" in filename:
                experiment = f"{experiment}_covid"
        elif (target_mip == "CMIP") & (int(year_end) < 2015):
            experiment = "historical"

        elif target_mip == "AerChemMIP":
            experiment = f"ssp{filename.split('ssp')[-1][:3]}"
            if "lowNTCF" in filename:
                experiment = f"{experiment}_lowNTTCF"

        else:
            self.logger.warning("WARNING: unknown target mip", target_mip)
            experiment = "None"

        return experiment

    def download_from_model(
        self,
        project: str = "CMIP6",
        default_frequency: str = "mon",
        default_version: str = "latest",
        default_grid_label: str = "gn",
    ):
        """
        Function handling the download of all variables that are associated wtih a model's output
        Searches for all filles associated with the respected variables and experiment that the downloader wsa initialized with.

        A search connection is established and the search is iterativeley constraint to meet all specifications.
        Data is downloaded and stored in a seperate file for each year. The default format is netCDF4.
        Resulting hierachy:
            CMIPx
                model_id
                    ensemble_member
                        experiment
                            variable
                                nominal_resolution
                                    frequency
                                        year.nc
        If the constraints cannot be met, per default behaviour for the downloader to select first other available value.


        params:
            project (str): umbrella project id e.g. CMIPx
            default_frequency (str): default frequency to download
            default_version (str): data upload version, if 'latest', the newest version will get selected always
            defaul_grid_label (str): default gridding method in which the data is provided

        """

        # iterate over respective vars
        for v in self.model_vars:
            self.logger.info(f"Downloading data for variable: {v} \n \n ")
            # iterate over experiments
            for e in self.experiments:
                # check if experiment is availabe
                if e in SUPPORTED_EXPERIMENTS:
                    self.logger.info(f"Downloading data for experiment: {e}\n")
                    self.download_from_model_single_var(v, e)
                else:
                    self.logger.error(
                        f"Chosen experiment {e} not supported. All supported experiments: {SUPPORTED_EXPERIMENTS}. "
                        "\n Skipping. \n"
                    )

    def download_raw_input(
        self,
        project="input4mips",
        institution_id="PNNL-JGCRI",  # make sure that we have the correct data
        default_frequency="mon",
        default_version="latest",
        default_grid_label="gn",
    ):
        """
        Function handling the download of all variables that are associated wtih a model's input (input4mips).
        Searches for all filles associated with the respected variables that the downloader was initialized with.
        A search connection is established and the search is iterativeley constraint to meet all specifications.
        Data is downloaded and stored in a seperate file for each year. The default format is netCDF4.
        Resulting hierachy:
            input4mips
                experiment
                    variable
                        nominal_resolution
                            frequency
                                year.nc
        If the constraints cannot be met, the default behaviour for the downloader is to select first other available value.

        params:
            project (str): umbrella project, in this case "input4mips"
            institution_id (str): institution that provided the data
            default_frequency (str): default frequency to download
            default_version (str): data upload version, if 'latest', the newest version will get selected always
            defaul_grid_label (str): default gridding method in which the data is provided

        """
        for v in self.raw_vars:
            self.logger.info(f"Downloading data for variable: {v} \n \n ")
            self.download_raw_input_single_var(v)


if __name__ == "__main__":
    # which case we are handeling
    input4mips = True

    # determine if we are on a slurm cluster
    cluster = "none"
    if "SLURM_TMPDIR" in os.environ:
        cluster = "slurm"

    if cluster == "slurm":
        data_dir = f"{os.environ['SLURM_TMPDIR']}/causalpaca/data/"
    else:
        data_dir = DATA_DIR / "tmp"

    if input4mips:
        experiments = ["ssp126", "ssp245", "ssp370", "ssp585"]
        # Funktioniert: BC_em_anthro, BC_em_AIR_anthro
        variables = ["CO2_em_anthro", "CO2_em_AIR_anthro"]
        # variables = ["BC", "CO2", "CH4", "SO2"]
        # variables=[
        #     "BC_em_anthro", "BC_em_AIR_anthro", "BC", "BC_em_openburning",
        #     "CO2_em_anthro", "CO2_em_AIR_anthro", "CO2", "CO2_em_openburning",
        #     "CH4_em_anthro", "CH4_em_AIR_anthro", "CH4", "CH4_em_openburning",
        #     "SO2_em_anthro", "SO2_em_AIR_anthro", "SO2", "SO2_em_openburning",
        #     ]
        ensemble_members = None
        model = None
        downloader = Downloader(
            experiments=experiments, vars=variables, model=model, data_dir=data_dir, ensemble_members=ensemble_members
        )
        downloader.download_raw_input()

    else:
        variables = VARS
        experiments = SCENARIOS
        model = "CanESM5"
        # variables=["pr", "tas"]

        # experiments=["ssp126", "historical"]
        # model="NorESM2-LM"
        max_ensemble_members = 1
        ensemble_members = ["r1i1p1f1"]

        # determine if we are on a slurm cluster
        cluster = "none"
        if "SLURM_TMPDIR" in os.environ:
            cluster = "slurm"

        if cluster == "slurm":
            data_dir = f"{os.environ['SLURM_TMPDIR']}/causalpaca/data/"
        else:
            data_dir = DATA_DIR / "tmp"

        downloader = Downloader(
            experiments=experiments, vars=variables, model=model, data_dir=data_dir, ensemble_members=ensemble_members
        )
        downloader.download_from_model()
        # downloader.download_raw_input()
