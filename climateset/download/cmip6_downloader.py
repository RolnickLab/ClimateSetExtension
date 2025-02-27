from pyesgf.search import SearchConnection

from climateset.download.abstract_downloader import AbstractDownloader
from climateset.download.constants.esgf import CMIP6
from climateset.download.downloader_config import (
    CMIP6DownloaderConfig,
    create_cmip6_downloader_config_from_file,
)
from climateset.download.utils import (
    download_model_variable,
    get_upload_version,
    handle_base_search_constraints,
)
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

    def download_from_model_single_var(  # noqa: C901
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
        conn = SearchConnection(url=self.config.node_link, distrib=False)

        facets = (
            "project,experiment_id,source_id,variable,frequency,variant_label,variable, nominal_resolution, "
            "version, grid_label, experiment_id"
        )

        self.logger.info("Using download_from_model_single_var() function")

        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            source_id=model,
            variable=variable,
            facets=facets,
        )

        ctx = handle_base_search_constraints(ctx, default_frequency, default_grid_label)

        variants = list(ctx.facet_counts["variant_label"])

        if len(variants) < 1:
            self.logger.info(
                "No items were found for this request. Please check on the esgf server if the combination of your "
                "model/scenarios/variables exists."
            )
            raise ValueError(
                f"Downloader did not find any items on esgf for your request with: Project {project}, "
                f"Experiment {experiment}, Model {model}, Variable {variable}."
            )

        self.logger.info(f"Available variants : {variants}\n")
        self.logger.info(f"Length : {len(variants)}")

        # TODO refactor logic of if/else
        if not self.config.ensemble_members:
            if self.config.max_ensemble_members > len(variants):
                self.logger.info("Less ensemble members available than maximum number desired. Including all variants.")
                ensemble_member_final_list = variants
            else:
                self.logger.info(
                    f"{len(variants)} ensemble members available than desired (max {self.config.max_ensemble_members}. "
                    f"Choosing only the first {self.config.max_ensemble_members}.)."
                )
                ensemble_member_final_list = variants[: self.config.max_ensemble_members]
        else:
            self.logger.info(f"Desired list of ensemble members given: {self.config.ensemble_members}")
            ensemble_member_final_list = list(set(variants) & set(self.config.ensemble_members))
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
                project=CMIP6,
                model_id=model,
                search_results=results,
                variable=variable,
                base_path=self.config.data_dir,
            )


def cmip6_download_from_config(config):
    config_object = create_cmip6_downloader_config_from_file(config)
    downloader = CMIP6Downloader(config=config_object)
    downloader.download()
