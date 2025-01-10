from abstract_downloader import AbstractDownloader
from pyesgf.search import SearchConnection

from climateset.download.utils import (
    _handle_base_search_constraints,
    download_model_variable,
    get_upload_version,
)
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class CMIP6Downloader(AbstractDownloader):
    def __init__(self):
        self.logger = LOGGER

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
