from abstract_downloader import AbstractDownloader
from pyesgf.search import SearchConnection

from climateset.download.utils import (
    _handle_base_search_constraints,
    download_metadata_variable,
    download_raw_input_variable,
    get_upload_version,
)
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class Input4MipsDownloader(AbstractDownloader):
    def __init__(self):
        self.raw_vars = ""
        self.logger = LOGGER

    def download(self):
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
