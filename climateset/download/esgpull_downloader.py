import asyncio
import logging
import shutil
from pathlib import Path

from esgpull.models import Options, Query, Selection

from climateset.download.abstract_downloader import AbstractDownloader
from climateset.download.constraints import CMIP6Constraints, Input4MIPsConstraints
from climateset.download.utils import isolated_esgpull_context
from climateset.utils import create_logger

# Configure esgpull Selection to accept additional custom facets
Selection.configure("target_mip", "version", replace=False)

LOGGER = create_logger(__name__)


def _download_and_move_files(esg, files, dest_dir: Path, logger: logging.Logger):
    """Downloads tracked files natively via esgpull (asyncio) and moves them from the isolated cache to the final
    requested target directory."""
    if not files:
        logger.info("No files to download.")
        return []

    logger.info(f"Adding {len(files)} files to esgpull DB...")
    # Add tracked files to the isolated internal DB queue
    esg.db.add(*files)

    async def _run_download():
        logger.info("Starting esg.download async...")
        return await esg.download(files, show_progress=False)

    logger.info("Executing asyncio.run(_run_download())...")
    # Execute async native download
    downloaded, errors = asyncio.run(_run_download())
    logger.info(f"Download complete. Downloaded: {len(downloaded)}, Errors: {len(errors)}")

    if errors:
        for err in errors:
            logger.error(f"Download error: {err}")

    # Transfer from cache to strictly formatted project tree
    if isinstance(dest_dir, str):
        dest_dir = Path(dest_dir)

    dest_dir.mkdir(parents=True, exist_ok=True)

    moved_files = []
    data_cache_dir = esg.config.paths.data
    if data_cache_dir.exists():
        for nc_file in data_cache_dir.rglob("*.nc"):
            dest_file = dest_dir / nc_file.name
            logger.info(f"Moving {nc_file.name} to {dest_dir}")
            shutil.move(str(nc_file), str(dest_file))
            moved_files.append(dest_file)

    return moved_files


def _apply_facet_fallback(esg, query: Query, facet_name: str, preferred_value: str | None, logger: logging.Logger):
    """Query the available facets and fall back if preferred_value is not found."""
    hints = esg.context.hints(query, file=False, facets=[facet_name])
    if hints and facet_name in hints[0] and hints[0][facet_name]:
        available_facets = list(hints[0][facet_name].keys())
        logger.info(f"Available {facet_name}: {available_facets}")

        if preferred_value and preferred_value in available_facets:
            logger.info(f"Choosing {facet_name}: {preferred_value}")
            query.selection[facet_name] = [preferred_value]
        else:
            if preferred_value:
                logger.warning(f"Preferred {facet_name} '{preferred_value}' not available.")
            chosen = available_facets[0]
            logger.info(f"Choosing {facet_name} {chosen} instead.")
            query.selection[facet_name] = [chosen]
    else:
        logger.warning(f"No {facet_name} found.")


def _apply_version_fallback(esg, query: Query, preferred_version: str | None, logger: logging.Logger):
    if preferred_version == "latest" or preferred_version is None:
        # Use latest=True in options. Since Options is an Enum-backed mapped model,
        # we can recreate it preserving the existing distrib option
        is_distrib = query.options.distrib.name == "true"
        query.options = Options(distrib=is_distrib, latest=True)
        logger.info("Choosing latest version.")
    else:
        hints = esg.context.hints(query, file=False, facets=["version"])
        if hints and "version" in hints[0] and hints[0]["version"]:
            available_versions = list(hints[0]["version"].keys())
            if preferred_version in available_versions:
                query.selection["version"] = [preferred_version]
            else:
                logger.warning(
                    f"Preferred version {preferred_version} does not exist. Resuming with latest {available_versions[0]}"
                )
                query.selection["version"] = [available_versions[0]]


def _apply_variants_filter(
    esg,
    query: Query,
    max_ensemble_members: int,
    ensemble_members: list[str],
    logger: logging.Logger,
) -> list[str]:
    hints = esg.context.hints(query, file=False, facets=["variant_label"])
    if not hints or "variant_label" not in hints[0] or not hints[0]["variant_label"]:
        return []

    variants = list(hints[0]["variant_label"].keys())
    logger.info(f"Available variants : {variants}\nLength : {len(variants)}")

    if not ensemble_members:
        if max_ensemble_members > len(variants):
            logger.info("Less ensemble members available than maximum number desired. Including all variants.")
            return variants
        logger.info(
            f"{len(variants)} ensemble members available, desired (max {max_ensemble_members}). Choosing only the first {max_ensemble_members}."
        )
        return variants[:max_ensemble_members]

    logger.info(f"Desired list of ensemble members given: {ensemble_members}")
    ensemble_member_final_list = list(set(variants) & set(ensemble_members))
    return ensemble_member_final_list


class EsgpullDownloader(AbstractDownloader):
    def __init__(self, config=None, distrib: bool = False):
        self.config = config
        self.distrib = distrib
        self.logger = LOGGER

    def download(self):
        # Dispatch based on config type (or could be an abstract base)
        # Note: EsgpullDownloader executes searches via esgpull.
        # The actual download logic via esgpull is in task 4,
        # so for now `download` can just invoke search to satisfy the interface.
        pass

    def search_and_download_esgf_raw_single_var(
        self,
        variable: str,
        institution_id: str,
        project: str,
        default_grid_label: str,
        default_frequency: str,
        preferred_version: str,
        data_dir: Path | str,
    ):
        with isolated_esgpull_context(data_dir) as esg:
            initial_constraints = Input4MIPsConstraints(
                project=project, institution_id=institution_id, variable=variable
            ).to_esgpull_query()

            query = Query(selection=initial_constraints)
            query.options.distrib = self.distrib

            _apply_facet_fallback(esg, query, "grid_label", default_grid_label, self.logger)
            _apply_facet_fallback(esg, query, "frequency", default_frequency, self.logger)

            # Esgpull handles multi-values natively. Fetch targets if any.
            hints = esg.context.hints(query, file=False, facets=["target_mip"])
            if hints and "target_mip" in hints[0] and hints[0]["target_mip"]:
                target_mips = list(hints[0]["target_mip"].keys())
                self.logger.info(f"Available target mips: {target_mips}")
                query.selection["target_mip"] = target_mips

            _apply_version_fallback(esg, query, preferred_version, self.logger)

            files = esg.context.search(query, file=True)
            self.logger.info(f"Result len: {len(files)}")

            dest_dir = Path(data_dir) / f"{project}/raw_input_vars/{institution_id}/{variable}"
            return _download_and_move_files(esg, files, dest_dir, self.logger)

    def search_and_download_esgf_biomass_single_var(
        self,
        variable: str,
        variable_id: str,
        institution_id: str,
        project: str,
        default_grid_label: str,
        default_frequency: str,
        preferred_version: str,
        data_dir: Path | str,
    ):
        with isolated_esgpull_context(data_dir) as esg:
            initial_constraints = Input4MIPsConstraints(
                project=project,
                institution_id=institution_id,
                variable=variable,
                variable_id=variable_id,
            ).to_esgpull_query()

            query = Query(selection=initial_constraints)
            query.options.distrib = self.distrib

            _apply_facet_fallback(esg, query, "grid_label", default_grid_label, self.logger)
            _apply_facet_fallback(esg, query, "frequency", default_frequency, self.logger)
            _apply_version_fallback(esg, query, preferred_version, self.logger)

            files = esg.context.search(query, file=True)
            self.logger.info(f"Result len: {len(files)}")

            dest_dir = Path(data_dir) / f"{project}/meta_vars/{institution_id}/{variable}"
            return _download_and_move_files(esg, files, dest_dir, self.logger)

    def search_and_download_esgf_model_single_var(
        self,
        model: str,
        variable: str,
        experiment: str,
        project: str,
        default_grid_label: str,
        default_frequency: str,
        preferred_version: str,
        max_ensemble_members: int,
        ensemble_members: list[str],
        data_dir: Path | str,
    ):
        with isolated_esgpull_context(data_dir) as esg:
            cmip_constraints = CMIP6Constraints(
                project=project, experiment_id=experiment, source_id=model, variable=variable
            ).to_esgpull_query()

            query = Query(selection=cmip_constraints)
            query.options.distrib = self.distrib

            _apply_facet_fallback(esg, query, "frequency", default_frequency, self.logger)
            _apply_facet_fallback(esg, query, "grid_label", default_grid_label, self.logger)

            ensemble_member_final_list = _apply_variants_filter(
                esg, query, max_ensemble_members, ensemble_members, self.logger
            )
            if not ensemble_member_final_list:
                self.logger.info("No items were found for this request.")
                return None

            # Esgpull supports multi-value list queries seamlessly
            query.selection["variant_label"] = ensemble_member_final_list
            _apply_version_fallback(esg, query, preferred_version, self.logger)

            files = esg.context.search(query, file=True)
            self.logger.info(f"Result len {len(files)}")

            dest_dir = Path(data_dir) / f"{project}/{model}/{variable}"
            return _download_and_move_files(esg, files, dest_dir, self.logger)
