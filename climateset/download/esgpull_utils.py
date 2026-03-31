import asyncio
import contextlib
import logging
import shutil
import uuid
from pathlib import Path
from typing import Generator

from esgpull import Esgpull, Query
from esgpull.models import Options, Selection

from climateset.download.constraints import CMIP6Constraints, Input4MIPsConstraints
from climateset.download.utils import _handle_ensemble_members
from climateset.utils import create_logger

# Configure esgpull Selection to accept additional custom facets
Selection.configure("target_mip", "version", replace=False)

LOGGER = create_logger(__name__)


@contextlib.contextmanager
def isolated_esgpull_context(raw_data_path: Path | str) -> Generator[Esgpull, None, None]:
    """
    Context manager that creates a unique, isolated execution environment for esgpull to avoid file lock collisions and
    pollution of the user's $HOME directory.

    Args:
        raw_data_path: The base path for RAW_DATA where .esgpull_jobs will be created.

    Yields:
        Esgpull: An isolated instance of Esgpull.
    """
    if isinstance(raw_data_path, str):
        raw_data_path = Path(raw_data_path)

    # Create a unique, isolated directory for this esgpull instance
    # using a UUID to prevent collisions between parallel jobs.
    unique_id = uuid.uuid4().hex
    esgpull_jobs_dir = raw_data_path / ".esgpull_jobs"
    isolated_path = esgpull_jobs_dir / unique_id

    # Ensure the parent directory exists
    esgpull_jobs_dir.mkdir(parents=True, exist_ok=True)

    esg = None
    try:
        esg = Esgpull(path=isolated_path, install=True)
        yield esg
    finally:
        # Tear down and safely purge the isolation folder and its SQLite DB
        if isolated_path.exists():
            shutil.rmtree(isolated_path, ignore_errors=True)


def esgpull_search_and_download_esgf_raw_single_var(
    esg: Esgpull,
    variable: str,
    institution_id: str,
    project: str,
    default_grid_label: str,
    default_frequency: str,
    preferred_version: str,
    data_dir: Path,
    distrib: bool = False,
    logger: logging.Logger = LOGGER,
):
    initial_constraints = Input4MIPsConstraints(
        project=project, institution_id=institution_id, variable=variable
    ).to_esgpull_query()

    query = Query(selection=initial_constraints)
    query.options.distrib = distrib

    _apply_facet_fallback(esg, query, "grid_label", default_grid_label, logger)
    _apply_facet_fallback(esg, query, "frequency", default_frequency, logger)

    # Esgpull handles multi-values natively. Fetch targets if any.
    hints = esg.context.hints(query, file=False, facets=["target_mip"])
    if hints and "target_mip" in hints[0] and hints[0]["target_mip"]:
        target_mips = list(hints[0]["target_mip"].keys())
        logger.info(f"Available target mips: {target_mips}")
        query.selection["target_mip"] = target_mips

    _apply_version_fallback(esg, query, preferred_version, logger)

    files = esg.context.search(query, file=True)
    logger.info(f"Result len: {len(files)}")

    dest_dir = Path(data_dir) / f"{project}/raw_input_vars/{institution_id}/{variable}"
    return _download_and_move_files(esg, files, dest_dir, logger)


def esgpull_search_and_download_esgf_biomass_single_var(
    esg: Esgpull,
    variable: str,
    variable_id: str,
    institution_id: str,
    project: str,
    default_grid_label: str,
    default_frequency: str,
    preferred_version: str,
    data_dir: Path,
    distrib: bool = False,
    logger: logging.Logger = LOGGER,
):
    initial_constraints = Input4MIPsConstraints(
        project=project,
        institution_id=institution_id,
        variable=variable,
        variable_id=variable_id,
    ).to_esgpull_query()

    query = Query(selection=initial_constraints)
    query.options.distrib = distrib

    _apply_facet_fallback(esg, query, "grid_label", default_grid_label, logger)
    _apply_facet_fallback(esg, query, "frequency", default_frequency, logger)
    _apply_version_fallback(esg, query, preferred_version, logger)

    files = esg.context.search(query, file=True)
    logger.info(f"Result len: {len(files)}")

    dest_dir = Path(data_dir) / f"{project}/meta_vars/{institution_id}/{variable}"
    logger.info(f"Destination folder: [{dest_dir}]")
    return _download_and_move_files(esg, files, dest_dir, logger)


def esgpull_search_and_download_esgf_model_single_var(
    esg: Esgpull,
    model: str,
    variable: str,
    experiment: str,
    project: str,
    default_grid_label: str,
    default_frequency: str,
    preferred_version: str,
    max_ensemble_members: int,
    ensemble_members: list[str],
    data_dir: Path,
    distrib: bool = False,
    logger: logging.Logger = LOGGER,
):
    cmip_constraints = CMIP6Constraints(
        project=project, experiment_id=experiment, source_id=model, variable=variable
    ).to_esgpull_query()

    query = Query(selection=cmip_constraints)
    query.options.distrib = distrib

    _apply_facet_fallback(esg, query, "frequency", default_frequency, logger)
    _apply_facet_fallback(esg, query, "grid_label", default_grid_label, logger)

    ensemble_member_final_list = _apply_variants_filter(esg, query, max_ensemble_members, ensemble_members, logger)
    if not ensemble_member_final_list:
        logger.info("No items were found for this request.")
        return None

    # Esgpull supports multi-value list queries seamlessly
    query.selection["variant_label"] = ensemble_member_final_list
    _apply_version_fallback(esg, query, preferred_version, logger)

    files = esg.context.search(query, file=True)
    logger.info(f"Result len {len(files)}")

    dest_dir = Path(data_dir) / f"{project}/{model}/{variable}"
    return _download_and_move_files(esg, files, dest_dir, logger)


def _download_and_move_files(esg, files, dest_dir: Path, logger: logging.Logger = LOGGER):
    """Downloads tracked files natively via esgpull (asyncio) and moves them from the isolated cache to the final
    requested target directory."""
    if not files:
        logger.info("No files to download.")
        return []

    unique_files = {}
    for f in files:
        if f.file_id not in unique_files:
            unique_files[f.file_id] = f
    files_to_add = list(unique_files.values())

    logger.info(f"Adding {len(files_to_add)} files to esgpull DB...")
    # Add tracked files to the isolated internal DB queue
    esg.db.add(*files_to_add)

    async def _run_download():
        logger.info("Starting esg.download async...")
        return await esg.download(files_to_add, show_progress=False)

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


def _apply_facet_fallback(
    esg, query: Query, facet_name: str, preferred_value: str | None, logger: logging.Logger = LOGGER
):
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


def _apply_version_fallback(esg, query: Query, preferred_version: str | None, logger: logging.Logger = LOGGER):
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
                    f"Preferred version {preferred_version} does not exist. "
                    f"Resuming with latest {available_versions[0]}"
                )
                query.selection["version"] = [available_versions[0]]


def _apply_variants_filter(
    esg,
    query: Query,
    max_ensemble_members: int,
    ensemble_members: list[str],
    logger: logging.Logger = LOGGER,
) -> list[str]:
    hints = esg.context.hints(query, file=False, facets=["variant_label"])
    if not hints or "variant_label" not in hints[0] or not hints[0]["variant_label"]:
        return []

    variants = list(hints[0]["variant_label"].keys())
    logger.info(f"Available variants : {variants}\nLength : {len(variants)}")

    ensemble_members_list = _handle_ensemble_members(
        variants=variants, ensemble_members=ensemble_members, max_ensemble_members=max_ensemble_members, logger=logger
    )
    return ensemble_members_list
