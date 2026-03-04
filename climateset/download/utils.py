import contextlib
import logging
import re
import shutil
import subprocess
import time
import uuid
from pathlib import Path
from typing import Generator

import xarray as xr
from esgpull import Esgpull

from climateset import RAW_DATA
from climateset.download.client import SearchClient, SearchSession
from climateset.download.constraints import (
    BaseSearchConstraints,
    CMIP6Constraints,
    Input4MIPsConstraints,
)
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


def extract_target_mip_exp_name(filename: str, target_mip: str, logger: logging.Logger = LOGGER):
    """
    Helper function extracting the target experiment name from a given file name and the target's umbrella MIP.

    supported target mips: "CMIP" "ScenarioMIP", "DAMIP", "AerChemMIP"

    Args:
        filename : name of the download url to extract the information from
        target_mip : name of the umbrella MIP
        logger: Logger instance
    """
    year_end = filename.split("_")[-1].split("-")[1].split(".")[0][:4]
    if target_mip in ["ScenarioMIP", "DAMIP"]:
        experiment = f"ssp{filename.split('ssp')[-1][:3]}"
        if "covid" in filename:
            experiment = f"{experiment}_covid"
    elif target_mip == "CMIP":
        cutoff_year_for_historical = 2015
        if int(year_end) > cutoff_year_for_historical:
            logger.info(f"TARGET MIP : {filename}")
            experiment = f"ssp{filename.split('ssp')[-1][:3]}"
        else:
            experiment = "historical"

    elif target_mip == "AerChemMIP":
        experiment = f"ssp{filename.split('ssp')[-1][:3]}"
        if "lowNTCF" in filename:
            experiment = f"{experiment}_lowNTTCF"

    else:
        logger.info(f"WARNING: unknown target mip : {target_mip}")
        experiment = "None"

    return experiment


def get_nominal_resolution(session, logger: logging.Logger = LOGGER):
    """
    Get nominal resolution from the search session's available facets.

    Args:
        session: climateset.download.client.SearchSession
        logger: Logger instance

    Returns:
        str: Selected nominal resolution
    """
    nominal_resolution = ""
    nominal_resolution_list = session.get_available_facets("nominal_resolution")

    if nominal_resolution_list:
        logger.info(f"Available nominal resolution : {nominal_resolution_list}")
    else:
        logger.warning("No nominal resolution")
        return nominal_resolution

    if len(nominal_resolution_list) > 1:
        logger.warning("Multiple nominal resolutions exist, will try to get smallest resolution.")

    nominal_resolution = nominal_resolution_list[0]
    logger.info(f"Choosing nominal resolution : {nominal_resolution}")
    return nominal_resolution


def infer_nominal_resolution(ds: xr.Dataset, nominal_resolution: str, logger: logging.Logger = LOGGER) -> str | int:
    """
    This method checks if there really is not nominal resolution by trying to compute it from the longitude increment.

    In principle lon and lat should be the same, however, this is just an approximation
    same approximation used by climate modeling centers information is just for
    informing the structure, resolution will be checked in preprocessing.

    Args:
        ds:
        nominal_resolution:
        logger:

    Returns:
    """
    nom_res = nominal_resolution
    try:
        degree = abs(ds.lon[0].item() - ds.lon[1].item())
        nom_res = int(degree * 100)
        logger.info(f"Inferring nominal resolution: {nom_res}")
    except Exception as error:  # pylint: disable=W0718
        logger.warning(f"Caught the following exception but continuing : {error}")
    return nom_res


def filter_download_script(wget_script_content, start_year, end_year):
    """
    Function to modify wget download script from ESGF. Only download files which contain dates between start_year and
    end_year. This is mostly useful for toy dataset creation.

    Args:
        wget_script_content (str): wget script from ESGF
        start_year (str): year parsed from config file
        end_year (str): year parsed from config file
    """
    lines = wget_script_content.split("\n")
    modified_script = []
    in_section = False
    finished = False
    for line in lines:
        if in_section and not finished:
            if re.match(r"^EOF", line):
                in_section = False
                finished = True
                modified_script.append(line)
            else:
                result = re.search(r"(\d{4})(\d{2})-(\d{4})(\d{2})\.nc", line)
                file_start = result.group(1)
                file_end = result.group(3)
                if int(file_end) >= int(start_year) and int(file_start) <= int(end_year):
                    modified_script.append(line)
        else:
            modified_script.append(line)
            if re.match(r"download_files=\"", line):
                in_section = True

    return "\n".join(modified_script)


def _download_result(result, download_path, logger: logging.Logger = LOGGER):
    max_retries = 3
    delay = 1
    for attempt in range(1, max_retries + 1):
        try:
            file_context = result.file_context()
            wget_script_content = file_context.get_download_script()
            subprocess.run(
                ["bash", "-c", wget_script_content, "download", "-s"], shell=False, cwd=download_path, check=False
            )
            break
        except Exception as e:  # pylint: disable=W0718
            logger.error(f"Attempt {attempt} failed: {e}")
            if attempt < max_retries:
                time.sleep(delay)
            else:
                raise e


def _download_process(temp_download_path, search_results, logger: logging.Logger = LOGGER):
    temp_download_path.mkdir(parents=True, exist_ok=True)
    for result in search_results:
        _download_result(result=result, download_path=temp_download_path, logger=logger)


def download_raw_input_variable(project, institution_id, search_results, variable, base_path: str | Path = RAW_DATA):
    if isinstance(base_path, str):
        base_path = Path(base_path)
    temp_download_path = base_path / f"{project}/raw_input_vars/{institution_id}/{variable}"
    _download_process(temp_download_path, search_results)


def download_model_variable(project, model_id, search_results, variable, base_path: str | Path = RAW_DATA):
    if isinstance(base_path, str):
        base_path = Path(base_path)
    temp_download_path = base_path / f"{project}/{model_id}/{variable}"
    _download_process(temp_download_path, search_results)


def download_metadata_variable(project, institution_id, search_results, variable, base_path: str | Path = RAW_DATA):
    if isinstance(base_path, str):
        base_path = Path(base_path)
    temp_download_path = base_path / f"{project}/meta_vars/{institution_id}/{variable}"
    _download_process(temp_download_path, search_results)


def get_grid_label(session, default_grid_label, logger=LOGGER):
    """
    Get grid label from the search session.

    Args:
        session: climateset.download.client.SearchSession
        default_grid_label: Default grid label to use if available
        logger: Logger instance

    Returns:
        str: Selected grid label
    """
    grid_label = ""
    grid_label_list = session.get_available_facets("grid_label")

    if grid_label_list:
        logger.info(f"Available grid labels : {grid_label_list}")
    else:
        logger.warning("No grid labels found")
        return grid_label

    if default_grid_label and default_grid_label in grid_label_list:
        logger.info(f"Choosing grid : {default_grid_label}")
        grid_label = default_grid_label
    else:
        logger.warning("Default grid label not available.")
        grid_label = grid_label_list[0]
        logger.info(f"Choosing grid {grid_label} instead.")
    return grid_label


def get_upload_version(session, preferred_version, logger=LOGGER):
    """
    Get upload version from the search session.

    Args:
        session: climateset.download.client.SearchSession
        preferred_version: Preferred version ('latest' or specific)
        logger: Logger instance

    Returns:
        str: Selected version
    """
    version = ""
    versions = session.get_available_facets("version")

    if not versions:
        logger.warning("No versions are available. Skipping.")
        return version

    logger.info(f"Available versions : {versions}")
    if preferred_version:
        if preferred_version == "latest":
            version = versions[0]
            logger.info(f"Choosing latest version: {version}")
        else:
            try:
                if preferred_version in versions:
                    version = preferred_version
                else:
                    raise KeyError
            except KeyError:
                logger.warning(f"Preferred version {preferred_version} does not exist.")
                version = versions[0]
                logger.info(f"Resuming with latest {version}:")
    return version


def get_frequency(session, default_frequency, logger=LOGGER):
    """
    Get frequency from the search session.

    Args:
        session: climateset.download.client.SearchSession
        default_frequency: Default frequency to use
        logger: Logger instance

    Returns:
        str: Selected frequency
    """
    frequency = ""
    frequency_list = session.get_available_facets("frequency")

    if frequency_list:
        logger.info(f"Available frequencies : {frequency_list}")
    else:
        logger.warning("No frequencies are available. Skipping")
        return frequency

    if default_frequency and default_frequency in frequency_list:
        frequency = default_frequency
        logger.info(f"Choosing default frequency : {frequency}")
    else:
        frequency = frequency_list[0]
        logger.info(f"Default frequency not available, choosing first available one instead: {frequency}")
    return frequency


def handle_yaml_config_path(config_file_name, config_path):
    if isinstance(config_path, str):
        config_path = Path(config_path)
    if not config_file_name.endswith(".yaml"):
        config_file_name = f"{config_file_name}.yaml"
    config_full_path = config_path / config_file_name
    return config_full_path


def match_key_in_list(input_key: str, key_list: list[str]) -> str | None:
    for key in key_list:
        if input_key.lower() == key.lower():
            return key
        if input_key.upper() == key.upper():
            return key
    return None


def search_and_download_esgf_raw_single_var(
    variable: str,
    institution_id: str,
    project: str,
    default_grid_label: str,
    default_frequency: str,
    preferred_version: str,
    data_dir,
    logger=LOGGER,
):
    """
    Search and download raw input variables using SearchClient.

    Follows iterative constraint strategy.
    """
    # Use distrib=False to match original behavior and avoid potential distributed search issues/warnings
    with SearchClient(distrib=False) as client:
        try:
            session = client.new_session()

            initial_constraints = Input4MIPsConstraints(
                project=project, institution_id=institution_id, variable=variable
            )
            session.constrain(initial_constraints)

            # 1. Grid Label
            grid_label = get_grid_label(session, default_grid_label, logger)
            if grid_label:
                session.constrain(BaseSearchConstraints(grid_label=grid_label))

            # 2. Nominal Resolution
            nominal_resolution = get_nominal_resolution(session, logger)
            if nominal_resolution:
                session.constrain(BaseSearchConstraints(nominal_resolution=nominal_resolution))

            # 3. Frequency
            frequency = get_frequency(session, default_frequency, logger)
            if frequency:
                session.constrain(BaseSearchConstraints(frequency=frequency))

            # 4. Target MIP (Iterative)
            mips_targets = session.get_available_facets("target_mip")
            logger.info(f"Available target mips: {mips_targets}")

            results_list = []

            base_constraints = [
                initial_constraints,
                BaseSearchConstraints(grid_label=grid_label) if grid_label else BaseSearchConstraints(),
                (
                    BaseSearchConstraints(nominal_resolution=nominal_resolution)
                    if nominal_resolution
                    else BaseSearchConstraints()
                ),
                BaseSearchConstraints(frequency=frequency) if frequency else BaseSearchConstraints(),
            ]

            for target in mips_targets:
                logger.debug(f"Checking target mip: {target}")
                sub_session = client.new_session()
                # Replay base constraints
                for c in base_constraints:
                    sub_session.constrain(c)

                sub_session.constrain(Input4MIPsConstraints(target_mip=target))

                version = get_upload_version(sub_session, preferred_version, logger)
                if version:
                    sub_session.constrain(BaseSearchConstraints(version=version))

                results = sub_session.search()
                logger.info(f"Result len for target {target}: {len(results)}")
                if results:
                    results_list.append(results)

            if results_list:
                for r in results_list:
                    download_raw_input_variable(
                        project=project,
                        institution_id=institution_id,
                        search_results=r,
                        variable=variable,
                        base_path=data_dir,
                    )
                return results_list
            logger.error("Could not find anything for configured nodes")

        except Exception as e:
            logger.error(f"Error during search/download: {e}")
            raise e

    raise RuntimeError("Could not find anything for all urls")


def search_and_download_esgf_biomass_single_var(
    variable: str,
    variable_id: str,
    institution_id: str,
    project: str,
    default_grid_label: str,
    default_frequency: str,
    preferred_version: str,
    base_path: Path,
    logger=LOGGER,
):
    with SearchClient(distrib=False) as client:
        try:
            session = client.new_session()

            initial_constraints = Input4MIPsConstraints(
                project=project, institution_id=institution_id, variable=variable, variable_id=variable_id
            )
            session.constrain(initial_constraints)

            # 1. Grid Label
            grid_label = get_grid_label(session, default_grid_label, logger)
            if grid_label:
                session.constrain(BaseSearchConstraints(grid_label=grid_label))

            # 2. Frequency
            frequency = get_frequency(session, default_frequency, logger)
            if frequency:
                session.constrain(BaseSearchConstraints(frequency=frequency))

            # 3. Version
            version = get_upload_version(session, preferred_version, logger)
            if version:
                session.constrain(BaseSearchConstraints(version=version))

            results = session.search()
            logger.info(f"Result len  {len(results)}")

            if results:
                try:
                    logger.info(results[0].file_context())
                except Exception:  # pylint: disable=broad-exception-caught
                    pass

                download_metadata_variable(
                    project=project,
                    institution_id=institution_id,
                    search_results=results,
                    variable=variable,
                    base_path=base_path,
                )
                return results
            logger.error("Could not find anything for configured nodes")

        except Exception as e:  # pylint: disable=broad-exception-caught
            logger.error(f"Error: {e}")

    raise RuntimeError("Could not find anything for all urls")


def _get_variants_and_filter(
    session: SearchSession, max_ensemble_members: int, ensemble_members: list[str], logger: logging.Logger
) -> list[str]:
    """Helper to retrieve and filter variant labels."""
    variants = session.get_available_facets("variant_label")

    if len(variants) < 1:
        # Note: Previous code raised ValueError here but logging info first
        return []

    logger.info(f"Available variants : {variants}\n")
    logger.info(f"Length : {len(variants)}")

    if not ensemble_members:
        if max_ensemble_members > len(variants):
            logger.info("Less ensemble members available than maximum number desired. Including all variants.")
            return variants
        logger.info(
            f"{len(variants)} ensemble members available, desired (max {max_ensemble_members}). "
            f"Choosing only the first {max_ensemble_members}."
        )
        return variants[:max_ensemble_members]

    logger.info(f"Desired list of ensemble members given: {ensemble_members}")
    ensemble_member_final_list = list(set(variants) & set(ensemble_members))
    return ensemble_member_final_list


def search_and_download_esgf_model_single_var(
    model: str,
    variable: str,
    experiment: str,
    project: str,
    default_grid_label: str,
    default_frequency: str,
    preferred_version: str,
    max_ensemble_members: int,
    ensemble_members: list[str],
    base_path: Path,
    logger=LOGGER,
):
    logger.info("Using download_from_model_single_var() function")

    with SearchClient() as client:
        try:
            session = client.new_session()

            cmip_constraints = CMIP6Constraints(
                project=project, experiment_id=experiment, source_id=model, variable=variable
            )
            session.constrain(cmip_constraints)

            # 1. Frequency
            frequency = get_frequency(session, default_frequency, logger)
            if frequency:
                session.constrain(BaseSearchConstraints(frequency=frequency))

            # 2. Grid Label
            grid_label = get_grid_label(session, default_grid_label, logger)
            if grid_label:
                session.constrain(BaseSearchConstraints(grid_label=grid_label))

            # 3. Variants (Ensemble Members)
            ensemble_member_final_list = _get_variants_and_filter(
                session, max_ensemble_members, ensemble_members, logger
            )

            if not ensemble_member_final_list:
                logger.info(
                    "No items were found for this request. Please check on the esgf server if the combination of your "
                    "model/scenarios/variables exists."
                )
                if not ensemble_members and len(session.get_available_facets("variant_label")) < 1:
                    # Replicate original ValueError text
                    raise ValueError(
                        f"Downloader did not find any items on esgf for your request with: Project {project}, "
                        f"Experiment {experiment}, Model {model}, Variable {variable}."
                    )
                logger.info("WARNING: no overlap between available and desired ensemble members!")
                logger.info("Skipping.")
                return None

            results_list = []

            base_constraints_list = [
                cmip_constraints,
                BaseSearchConstraints(frequency=frequency) if frequency else BaseSearchConstraints(),
                BaseSearchConstraints(grid_label=grid_label) if grid_label else BaseSearchConstraints(),
            ]

            for ensemble_member in ensemble_member_final_list:
                logger.info(f"Ensembles member: {ensemble_member}")

                sub_session = client.new_session()
                for c in base_constraints_list:
                    sub_session.constrain(c)

                sub_session.constrain(CMIP6Constraints(variant_label=ensemble_member))

                version = get_upload_version(sub_session, preferred_version, logger)
                if version:
                    sub_session.constrain(BaseSearchConstraints(version=version))

                results = sub_session.search()
                if results:
                    results_list.append(results)
                logger.info(f"Result len {len(results)}")

            logger.info(results_list)
            if results_list:
                for results in results_list:
                    download_model_variable(
                        project=project,
                        model_id=model,
                        search_results=results,
                        variable=variable,
                        base_path=base_path,
                    )
                return results_list
            logger.error("Could not find anything for configured nodes")

        except ValueError:
            raise
        except Exception as e:  # pylint: disable=broad-exception-caught
            logger.error(f"Error: {e}")

    raise RuntimeError("Could not find anything for all urls")


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
