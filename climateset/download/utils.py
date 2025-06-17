import logging
import re
import subprocess
import time
from pathlib import Path

import xarray as xr
from pyesgf.search import SearchConnection
from pyesgf.search.context import DatasetSearchContext

from climateset import RAW_DATA
from climateset.download.constants import NODE_LINK_URLS
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


def get_nominal_resolution(context, logger: logging.Logger = LOGGER):
    """

    Args:
        context:
        logger:

    Returns:

    """
    nominal_resolution = ""
    nominal_resolution_list = []
    if "nominal_resolution" in context.facet_counts:
        nominal_resolution_list = list(context.facet_counts["nominal_resolution"])
        logger.info(f"Available nominal resolution : {nominal_resolution_list}")
    if not nominal_resolution_list:
        logger.warning("No nominal resolution")
        return nominal_resolution
    if len(nominal_resolution_list) > 1:
        logger.warning("Multiple nominal resolutions exist, will try to get smallest resolution.")
    nominal_resolution = nominal_resolution_list[0]
    logger.info(f"Choosing nominal resolution : {nominal_resolution}")
    return nominal_resolution


def infer_nominal_resolution(ds: xr.Dataset, nominal_resolution: str, logger: logging.Logger = LOGGER) -> str:
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


def get_grid_label(context, default_grid_label, logger=LOGGER):
    grid_label = ""
    grid_label_list = []
    if "grid_label" in context.facet_counts:
        grid_label_list = list(context.facet_counts["grid_label"])
        logger.info(f"Available grid labels : {grid_label_list}")
    if not grid_label_list:
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


def get_upload_version(context, preferred_version, logger=LOGGER):
    version = ""
    versions = []
    if "version" in context.facet_counts:
        versions = list(context.facet_counts["version"])
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
                version = versions[preferred_version]
            except KeyError:
                logger.warning(f"Preferred version {preferred_version} does not exist.")
                version = versions[0]
                logger.info(f"Resuming with latest {version}:")
    return version


def get_frequency(context, default_frequency, logger=LOGGER):
    frequency = ""
    frequency_list = []
    if "frequency" in context.facet_counts:
        frequency_list = list(context.facet_counts["frequency"])
        logger.info(f"Available frequencies : {frequency_list}")
    if not frequency_list:
        logger.warning("No frequencies are available. Skipping")
        return frequency
    if default_frequency and default_frequency in frequency_list:
        frequency = default_frequency
        logger.info(f"Choosing default frequency : {frequency}")
    else:
        frequency = frequency_list[0]
        logger.info(f"Default frequency not available, choosing first available one instead: {frequency}")
    return frequency


def handle_base_search_constraints(ctx, default_frequency, default_grid_label):
    grid_label = get_grid_label(context=ctx, default_grid_label=default_grid_label)
    if grid_label:
        ctx = ctx.constrain(grid_label=grid_label)
    nominal_resolution = get_nominal_resolution(context=ctx)
    if nominal_resolution:
        ctx = ctx.constrain(nominal_resolution=nominal_resolution)
    frequency = get_frequency(context=ctx, default_frequency=default_frequency)
    if frequency:
        ctx = ctx.constrain(frequency=frequency)
    return ctx


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


def get_base_search_context(
    url: str = None,
    facets: str = None,
    variable: str = None,
    variable_id: str = None,
    institution_id: str = None,
    project: str = None,
    experiment_id: str = None,
    source_id: str = None,
    default_grid_label: str = None,
    default_frequency: str = None,
) -> DatasetSearchContext:
    conn = SearchConnection(url=url, distrib=False)
    ctx = conn.new_context(
        project=project,
        variable=variable,
        variable_id=variable_id,
        institution_id=institution_id,
        experiment_id=experiment_id,
        source_id=source_id,
        facets=facets,
    )
    ctx = handle_base_search_constraints(ctx, default_frequency, default_grid_label)
    return ctx


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
    facets = "project,frequency,variable,nominal_resolution,version,target_mip,grid_label"
    for url in NODE_LINK_URLS:
        results_list = []
        try:

            ctx = get_base_search_context(
                url=url,
                project=project,
                institution_id=institution_id,
                variable=variable,
                facets=facets,
                default_grid_label=default_grid_label,
                default_frequency=default_frequency,
            )

            mips_targets = list(ctx.facet_counts["target_mip"])
            logger.info(f"Available target mips: {mips_targets}")

            for target in mips_targets:
                ctx_target = ctx.constrain(target_mip=target)
                version = get_upload_version(context=ctx_target, preferred_version=preferred_version)
                if version:
                    ctx_target = ctx_target.constrain(version=version)

                results = ctx_target.search()
                logger.info(f"Result len  {len(results)}")
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
            logger.error(f"Could not find anything for {url}")
        except Exception as e:
            logger.error(f"Error: {e}")

    raise Exception(f"Could not find anything for all urls: {NODE_LINK_URLS}")


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
    facets = "nominal_resolution,version"
    for url in NODE_LINK_URLS:
        try:
            ctx = get_base_search_context(
                url=url,
                facets=facets,
                variable=variable,
                variable_id=variable_id,
                institution_id=institution_id,
                project=project,
                default_grid_label=default_grid_label,
                default_frequency=default_frequency,
            )

            version = get_upload_version(context=ctx, preferred_version=preferred_version)
            if version:
                ctx = ctx.constrain(version=version)

            results = ctx.search()
            logger.info(f"Result len  {len(results)}")

            result_list = [r.file_context().search() for r in results]
            logger.info(f"List of results :\n{result_list}")
            if results:
                logger.info(results[0].file_context())
                download_metadata_variable(
                    project=project,
                    institution_id=institution_id,
                    search_results=results,
                    variable=variable,
                    base_path=base_path,
                )
                return results
            logger.error(f"Could not find anything for {url}")
        except Exception as e:
            logger.error(f"Error: {e}")

    raise Exception(f"Could not find anything for all urls: {NODE_LINK_URLS}")


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
    facets = (
        "project,experiment_id,source_id,variable,frequency,variant_label,variable, nominal_resolution, "
        "version, grid_label, experiment_id"
    )

    for url in NODE_LINK_URLS:
        results_list = []
        try:
            logger.info("Using download_from_model_single_var() function")

            ctx = get_base_search_context(
                url=url,
                facets=facets,
                variable=variable,
                experiment_id=experiment,
                source_id=model,
                default_frequency=default_frequency,
                default_grid_label=default_grid_label,
            )

            logger.info(ctx)

            variants = list(ctx.facet_counts["variant_label"])

            if len(variants) < 1:
                logger.info(
                    "No items were found for this request. Please check on the esgf server if the combination of your "
                    "model/scenarios/variables exists."
                )
                raise ValueError(
                    f"Downloader did not find any items on esgf for your request with: Project {project}, "
                    f"Experiment {experiment}, Model {model}, Variable {variable}."
                )

            logger.info(f"Available variants : {variants}\n")
            logger.info(f"Length : {len(variants)}")

            # TODO refactor logic of if/else
            if not ensemble_members:
                if max_ensemble_members > len(variants):
                    logger.info("Less ensemble members available than maximum number desired. Including all variants.")
                    ensemble_member_final_list = variants
                else:
                    logger.info(
                        f"{len(variants)} ensemble members available than desired (max {max_ensemble_members}. "
                        f"Choosing only the first {max_ensemble_members}.)."
                    )
                    ensemble_member_final_list = variants[:max_ensemble_members]
            else:
                logger.info(f"Desired list of ensemble members given: {ensemble_members}")
                ensemble_member_final_list = list(set(variants) & set(ensemble_members))
                if len(ensemble_member_final_list) == 0:
                    logger.info("WARNING: no overlap between available and desired ensemble members!")
                    logger.info("Skipping.")
                    return None

            for ensemble_member in ensemble_member_final_list:
                logger.info(f"Ensembles member: {ensemble_member}")
                ctx_ensemble = ctx.constrain(variant_label=ensemble_member)
                logger.info(ctx_ensemble)

                version = get_upload_version(context=ctx, preferred_version=preferred_version)
                if version:
                    ctx_ensemble = ctx_ensemble.constrain(version=version)

                results = ctx_ensemble.search()
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
            logger.error(f"Could not find anything for {url}")
        except Exception as e:
            logger.error(f"Error: {e}")

    raise Exception(f"Could not find anything for all urls: {NODE_LINK_URLS}")
