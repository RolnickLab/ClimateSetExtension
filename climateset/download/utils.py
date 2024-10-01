import logging
import re
import subprocess
import time

import xarray as xr

from climateset import RAW_DATA
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
    # logger.info(f'years from {year_from} to {year_end}')

    if (target_mip == "ScenarioMIP") or (target_mip == "DAMIP"):
        # extract ssp experiment from file name
        experiment = f"ssp{filename.split('ssp')[-1][:3]}"
        if "covid" in filename:
            experiment = f"{experiment}_covid"
    elif target_mip == "CMIP":
        if int(year_end) > 2015:
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


def get_nominal_resolution(ctx, logger: logging.Logger = LOGGER):
    """

    Args:
        ctx:
        logger:

    Returns:

    """
    nominal_resolution = ""
    nominal_resolution_list = []
    if "nominal_resolution" in ctx.facet_counts.keys():
        nominal_resolution_list = list(ctx.facet_counts["nominal_resolution"].keys())
        logger.info(f"Available nominal resolution : {nominal_resolution_list}")
    if not nominal_resolution_list:
        logger.warning("No nominal resolution")
        return nominal_resolution
    # deal with multiple nominal resolutions, taking smallest one as default
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
    except Exception as error:
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
                continue
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
            subprocess.run(["bash", "-c", wget_script_content, "download", "-s"], shell=False, cwd=download_path)
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


def download_raw_input_variable(institution_id, search_results, variable):
    temp_download_path = RAW_DATA / f"raw_input_vars/{institution_id}/{variable}"
    _download_process(temp_download_path, search_results)


def download_model_variable(model_id, search_results, variable):
    temp_download_path = RAW_DATA / f"model_vars/{model_id}/{variable}"
    _download_process(temp_download_path, search_results)


def download_metadata_variable(institution_id, search_results, variable):
    temp_download_path = RAW_DATA / f"meta_vars/{institution_id}/{variable}"
    _download_process(temp_download_path, search_results)


def get_grid_label(ctx, default_grid_label, logger=LOGGER):
    grid_label = ""
    grid_label_list = []
    if "grid_label" in ctx.facet_counts.keys():
        grid_label_list = list(ctx.facet_counts["grid_label"].keys())
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


def get_max_ensemble_member_number(df_model_source, experiments, model, logger=LOGGER):
    if model is not None:
        if model not in df_model_source["source_id"].tolist():
            logger.info(f"Model {model} not supported.")
            raise AttributeError
        # extract member information
        model_id = df_model_source.index[df_model_source["source_id"] == model].values
        # get ensemble members per scenario
        max_ensemble_members_list = df_model_source["num_ensemble_members"][model_id].values.tolist()[0].split(" ")
        scenarios = df_model_source["scenarios"][model_id].values.tolist()[0].split(" ")
        # create lookup
        max_ensemble_members_lookup = {}
        for s, m in zip(scenarios, max_ensemble_members_list):
            max_ensemble_members_lookup[s] = int(m)
        max_possible_member_number = min(
            [max_ensemble_members_lookup[e] for e in experiments if e != "historical"]
        )  # TODO fix historical
    return max_possible_member_number


def get_upload_version(context, preferred_version, logger=LOGGER):
    version = ""
    versions = []
    if "version" in context.facet_counts.keys():
        versions = list(context.facet_counts["version"].keys())
    if not versions:
        logger.warning("No versions are available. Skipping.")
        return version
    logger.info(f"Available versions : {versions}")
    if preferred_version:
        # deal with different versions
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


def get_frequency(ctx, default_frequency, logger=LOGGER):
    frequency = ""
    frequency_list = []
    if "frequency" in ctx.facet_counts.keys():
        frequency_list = list(ctx.facet_counts["frequency"].keys())
        logger.warning(f"Available frequencies : {frequency_list}")
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
