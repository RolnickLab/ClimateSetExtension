import logging
import re
import subprocess

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
    nominal_resolutions = list(ctx.facet_counts["nominal_resolution"].keys())
    logger.info(f"Available nominal resolution : {nominal_resolutions}")
    # deal with multiple nominal resolutions, taking smallest one as default
    if len(nominal_resolutions) > 1:
        logger.info(
            "Multiple nominal resolutions exist, choosing smallest_nominal resolution (trying), " "please do a check up"
        )
    nominal_resolution = nominal_resolutions[0]
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


# TODO add retry + logger so failure can be tracked
def raw_download_process(model_or_institution_id, search_results, variable):
    temp_download_path = RAW_DATA / f"{model_or_institution_id}/{variable}"
    temp_download_path.mkdir(parents=True, exist_ok=True)
    for result in search_results:
        file_context = result.file_context()
        wget_script_content = file_context.get_download_script()

        # Optionally filter file list for download
        # if self.start_year is not None and self.end_year is not None:
        # wget_script_content = filter_download_script(wget_script_content, self.start_year, self.end_year)

        subprocess.run(["bash", "-c", wget_script_content, "download", "-s"], shell=False, cwd=temp_download_path)


def get_grid_label(ctx, default_grid_label, logger=LOGGER):
    grid_labels = list(ctx.facet_counts["grid_label"].keys())
    logger.info(f"Available grid labels : {grid_labels}")
    if default_grid_label is not None:
        if default_grid_label in grid_labels:
            logger.info(f"Choosing grid : {default_grid_label}")
            grid_label = default_grid_label
        else:
            logger.info("Default grid label not available.")
            grid_label = grid_labels[0]
            logger.info(f"Choosing grid {grid_label} instead.")
    else:
        grid_label = ""
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


def get_upload_version(ctx, preferred_version, logger=LOGGER):
    version = ""
    versions = list(ctx.facet_counts["version"].keys())
    if not versions:
        logger.info("No versions are available. Skipping.")
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
                logger.info(f"Preferred version {preferred_version} does not exist.")
                version = versions[0]
                logger.info(f"Resuming with latest {version}:")
    return version


def get_frequency(ctx, default_frequency, logger=LOGGER):
    if default_frequency is not None:
        # choose default frequency if wanted
        frequencies = list(ctx.facet_counts["frequency"].keys())
        logger.info(f"Available frequencies : {frequencies}")

        if default_frequency in frequencies:
            frequency = default_frequency
            logger.info(f"Choosing default frequency : {frequencies}")
        else:
            frequency = frequencies[0]
            logger.info(
                "Default frequency not available, choosing first available one instead: ",
                frequency,
            )
    else:
        frequency = ""
    return frequency
