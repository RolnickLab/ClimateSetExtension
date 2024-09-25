import logging
import re

import xarray as xr

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
