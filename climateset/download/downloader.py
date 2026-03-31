import logging
import pathlib

from climateset.download.cmip6_downloader import cmip6_download_from_config
from climateset.download.constants.esgf import CMIP6, INPUT4MIPS
from climateset.download.downloader_config import AVAILABLE_CONFIGS
from climateset.download.input4mips_downloader import input4mips_download_from_config
from climateset.download.utils import match_key_in_list
from climateset.utils import create_logger, get_yaml_config

LOGGER = create_logger(__name__)


def download_from_config_file(config_file: str | pathlib.Path, logger: logging.Logger = LOGGER):
    """
    This function downloads variables automatically from input config file
    Args:
        config_file: Path to a configuration yaml file
        logger: Logging instance
    """
    if isinstance(config_file, str):
        config_file = pathlib.Path(config_file)
    config_dict = get_yaml_config(config_file)

    downloader_factory = {
        INPUT4MIPS: input4mips_download_from_config,
        CMIP6: cmip6_download_from_config,
    }

    for config_key in config_dict:
        verified_key = match_key_in_list(input_key=config_key, key_list=AVAILABLE_CONFIGS)
        if verified_key:
            downloader_factory[verified_key](config=config_file)
        else:
            logger.error(
                f"Input project [{config_key}] from [{config_file}] was not found in available projects. "
                "Removing it from download list"
            )
