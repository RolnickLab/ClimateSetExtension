import pathlib
from typing import Union

from climateset.download.cmip6_downloader import CMIP6Downloader
from climateset.download.constants.esgf import CMIP6, INPUT4MIPS
from climateset.download.downloader_config import (
    AVAILABLE_CONFIGS,
    create_cmip6_downloader_config_from_file,
    create_input4mips_downloader_config_from_file,
)
from climateset.download.input4mips_downloader import Input4MipsDownloader
from climateset.download.utils import match_key_in_list
from climateset.utils import create_logger, get_yaml_config

LOGGER = create_logger(__name__)


def download_from_config_file(config_file: Union[str, pathlib.Path]):
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
        INPUT4MIPS: {"configs": create_input4mips_downloader_config_from_file, "downloader": Input4MipsDownloader},
        CMIP6: {"configs": create_cmip6_downloader_config_from_file, "downloader": CMIP6Downloader},
    }

    verified_config_keys = []
    for config_key in config_dict:
        verified_key = match_key_in_list(input_key=config_key, key_list=AVAILABLE_CONFIGS)
        if verified_key:
            verified_config_keys.append(verified_key)

    for config_key in verified_config_keys:
        configs = downloader_factory[config_key]["configs"](config_file=config_file)
        downloader = downloader_factory[config_key]["downloader"](config=configs)
        downloader.download()
