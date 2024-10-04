import logging
import pathlib
import sys
from typing import Union

import yaml

from climateset import CONFIGS


def create_logger(logger_name: str) -> logging.Logger:
    """
    Creates a logger object using input name parameter that outputs to stdout.

    Args:
        logger_name (str) :Name of logger

    Returns:
        logging.Logger:
        Created logger object
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        fmt="[%(asctime)s] %(levelname)-10.10s [%(threadName)s][%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


LOGGER = create_logger(__name__)


def get_keys_from_value(d, val, logger=LOGGER):
    keys = [k for k, v in d.items() if val in v]
    if keys:
        return keys[0]
    logger.warning(f"WARNING: source not found vor var {val}")
    return None


def get_mip(experiment: str):
    """Return name of MIP group given the specific experiment name."""
    if experiment == "ssp245-covid":
        return "DAMIP"
    if experiment == "ssp370-lowNTCF":
        return "AerChemMIP"
    if experiment.startswith("ssp"):
        return "ScenarioMIP"
    if experiment.startswith("hist-"):
        return "DAMIP"
    return "CMIP"


def get_yaml_config(yaml_config_file: Union[str, pathlib.Path], logger: logging.Logger = LOGGER) -> dict:
    """
    Reads a YAML configuration file and returns its contents as a dictionary.

    This function searches for the specified YAML file in the `config/`
    directory. If the file is found, its contents are parsed and returned as a
    dictionary.

    Args:
        yaml_config_file: The path to the YAML config file. If the file is
            located in the `config/` directory, you can provide the file's
            name without the extension.
        logger: The logger instance to handle messaging. Defaults to the
            global `LOGGER`.

    Returns:
        A dictionary containing the parsed YAML configuration values.

    Raises:
        FileNotFoundError: If the specified YAML file is not found.
        YAMLError: If there's an error parsing the YAML file.

    Examples:
        # For a file named 'app_config.yml' in the 'config/' folder:
        params = get_yaml_config('app_config')
    """
    if isinstance(yaml_config_file, str):
        yaml_config_file = pathlib.Path(yaml_config_file)
    potential_paths = [
        pathlib.Path(yaml_config_file),
        CONFIGS / yaml_config_file,
        CONFIGS / f"{yaml_config_file}.yaml",
        CONFIGS / f"{yaml_config_file}.yml",
    ]

    config_filepath = None
    for path in potential_paths:
        if path.exists():
            config_filepath = path
            logger.info(f"Yaml config file [{str(path)}] found.")
            break

    params = {}
    if not config_filepath:
        logger.error(f"Yaml config file [{yaml_config_file}] was not found.")
        return params

    try:
        with config_filepath.open("r", encoding="UTF-8") as file:
            logger.info(f"Loading YAML config file [{config_filepath}].")
            return yaml.safe_load(file)
    except yaml.YAMLError as e:
        logger.warning(f"Error loading YAML file [{config_filepath}]: {e}")
        return {}
