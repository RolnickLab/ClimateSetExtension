import logging
import sys


def create_logger(logger_name: str) -> logging.Logger:
    """Creates a logger object using input name parameter that outputs to stdout.

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
    return logger


LOGGER = create_logger(__name__)


def get_keys_from_value(d, val, logger=LOGGER):
    keys = [k for k, v in d.items() if val in v]
    if keys:
        return keys[0]
    logger.warning(f"WARNING: source not found vor var {val}")
    return None


def get_mip(experiment: str):
    """
    Return name of MIP group given the specific experiment name.
    """
    if experiment == "ssp245-covid":
        return "DAMIP"
    if experiment == "ssp370-lowNTCF":
        return "AerChemMIP"
    if experiment.startswith("ssp"):
        return "ScenarioMIP"
    if experiment.startswith("hist-"):
        return "DAMIP"
    return "CMIP"
