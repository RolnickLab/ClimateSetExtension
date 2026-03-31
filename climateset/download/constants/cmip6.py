# TODO remove raw variables from here
# pylint: disable=C0103
from dataclasses import dataclass
from typing import Final

from climateset.utils import get_yaml_config


@dataclass(frozen=True)
class Cmip6Constants:
    """
    Dataclass to represent CMIP6 constants that are used by the download module.

    Attributes:
        NODE_LINK : Where the data can be accessed
        MODEL_SOURCES : Identifiers for supported climate models
        VAR_SOURCE_LOOKUP : model and raw variables
        SUPPORTED_EXPERIMENTS : experiments of climate models (runs) that are supported
    """

    NODE_LINK: Final[str]
    MODEL_SOURCES: Final[tuple[str, ...]]
    VAR_SOURCE_LOOKUP: Final[tuple[str, ...]]
    SUPPORTED_EXPERIMENTS: Final[tuple[str, ...]]


_data = get_yaml_config("downloader/constants/cmip6.yaml")

CMIP6_CONSTANTS = Cmip6Constants(
    NODE_LINK=_data["node_link"],
    MODEL_SOURCES=tuple(_data["model_sources"]),
    SUPPORTED_EXPERIMENTS=tuple(_data["supported_experiments"]),
    VAR_SOURCE_LOOKUP=tuple(_data["var_source_lookup"]),
)
