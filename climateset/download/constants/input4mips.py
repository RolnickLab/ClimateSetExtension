# TODO add VAR_SOURCE_LOOKUP with raw variables
# TODO add supported experiments
# TODO do we really need emission endings, meta_endings_prc, meta_endings_shar?? how is this used so far?
from dataclasses import dataclass
from typing import Final

from climateset.utils import get_yaml_config


@dataclass(frozen=True)
class Input4MIPSConstants:
    NODE_LINK: Final[str]
    EMISSIONS_ENDINGS: Final[tuple[str, ...]]
    META_ENDINGS_PRC: Final[tuple[str, ...]]
    META_ENDINGS_SHAR: Final[tuple[str, ...]]
    MIP_ERA: Final[str]
    TARGET_MIP: Final[str]
    SUPPORTED_EXPERIMENTS: Final[tuple[str, ...]]
    VAR_SOURCE_LOOKUP: Final[tuple[str, ...]]


_data = get_yaml_config("downloader/constants/imput4MIPs.yaml")

INPUT4MIPS_CONSTANTS = Input4MIPSConstants(
    NODE_LINK=_data["node_link"],
    EMISSIONS_ENDINGS=tuple(_data["emissions_endings"]),
    META_ENDINGS_PRC=tuple(_data["meta_endings_prc"]),
    META_ENDINGS_SHAR=tuple(_data["meta_endings_shar"]),
    MIP_ERA=_data["mip_era"],
    TARGET_MIP=_data["target_mip"],
    SUPPORTED_EXPERIMENTS=tuple(_data["supported_experiments"]),
    VAR_SOURCE_LOOKUP=tuple(_data["var_source_lookup"]),
)
