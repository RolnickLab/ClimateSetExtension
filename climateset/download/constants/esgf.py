from .cmip6 import Cmip6Constants
from .cmip6plus import Cmip6plusConstants
from .input4mips import Input4mipsConstants

CMIP6 = "CMIP6"
CMIP6PLUS = "CMIP6Plus"
INPUT4MIPS = "input4MIPs"

ESGF_PROJECTS = frozenset([CMIP6, CMIP6PLUS, INPUT4MIPS])

# constant classes for esgf projects implemented here
# add your own esgf project for downloading to download/constants/ and add the constant class to the dict and lists here
ESGF_PROJECTS_CONSTANTS = {
    CMIP6: Cmip6Constants,
    CMIP6PLUS: Cmip6plusConstants,
    INPUT4MIPS: Input4mipsConstants,
}

# datasets that provide inputs to climate models
ESGF_RAW_INPUT_LIST = [INPUT4MIPS]

# datasets that provide outputs from climate models
ESGF_MODEL_OUTPUT_LIST = [CMIP6, CMIP6PLUS]
