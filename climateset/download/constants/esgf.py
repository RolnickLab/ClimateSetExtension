from .cmip6 import Cmip6Constants
from .cmip6plus import Cmip6plusConstants
from .input4mips import Input4mipsConstants

# constant classes for esgf projects implemented here
# add your own esgf project for downloading to download/constants/ and add the constant class to the dict and lists here
ESGF_PROJECTS = {
    "CMIP6": Cmip6Constants,
    "CMIP6Plus": Cmip6plusConstants,
    "input4MIPs": Input4mipsConstants,
}

# datasets that provide inputs to climate models
ESGF_RAW_INPUT_LIST = ["input4MIPs"]

# datasets that provide outputs from climate models
ESGF_MODEL_OUTPUT_LIST = ["CMIP6", "CMIP6Plus"]
