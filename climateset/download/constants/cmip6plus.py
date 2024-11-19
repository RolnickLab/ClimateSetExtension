# TODO remove raw variables from here
class Cmip6plusConstants:
    """
    Attributes:
        NODE_LINK (str): Where the data can be accessed
        MODEL_SOURCES (List<str>): Identifiers for supported climate models
        VAR_SOURCE_LOOKUP (Dict<str, List<str>>): model and raw variables
        SUPPORTED_EXPERIMENTS (list<str>): experiments of climate models (runs) that are supported
    """

    NODE_LINK = "http://esgf-data2.llnl.gov"

    MODEL_SOURCES = [
        "HasGEM3-GC31-LL",
    ]

    VAR_SOURCE_LOOKUP = [
        "areacella",
        "mrsofc",
    ]

    SUPPORTED_EXPERIMENTS = [
        "hist-lu",
        "hist-piAer",
        "hist-piVolc",
    ]
