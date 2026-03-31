class OmipConstants:
    """
    Attributes:
        NODE_LINK (str): Where the data can be accessed
        MODEL_SOURCES (List<str>): Identifiers for supported climate models
        VAR_SOURCE_LOOKUP (Dict<str, List<str>>): model and raw variables
        SUPPORTED_EXPERIMENTS (list<str>): experiments of climate models (runs) that are supported
    """

    NODE_LINK = "http://esgf-data2.llnl.gov"

    MODEL_SOURCES = [
        "NorESM2-LM",
    ]

    VAR_SOURCE_LOOKUP = [
        "omldamax",
    ]

    SUPPORTED_EXPERIMENTS = [
        "omip1",
    ]
