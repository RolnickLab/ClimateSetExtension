NODE_LINK = "http://esgf-data2.llnl.gov"

MODEL_SOURCES = [
    "HasGEM3-GC31-LL",
]

VAR_SOURCE_LOOKUP = {
    "model": [
        "areacella",
        "mrsofc",
    ],
    "raw": [
        "areacella",
        "mrsofc",
    ],
}

SUPPORTED_EXPERIMENTS = [
    "hist-lu",
    "hist-piAer",
    "hist-piVolc",
]

GRIDDING_HIERACHY = ["gn"]

RES_TO_CHUNKSIZE = {"year": 1, "mon": 12, "6hr": 1460, "3hr": 2920, "day": 364}
