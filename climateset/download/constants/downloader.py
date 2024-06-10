### DOWNLOADER PARAMS ##########################################################

# list of years that are considered for the data
YEARS = [0]

CO2 = ["CO2", "CO2_em_anthro", "CO2_em_openburning", "CO2_em_AIR_anthro"]
BC = ["BC", "BC_em_anthro", "BC_em_openburning", "BC_em_AIR_anthro"]
CH4 = ["CH4", "CH4_em_anthro", "CH4_em_openburning", "CH4_em_AIR_anthro"]
SO2 = ["SO2", "SO2_em_anthro", "SO2_em_openburning", "SO2_em_AIR_anthro"]

IN_VARS = CO2 + BC + CH4 + SO2
OUT_VARS = ["pr", "tas"]

VARS = IN_VARS + OUT_VARS

# scenarios
SCENARIOS = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
ADDITIONAL_SCENARIOS = ["hist-aer", "hist-GHG", "piControl", "ssp370-lowNTCF"]

# model
MODELS = ["nan"]

# number of esemble members to be considered
NUM_ENSEMBLE = 1

# which type of grid
GRID = "grid"
