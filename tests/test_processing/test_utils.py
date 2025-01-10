from climateset import TEST_DIR
from climateset.processing.raw.input4mips.processing import rename_biomassburning_file


def test_rename_biomassburning_files():
    test_file = (
        TEST_DIR / "resources/CO2-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-1_gn_185001-201512.nc"
    )
    result = rename_biomassburning_file(test_file)
    print(result)
    assert False
