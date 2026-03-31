import shutil
from unittest.mock import ANY, call, patch

import pytest

from climateset import TEST_DIR
from climateset.download.cmip6_downloader import CMIP6Downloader
from climateset.download.constants.esgf import CMIP6, INPUT4MIPS
from climateset.download.downloader import download_from_config_file
from climateset.download.downloader_config import (
    create_cmip6_downloader_config_from_file,
    create_input4mips_downloader_config_from_file,
)
from climateset.download.input4mips_downloader import Input4MipsDownloader

MINIMAL_DATASET_CONFIG_PATH = TEST_DIR / "resources/test_minimal_dataset.yaml"
TEST_TMP_DIR = TEST_DIR / "resources/.tmp"
MAX_ENSEMBLE_MEMBERS = 10

DOWNLOAD_RAW_INPUT_SINGLE_VAR = (
    "climateset.download.input4mips_downloader.Input4MipsDownloader.download_raw_input_single_var"
)
DOWNLOAD_META_HISTORIC_SINGLE_VAR = (
    "climateset.download.input4mips_downloader.Input4MipsDownloader.download_meta_historic_biomassburning_single_var"
)
DOWNLOAD_MODEL_SINGLE_VAR = "climateset.download.cmip6_downloader.CMIP6Downloader.download_from_model_single_var"
SUBPROCESS_RUN = "subprocess.run"

DOWNLOAD_RAW_INPUT_SINGLE_VAR_V2 = (
    "climateset.download.input4mips_downloader.Input4MipsDownloaderV2.download_raw_input_single_var"
)
DOWNLOAD_META_HISTORIC_SINGLE_VAR_V2 = (
    "climateset.download.input4mips_downloader.Input4MipsDownloaderV2.download_meta_historic_biomassburning_single_var"
)
DOWNLOAD_MODEL_SINGLE_VAR_V2 = "climateset.download.cmip6_downloader.CMIP6DownloaderV2.download_from_model_single_var"

EXPECTED_MINIMAL_RAW_INPUT_CALLS = [
    call(variable="CO2_em_anthro", institution_id="PNNL-JGCRI"),
    call(variable="CO2_em_AIR_anthro", institution_id="PNNL-JGCRI"),
    call(variable="CH4_em_openburning", institution_id="IAMC"),
    call(variable="CH4_em_anthro", institution_id="PNNL-JGCRI"),
    call(variable="CH4_em_AIR_anthro", institution_id="PNNL-JGCRI"),
    call(variable="CO2", institution_id="VUA"),
    call(variable="CH4", institution_id="VUA"),
    call(variable="CH4_openburning_share", institution_id="IAMC"),
]
RAW_INPUT_NUM_OF_CALLS = 8

EXPECTED_MINIMAL_META_HISTORIC_CALLS = [
    call(variable="CH4_percentage_AGRI", institution_id="VUA"),
    call(variable="CH4_percentage_BORF", institution_id="VUA"),
    call(variable="CH4_percentage_DEFO", institution_id="VUA"),
    call(variable="CH4_percentage_PEAT", institution_id="VUA"),
    call(variable="CH4_percentage_SAVA", institution_id="VUA"),
    call(variable="CH4_percentage_TEMF", institution_id="VUA"),
]
META_HISTORIC_NUM_OF_CALLS = 6

EXPECTED_MINIMAL_MODEL_CALLS = [
    call(model="NorESM2-LM", project="CMIP6", variable="tas", experiment="ssp126"),
]
MODEL_SINGLE_NUM_OF_CALLS = 1

EXPECTED_MINIMAL_RAW_INPUT_CALLS_V2 = [
    call(esg=ANY, variable="CO2_em_anthro", institution_id="PNNL-JGCRI"),
    call(esg=ANY, variable="CO2_em_AIR_anthro", institution_id="PNNL-JGCRI"),
    call(esg=ANY, variable="CH4_em_openburning", institution_id="IAMC"),
    call(esg=ANY, variable="CH4_em_anthro", institution_id="PNNL-JGCRI"),
    call(esg=ANY, variable="CH4_em_AIR_anthro", institution_id="PNNL-JGCRI"),
    call(esg=ANY, variable="CO2", institution_id="VUA"),
    call(esg=ANY, variable="CH4", institution_id="VUA"),
    call(esg=ANY, variable="CH4_openburning_share", institution_id="IAMC"),
]

EXPECTED_MINIMAL_META_HISTORIC_CALLS_V2 = [
    call(esg=ANY, variable="CH4_percentage_AGRI", institution_id="VUA"),
    call(esg=ANY, variable="CH4_percentage_BORF", institution_id="VUA"),
    call(esg=ANY, variable="CH4_percentage_DEFO", institution_id="VUA"),
    call(esg=ANY, variable="CH4_percentage_PEAT", institution_id="VUA"),
    call(esg=ANY, variable="CH4_percentage_SAVA", institution_id="VUA"),
    call(esg=ANY, variable="CH4_percentage_TEMF", institution_id="VUA"),
]

EXPECTED_MINIMAL_MODEL_CALLS_V2 = [
    call(esg=ANY, model="NorESM2-LM", project="CMIP6", variable="tas", experiment="ssp126"),
]


def delete_tmp_dir():
    shutil.rmtree(TEST_TMP_DIR, ignore_errors=True)


@pytest.fixture()
def mock_raw_input_single_var():
    with patch(DOWNLOAD_RAW_INPUT_SINGLE_VAR) as mock_function:
        yield mock_function


@pytest.fixture()
def mock_meta_historic_single_var():
    with patch(DOWNLOAD_META_HISTORIC_SINGLE_VAR) as mock_function:
        yield mock_function


@pytest.fixture()
def mock_model_single_var():
    with patch(DOWNLOAD_MODEL_SINGLE_VAR) as mock_function:
        yield mock_function


@pytest.fixture()
def mock_subprocess_run():
    with patch(SUBPROCESS_RUN) as mock_function:
        yield mock_function


@pytest.fixture
def input4mips_downloader_object():
    config_source = MINIMAL_DATASET_CONFIG_PATH
    config = create_input4mips_downloader_config_from_file(config_source)
    config.data_dir = TEST_TMP_DIR
    yield Input4MipsDownloader(config=config)
    delete_tmp_dir()


@pytest.fixture
def cmip6_downloader_object():
    config_source = MINIMAL_DATASET_CONFIG_PATH
    config = create_cmip6_downloader_config_from_file(config_source)
    config.data_dir = TEST_TMP_DIR
    yield CMIP6Downloader(config=config)
    delete_tmp_dir()


def test_downloader_init(input4mips_downloader_object, cmip6_downloader_object):
    assert isinstance(input4mips_downloader_object, Input4MipsDownloader)
    assert isinstance(cmip6_downloader_object, CMIP6Downloader)


def test_downloader_base_params(input4mips_downloader_object, cmip6_downloader_object):
    assert input4mips_downloader_object.config.project == INPUT4MIPS
    assert input4mips_downloader_object.config.experiments == ["historical", "ssp126"]
    assert cmip6_downloader_object.config.project == CMIP6
    assert cmip6_downloader_object.config.models == ["NorESM2-LM"]
    assert cmip6_downloader_object.config.experiments == ["ssp126"]


def test_downloader_max_possible_member_number(cmip6_downloader_object):
    assert cmip6_downloader_object.config.max_ensemble_members == MAX_ENSEMBLE_MEMBERS


def test_downloader_variables(input4mips_downloader_object, cmip6_downloader_object):
    assert cmip6_downloader_object.config.variables == ["tas"]
    assert input4mips_downloader_object.config.variables == [
        "CO2_em_anthro",
        "CO2_em_AIR_anthro",
        "CH4_em_openburning",
        "CH4_em_anthro",
        "CH4_em_AIR_anthro",
    ]
    assert input4mips_downloader_object.config.biomass_vars == ["CO2", "CH4"]
    assert input4mips_downloader_object.config.meta_vars_percentage == [
        "CH4_percentage_AGRI",
        "CH4_percentage_BORF",
        "CH4_percentage_DEFO",
        "CH4_percentage_PEAT",
        "CH4_percentage_SAVA",
        "CH4_percentage_TEMF",
    ]
    assert input4mips_downloader_object.config.meta_vars_share == ["CH4_openburning_share"]


def test_downloader_model_params(cmip6_downloader_object):
    # TODO refactor this test for new Node list
    assert cmip6_downloader_object.config.node_link in "https://esgf-node.llnl.gov/esg-search/"


def test_download_raw_input(input4mips_downloader_object, mock_raw_input_single_var, mock_meta_historic_single_var):
    input4mips_downloader_object.download()
    assert mock_raw_input_single_var.call_args_list == EXPECTED_MINIMAL_RAW_INPUT_CALLS
    assert mock_raw_input_single_var.call_count == RAW_INPUT_NUM_OF_CALLS
    assert mock_meta_historic_single_var.call_args_list == EXPECTED_MINIMAL_META_HISTORIC_CALLS
    assert mock_meta_historic_single_var.call_count == META_HISTORIC_NUM_OF_CALLS


def test_download_from_model(cmip6_downloader_object, mock_model_single_var):
    cmip6_downloader_object.download()
    assert mock_model_single_var.call_args_list == EXPECTED_MINIMAL_MODEL_CALLS
    assert mock_model_single_var.call_count == 1


@patch(DOWNLOAD_MODEL_SINGLE_VAR_V2)
@patch(DOWNLOAD_META_HISTORIC_SINGLE_VAR_V2)
@patch(DOWNLOAD_RAW_INPUT_SINGLE_VAR_V2)
def test_download_from_config_file(
    mock_raw_input_single_var_v2, mock_meta_historic_single_var_v2, mock_model_single_var_v2
):
    download_from_config_file(config_file=MINIMAL_DATASET_CONFIG_PATH)

    assert mock_raw_input_single_var_v2.call_args_list == EXPECTED_MINIMAL_RAW_INPUT_CALLS_V2
    assert mock_raw_input_single_var_v2.call_count == RAW_INPUT_NUM_OF_CALLS
    assert mock_meta_historic_single_var_v2.call_args_list == EXPECTED_MINIMAL_META_HISTORIC_CALLS_V2
    assert mock_meta_historic_single_var_v2.call_count == META_HISTORIC_NUM_OF_CALLS
    assert mock_model_single_var_v2.call_args_list == EXPECTED_MINIMAL_MODEL_CALLS_V2
    assert mock_model_single_var_v2.call_count == MODEL_SINGLE_NUM_OF_CALLS


def _assert_content_is_in_wget_script(mock_call, string_content):
    # This seems very obscure, so here's what's happening,
    # With the provided inputs, there should be only 1 call.
    # We then access the call's arguments. We are interested in
    # the content of the wget script that is generated, and we
    # want to make sure that for the same inputs, we get the same files
    call_list = mock_call.call_args_list
    first_and_only_call = call_list[0]
    call_arguments = first_and_only_call.args[0]
    wget_script_content = call_arguments[2]
    print(string_content)
    print(wget_script_content)
    assert string_content in wget_script_content


def test_download_raw_input_single_var(input4mips_downloader_object, mock_subprocess_run):
    download_subprocess = mock_subprocess_run
    input4mips_downloader_object.download_raw_input_single_var(variable="CO2_em_anthro", institution_id="PNNL-JGCRI")

    # These are partial file strings. Since we download multiple variables at the same time, it's very complicated
    # to specify versions for each without becoming cubbersome. Therefore, this test just looks for the file parts
    # That don't change once a new version gets published (which made this test crash and required updating this
    # variable
    expected_files = [
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-",
        "_gn_175001-179912.nc",
        "_gn_180001-184912.nc",
        "_gn_185001-189912.nc",
        "_gn_190001-194912.nc",
        "_gn_195001-199912.nc",
    ]
    download_subprocess.assert_called_once()
    for f in expected_files:
        _assert_content_is_in_wget_script(download_subprocess, f)


@pytest.mark.xfail
def test_download_meta_historic_biomassburning_single_var(input4mips_downloader_object, mock_subprocess_run):
    input4mips_downloader_object.download_meta_historic_biomassburning_single_var(
        variable="CH4_percentage_AGRI", institution_id="VUA"
    )

    expected_files = [
        "CH4-percentage-AGRI-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-2_gn_175001-201512.nc"
    ]
    mock_subprocess_run.assert_called_once()
    for f in expected_files:
        _assert_content_is_in_wget_script(mock_call=mock_subprocess_run, string_content=f)


def test_download_from_model_single_var(cmip6_downloader_object, mock_subprocess_run):
    cmip6_downloader_object.download()

    expected_files = [
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_201501-202012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_202101-203012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_203101-204012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_204101-205012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_205101-206012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_206101-207012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_207101-208012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_208101-209012.nc",
        "tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_209101-210012.nc",
    ]
    mock_subprocess_run.assert_called()
    for f in expected_files:
        _assert_content_is_in_wget_script(mock_call=mock_subprocess_run, string_content=f)
