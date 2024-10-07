from unittest.mock import call, patch

import pytest

from climateset import TEST_DIR
from climateset.download.downloader import Downloader, download_from_config_file
from climateset.utils import get_yaml_config

MINIMAL_DATASET_CONFIG_PATH = TEST_DIR / "resources/test_minimal_dataset.yaml"

DOWNLOAD_RAW_INPUT_SINGLE_VAR = "climateset.download.downloader.Downloader.download_raw_input_single_var"
DOWNLOAD_META_HISTORIC_SINGLE_VAR = (
    "climateset.download.downloader.Downloader.download_meta_historic_biomassburning_single_var"
)
DOWNLOAD_FROM_MODEL_SINGLE_VAR = "climateset.download.downloader.Downloader.download_from_model_single_var"
SUBPROCESS_RUN = "subprocess.run"

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

EXPECTED_MINIMAL_META_HISTORIC_CALLS = [
    call(variable="CH4_percentage_AGRI", institution_id="VUA"),
    call(variable="CH4_percentage_BORF", institution_id="VUA"),
    call(variable="CH4_percentage_DEFO", institution_id="VUA"),
    call(variable="CH4_percentage_PEAT", institution_id="VUA"),
    call(variable="CH4_percentage_SAVA", institution_id="VUA"),
    call(variable="CH4_percentage_TEMF", institution_id="VUA"),
]

EXPECTED_MINIMAL_MODEL_CALLS = [
    call(variable="tas", experiment="historical"),
    call(variable="tas", experiment="ssp126"),
]


@pytest.fixture()
def mock_raw_input_single_var():
    with patch(DOWNLOAD_RAW_INPUT_SINGLE_VAR) as mock_raw_input:
        yield mock_raw_input


@pytest.fixture()
def mock_meta_historic_single_var():
    with patch(DOWNLOAD_META_HISTORIC_SINGLE_VAR) as mock_meta_input:
        yield mock_meta_input


@pytest.fixture()
def mock_model_single_var():
    with patch(DOWNLOAD_FROM_MODEL_SINGLE_VAR) as mock_model:
        yield mock_model


@pytest.fixture()
def mock_subprocess_run():
    with patch(SUBPROCESS_RUN) as mock_subprocess_run:
        yield mock_subprocess_run


@pytest.fixture
def simple_downloader_object():
    config_source = MINIMAL_DATASET_CONFIG_PATH
    config = get_yaml_config(config_source)
    model = config["models"][0]
    downloader_kwargs = config["downloader_kwargs"]
    return Downloader(model=model, **downloader_kwargs)


def test_downloader_init(simple_downloader_object):
    assert isinstance(simple_downloader_object, Downloader)


def test_downloader_base_params(simple_downloader_object):
    assert simple_downloader_object.model == "NorESM2-LM"
    assert simple_downloader_object.experiments == ["historical", "ssp126"]


def test_downloader_max_possible_member_number(simple_downloader_object):
    assert simple_downloader_object.max_ensemble_members == 1


def test_downloader_variables(simple_downloader_object):
    assert simple_downloader_object.raw_vars == [
        "CO2_em_anthro",
        "CO2_em_AIR_anthro",
        "CH4_em_openburning",
        "CH4_em_anthro",
        "CH4_em_AIR_anthro",
    ]
    assert simple_downloader_object.biomass_vars == ["CO2", "CH4"]
    assert simple_downloader_object.model_vars == ["tas"]
    assert simple_downloader_object.meta_vars_percentage == [
        "CH4_percentage_AGRI",
        "CH4_percentage_BORF",
        "CH4_percentage_DEFO",
        "CH4_percentage_PEAT",
        "CH4_percentage_SAVA",
        "CH4_percentage_TEMF",
    ]
    assert simple_downloader_object.meta_vars_share == ["CH4_openburning_share"]


def test_downloader_model_params(simple_downloader_object):
    assert simple_downloader_object.model_node_link == "https://esgf-data.dkrz.de/esg-search"
    assert simple_downloader_object.model_source_center == "NCC"


def test_download_raw_input(simple_downloader_object, mock_raw_input_single_var, mock_meta_historic_single_var):
    mock_raw_input = mock_raw_input_single_var
    mock_meta_input = mock_meta_historic_single_var
    simple_downloader_object.download_raw_input()
    assert mock_raw_input.call_args_list == EXPECTED_MINIMAL_RAW_INPUT_CALLS
    assert mock_meta_input.call_args_list == EXPECTED_MINIMAL_META_HISTORIC_CALLS


def test_download_from_model(simple_downloader_object, mock_model_single_var):
    mock_model = mock_model_single_var
    simple_downloader_object.download_from_model()
    assert mock_model.call_args_list == EXPECTED_MINIMAL_MODEL_CALLS


def test_download_from_config_file(
    simple_downloader_object, mock_raw_input_single_var, mock_meta_historic_single_var, mock_model_single_var
):
    mock_raw_input = mock_raw_input_single_var
    mock_meta_input = mock_meta_historic_single_var
    mock_model = mock_model_single_var

    download_from_config_file(config=MINIMAL_DATASET_CONFIG_PATH)

    assert mock_raw_input.call_args_list == EXPECTED_MINIMAL_RAW_INPUT_CALLS
    assert mock_meta_input.call_args_list == EXPECTED_MINIMAL_META_HISTORIC_CALLS
    assert mock_model.call_args_list == EXPECTED_MINIMAL_MODEL_CALLS


def _assert_content_is_in_wget_script(download_subprocess, f):
    # This seems very obscure, so here's what's happening,
    # With the provided inputs, there should be only 1 call.
    # We then access the call's arguments. We are interested in
    # the content of the wget script that is generated, and we
    # want to make sure that for there inputs, we get the same files
    call_list = download_subprocess.call_args_list
    first_and_only_call = call_list[0]
    call_arguments = first_and_only_call.args[0]
    wget_script_content = call_arguments[2]
    assert f in wget_script_content


def test_download_raw_input_single_var(simple_downloader_object, mock_subprocess_run):
    download_subprocess = mock_subprocess_run
    simple_downloader_object.download_raw_input_single_var(variable="CO2_em_anthro", institution_id="PNNL-JGCRI")

    expected_files = [
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_175001-179912.nc",
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_180001-184912.nc",
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_185001-185012.nc",
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_185101-189912.nc",
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_190001-194912.nc",
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_195001-199912.nc",
        "CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_200001-201412.nc",
    ]
    download_subprocess.assert_called_once()
    for f in expected_files:
        _assert_content_is_in_wget_script(download_subprocess, f)


def test_download_meta_historic_biomassburning_single_var(simple_downloader_object, mock_subprocess_run):
    download_subprocess = mock_subprocess_run
    simple_downloader_object.download_meta_historic_biomassburning_single_var(
        variable="CH4_percentage_AGRI", institution_id="VUA"
    )

    expected_files = [
        "CH4-percentage-AGRI-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-2_gn_175001-201512.nc"
    ]
    download_subprocess.assert_called_once()
    for f in expected_files:
        _assert_content_is_in_wget_script(download_subprocess, f)


def test_download_from_model_single_var(simple_downloader_object, mock_subprocess_run):
    download_subprocess = mock_subprocess_run
    simple_downloader_object.download_from_model_single_var(variable="tas", experiment="ssp126")

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
    download_subprocess.assert_called_once()
    for f in expected_files:
        _assert_content_is_in_wget_script(download_subprocess, f)
