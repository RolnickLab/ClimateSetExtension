import pytest

from climateset import TEST_DIR
from climateset.download.downloader import Downloader
from climateset.utils import get_yaml_config


@pytest.fixture
def simple_downloader_object():
    config_source = TEST_DIR / "resources/test_minimal_dataset.yaml"
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
