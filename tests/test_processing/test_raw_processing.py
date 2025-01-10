from unittest.mock import patch

import pytest
from xarray import open_dataset

from climateset import DATA_DIR
from climateset.processing.raw.checker import BasicDirectoryChecker
from climateset.processing.raw.input4mips.input4mips_processing import (
    AVAILABLE_INPUT4MIPS_PROCESSING_STEPS,
    Input4MipsEmissionProcesser,
)


@pytest.fixture()
def mock_xarray_load_dataset():
    with patch("xarray.load_dataset", new=open_dataset) as mock_function:
        yield mock_function


@pytest.fixture
def simple_emission_processing_object():
    data_dir = DATA_DIR / "raw/raw_input_vars/PNNL-JGCRI/CO2_em_anthro"
    return Input4MipsEmissionProcesser(input_directory=data_dir)


def test_emission_processing_init(simple_emission_processing_object):
    assert simple_emission_processing_object.input_directory == DATA_DIR / "raw/raw_input_vars/PNNL-JGCRI/CO2_em_anthro"
    assert isinstance(simple_emission_processing_object.checker, BasicDirectoryChecker)


def test_emission_processing_basic_directory_check(simple_emission_processing_object, mock_xarray_load_dataset):
    check_results = simple_emission_processing_object.checker.check_directory()
    assert check_results


def test_emission_processing_check_available_steps(simple_emission_processing_object):
    steps = simple_emission_processing_object.list_available_steps()
    is_step_available = [step in AVAILABLE_INPUT4MIPS_PROCESSING_STEPS for step in steps]
    assert all(is_step_available)


def test_emission_processing_add_processing_steps(simple_emission_processing_object):
    step_list = list(AVAILABLE_INPUT4MIPS_PROCESSING_STEPS)
    simple_emission_processing_object.add_processing_step(step_list)
    is_step_added = [step in step_list for step in simple_emission_processing_object.processing_steps]
    assert all(is_step_added)


def test_emission_processing_add_individual_steps(simple_emission_processing_object):
    simple_emission_processing_object.add_correct_names_step()
    simple_emission_processing_object.add_reorder_ssp_co2_dimensions_step()
    simple_emission_processing_object.add_create_fire_files()
    simple_emission_processing_object.add_correct_units_step()
    simple_emission_processing_object.add_correct_calendar_step()
    simple_emission_processing_object.add_correct_time_axis_step()
    simple_emission_processing_object.add_sum_levels_step()
    simple_emission_processing_object.add_sum_sectors_step()
    simple_emission_processing_object.add_create_totals_step()
    simple_emission_processing_object.add_merge_ghg_step()
    for step in simple_emission_processing_object.processing_steps:
        assert step in list(AVAILABLE_INPUT4MIPS_PROCESSING_STEPS)
