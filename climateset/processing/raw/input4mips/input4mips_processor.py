import re
from pathlib import Path
from typing import Union

from climateset import PROCESSED_DATA
from climateset.processing.abstract_processor_step import AbstractProcessorStep
from climateset.processing.raw.abstract_raw_processor import AbstractRawProcessor
from climateset.processing.raw.checker import AbstractDirectoryChecker
from climateset.processing.raw.input4mips.processing import (
    CDO_SET_CALENDAR,
    CDO_SET_DAY,
    CDO_SET_TIME,
    CDO_VERTICAL_SUM,
    cdo_process_directory,
    correct_units_all_files_in_directory,
    create_anthro_fire_directory,
    create_emissions_totals,
    merge_ghg_directory,
    rename_biomassburning_files_in_directory,
    sum_sectors_directory,
)
from climateset.utils import create_logger

#
# Processing steps
#

# private step
_CDO_PROCESSING = "_cdo_processing"

# user available steps
CORRECT_NAMES = "correct_names"
REORDER_SSP_CO2 = "reorder_ssp_co2"
CREATE_FIRE_FILES = "create_fire_files"
CORRECT_UNITS = "correct_units"
CORRECT_CALENDAR = "correct_calendar"
CORRECT_TIME_AXIS = "correct_time_axis"
SUM_LEVELS = "sum_levels"
SUM_SECTORS = "sum_sectors"
CREATE_TOTALS = "create_totals"
MERGE_GHG = "merge_ghg"

AVAILABLE_INPUT4MIPS_PROCESSING_STEPS = frozenset(
    [
        CORRECT_NAMES,
        REORDER_SSP_CO2,
        CREATE_FIRE_FILES,
        CORRECT_UNITS,
        CORRECT_CALENDAR,
        CORRECT_TIME_AXIS,
        SUM_LEVELS,
        SUM_SECTORS,
        CREATE_TOTALS,
        MERGE_GHG,
    ]
)
ORDER = "order"
PROCESSOR_STEP = "processor_step"

LOGGER = create_logger(__name__)


class Input4MipsEmissionProcessor(AbstractRawProcessor):
    """
    Can be called to summarize emission over sectors.

    Only applied to Input4MIPs data.
    """

    def __init__(
        self,
        input_directory: Union[str, Path],
        working_directory: Union[str, Path] = PROCESSED_DATA / "input4mips",
        checker: AbstractDirectoryChecker = None,
        processing_steps=None,
        spat_load: str = "map",
        temp_load: str = "mon",
        ghg_list: list = None,
        ghg_sum_name: str = None,
        meta_data: Path = None,
        fire_cases: bool = True,
        force_sum_sectors: list = None,
        input_freq: str = "mon",
        input_res: str = "50_km",
    ):
        """Init emission / ghg processer
        Args:
            verify_units (bool): Makes sure that all ghg conform to the same units.
                Default false because slow.
            verify_resolution (bool): Checks if the res in the file name is the same like in the file.
                Default false because slow.
            verify_lonlat (bool): Checks if the longitude, latitude values are as expected.
                Default false because slow.
            verify_num_of_files (bool): Checks how many files are in the final dir.
                raises a warning if more or less than one.
            verify_variables (bool): Checks if the variables are also the expected vars.
            fix_co2 (bool): Fixes co2 data if necessary. Needed for ssp scenarios!
                Default false, because we are providing raw data where it's already fixed.
            correct_names (bool): renames all files that are expected
                to be biomassburning files. Attention! apply only to raw data,
                will mess things up if applied later!
            correct_units (bool): Default False, because it takes long. Corrects
                the units according to meta data saved in the raw param file.
            correct_calendar (bool): Makes sure that the right calendar is used.
            correct_time_axis (bool): Makes sure that each file starts at the first
                of each time-unit and shifts if not.
            sum_levels (bool): Makes sure that all emissions are sumed over the
                different levels and sectors. For AIR files, the emissions are
                summed over 25 different height levels. For ANTHRO files, the
                emissions are summed over 8 different sectors. Openburning and
                biomassburning have only one level (sectors are diff in external
                percentage files, see 'fire_cases').
            sum_sectors (bool): Makes sure that ghg are summarized over diff sectors.
                Historical data is treated differently. Different versions
                are stored in the output_dir for historical data, since different
                climate models need different versions here. There is a param/config file
                mapping the right climate models to the right sum_sector types.
            create_totals (bool): For some tasks you might not want to use
                emission files, but arrays of totals over a specific time frame.
                With this you can creat either map or total data; monthly or yearly.
                See spat_load and temp_load for that.
            spat_load (str): Can be "map" or "tot". Map is a long / lat map.
                Total is an array with the totals summed over the spatial axis.
            temp_load (str): Can be "mon" or "year". Decides if the loaded data
                is monthly or yearly averaged. Can be extended later on.
            merge_ghg (list or bool): If True:
                Simply merges all ghg that can be found in the sub_dir (e.g.).
                Instead one can hand over a list, which means that those ghg
                specifically will be merged into a new file.
            ghg_list (list): The greenhouse gases (ghg) that should be merged in
                merge_ghg. Should be the folder names used in the respective dir.
            meta_data (Path): Default None. Only necessary when sum_sectors is
                used and anthro fire files should be generated.
            fire_cases (bool): If the different fire cases should be considered
                during sum_sectors.
            force_sum_sectors (list): Which GHG are forced to be summed (sum_sectors) even if not all
                sectors (AIR, anthro, biomassburning/openburning) is available.
                Default is ["CO2"], since CO2 is the only GHG that does not have
                biomassburning/openburning files as the rest does.
            input_freq (str): Default "mon". Only used for create_totals_subdir.
            input_res (str): Default "50 km". Only used for create_totals_subdir.
            create_fire_files (bool): Default false, because it only needs to be
                created one time.
        """
        # init abstract class with checks
        super().__init__(
            input_directory=input_directory,
            working_directory=working_directory,
            checker=checker,
            processing_steps=processing_steps,
        )
        self.cdo_operators = []
        self.spat_load = spat_load
        self.temp_load = temp_load

        # ghg list that should be merged
        self.ghg_list = ghg_list
        if not self.ghg_list:
            self.ghg_list = ["BC_sum", "CH4_sum", "SO2_sum"]
        self.ghg_sum_name = ghg_sum_name
        self.meta_data = meta_data
        self.fire_cases = fire_cases
        self.force_sum_sectors = force_sum_sectors
        if not self.force_sum_sectors:
            self.force_sum_sectors = ["CO2"]
        self.input_freq = input_freq
        self.input_res = input_res

        self.available_steps = {
            CORRECT_NAMES: {ORDER: 1, PROCESSOR_STEP: RenameBiomassBurningFilesStep()},
            REORDER_SSP_CO2: {ORDER: 2, PROCESSOR_STEP: ReorderSSPCO2DimensionsStep()},
            CREATE_FIRE_FILES: {
                ORDER: 3,
                PROCESSOR_STEP: CreateAnthroFireFilesStep(
                    metadata_directory=self.meta_data, cdo_operators=self.cdo_operators
                ),
            },
            CORRECT_CALENDAR: {ORDER: 4, PROCESSOR_STEP: self._cdo_add_correct_calendar},
            CORRECT_TIME_AXIS: {ORDER: 5, PROCESSOR_STEP: self._cdo_add_correct_time_axis},
            SUM_LEVELS: {ORDER: 6, PROCESSOR_STEP: self._cdo_add_sum_levels},
            _CDO_PROCESSING: {ORDER: 7, PROCESSOR_STEP: CdoEmissionsProcessingStep(cdo_operators=self.cdo_operators)},
            CORRECT_UNITS: {ORDER: 8, PROCESSOR_STEP: CorrectUnitsStep(desired_units=self.desired_units)},
            SUM_SECTORS: {
                ORDER: 9,
                PROCESSOR_STEP: SumSectors(working_dir=self.working_directory, meta_data=self.meta_data),
            },
            CREATE_TOTALS: {
                ORDER: 10,
                PROCESSOR_STEP: CreateEmissionsTotals(
                    working_dir=self.working_directory,
                    spat=self.spat_load,
                    temp=self.temp_load,
                    input_res=self.input_res,
                    input_freq=self.input_freq,
                ),
            },
            MERGE_GHG: {
                ORDER: 11,
                PROCESSOR_STEP: MergeGHG(
                    working_dir=self.working_directory, ghg_list=self.ghg_list, ghg_sum_name=self.ghg_sum_name
                ),
            },
        }

        # TODO CO2 baseline - where??
        # CO2 baseline: do this outside, because it means that the data looks
        # different for each climate model. best case you just store that somewhere
        # in a config, store that well and apply it in the dataloader for each model separatetly?
        # self.substract_co2_baseline = substract_co2_baseline

    def add_correct_names_step(self):
        self.processing_steps.append(CORRECT_NAMES)

    def add_reorder_ssp_co2_dimensions_step(self):
        self.processing_steps.append(REORDER_SSP_CO2)

    def add_create_fire_files(self):
        self.processing_steps.append(CREATE_FIRE_FILES)

    def add_correct_units_step(self):
        self.processing_steps.append(CORRECT_UNITS)

    def add_correct_calendar_step(self):
        self.processing_steps.append(CORRECT_CALENDAR)

    def add_correct_time_axis_step(self):
        self.processing_steps.append(CORRECT_TIME_AXIS)

    def add_sum_levels_step(self):
        self.processing_steps.append(SUM_LEVELS)

    def add_sum_sectors_step(self):
        self.processing_steps.append(SUM_SECTORS)

    def add_create_totals_step(self):
        self.processing_steps.append(CREATE_TOTALS)

    def add_merge_ghg_step(self):
        self.processing_steps.append(MERGE_GHG)

    def set_load_types(self, spat_load: str = "map", temp_load: str = "mon"):
        """
        Setting types for the loading processor. Can be set to update the types, so the same processor instance can be
        applied afterward.

        Args:
            spat_load (str): Can be "map" or "tot". Map is a long / lat map.
                Total is an array with the totals summed over the spatial axis.
            temp_load (str): Can be "mon" or "year". Decides if the loaded data
                is monthly or yearly averaged. Can be extended later on.
        """
        self.spat_load = spat_load
        self.temp_load = temp_load

    def type_class_meta(self) -> str:
        """Returns the name tag of the subclass."""
        return "input4mips"

    # TODO share this in process class with input4mips / cmip6
    def file_belongs_to_type(self, input_file: Path) -> bool:
        """
        Check if the name input4mip(s) appears in the path name.

        Args:
            input_file (Path): the file that should be checked for input4mips
        Returns:
            True if name indicates it's a cmip6 file, False if not.
        """
        return bool(re.search("input4mip", str(input_file), re.IGNORECASE))

    def process_directory(self):
        """
        Applying all the stuff as decided in the init function.

        Most of the functions build upon each other, i.e. you cannot apply store_totals if you haven't used merge_ghg
        before. To keep things modular the user can still decide which ones to apply (e.g. because some things may have
        been applied earlier / are not necessary). Just be aware that things break if you are not making sure that the
        order is followed.
        """
        special_processes = [
            process for process in self.processing_steps if process in [CORRECT_CALENDAR, CORRECT_TIME_AXIS, SUM_LEVELS]
        ]
        if special_processes:
            for process in special_processes:
                self.available_steps[process][PROCESSOR_STEP]()
                self.processing_steps.remove(process)
            self.processing_steps.append(_CDO_PROCESSING)

        process_list = [self.available_steps[step] for step in self.processing_steps]
        ordered_processes = sorted(process_list, key=lambda step: step[ORDER])

        current_processing_directory = self.input_directory
        for process_dict in ordered_processes:
            processing_step: AbstractProcessorStep = process_dict[PROCESSOR_STEP]
            processing_step.execute(current_processing_directory)
            current_processing_directory = processing_step.get_results_directory()

        return current_processing_directory

    def _cdo_add_correct_calendar(self):
        self.cdo_operators.append(f"{CDO_SET_CALENDAR},{self.calendar}")

    def _cdo_add_correct_time_axis(self):
        # TODO add check for mon data only
        self.cdo_operators.append(CDO_SET_DAY)
        self.cdo_operators.append(CDO_SET_TIME)

    def _cdo_add_sum_levels(self):
        self.cdo_operators.append(CDO_VERTICAL_SUM)

    def set_meta_data(self, meta_data: Path):
        """Set meta data path if not done during initialization."""
        self.meta_data = meta_data


class RenameBiomassBurningFilesStep(AbstractProcessorStep):
    def execute(self, input_directory):
        """Renames files (biomassburning) if needed."""
        output_dir = rename_biomassburning_files_in_directory(input_directory)

        LOGGER.info(f"... Finished renaming biomassburning files in {input_directory}")
        self.results_directory = output_dir


class ReorderSSPCO2DimensionsStep(AbstractProcessorStep):
    def execute(self, input_directory):
        output_dir = rename_biomassburning_files_in_directory(input_directory)
        LOGGER.info(f"... Finished renaming biomassburning files in {input_directory}")
        self.results_directory = output_dir


class CreateAnthroFireFilesStep(AbstractProcessorStep):
    def __init__(self, metadata_directory, cdo_operators, overwrite=False):
        super().__init__()
        self.metadata_directory = metadata_directory
        self.cdo_operators = cdo_operators
        self.overwrite = overwrite

    def execute(self, input_directory):
        self.results_directory = create_anthro_fire_directory(
            input_directory,
            meta_data=self.metadata_directory,
            cdo_operators=self.cdo_operators,
            overwrite=self.overwrite,
        )


class CdoEmissionsProcessingStep(AbstractProcessorStep):
    def __init__(self, cdo_operators, overwrite=False):
        super().__init__()
        self.cdo_operators = cdo_operators
        self.overwrite = overwrite

    def execute(self, input_directory):
        self.results_directory = cdo_process_directory(
            input_directory=input_directory, cdo_operators=self.cdo_operators, overwrite=self.overwrite
        )


class CorrectUnitsStep(AbstractProcessorStep):
    def __init__(self, desired_units):
        super().__init__()
        self.desired_units = desired_units

    def execute(self, input_directory):
        """
        Correcting the ghg units (if possible).

        Attention, this module
        is slow since it is not handled with cdo. Do only if necessary.
        Args:
            input_dir (Path): All files in that dir will be adapted
        """
        output_dir = correct_units_all_files_in_directory(input_dir=input_directory, desired_units=self.desired_units)
        self.results_directory = output_dir


class SumSectors(AbstractProcessorStep):
    def __init__(self, working_dir, meta_data):
        super().__init__()
        self.meta_data = meta_data
        self.working_dir = working_dir

    def execute(self, input_directory):
        self.results_directory = sum_sectors_directory(
            input_dir=input_directory, output_dir=self.working_dir, meta_data=self.meta_data
        )


class CreateEmissionsTotals(AbstractProcessorStep):
    def __init__(
        self,
        working_dir,
        spat: str = "map",
        temp: str = "mon",
        input_freq: str = "mon",
        input_res: str = "50_km",
        overwrite: bool = True,
    ):
        super().__init__()
        self.working_dir = working_dir
        self.spat = spat
        self.temp = temp
        self.input_freq = input_freq
        self.input_res = input_res
        self.overwrite = overwrite

    def execute(self, input_directory):
        self.results_directory = create_emissions_totals(
            input_dir=input_directory,
            output_dir=self.working_dir,
            spat=self.spat,
            temp=self.temp,
            input_freq=self.input_freq,
            input_res=self.input_res,
            overwrite=self.overwrite,
        )


class MergeGHG(AbstractProcessorStep):
    def __init__(self, working_dir, ghg_list, ghg_sum_name, overwrite=False):
        super().__init__()
        self.working_dir = working_dir
        self.ghg_list = ghg_list
        self.ghg_sum_name = ghg_sum_name
        self.overwrite = overwrite

    def execute(self, input_directory):
        self.results_directory = merge_ghg_directory(
            input_dir=input_directory,
            output_dir=self.working_dir,
            ghg_list=self.ghg_list,
            ghg_sum_name=self.ghg_sum_name,
            overwrite=self.overwrite,
        )
