import os
import re
import subprocess
import warnings
from copy import copy
from pathlib import Path
from typing import Union

import pint
import xarray as xr
from tqdm import tqdm

from climateset.processing.raw.abstract_raw_processing import AbstractRawProcesser
from climateset.processing.raw.checker import AbstractDirectoryChecker
from climateset.processing.raw.utils import create_generic_output_path
from climateset.utils import create_logger

#
# Processing steps
#

# private step
_CDO_PROCESSING = "_cdo_processing"

# user available steps
CORRECT_NAMES = "correct_names"
FIX_CO2 = "fix_co2"
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
        FIX_CO2,
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

SET_DAY = "-setday,1"
SET_TIME = "-settime,00:00:00"

LOGGER = create_logger(__name__)


class Input4MipsEmissionProcesser(AbstractRawProcesser):
    """
    Can be called to summarize emission over sectors.

    Only applied to Input4MIPs data.
    """

    def __init__(
        self,
        input_directory: Union[str, Path],
        working_directory: Union[str, Path],
        checker: AbstractDirectoryChecker = None,
        processing_steps=None,
        spat_load: str = "map",
        temp_load: str = "mon",
        ghg_list: list = None,
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

        self.available_steps = {
            CORRECT_NAMES: {"order": 1, "function": self.correct_names_subdir},
            FIX_CO2: {"order": 2, "function": self.fix_co2_subdir},
            CREATE_FIRE_FILES: {"order": 3, "function": self.create_anthro_fire_subdir},
            CORRECT_UNITS: {"order": 4, "function": self.correct_units_subdir},
            CORRECT_CALENDAR: {"order": 5, "function": self._correct_calendar},
            CORRECT_TIME_AXIS: {"order": 6, "function": self._correct_time_axis},
            SUM_LEVELS: {"order": 7, "function": self._sum_levels},
            _CDO_PROCESSING: {"order": 8, "function": self.cdo_preprocess_subdir},
            SUM_SECTORS: {"order": 9, "function": self.sum_sectors_subdir},
            CREATE_TOTALS: {"order": 10, "function": self.create_totals_subdir},
            MERGE_GHG: {"order": 11, "function": self.merge_ghg_subdir},
        }

        self.cdo_operators = []
        self.spat_load = spat_load
        self.temp_load = temp_load

        # ghg list that should be merged
        self.ghg_list = ghg_list
        if not self.ghg_list:
            self.ghg_list = ["BC_sum", "CH4_sum", "SO2_sum"]
        self.meta_data = meta_data
        self.fire_cases = fire_cases
        self.force_sum_sectors = force_sum_sectors
        if not self.force_sum_sectors:
            self.force_sum_sectors = ["CO2"]
        self.input_freq = input_freq
        self.input_res = input_res

        # TODO CO2 baseline - where??
        # CO2 baseline: do this outside, because it means that the data looks
        # different for each climate model. best case you just store that somewhere
        # in a config, store that well and apply it in the dataloader for each model separatetly?
        # self.substract_co2_baseline = substract_co2_baseline

    def add_correct_names_step(self):
        self.processing_steps.append(CORRECT_NAMES)

    def add_fix_co2_step(self):
        self.processing_steps.append(FIX_CO2)

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
        Setting types for the loading processer. Can be set to update the types, so the same processer instance can be
        applied afterwards.

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

    def preprocess_subdir(
        self,
        input_dir: Path,
        cleaned_dir: Path,
        processed_dir: Path,
        load_dir: Path,
        overwrite: bool = False,
        silent: bool = False,
        sum_sec_input_res: str = "50_km",
        sum_sec_input_freq: str = "mon",
    ):
        """
        Applying all the stuff as decided in the init function. Most of the functions build upon each other, i.e. you
        cannot apply store_totals if you haven't used merge_ghg before. To keep things modular the user can still decide
        which ones to apply (e.g. because some things may have been applied earlier / are not necessary). Just be aware
        that things break if you are not making sure that the order is followed.

        Args:
            input_dir (Path): To which directory the processing should be applied.
            cleaned_dir (Path): Where cleaned data should be stored. This is used
                for the "preprocess" steps.
            processed_dir (Path): Where processed data is stored. This is used
                for sum_sectors. Simple preprocessing is directly applied on
                raw data here.
            load_dir (Path): Where data is stored that is strongly modified and
                ready to be loaded. This is used for merge_ghg and store_totals.
            overwrite (bool): If the data should be overwritten if it already
                exists. Default False.
            silent (bool): If this should be processed silently.
            sum_sec_input_res:
            sum_sec_input_freq:
        """
        if any([process in self.processing_steps for process in [CORRECT_CALENDAR, CORRECT_TIME_AXIS, SUM_LEVELS]]):
            self.processing_steps.append(_CDO_PROCESSING)

        process_list = [self.available_steps[step] for step in self.available_steps]
        ordered_processes = sorted(process_list, key=lambda step: step["order"])

        for _, process in ordered_processes:
            # process["function"]()
            print(process)

    def fix_co2_subdir(self, input_dir: Path):
        """
        Fix CO2 files.

        Args:
            input_dir (Path): All ssp CO2 files in that dir will be adapted
        """
        print(f"Starting to fix co2 files in {input_dir}.")
        total_files = len(list(input_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
            if len(files) > 0:
                for file in files:
                    if ("ssp" in file) and ("CO2" in file):
                        self._fix_co2(Path(path) / file)

        print("...Finished fixing CO2 files.")

    def _fix_co2(self, input_file: Path):
        """
        Fix single co2 file.

        CDO can be applied again afterwards.
        Memory issues might arise (core dumped). Not sure why. Runs on compute canada though.
        Args:
            input_file (Path): The co2 file that should be fixed.
        """
        # ncpdq --ovr -a time,level,lon,lat input.nc output.nc
        commands = ["ncpdq", "--ovr", "-a", "time,level,lon,lat", input_file, input_file]
        subprocess.call(commands)

    def correct_names_subdir(self, input_dir: Path):
        """
        Renames files (biomassburning) if needed.

        Args:
            input_dir (Path): All files in that dir will be adapted
        """
        print(f"Starting to rename biomassburning files in {input_dir}.")
        total_files = len(list(input_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
            if len(files) > 0:
                for file in files:
                    self._rename_biomassburning_files(Path(path) / file)

        # remove empty subdirs
        for p in input_dir.glob("**/*"):
            if p.is_dir() and len(list(p.iterdir())) == 0:
                os.removedirs(p)

        print(f"... Finished renaming biomassburning files in {input_dir}")

    # TODO test
    def correct_units_subdir(self, input_dir: Path):
        """
        Correcting the ghg units (if possible).

        Attention, this module
        is slow since it is not handled with cdo. Do only if necessary.
        Args:
            input_dir (Path): All files in that dir will be adapted
        """
        print(f"Starting to correct all units (if needed) in {input_dir}.")
        total_files = len(list(input_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
            if len(files) > 0:
                for file in files:
                    self._correct_units(Path(path) / file)

        print(f"... Finished correcting the units (if needed) in {input_dir}")

    def _rename_biomassburning_files(self, input_file: Path):
        """
        Rename file for biomassburning in case it is in the naked version (from the downloader). Originally, the files
        are just called 'GHG' instead of GHG_em_biomassburning. Overwrites the files directly. Attention,

        apply this only directly to the raw data!! The Input4mips preprocesser
        renames variables later to the pure / naked GHG version - do not rename
        them!
        Args:
            input_file (Path): File that should be renamed
        """
        # fetch the ghg
        ghg = str(input_file.parents[0]).split("/")[-4]
        # check if this is the naked case
        if "_" not in ghg:
            new_ghg = f"{ghg}_em_biomassburning"
            new_path_name = Path(str(input_file).replace(ghg, new_ghg))
            new_path_name.parent.mkdir(parents=True, exist_ok=True)
            # use cdo to rename the var internally
            commands = [
                "cdo",
                "-s",
                "-L",
                f"-chname,{ghg},{ghg}_em_biomassburning",
                f"-setattribute,variable_id={ghg}_em_biomassburning",
                input_file,
                new_path_name,
            ]
            subprocess.call(commands)

            # remove the old file
            input_file.unlink()

    # TODO test
    def _correct_units(self, input_file: Path):
        """
        Corrects units of a single file.

        Be aware this is slow if applied
        to a large directory. Overwrites the old file
        Args:
            input_file (Path): The file whose units should be adapted.
        """
        ureg = pint.UnitRegistry()
        ds = xr.load_dataset(input_file)

        # only use the data variables (ignore bnds vars)
        for var in ds.data_vars:
            if "bnds" not in str(var):
                # check the unit
                found_unit = ds[var].units.replace("-", "^-")
                desired_unit = self.desired_units[var].replace("-", "^-")

                # exit if it's the right one
                if found_unit == desired_unit:
                    return

                else:
                    try:
                        ureg(found_unit)
                        ureg(desired_unit)
                    except pint.errors.UndefinedUnitError:
                        print(
                            "Unit is not defined for pint. Check here: https://github.com/hgrecco/pint/blob/master/pint/default_en.txt"
                        )
                        return
                    # update ds
                    # this function works for any kind of unit transformations (if they can be transformed)
                    ds.update(
                        {
                            var: xr.apply_ufunc(
                                lambda x: ureg.Quantity(x, found_unit).to(desired_unit).magnitude, ds[var]
                            )
                        }
                    )
                    # update attrs
                    ds[var].attrs["units"] = desired_unit
                    # overwrite old file
                    ds.to_netcdf(input_file, mode="w", format="NETCDF4", engine="netcdf4")

    def merge_ghg_subdir(
        self, input_dir: Path, output_dir: Path, ghg_list: list, ghg_sum_name: str = "ALL", overwrite: bool = False
    ):
        """
        Merging GHG emission files together.

        Can only be applied after
        'create_totals_subdir'!
        Args:
            input_dir (Path): Data dir that should be processed
            output_dir (Path): Where the merged data should be stored
            ghg_sum_name (str): How the dir and files where the merged data is
                stored should be tagged.
            overwrite (bool): If files should be overwritten when they already exist
        """
        print(f"Start creating merged GHG emission files for {ghg_list} in {input_dir}...")
        # create a dict with the ghg files that should be merged
        file_dict = {}
        # total_files = len(list(input_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(input_dir)):
            if len(files) > 0:
                for file in files:
                    full_path = str(Path(path) / file)
                    # do this only for the ghg listed above:
                    for ghg in ghg_list:
                        if ghg in full_path:
                            ghg_file_str = full_path.replace(ghg, ghg_sum_name)
                            # add ghg group as key if they don't exist yet
                            if not (ghg_file_str in file_dict):
                                file_dict[ghg_file_str] = [ghg]
                            else:
                                file_dict[ghg_file_str].append(ghg)

        # create the needed dirs
        for file, ghgs in file_dict.items():
            # create output path
            second_path_part = file.split("input4mips/")[-1]
            out_path = Path(output_dir / "input4mips" / second_path_part)
            out_path.parent.mkdir(parents=True, exist_ok=True)

            # create input files
            input_files = [file.replace(ghg_sum_name, ghg) for ghg in ghgs]

            if (not out_path.is_file()) or overwrite:
                self._merge_ghg_files(
                    input_files,
                    out_path,
                )

        print(f"... Finished merging GHG {ghg_list} in {output_dir}.")

    def _merge_ghg_files(self, input_files: list, output_file: Path, threads: int = 1):
        """
        Merging GHG emission files together.

        Can only be applied after
        'create_totals_subdir'!
        Args:
            input_files (list): List of files (Path) that should be merged.
            output_file (Path): Where the merged data should be stored
            threads (int): How many threads should be used.
        """
        if output_file.is_file():
            output_file.unlink()

        commands = ["cdo", "-w", "-s", "-L", "-P", str(threads), "-merge", *input_files, output_file]

        subprocess.call(commands)

    def create_totals_subdir(
        self,
        input_dir: Path,
        output_dir: Path,
        spat: str = "map",
        temp: str = "mon",
        input_freq: str = "mon",
        input_res: str = "50_km",
        overwrite: bool = True,
    ):
        """
        Creating the data that can actually be loaded by pytorch.

        Args:
            input_dir (Path): Data dir that should be processed
            output_dir (Path): Where the merged data should be stored
            spat (str): Can be "map" or "total". Map stores ghg emission maps
                (longitude / latitude). Total sums the values over the spatial axes.
            temp (str): Can be "mon" or "year". In case of year the data is summed
                over years. Assumes that monthly data is used!
            input_freq (str): The frequency of the input data (can only process)
                one type after another, so this must be fixed. Default "mon".
            input_res (str): The nominal resolution of the input data (can only process)
                one type after another, so this must be fixed. Default "50_km".
            overwrite (bool): If files should be overwritten when they already exist
        """
        print(f"Start creating loadable emission data {input_dir}....")
        if spat not in ["map", "total"]:
            raise ValueError("Parameter spat must be either 'map' or 'total'.")
        if temp not in ["mon", "year"]:
            raise ValueError("Parameter temp must be either 'mon' or 'year'.")

        # Match all the files that belong together and contain mon (only difference is the year ending!)
        # TODO move this function into an external one
        file_dict = {}
        # total_files = len(list(input_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(input_dir)):
            if len(files) > 0:
                for file in files:
                    full_path = str(Path(path) / file)
                    if (input_res in full_path) and (input_freq in full_path):
                        year = int(file.split(".")[0].split("_")[-1])
                        year_file_str = full_path.replace(str(year), "YEAR")
                        # add ghg group as key if they don't exist yet
                        if year_file_str not in file_dict:
                            file_dict[year_file_str] = [year]
                        else:
                            file_dict[year_file_str].append(year)

        # create the needed dirs
        for file, years in file_dict.items():
            # create output path
            second_path_part = file.split("input4mips/")[-1]
            second_path_part = second_path_part.replace(input_res, f"{spat}_{input_res}")
            second_path_part = second_path_part.replace(input_freq, f"{temp}")
            second_path_part = second_path_part.replace("YEAR", f"{min(years)}-{max(years)}")
            out_path = Path(output_dir / "input4mips" / second_path_part)
            out_path.parent.mkdir(parents=True, exist_ok=True)

            # create input files
            input_files = [file.replace("YEAR", str(y)) for y in years]

            if (not out_path.is_file()) or overwrite:
                self._create_totals(
                    input_files,
                    out_path,
                    spat,
                    temp,
                )

        print(f"... Finished creating loadable emission data {output_dir}....")

    def _create_totals(self, input_files: list, output_path: Path, spat: str, temp: str, threads: int = 1):
        """
        Creating totals for a given list of files for different spatial and temporal settings.

        Args:
            input_files (list): List of paths that are used to create totals
            output_path (Path): Where the merged totals should be stored.
            spat (str): Can be "map" or "total". Map stores ghg emission maps
                (longitude / latitude). Total sums the values over the spatial axes.
            temp (str): Can be "mon" or "year". In case of year the data is summed
                over years. Assumes that monthly data is used!
            threads (int): number of threads
        """
        # remove the file if it already exists
        if output_path.is_file():
            output_path.unlink()

        # first basic command part
        commands = [
            "cdo",
            "-w",
            "-s",
            "-L",
            "-P",
            str(threads),
        ]

        # add depending on cases
        if spat == "total":
            # commands.append("-fldsum")
            # ATTTENTION this assumes that sum sectors has been applied before!!
            ghg = str(input_files[0]).split("/")[-1].split("_")[2]
            commands.append("expr,sum=fldsum({0});min=fldmin({0});max=fldmax({0});mean=fldmean({0});".format(ghg))

        if temp == "year":
            commands.append("-yearsum")

        # add relevant commands for all cases
        commands.extend(["-mergetime", *input_files, output_path])

        subprocess.call(commands)

    def _correct_calendar(self):
        self.cdo_operators.append(f"-setcalendar,{self.calendar}")

    def _correct_time_axis(self):
        # TODO add check for mon data only
        self.cdo_operators.append(SET_DAY)
        self.cdo_operators.append(SET_TIME)

    def _sum_levels(self):
        self.cdo_operators.append("-vertsum")

    # TODO if a year / file is missing -> write a warning!
    # TODO LATER add feature to check lon / lat consistency for input4mips
    def cdo_preprocess_subdir(self, sub_dir: Path, output_dir: Path, overwrite: bool = False):
        """
        Apply desired emission processing steps to a subdir and store it in output_dir.

        Can only be applied on raw input4mpis data, not on those
        that have already been processed (i.e. merged / summed up).
        Args:
            sub_dir (Path): Where the preprocessing should be applied
            output_dir (Path): Where the results should be stored
            overwrite (bool): If the data should be overwritten if it already exists
        """
        # 1. apply the "preprocessing" - all the functions that can
        # be applied to a single file
        print(f"Start preprocessing of emission files {sub_dir}....")
        total_files = len(list(sub_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(sub_dir), total=total_files):
            if len(files) > 0:
                for file in files:
                    # create output dir
                    input_file = Path(path) / file
                    if self.file_belongs_to_type(input_file):
                        output_file = create_generic_output_path(output_dir, path, file)
                        # skip if file already exists and we dont wanna overwrite it
                        if (not output_file.is_file()) or overwrite:
                            print(f"\nProcessing the following file: {input_file}")
                            self._cdo_preprocess(input_file, output_file)

        print(f"...Finished preprocessing \nof {sub_dir} and saved it at {output_dir}.")

    def _cdo_preprocess(self, input_file: Path, output_file: Path, threads: int = 1, cdo_operators=None):
        """Applies all the emission processing tasks defined during init on
        a given netcdf file - but only those ones that can operate on a single
        file.
        Args:
            input_file (Path): The file that should be processed.
            output_file (Path): Where the resulting file is stored.
            threads (int): threads for cdo
        """
        # RENAME variables in case of biomassburning
        # needs to be done beforehand becauses messes too much stuff up
        # is done immediately in place
        if "em_biomassburning" in str(input_file):
            ds = xr.load_dataset(input_file)
            var = ds.variable_id
            old_name = str(input_file)
            renamed_file = Path(old_name.replace(".nc", "_renamed.nc"))
            if "em_biomassburning" not in var:
                commands = [
                    "cdo",
                    "-s",
                    "-L",
                    f"-chname,{var},{var}_em_biomassburning",
                    f"-setattribute,variable_id={var}_em_biomassburning",
                    input_file,
                    renamed_file,
                ]
                subprocess.call(commands)
                # rename files & reassign the input_file
                input_file.unlink()
                input_file = renamed_file.rename(old_name)

        chained_operators = copy(self.cdo_operators)
        # cdo vertsum infile outfile [CHECKED]
        # AIR: sum over height levels
        # anthro: sum over different sectors

        # ATTENTION only relevant for "monthly" data!
        # cdo -setday,1 -settime,00:00:00
        # TODO add check for mon data only
        if "mon" not in "":
            chained_operators.remove(SET_DAY)
            chained_operators.remove(SET_TIME)
            LOGGER.warning("Tried to set the time axis on non-monthly file; removing related cdo operators for file.")

        # run commands
        commands = ["cdo", "-s", "-w", "-L", "-P", str(threads), *chained_operators, input_file, output_file]
        subprocess.call(commands)

        # unit_key_list = [
        #     "units", "unit",
        #     "variable_units", "variable_unit",
        #     "var_units", "var_unit"
        # ]

        # for k in unit_key_list:
        #     try:
        #         # TODO try to get the units
        #         pass
        #         # if successful: break loop
        #     except AttributeError:
        #         # TODO continue
        #         continue
        # if unit is None:
        #     raise KeyError("Unit could not be found in file {}".format(input_file))
        # else:
        #     # TODO check if the unit is the one saved in the jason file
        #     calendar = self.meta_raw_dict["calendar"][var]
        #     pass
        # simply raises an error in case we encounter units that we dont expect

    # TODO break this into two functions
    def sum_sectors_subdir(
        self,
        sub_dir: Path,
        output_dir: Path,
        force_sum_sectors: list = None,
        fire_cases: bool = False,
        overwrite: bool = True,
        silent: bool = False,
        input_res: str = "50_km",
        input_freq: str = "mon",
    ):
        """
        Args:
            force_sum_sectors (str): Default is empty, i.e. no merge is forced.
                Put the name of GHG here that should be forced to merge.
            fire_cases (bool): Default False means that simply all the data  is
                summed up and that's it. In case this is set to true, three different
                things are calculated: anthro-fires, no-fires, all-fires. Can
                only be done if "em_openburning" or "em_biomassburning" files
                are included.
            overwrite (bool): if the files should be overwritten if they already exist
            silent (bool): If this should be processed silently.
        """
        print(f"Start summing up different types of emission files in {sub_dir}...")
        if not force_sum_sectors:
            force_sum_sectors = []
        # params of the files to sum over
        # TODO make these params in a config
        # TODO find out which files are skipped
        # res = "50_km"
        # freq = "mon"
        # years = []
        # dicts needed to match files together
        # "example_path_SECTOR_sth": [None, None, None]
        ssp_file_dict = {}
        historical_file_dict = {}
        sector_idx = {"em_AIR_anthro": 0, "em_anthro": 1, "em_openburning": 2, "em_biomassburning": 2}

        # match all files that have the exact same name and end on
        total_files = len(list(sub_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(sub_dir), total=total_files):
            if len(files) > 0:
                for file in files:
                    full_path = str(Path(path) / file)
                    # only process the files with the right res and frequency
                    if (input_res in full_path) and (input_freq in full_path):
                        # sector descriptions
                        # check if this is historical or ssp and match
                        if "ssp" in file:
                            sector = re.search(r"(em_AIR_anthro|em_anthro|em_openburning)", full_path).group(0)
                            sector_file_str = full_path.replace(sector, "SECTOR")
                            # add ghg group as key if they don't exist yet
                            if sector_file_str not in ssp_file_dict:
                                ssp_file_dict[sector_file_str] = [None, None, None]
                            # add the found sector at the right sector index
                            # ssp_file_dict[sector_file_str][sector_idx[sector]] = full_path
                            ssp_file_dict[sector_file_str][sector_idx[sector]] = sector

                        elif "historical" in file:
                            sector = re.search(r"(em_AIR_anthro|em_anthro|em_biomassburning)", full_path).group(0)
                            sector_file_str = full_path.replace(sector, "SECTOR")
                            # add ghg group as key if they don't exist yet
                            if sector_file_str not in historical_file_dict:
                                historical_file_dict[sector_file_str] = [None, None, None]
                            # add the sector at the right sector index
                            # history_file_dict[sector_file_str][sector_idx[sector]] = full_path
                            historical_file_dict[sector_file_str][sector_idx[sector]] = sector

                        else:
                            raise ValueError(
                                "Input4mips files must include the scenario ('ssp' or 'historical') in their name."
                            )

        self._sum_sector_cases(
            file_dict=ssp_file_dict,
            output_dir=output_dir,
            fire_cases=fire_cases,
            force_sum_sectors=force_sum_sectors,
            overwrite=False,
            silent=silent,
        )

        self._sum_sector_cases(
            file_dict=historical_file_dict,
            output_dir=output_dir,
            fire_cases=fire_cases,
            force_sum_sectors=force_sum_sectors,
            overwrite=False,
            silent=silent,
        )

        print(f"...Finished summing up ghg emissions and stored results in {output_dir}")

    def _sum_sector_cases(
        self,
        file_dict: dict,
        output_dir: Path,
        fire_cases: bool,
        force_sum_sectors: list = None,
        overwrite: bool = False,
        silent: bool = False,
    ):
        """
        Calls sum sector for different cases.

        Args:
            file_dict (dict):  "example_path_SECTOR_sth": [file_path, file_path, file_path]
            output_dir (Path): output directory, where to save the diff cases.
            fire_cases (bool): if the diff fire cases should be considered
            force_sum_sectors (str): Default is empty, i.e. no merge is forced.
                Put the name of GHG here that should be forced to merge.
            overwrite (bool): if the files should be overwritten if they already exist
            silent (bool): If this should be processed silently.
        """
        if not force_sum_sectors:
            force_sum_sectors = []

        for file, sectors in tqdm(file_dict.items()):
            num_nones = sum([sec is None for sec in sectors])

            # create file lists and names for the case where everything is added
            ghg = file.split("input4mips/")[-1].split("/")[1].split("_")[0]
            input_files = [Path(file.replace("SECTOR", sec)) for sec in sectors if sec is not None]
            sector_list = [f"{ghg}_{sec}" for sec in sectors if sec is not None]
            out_path = Path(output_dir / "input4mips" / file.split("input4mips/")[-1].replace("SECTOR", "sum"))
            out_path.parent.mkdir(parents=True, exist_ok=True)

            if num_nones == 0:
                # normal case
                if fire_cases:
                    self.sum_sectors_fires(input_files, sector_list, out_path, overwrite=overwrite)
                else:
                    self._sum_sectors(input_files, sector_list, out_path, overwrite=overwrite)

            elif (num_nones == 1) and (ghg in force_sum_sectors):
                # just the simple summing if desired + warning
                if not silent:
                    warnings.warn(f"Merged files for {file}, however some sectors are missing", stacklevel=2)
                self._sum_sectors(input_files, sector_list, out_path, overwrite=overwrite)

            elif (num_nones == 1) and ghg not in force_sum_sectors:
                if not silent:
                    warnings.warn(
                        "Skipping case {}. Not all three sectors are available. You can force merging with 'force_sum_sectors=[GHG1, GHG2]'.".format(
                            file
                        ),
                        stacklevel=2,
                    )

            elif num_nones == 2:
                # warning, skipping
                if not silent:
                    warnings.warn(
                        f"Skipping case {file}. At least two files must exist to be summed together.",
                        stacklevel=2,
                    )
            else:
                raise RuntimeError("'num_nones' must be between 0 and 3.")

    def set_meta_data(self, meta_data: Path):
        """Set meta data path if not done during initialization."""
        self.meta_data = meta_data

    def create_anthro_fire_subdir(self, input_dir: Path, overwrite: bool = False):
        """
        Creating anthropogenic fire files for a subdir.

        Args:
            input_dir (Path): Input dir
            overwrite (bool): If data should be overwritten. Default False.
        """
        # loop through input dir, get biomassburning and openburning files
        print(f"Starting to create anthropogenic fire files for files in {input_dir}.")
        fire_output_dir = self.meta_data / "anthro-fire-data"
        total_files = len(list(input_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
            if len(files) > 0:
                for file in files:
                    input_file = Path(path) / file
                    # must match the right nominal and temporal resolution
                    if ("biomassburning" in file) and ("25_km" in file) and ("mon" in file):
                        type = "biomassburning"
                    elif ("openburning" in file) and ("50_km" in file) and ("mon" in file):
                        type = "openburning"
                    else:
                        continue
                    scenario = file.split("_")[1]
                    ghg = file.split("_")[2]
                    fire_name = file.replace("_em_", "_anthro_")
                    output_file = Path(fire_output_dir / type / scenario / ghg / fire_name)
                    output_file.parent.mkdir(parents=True, exist_ok=True)
                    if (not output_file.is_file()) or overwrite:
                        self._create_anthro_fire_file(input_file, type, output_file)

        print(f"... Finished creating anthro fire files and stored them in {fire_output_dir}")

    def _create_anthro_fire_file(self, input_file: Path, type: str, output_file: Path, threads: int = 1):
        """
        Create a new openburning or biomassburning file from a given file.

        Returns the path where the new file lifes.
        Args:
            input_file (Path):
            type (str): either 'biomassburning' or 'openburning'
            output_file (Path): Where to store the anthro fire file
            threads (int): threads for cdo
        """
        if self.meta_data is None:
            print("You need to provide meta data for this function. Use set_meta_data func for that.")
        file_parts = str(input_file.stem).split("_")
        scenario, ghg, year = file_parts[1], file_parts[2], file_parts[-1]
        meta_dir = self.meta_data / type

        # prepare additional operators (calendar and time-axis) for later compatability
        chained_operators = []
        # correct calendar and time axis right away if necessary
        if self.correct_time_axis:
            chained_operators.append("-setday,1")
            chained_operators.append("-settime,00:00:00")
        if self.correct_calendar:
            # default calendar: 365_day
            chained_operators.append(f"-setcalendar,{self.calendar}")

        if type == "openburning":
            # find a file that contains openburning, scenario and ghg in the meta_data / "future-openburning"
            meta_dir = self.meta_data / "future-openburning"
            matched_files = [
                f for f in meta_dir.rglob("*.nc") if (scenario in str(f)) and (type in str(f)) and (ghg in str(f))
            ]
            # take first file that matches that
            raw_percentage_file = matched_files[0]
            tmp_perc_file = input_file.parents[0] / f"tmp_percentage_{year}.nc"

            # 1. drop the years you don't need from the meta file. -selyear,2015
            # 2. drop the levels you dont need - sellevel,0,1
            # 3. sum the remaining levels / sectors -vertsum
            # -> cdo -vertsum -sellevel,0,1 -selyear,year meta_input.nc 2015_share_output.nc
            commands = [
                "cdo",
                "-w",
                "-s",
                "-L",
                "-P",
                str(threads),
                "-vertsum",
                "-sellevel,0,1",  # 0 is agriculture, 1 is deforestation
                f"-selyear,{year}",
                raw_percentage_file,
                tmp_perc_file,
            ]
            subprocess.call(commands)

            # multiply -> cdo mul fire_file.nc percentage_file new_anthro_fire_file.nc
            # [IMPORTANT first file must be input-file -> gets meta data from first input file!]
            commands = [
                "cdo",
                "-w",
                "-s",
                "-L",
                "-P",
                str(threads),
                *chained_operators,
                "-mul",
                input_file,
                tmp_perc_file,
                output_file,
            ]
            subprocess.call(commands)

            # remove 2015_share_output file to save space (unlink)
            tmp_perc_file.unlink()

        elif type == "biomassburning":
            meta_dir = self.meta_data / "historic-biomassburning"
            matched_files = [
                f for f in meta_dir.rglob("*.nc") if (type in str(f)) and (ghg in str(f)) and (year in str(f))
            ]

            # find the agri / defo files (written this way for later extension)
            perc_dict = {
                "AGRI": None,
                "DEFO": None,
            }
            for f in matched_files:
                for key in perc_dict.keys():
                    if key in str(f):
                        perc_dict[key] = f
            if (perc_dict["AGRI"] is None) or (perc_dict["DEFO"] is None):
                warnings.warn("The AGRI or DEFO file does not exist, see dictionary: \n", stacklevel=2)
                print(perc_dict)
                print("... skipping this file.")
                return

            # cdo 'expr,percentage_AGRI=percentage_AGRI/100' add percentage_AGRI_year.nc percentage_DEFO_year.nc percentage_ANTHRO_year.nc
            tmp_anthro_perc_file = input_file.parents[0] / f"tmp_percentage_anthro_{year}.nc"
            commands = [
                "cdo",
                "-w",
                "-s",
                "-L",
                "-P",
                str(threads),
                "expr,percentage_AGRI=percentage_AGRI/100;",
                "-add",
                perc_dict["AGRI"],
                perc_dict["DEFO"],
                tmp_anthro_perc_file,
            ]
            subprocess.call(commands)

            # cdo mul input_file.nc percentage_ANTHRO_year.nc output.nc
            commands = [
                "cdo",
                "-w",
                "-s",
                "-L",
                "-P",
                str(threads),
                *chained_operators,
                "-mul",
                input_file,
                tmp_anthro_perc_file,
                output_file,
            ]
            subprocess.call(commands)

            # delete temp files
            tmp_anthro_perc_file.unlink()

            # BIOMASSBURNING
            # TODO get the right biomassburning files: year, GHG, biomassburning and AGRI + DEFO
            # 1. cdo 'expr,percentage_AGRI=percentage_AGRI/100' -selyear,year meta_input.nc percentage_AGRI_year.nc
            # 2. cdo 'expr,percentage_DEFO=percentage_DEFO/100' -selyear,year meta_input.nc percentage_DEFO_year.nc
            # 3. cdo add percentage_AGRI_year.nc percentage_DEFO_year.nc percentage_ANTHRO_year.nc
            # 4. cdo mul input_file.nc percentage_ANTHRO_year.nc output.nc
            # 5 remove old files: all the percentage files
        else:
            raise ValueError("Type must be either 'openburing' or 'biomassburning'.")

    def sum_sectors_fires(self, input_files: list, sectors: list, output_path: Path, overwrite: bool = True):
        """
        Args:
        input_files (list): List of paths. The paths contain different
            sectors and the information is summed up and stored in the
            desired output_file.
        sectors (list): List of strings. The different files that are stored
            in input_files. The variable names of input_files.
        output_path (Path): File where the resulting emissions are stored.
            The output_path is adapted in case different fire cases are created.
        overwrite (bool): Default is True, i.e. data will be overwritten if it already exists

        """
        # check where (and if) biomassburning or openburning are contained
        fire_idx = None
        for i, var in enumerate(sectors):
            if ("openburning" in var) or ("biomassburning" in var):
                fire_idx = i
        if fire_idx is None:
            raise ValueError(
                "Fire cases can only be created if openburning or biomassburning files are given. (Must be contained in the file name as well)."
            )

        # Get ghg and fire type
        ghg = sectors[0].split("_")[0]
        # fire_type = sectors[fire_idx].split("_")[-1]

        # CASE 1: ALL FIRES
        all_fire_output_path = Path(str(output_path).replace(f"_{ghg}_", f"_{ghg}_all-fires_"))
        all_fire_files = input_files
        all_fire_sectors = sectors
        self._sum_sectors(all_fire_files, all_fire_sectors, all_fire_output_path, overwrite=overwrite)

        # CASE 2: ANTHRO FIRES
        anthro_fire_output_path = Path(str(output_path).replace(f"_{ghg}_", f"_{ghg}_anthro-fires_"))
        anthro_fire_sectors = sectors
        # get the right anthro fire file from meta
        # 1. get the current fire file
        all_fire_file = input_files[fire_idx]
        anthro_file_name = str(all_fire_file.name).replace("_em_", "_anthro_")
        potential_anthro_files = [f for f in self.meta_data.rglob(anthro_file_name)]
        if len(potential_anthro_files) < 1:
            print(f"No anthro fire file found for {all_fire_file.name}.")
            pass
        else:
            if len(potential_anthro_files) > 1:
                print(f"Several fire files found for {all_fire_file.name}. Choosing only the first one.")
            anthro_fire_files = input_files.copy()
            anthro_fire_files[fire_idx] = potential_anthro_files[0]
            self._sum_sectors(anthro_fire_files, anthro_fire_sectors, anthro_fire_output_path, overwrite=overwrite)

        # CASE 3: NO FIRES
        # drop the openburning or biomassburning sectors and files
        no_fire_output_path = Path(str(output_path).replace(f"_{ghg}_", f"_{ghg}_no-fires_"))
        no_fire_files = input_files.copy()
        no_fire_sectors = sectors.copy()
        no_fire_files.pop(fire_idx)
        no_fire_sectors.pop(fire_idx)
        self._sum_sectors(no_fire_files, no_fire_sectors, no_fire_output_path, overwrite=overwrite)

    def _sum_sectors(
        self, input_files: list, sectors: list, output_path: Path, overwrite: bool = True, threads: str = 1
    ):
        """
        Args:
            input_files (list): List of paths. The paths contain different
                sectors and the information is summed up and stored in the
                desired output_file.
            sectors (list): List of strings. The different files that are stored
                in input_files. The variable names of input_files.
            output_path (Path): File where the resulting emissions are stored.
                The output_path is adapted in case different fire cases are created.
            overwrite (bool): True, i.e. if the file in output_path already
                exists, it is overwritten. Set to False if you don't wanna overwrite
                these files.
            threads (str): Default 1.

        """
        # exit this if the file already exists and we dont wanna overwrite them
        if output_path.is_file() and (not overwrite):
            return

        ghg = sectors[0].split("_")[0]
        # math sum expression
        sum_expr = f"{ghg}={sectors[0]}"
        for sec in sectors[1:]:
            sum_expr = f"{sum_expr}+{sec}"

        # cdo expr,'var_new=var1+var2;' -merge infile1 infile2 outfile # works only on cmd
        # cdo 'expr,var_new=var1+var2;' -merge infile1 infile2 outfile # works for both
        commands = [
            "cdo",
            "-w",
            "-s",
            "-L",
            "-P",
            str(threads),
            f"expr,{sum_expr};",
            "-merge",
            *input_files,
            output_path,
        ]
        subprocess.call(commands)
