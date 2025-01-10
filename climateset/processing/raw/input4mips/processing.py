import os
import re
import subprocess
from copy import copy
from pathlib import Path

import pint
import xarray as xr
from tqdm import tqdm

from climateset.processing.raw.utils import create_generic_output_path
from climateset.utils import create_logger

LOGGER = create_logger(__name__)

CDO_SET_CALENDAR = "-setcalendar"
CDO_SET_DAY = "-setday,1"
CDO_SET_TIME = "-settime,00:00:00"
CDO_VERTICAL_SUM = "-vertsum"


def rename_biomassburning_file(input_file: Path):
    """
    Rename file for biomassburning in case it is in the naked version (from the downloader). Originally, the files are
    just called 'GHG' instead of GHG_em_biomassburning. Overwrites the files directly. Attention,

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


def rename_biomassburning_files_in_directory(input_dir):
    LOGGER.info(f"Starting to rename biomassburning files in {input_dir}.")
    total_files = len(list(input_dir.rglob("*.nc")))
    for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
        if len(files) > 0:
            for file in files:
                rename_biomassburning_file(Path(path) / file)
    # remove empty subdirs
    for p in input_dir.glob("**/*"):
        if p.is_dir() and len(list(p.iterdir())) == 0:
            os.removedirs(p)
    return input_dir


def reorder_ssp_co2_file_dimensions(input_file: Path):
    """
    Fix single co2 file.

    CDO can be applied again afterwards.
    Memory issues might arise (core dumped). Not sure why. Runs on compute canada though.

    This function is essentially a wrapper around the following CDO command:

    cdo ncpdq --ovr -a time,level,lon,lat input.nc output.nc


    Args:
        input_file (Path): The co2 file that should be fixed.
    """
    # ncpdq --ovr -a time,level,lon,lat input.nc output.nc
    commands = ["ncpdq", "--ovr", "-a", "time,level,lon,lat", input_file, input_file]
    subprocess.call(commands)


def reorder_all_ssp_co2_files_dimensions(input_dir: Path):
    LOGGER.info(f"Starting to fix co2 files in {input_dir}.")
    total_files = len(list(input_dir.rglob("*.nc")))
    for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
        if len(files) > 0:
            for file in files:
                if ("ssp" in file) and ("CO2" in file):
                    file_path = Path(path) / file
                    reorder_ssp_co2_file_dimensions(file_path)

    LOGGER.info("...Finished fixing CO2 files.")
    return input_dir


def create_anthro_fire_file(
    input_file: Path, metadata_directory: Path, filetype: str, cdo_operators, output_file: Path, threads: int = 1
):
    """
    Create a new openburning or biomassburning file from a given file.

    Returns the path where the new file lifes.
    Args:
        input_file (Path):
        type (str): either 'biomassburning' or 'openburning'
        output_file (Path): Where to store the anthro fire file
        threads (int): threads for cdo
    """
    if metadata_directory is None:
        LOGGER.info("You need to provide meta data for this function. Use set_meta_data func for that.")
    file_parts = str(input_file.stem).split("_")
    scenario, ghg, year = file_parts[1], file_parts[2], file_parts[-1]
    meta_dir = metadata_directory / filetype

    # prepare additional operators (calendar and time-axis) for later compatability
    chained_operators = copy(cdo_operators)

    if filetype == "openburning":
        # find a file that contains openburning, scenario and ghg in the meta_data / "future-openburning"
        meta_dir = metadata_directory / "future-openburning"
        matched_files = [
            f for f in meta_dir.rglob("*.nc") if (scenario in str(f)) and (filetype in str(f)) and (ghg in str(f))
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
            CDO_VERTICAL_SUM,
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

    elif filetype == "biomassburning":
        meta_dir = metadata_directory / "historic-biomassburning"
        matched_files = [
            f for f in meta_dir.rglob("*.nc") if (filetype in str(f)) and (ghg in str(f)) and (year in str(f))
        ]

        # find the agri / defo files (written this way for later extension)
        perc_dict = {
            "AGRI": None,
            "DEFO": None,
        }
        for f in matched_files:
            for key in perc_dict:
                if key in str(f):
                    perc_dict[key] = f
        if (perc_dict["AGRI"] is None) or (perc_dict["DEFO"] is None):
            LOGGER.warning("The AGRI or DEFO file does not exist, see dictionary: \n", stacklevel=2)
            LOGGER.info(perc_dict)
            LOGGER.info("... skipping this file.")
            return

        # cdo 'expr,percentage_AGRI=percentage_AGRI/100' add percentage_AGRI_year.nc
        #     percentage_DEFO_year.nc percentage_ANTHRO_year.nc
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


# TODO fix output directory
def create_anthro_fire_directory(
    input_directory: Path, meta_data: Path, cdo_operators, overwrite: bool = False
) -> Path:
    """
    Creating anthropogenic fire files for a subdir.

    Args:
        input_directory (Path): Input dir
    """
    # loop through input dir, get biomassburning and openburning files
    LOGGER.info(f"Starting to create anthropogenic fire files for files in {input_directory}.")
    fire_output_dir = meta_data / "anthro-fire-data"
    total_files = len(list(input_directory.rglob("*.nc")))
    for path, subdirs, files in tqdm(os.walk(input_directory), total=total_files):
        if len(files) > 0:
            for file in files:
                input_file = Path(path) / file
                # must match the right nominal and temporal resolution
                if ("biomassburning" in file) and ("25_km" in file) and ("mon" in file):
                    filetype = "biomassburning"
                elif ("openburning" in file) and ("50_km" in file) and ("mon" in file):
                    filetype = "openburning"
                else:
                    continue
                scenario = file.split("_")[1]
                ghg = file.split("_")[2]
                fire_name = file.replace("_em_", "_anthro_")
                output_file = Path(fire_output_dir / filetype / scenario / ghg / fire_name)
                output_file.parent.mkdir(parents=True, exist_ok=True)
                if (not output_file.is_file()) or overwrite:
                    create_anthro_fire_file(
                        input_file=input_file,
                        metadata_directory=meta_data,
                        filetype=filetype,
                        cdo_operators=cdo_operators,
                        output_file=output_file,
                    )

    LOGGER.info(f"... Finished creating anthro fire files and stored them in {fire_output_dir}")
    return input_directory


def cdo_process_file(input_file: Path, output_file: Path, threads: int = 1, cdo_operators=None):
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

    chained_operators = copy(cdo_operators)
    # cdo vertsum infile outfile [CHECKED]
    # AIR: sum over height levels
    # anthro: sum over different sectors

    # ATTENTION only relevant for "monthly" data!
    # cdo -setday,1 -settime,00:00:00
    # TODO add check for mon data only
    if "mon" not in "":
        chained_operators.remove(CDO_SET_DAY)
        chained_operators.remove(CDO_SET_TIME)
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


# TODO fix output directory
def cdo_process_directory(input_directory, cdo_operators, output_directory=None, overwrite=False):
    """
    Apply desired emission processing steps to a subdir and store it in output_dir.

    Can only be applied on raw input4mpis data, not on those
    that have already been processed (i.e. merged / summed up).
    Args:
        input_directory (Path): Where the preprocessing should be applied
        output_dir (Path): Where the results should be stored
        overwrite (bool): If the data should be overwritten if it already exists
    """
    # TODO if a year / file is missing -> write a warning!
    # TODO LATER add feature to check lon / lat consistency for input4mips
    # 1. apply the "preprocessing" - all the functions that can
    # be applied to a single file
    LOGGER.info(f"Start preprocessing of emission files {input_directory}....")
    total_files = len(list(input_directory.rglob("*.nc")))
    for path, subdirs, files in tqdm(os.walk(input_directory), total=total_files):
        if len(files) > 0:
            for file in files:
                # create output dir
                input_file = Path(path) / file
                # TODO fix filetype check
                # if self.file_belongs_to_type(input_file):
                output_file = create_generic_output_path(output_directory, path, file)
                # skip if file already exists and we dont wanna overwrite it
                if (not output_file.is_file()) or overwrite:
                    LOGGER.info(f"\nProcessing the following file: {input_file}")
                    cdo_process_file(input_file, output_file, cdo_operators=cdo_operators)

    LOGGER.info(f"...Finished preprocessing \nof {input_directory} and saved it at {output_directory}.")
    return output_directory


def merge_ghg_files(input_files: list, output_file: Path, threads: int = 1):
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


def create_totals(input_files: list, output_path: Path, spat: str, temp: str, threads: int = 1):
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
        commands.append(f"expr,sum=fldsum({ghg});min=fldmin({ghg});max=fldmax({ghg});mean=fldmean({ghg});")

    if temp == "year":
        commands.append("-yearsum")

    # add relevant commands for all cases
    commands.extend(["-mergetime", *input_files, output_path])

    subprocess.call(commands)


def merge_ghg_directory(
    input_dir: Path, output_dir: Path, ghg_list: list, ghg_sum_name: str = "ALL", overwrite: bool = False
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
    LOGGER.info(f"Start creating merged GHG emission files for {ghg_list} in {input_dir}...")
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
                        if ghg_file_str not in file_dict:
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
            merge_ghg_files(
                input_files,
                out_path,
            )

    LOGGER.info(f"... Finished merging GHG {ghg_list} in {output_dir}.")
    return output_dir


def sum_sectors(input_files: list, sectors: list, output_path: Path, overwrite: bool = True, threads: str = 1):
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


def correct_units_in_file(input_file: Path, desired_units: dict):
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
            desired_unit = desired_units[var].replace("-", "^-")

            # exit if it's the right one
            if found_unit == desired_unit:
                return

            try:
                ureg(found_unit)
                ureg(desired_unit)
            except pint.errors.UndefinedUnitError:
                LOGGER.info(
                    "Unit is not defined for pint. Check here: "
                    "https://github.com/hgrecco/pint/blob/master/pint/default_en.txt"
                )
                return
            # update ds
            # this function works for any kind of unit transformations (if they can be transformed)
            ds.update({var: xr.apply_ufunc(lambda x: ureg.Quantity(x, found_unit).to(desired_unit).magnitude, ds[var])})
            # update attrs
            ds[var].attrs["units"] = desired_unit
            # overwrite old file
            ds.to_netcdf(input_file, mode="w", format="NETCDF4", engine="netcdf4")


def correct_units_all_files_in_directory(input_dir: Path, desired_units: dict):
    LOGGER.info(f"Starting to correct all units (if needed) in {input_dir}.")
    total_files = len(list(input_dir.rglob("*.nc")))
    for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
        if len(files) > 0:
            for file in files:
                correct_units_in_file(Path(path) / file, desired_units)

    LOGGER.info(f"... Finished correcting the units (if needed) in {input_dir}")
    return input_dir


def sum_sectors_fires(input_files: list, sectors: list, output_path: Path, meta_data: Path, overwrite: bool = True):
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
            "Fire cases can only be created if openburning or biomassburning files are given. "
            "(Must be contained in the file name as well)."
        )

    # Get ghg and fire type
    ghg = sectors[0].split("_")[0]
    # fire_type = sectors[fire_idx].split("_")[-1]

    # CASE 1: ALL FIRES
    all_fire_output_path = Path(str(output_path).replace(f"_{ghg}_", f"_{ghg}_all-fires_"))
    all_fire_files = input_files
    all_fire_sectors = sectors
    sum_sectors(all_fire_files, all_fire_sectors, all_fire_output_path, overwrite=overwrite)

    # CASE 2: ANTHRO FIRES
    anthro_fire_output_path = Path(str(output_path).replace(f"_{ghg}_", f"_{ghg}_anthro-fires_"))
    anthro_fire_sectors = sectors
    # get the right anthro fire file from meta
    # 1. get the current fire file
    all_fire_file = input_files[fire_idx]
    anthro_file_name = str(all_fire_file.name).replace("_em_", "_anthro_")
    potential_anthro_files = [f for f in meta_data.rglob(anthro_file_name)]
    if len(potential_anthro_files) < 1:
        LOGGER.info(f"No anthro fire file found for {all_fire_file.name}.")
        pass
    else:
        if len(potential_anthro_files) > 1:
            LOGGER.info(f"Several fire files found for {all_fire_file.name}. Choosing only the first one.")
        anthro_fire_files = input_files.copy()
        anthro_fire_files[fire_idx] = potential_anthro_files[0]
        sum_sectors(anthro_fire_files, anthro_fire_sectors, anthro_fire_output_path, overwrite=overwrite)

    # CASE 3: NO FIRES
    # drop the openburning or biomassburning sectors and files
    no_fire_output_path = Path(str(output_path).replace(f"_{ghg}_", f"_{ghg}_no-fires_"))
    no_fire_files = input_files.copy()
    no_fire_sectors = sectors.copy()
    no_fire_files.pop(fire_idx)
    no_fire_sectors.pop(fire_idx)
    sum_sectors(no_fire_files, no_fire_sectors, no_fire_output_path, overwrite=overwrite)


def sum_sector_cases(
    file_dict: dict,
    output_dir: Path,
    meta_data: Path,
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
                sum_sectors_fires(input_files, sector_list, out_path, meta_data=meta_data, overwrite=overwrite)
            else:
                sum_sectors(input_files, sector_list, out_path, overwrite=overwrite)

        elif (num_nones == 1) and (ghg in force_sum_sectors):
            # just the simple summing if desired + warning
            if not silent:
                LOGGER.warning(f"Merged files for {file}, however some sectors are missing", stacklevel=2)
            sum_sectors(input_files, sector_list, out_path, overwrite=overwrite)

        elif (num_nones == 1) and ghg not in force_sum_sectors:
            if not silent:
                LOGGER.warning(
                    f"Skipping case {file}. Not all three sectors are available. "
                    f"You can force merging with 'force_sum_sectors=[GHG1, GHG2]'."
                )

        elif num_nones == 2:
            # warning, skipping
            if not silent:
                LOGGER.warning(
                    f"Skipping case {file}. At least two files must exist to be summed together.",
                    stacklevel=2,
                )
        else:
            raise RuntimeError("'num_nones' must be between 0 and 3.")


# TODO fix output dir and overwrite
# TODO break this into two functions
def sum_sectors_directory(
    input_dir: Path,
    output_dir: Path,
    meta_data: Path,
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
    LOGGER.info(f"Start summing up different types of emission files in {input_dir}...")
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
    total_files = len(list(input_dir.rglob("*.nc")))
    for path, subdirs, files in tqdm(os.walk(input_dir), total=total_files):
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

    sum_sector_cases(
        file_dict=ssp_file_dict,
        output_dir=output_dir,
        fire_cases=fire_cases,
        meta_data=meta_data,
        force_sum_sectors=force_sum_sectors,
        overwrite=False,
        silent=silent,
    )

    sum_sector_cases(
        file_dict=historical_file_dict,
        output_dir=output_dir,
        fire_cases=fire_cases,
        meta_data=meta_data,
        force_sum_sectors=force_sum_sectors,
        overwrite=False,
        silent=silent,
    )

    LOGGER.info(f"...Finished summing up ghg emissions and stored results in {output_dir}")
    return output_dir


def create_emissions_totals(
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
    LOGGER.info(f"Start creating loadable emission data {input_dir}....")
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
            create_totals(
                input_files,
                out_path,
                spat,
                temp,
            )

    LOGGER.info(f"... Finished creating loadable emission data {output_dir}....")
    return output_dir
