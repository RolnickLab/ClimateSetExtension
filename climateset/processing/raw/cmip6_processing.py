import os
import re
import subprocess
import warnings
from pathlib import Path

import xarray as xr
import xmip
import xmip.preprocessing as xmip_preprocessing
from tqdm import tqdm

from climateset.processing.raw.abstract_raw_processing import AbstractRawProcesser
from climateset.processing.raw.utils import create_generic_output_path


# Attention: we need to write _desired_units for our all our data
# Attention: Before you apply, make sure this is only applied to CMIP6 data!
class ClimateModelProcesser(AbstractRawProcesser):
    """
    Can be called to apply xmip preprocessing.

    Works only on Cmip6 Data!
    """

    def __init__(
        self,
        rename_ds: bool = True,
        correct_lonlat: bool = True,
        correct_units: bool = True,
        fix_bounds: bool = True,
        verify_units: bool = False,
        verify_resolution: bool = False,
        verify_lonlat: bool = False,
        verify_num_of_files: bool = False,
        verify_variables: bool = False,
        correct_time_axis: bool = True,
        correct_calendar: bool = True,
        sum_levels: bool = True,
        **kwargs,
    ):
        """Init function for xmip processer
        Args:
            rename_ds (bool): Make naming in cmip6 consistent
            correct_lonlat (bool): Make lonlat naming consistent
            correct_units (bool): Make units consistent
            fix_bounds (bool): Fix bounds and vertex
            verify_units (bool): Makes sure that all ghg conform to the same units
            verify_resolution (bool): Checks if the res in the file name is the same like in the file
            verify_lonlat (bool): Checks if the longitude, latitude values are as expected
            verify_num_of_files (bool): Checks how many files are in the final dir.
                raises a warning if more or less than one.
            verify_variables (bool): Checks if the variables are also the expected vars.
            correct_calendar (bool): Makes sure that the right calendar is used.
            correct_time_axis (bool): Makes sure that each file starts at the first
                of each time-unit and shifts if not.
            sum_levels (bool): Makes sure that all emissions are sumed over the
                different levels.
        """
        # init abstract class with checks
        super().__init__(
            verify_units=verify_units,
            verify_resolution=verify_resolution,
            verify_lonlat=verify_lonlat,
            verify_num_of_files=verify_num_of_files,
            verify_variables=verify_variables,
        )

        # xmip processing steps
        self.rename_ds = rename_ds
        self.correct_lonlat = correct_lonlat
        self.correct_units = correct_units
        self.fix_bounds = fix_bounds

        # climatevers processing steps
        self.correct_calendar = correct_calendar
        self.correct_time_axis = correct_time_axis
        self.sum_levels = sum_levels

        # update the desired unit dicts for unit correction of xmip
        for var, unit in self.desired_units.items():
            self.update_desired_units_dict(var, unit)

    def type_class_meta(self) -> str:
        """Returns the name tag of the subclass."""
        return "cmip6"

    # TODO share this in process class with input4mips / cmip6
    def file_belongs_to_type(self, input_file: Path) -> bool:
        """
        Check if the name cmip6 appears in the path name.

        Args:
            input_file (Path): the file that should be checked for cmip6
        Returns:
            True if name indicates it's a cmip6 file, False if not.
        """
        return bool(re.search("cmip6", str(input_file), re.IGNORECASE))

    # Is used internally and can be updated externally
    def update_desired_units_dict(self, var: str, unit: str):
        """
        Must be applied before running xmip processing.

        Updates the global
        dict of xmip to include which units you want to use throughout all models.
        Args:
            var (str): the variable that shall be added as key
            unit (str): the desired unit that all models should use as value
        """
        # adapt global variable from xmip
        xmip.preprocessing._desired_units[var] = unit

    # share this with input4mips?
    def preprocess_subdir(
        self,
        sub_dir: Path,
        output_dir: Path,
        xmip: bool = True,
        climateverse: bool = True,
        model: str = "",
        overwrite: bool = False,
        threads: int = 1,
    ):
        """
        Applies the preprocessing to a complete directory.

        Args:
            sub_dir (Path): directory on which xmip processing should be applied
            output_dir (Path): where the preprocessed files should be stored
            xmip (bool): if xmip processing should be applied
            climateverse (bool): if additional processing from us should be applied
            model (str): Default is "", i.e. all models are processed the
                same, all at once. Set this to a specific model if you want instead.
            overwrite (bool): If data can be overwritten or not
        """
        print(f"Start Cmip6 preprocessing of {sub_dir}.")

        # loop through sub_dir to remap all files in here
        # total_files = len(list(sub_dir.rglob("*.nc")))
        for path, subdirs, files in tqdm(os.walk(sub_dir)):
            if len(files) > 0:
                for file in files:
                    # create output dir
                    input_file = Path(path) / file

                    if self.file_belongs_to_type(input_file) and (model in str(input_file)):
                        # create output file
                        if sub_dir == output_dir:
                            output_file = input_file
                        else:
                            output_file = create_generic_output_path(output_dir, path, file)

                        # process
                        if (not output_file.is_file()) or overwrite:
                            if climateverse and xmip:
                                self.climateverse_process(input_file, output_file, threads=threads)
                                self.xmip_process(output_file, output_file)
                            elif climateverse:
                                self.climateverse_process(input_file, output_file, threads=threads)
                            elif xmip:
                                self.xmip_process(input_file, output_file)
                            else:
                                warnings.warn("No processing requested.")

        print(f"Finished Cmip6 preprocessing of {sub_dir} and saved it at {output_dir}.")

    def xmip_process(self, input_file: Path, output_file: Path):
        """
        Applies the xmip processing to a specific file.

        Make sure this is
        a CMIP6 file beforehand.
        Args:
            input_file (Path): Input netcdf cmip6 file that shall be processed.
            output_file (Path): Where the processed dataset shall be stored.
        """
        # load file as dataset
        ds = xr.load_dataset(input_file)

        # call preprocessing
        if self.rename_ds:
            ds = xmip_preprocessing.rename_cmip6(ds)
        if self.correct_lonlat:
            ds = xmip_preprocessing.promote_empty_dims(ds)
            ds = xmip_preprocessing.correct_coordinates(ds)
            ds = xmip_preprocessing.broadcast_lonlat(ds)
            ds = xmip_preprocessing.correct_lon(ds)
        if self.correct_units:
            ds = xmip_preprocessing.correct_units(ds)
        if self.fix_bounds:
            ds = xmip_preprocessing.parse_lon_lat_bounds(ds)
            ds = xmip_preprocessing.sort_vertex_order(ds)
            ds = xmip_preprocessing.maybe_convert_bounds_to_vertex(ds)
            ds = xmip_preprocessing.maybe_convert_vertex_to_bounds(ds)
            ds = xmip_preprocessing.fix_metadata(ds)
            ds = ds.drop_vars(xmip_preprocessing._drop_coords, errors="ignore")

        # store at new location (not raw anymore)
        ds.to_netcdf(output_file, mode="w", format="NETCDF4", engine="netcdf4")

    def climateverse_process(self, input_file: Path, output_file: Path, threads: int = 1):
        """
        Applies additional processing steps that are not included in xmip. Needs to be applied BEFORE xmip processing
        because it's using cdo and xmip breaks cdo formatting (at least for me so far).

        Args:
            input_file (Path): The file that should be processed.
            output_file (Path): Where the resulting file is stored.
            threads (int): threads for cdo
        """
        chained_operators = []

        # cdo vertsum infile outfile
        if self.sum_levels:
            chained_operators.append("-vertsum")

        # ATTENTION only relevant for "monthly" data!
        # cdo -setday,1 -settime,00:00:00
        if self.correct_time_axis:
            if "mon" in str(input_file):
                chained_operators.append("-setday,1")
                chained_operators.append("-settime,00:00:00")
            else:
                warnings.warn("Skipping correcting time axis, since it should only be applied to monthly data.")

        if self.correct_calendar:
            # default calendar: 365_day
            chained_operators.append(f"-setcalendar,{self.calendar}")

        # run commands
        commands = ["cdo", "-s", "-w", "-L", "-P", str(threads), *chained_operators, input_file, output_file]
        subprocess.call(commands)
