import gc
import json
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Union

import numpy as np
import xarray as xr

from climateset import CONFIGS
from climateset.utils import create_logger, get_json_config

LOGGER = create_logger(__name__)

VERIFY_FILE = "verify_file_ok"
VERIFY_UNITS = "verify_units"
VERIFY_RESOLUTION = "verify_resolution"
VERIFY_LONLAT = "verify_lonlat"
VERIFY_VARIABLES = "verify_variables"

AVAILABLE_RAW_CHECKS = frozenset([VERIFY_UNITS, VERIFY_RESOLUTION, VERIFY_LONLAT, VERIFY_VARIABLES])


class AbstractDirectoryChecker(ABC):
    def __init__(self, directory: Union[str, Path], checker_steps: list = None):
        if isinstance(directory, str):
            directory = Path(directory)
        self.directory: Union[str, Path] = directory
        self.checker_steps: list = checker_steps
        if not self.checker_steps:
            self.checker_steps = []

        self.results = []
        self.errors = []

    @abstractmethod
    def check_directory(self):
        raise NotImplementedError

    def write_results(self, base_path: Path):
        if self.results:
            filename = base_path / "results.json"
            with open(filename, "w", encoding="utf-8") as json_file:
                json.dump(self.results, json_file, indent=2)
        if self.errors:
            filename = base_path / "errors.json"
            with open(filename, "w", encoding="utf-8") as json_file:
                json.dump(self.errors, json_file, indent=2)


class AbstractFileChecker(ABC):
    def __init__(self, input_file: Union[str, Path], checker_steps: list = None, available_checks: dict = None):
        self.input_file: Union[str, Path] = input_file
        self.checker_steps: list = checker_steps
        if not self.checker_steps:
            self.checker_steps = []
        self.available_steps = available_checks
        self.results = {}
        self.dataset: Optional[xr.Dataset] = None

    def check(self):
        self.check_file_ok()
        for step in self.checker_steps:
            self.available_steps[step]()
        return {"filename": self.input_file.name, "results": self.results}

    def check_file_ok(self):
        try:
            LOGGER.info(f"Trying to load file {self.input_file}")
            ds = xr.load_dataset(self.input_file)
            self.dataset = ds
            LOGGER.info(f"Successfully loaded file {self.input_file}")
            self.results["file_ok"] = True

        except ValueError:
            self.results["file_ok"] = False
            LOGGER.warning(f"The following file is corrupt:\n {self.input_file}")
            self.dataset = None

    def cleanup_memory(self):
        del self.dataset
        gc.collect()


class BasicDirectoryChecker(AbstractDirectoryChecker):
    def __init__(
        self,
        directory: Union[str, Path] = None,
        checker_steps: list = None,
    ):
        super().__init__(directory=directory, checker_steps=checker_steps)

    def check_directory(self) -> bool:
        """
        Checking all files in a sub dir.

        TAKES TIME.
        Args:

        Returns:
            bool: True if all checks were successful, False if not
        """
        for f in self.directory.glob("*.nc"):
            file_checker = BasicFileChecker(input_file=f, checker_steps=self.checker_steps)
            file_checker.check()
            if not all(file_checker.results.values()):
                self.errors.append(file_checker.results)
            self.results.append(file_checker.results)
            file_checker.cleanup_memory()

        if self.errors:
            LOGGER.warning("There were errors running while running the checks")
            LOGGER.warning(f"Errors : {self.errors}")
            return False

        LOGGER.info("All checks passed")
        return True


class BasicFileChecker(AbstractFileChecker):
    def __init__(
        self,
        input_file,
        checker_steps: list = None,
    ):
        super().__init__(
            input_file=input_file, checker_steps=checker_steps, available_checks={VERIFY_FILE: self.check_file_ok}
        )

    def check(self):
        """
        Logs warnings if there are inconsistent files.

        Checks always
        for corruptness, the rest is optional. See initialization of class.
        Args:
            input_file (Path): The path that should be checked
            log_file (Path): Where error / warnings should be logged additionally.
            mode (str): 'w' writes the logs into the file, 'a' appends it to
                an existing file.
        Returns:
            True if all checks were passed, False if not.
        """

        self.check_file_ok()


class RawDirectoryChecker(AbstractDirectoryChecker):
    def __init__(
        self,
        directory: Union[str, Path] = None,
        checker_steps: list = None,
        type_class=None,
        lonlat_resolution_to_verify: str = "50_km",
        processing_params_config: Union[str, Path] = CONFIGS / "processing" / "raw_processing_params.json",
    ):
        super().__init__(directory=directory, checker_steps=checker_steps)
        self.type_class = type_class
        self.lonlat_resolution_to_verify = lonlat_resolution_to_verify
        self.meta_raw_dict = get_json_config(processing_params_config)
        self.calendar = self.meta_raw_dict["calendar"]
        self.desired_units = self.meta_raw_dict["units"]

        lons = self.meta_raw_dict[self.lonlat_resolution_to_verify][f"{self.type_class}_lon"]
        lats = self.meta_raw_dict[self.lonlat_resolution_to_verify][f"{self.type_class}_lat"]
        self.lon_specs = (float(lons["min"]), float(lons["max"]), float(lons["step"]))  # (-179.75, 179.75, 0.5)
        self.lat_specs = (float(lats["min"]), float(lats["max"]), float(lats["step"]))  # (-89.75, 89.75, 0.5)

    def add_unit_check(self):
        self.checker_steps.append(VERIFY_UNITS)

    def add_resolution_check(self):
        self.checker_steps.append(VERIFY_RESOLUTION)

    def add_lonlat_check(self):
        self.checker_steps.append(VERIFY_LONLAT)

    def add_variables_check(self):
        self.checker_steps.append(VERIFY_VARIABLES)

    def check_directory(self) -> bool:
        """
        Checking all files in a sub dir.

        TAKES TIME.
        Returns:
            bool: True if all checks were successful, False if not
        """

        for f in self.directory.rglob("*.nc"):
            file_checker = RawFileChecker(
                input_file=f,
                type_class=self.type_class,
                checker_steps=self.checker_steps,
                lonlat_resolution_to_verify=None,
                calendar=self.calendar,
                units=self.desired_units,
                lon_specs=self.lon_specs,
                lat_specs=self.lat_specs,
                expected_units=self.meta_raw_dict["units"],
            )
            file_checker.check()
            if not all(file_checker.results.values()):
                self.errors.append(file_checker.results)
            self.results.append(file_checker.results)
            file_checker.cleanup_memory()

        if self.errors:
            LOGGER.warning("There were errors running while running the checks")
            LOGGER.warning(f"Errors : {self.errors}")
            return False

        LOGGER.info("All checks passed")
        return True


class RawFileChecker(AbstractFileChecker):
    def __init__(
        self,
        input_file,
        type_class=None,
        checker_steps: list = None,
        lonlat_resolution_to_verify=None,
        calendar=None,
        units=None,
        lon_specs=None,
        lat_specs=None,
        expected_units=None,
    ):
        super().__init__(
            input_file=input_file,
            checker_steps=checker_steps,
            available_checks={
                VERIFY_UNITS: self.check_units_ok,
                VERIFY_RESOLUTION: self.check_resolution,
                VERIFY_LONLAT: self.check_lonlat,
                VERIFY_VARIABLES: self.check_variables,
            },
        )
        self.type_class = type_class
        self.lonlat_resolution_to_verify = lonlat_resolution_to_verify
        self.calendar = calendar
        self.desired_units = units
        self.lon_specs = lon_specs
        self.lat_specs = lat_specs
        self.results = {}
        self.expected_units = expected_units

    def check(self):
        """
        Logs warnings if there are inconsistent files.

        Checks always
        for corruptness, the rest is optional. See initialization of class.
        Args:
            input_file (Path): The path that should be checked
            log_file (Path): Where error / warnings should be logged additionally.
            mode (str): 'w' writes the logs into the file, 'a' appends it to
                an existing file.
        Returns:
            True if all checks were passed, False if not.
        """

        self.check_file_ok()
        for step in self.checker_steps:
            self.available_steps[step]()

    def check_variables(self):
        try:
            actual_var = self.dataset.variable_id
            # cmip6 case
            if self.type_class == "cmip6":
                named_var = str(self.input_file.name).split("_")[4]
            # input4mips case
            else:
                named_var = str(self.input_file.parents[0]).split("/")[-4]
                if "sum" in named_var:
                    named_var = named_var.split("_")[0]

            if actual_var != named_var:
                self.results["var_ok"] = False
                LOGGER.warning(
                    f"The following file contains the variable {actual_var}, but the variable {named_var} "
                    f"was expected: \n{self.input_file}"
                )
            else:
                self.results["var_ok"] = True
        except AttributeError:
            self.results["var_ok"] = False
            LOGGER.warning(f"Dataset does not contain variable_id: \n {self.input_file}")

    def check_lonlat(self, ds, check_dict, input_file, passed_checks):
        lon_min, lon_max, lon_step = self.lon_specs  # (-179.75, 179.75, 0.5)
        lat_min, lat_max, lat_step = self.lat_specs  # (-89.75, 89.75, 0.5)
        try:
            found_lon_min, found_lon_max = ds.lon[0].item(), ds.lon[-1].item()
            found_lat_min, found_lat_max = ds.lat[0].item(), ds.lat[-1].item()
            found_lon_step = abs(ds.lon[0].item() - ds.lon[1].item())
            found_lat_step = abs(ds.lat[0].item() - ds.lat[1].item())

            # longitude checks
            if (lon_min == found_lon_min) and (lon_max == found_lon_max) and (np.isclose(lon_step, found_lon_step)):
                check_dict["lon_ok"] = True
            else:
                passed_checks = False
                check_dict["lon_ok"] = False
                LOGGER.warning(f"The following file has unexpected lon values: \n{input_file}")

            # latitude checks
            if (lat_min == found_lat_min) and (lat_max == found_lat_max) and (np.isclose(lat_step, found_lat_step)):
                check_dict["lat_ok"] = True
            else:
                passed_checks = False
                check_dict["lat_ok"] = False
                LOGGER.warning(f"The following file has unexpected lat values: \n{input_file}")

        except AttributeError:
            passed_checks = False
            check_dict["lon_ok"] = False
            check_dict["lat_ok"] = False
            LOGGER.warning(f"Dataset contains different namings for longitude / latitude. \n {input_file}")
        return passed_checks

    def check_resolution(self):
        try:
            actual_spat_res = self.dataset.attrs["nominal_resolution"]
            actual_spat_res = actual_spat_res.replace("_", " ")  # just in case
            named_spat_res = " ".join(str(self.input_file.name).split("_")[-5:-3])
            if actual_spat_res != named_spat_res:
                self.results["res_ok"] = False
                LOGGER.warning(
                    f"The following file has the nominal resolution {actual_spat_res}, but the resolution "
                    f"{named_spat_res} was expected: \n{self.input_file}"
                )
            else:
                self.results["res_ok"] = True

            actual_temp_res = self.dataset.attrs["frequency"]
            named_temp_res = str(self.input_file.name).split("_")[-3]
            if actual_temp_res != named_temp_res:
                self.dataset["freq_ok"] = False
                LOGGER.warning(
                    f"The following file has the temporal frequency {actual_temp_res}, but the frequency "
                    f"{named_temp_res} was expected: \n{self.input_file}"
                )
            else:
                self.dataset["freq_ok"] = True
        except AttributeError:
            self.dataset["res_ok"] = False
            self.dataset["freq_ok"] = False
            LOGGER.warning(f"The following file gives no access to nominal resolution or frequency: {self.input_file}")

    def check_units_ok(self):
        try:
            var = self.dataset.variable_id
            found_unit = None

            try:
                found_unit = self.dataset[var].units
            except AttributeError:
                try:
                    found_unit = self.dataset["variable_units"]
                except AttributeError:
                    raise KeyError(f"Unit could not be found in file {self.input_file}")

            expected_unit = self.expected_units[var]
            if found_unit != expected_unit:
                self.results["units_ok"] = False
                LOGGER.warning(
                    f"The following file has the unit {found_unit}, but the unit {expected_unit} "
                    f"was expected: \n{self.input_file}"
                )
            else:
                self.results["units_ok"] = True

        except AttributeError:
            self.results["units_ok"] = False
            LOGGER.warning(f"The following file does not allow us to call ds.variable_id:\n{self.input_file}")
