from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union

from climateset import CONFIGS
from climateset.processing.raw.checker import (
    AbstractDirectoryChecker,
    BasicDirectoryChecker,
)
from climateset.utils import create_logger, get_json_config

LOGGER = create_logger(__name__)


class AbstractRawProcesser(ABC):
    """Abstract class for raw processing."""

    def __init__(
        self,
        input_directory: Union[str, Path],
        working_directory: Union[str, Path],
        processing_steps: list = None,
        checker: AbstractDirectoryChecker = None,
        processing_parameters_config: Union[str, Path] = CONFIGS / "processing" / "raw_processing_params.json",
    ):
        """Init
        Args:

        """
        # type of class
        self.type_class = self.type_class_meta()
        self.input_directory = input_directory
        if isinstance(input_directory, str):
            self.input_directory = Path(input_directory)
        self.working_directory = working_directory
        if isinstance(working_directory, str):
            self.working_directory = Path(working_directory)
        self.processing_steps = processing_steps
        if not self.processing_steps:
            self.processing_steps = []
        self.available_steps = {}
        self.processing_params_config: Union[str, Path] = processing_parameters_config
        self.meta_raw_dict = get_json_config(self.processing_params_config)
        self.calendar = self.meta_raw_dict["calendar"]
        self.desired_units = self.meta_raw_dict["units"]

        self.checker = checker
        if not self.checker:
            self.checker: AbstractDirectoryChecker = BasicDirectoryChecker(directory=self.input_directory)

    @abstractmethod
    def type_class_meta(self) -> str:
        """Returns the name tag of the subclass."""

    @abstractmethod
    def preprocess_subdir(
        self,
        input_dir: Path,
        cleaned_dir: Path,
        processed_dir: Path,
        load_dir: Path,
        overwrite: bool,
        silent: bool,
        sum_sec_input_res: str,
        sum_sec_input_freq: str,
    ):
        """
        Preprocessing a subdir - must be implemented for each subclass.
        Args:
            input_dir:
            cleaned_dir:
            processed_dir:
            load_dir:
            overwrite:
            silent:
            sum_sec_input_res:
            sum_sec_input_freq:

        Returns:

        """

    @abstractmethod
    def file_belongs_to_type(self, input_file: Path) -> bool:
        """
        Checks if a file belongs to input4mips or cmip6 category.

        Args:
            input_file (Path): the file that should be checked
        Returns:
            True if it belong to the class type, False if not
        """

    def add_processing_step(self, steps: Union[str, list]):
        if isinstance(steps, str):
            steps = [steps]
        for step in steps:
            if step in self.available_steps:
                self.processing_steps.append(step)

    def list_available_steps(self):
        step_list = [step for step in self.available_steps.keys() if not step.startswith("_")]
        LOGGER.info(f"Available steps: {step_list}")
        return step_list
