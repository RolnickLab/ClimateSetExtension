from abc import ABC, abstractmethod
from pathlib import Path


class AbstractProcessorStep(ABC):
    def __init__(self):
        self.results_directory = None

    @abstractmethod
    def execute(self, input_directory):
        pass

    def get_results_directory(self):
        return self.results_directory


def process_steps(input_directory: Path, list_of_steps: list[AbstractProcessorStep]):
    current_input_directory = input_directory
    for step in list_of_steps:
        step.execute(current_input_directory)
        # current_input_directory = step.get_results_directory()

    return current_input_directory
