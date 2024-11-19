from abc import ABC, abstractmethod


class AbstractProcessorStep(ABC):
    def __init__(self):
        self.results_directory = None

    @abstractmethod
    def execute(self, input_directory):
        pass

    def get_results_directory(self):
        return self.results_directory
