from abc import ABC, abstractmethod


class AbstractDownloader(ABC):
    @abstractmethod
    def download(self):
        pass
