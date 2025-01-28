from abc import ABC, abstractmethod

from climateset.download.abstract_downloader_config import AbstractDownloaderConfig


class AbstractDownloader(ABC):
    def __init__(self, config: AbstractDownloaderConfig):
        self.config = config

    @abstractmethod
    def download(self):
        pass
