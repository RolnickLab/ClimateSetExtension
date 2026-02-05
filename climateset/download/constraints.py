from dataclasses import asdict, dataclass
from typing import Any, Optional


@dataclass(frozen=True)
class BaseSearchConstraints:
    """
    Immutable base constraints for ESGF searches.

    Attributes:
        project (str): The project name (e.g., "CMIP6").
        variable (str): The variable name (e.g., "tas").
        frequency (str): The frequency of the data (e.g., "mon").
        grid_label (str): The grid label (e.g., "gn").
        nominal_resolution (str): The nominal resolution (e.g., "100 km").
        version (str): The version of the dataset (e.g., "20190101").
    """

    project: Optional[str] = None
    variable: Optional[str] = None
    frequency: Optional[str] = None
    grid_label: Optional[str] = None
    nominal_resolution: Optional[str] = None
    version: Optional[str] = None

    def to_esgf_params(self) -> dict[str, Any]:
        """
        Convert to ESGF search parameters (kwargs), filtering out None values.

        Returns:
            dict[str, Any]: A dictionary identifying parameters suitable for esgf-pyclient.
        """
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass(frozen=True)
class Input4MIPsConstraints(BaseSearchConstraints):
    """
    Constraints specific to Input4MIPs searches.

    Attributes:
        institution_id (str): The institution ID.
        variable_id (str): The variable ID.
        target_mip (str): The target MIP.
    """

    institution_id: Optional[str] = None
    variable_id: Optional[str] = None
    target_mip: Optional[str] = None


@dataclass(frozen=True)
class CMIP6Constraints(BaseSearchConstraints):
    """
    Constraints specific to CMIP6 searches.

    Attributes:
        experiment_id (str): The experiment ID.
        source_id (str): The source model ID.
        variant_label (str): The variant label (ensemble member).
    """

    experiment_id: Optional[str] = None
    source_id: Optional[str] = None
    variant_label: Optional[str] = None
