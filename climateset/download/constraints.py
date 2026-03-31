from dataclasses import asdict, dataclass
from typing import Any

# Define a type alias for fields that can support esgpull multi-value lists
StrOrList = str | list[str] | None


@dataclass(frozen=True)
class BaseSearchConstraints:
    """
    Immutable base constraints for ESGF searches.

    Attributes:
        project (str | list[str] | None): The project name (e.g., "CMIP6").
        variable (str | list[str] | None): The variable name (e.g., "tas").
        frequency (str | list[str] | None): The frequency of the data (e.g., "mon").
        grid_label (str | list[str] | None): The grid label (e.g., "gn").
        nominal_resolution (str | list[str] | None): The nominal resolution (e.g., "100 km").
        version (str | list[str] | None): The version of the dataset (e.g., "20190101").
    """

    project: StrOrList = None
    variable: StrOrList = None
    frequency: StrOrList = None
    grid_label: StrOrList = None
    nominal_resolution: StrOrList = None
    version: StrOrList = None

    def to_esgf_params(self) -> dict[str, Any]:
        """
        Convert to ESGF search parameters (kwargs), filtering out None values.

        Returns:
            dict[str, Any]: A dictionary identifying parameters suitable for esgf-pyclient.
        """
        return {k: v for k, v in asdict(self).items() if v is not None}

    def to_esgpull_query(self) -> dict[str, Any]:
        """
        Convert constraints to parameters compatible with esgpull.models.Query selection. This explicitly handles multi-
        value lists, exclusion (!), and wildcard (*) capabilities provided by esgpull.

        Returns:
            dict[str, Any]: A dictionary suitable to be unpacked into an esgpull Query selection.
        """
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass(frozen=True)
class Input4MIPsConstraints(BaseSearchConstraints):
    """
    Constraints specific to Input4MIPs searches.

    Attributes:
        institution_id (str | list[str] | None): The institution ID.
        variable_id (str | list[str] | None): The variable ID.
        target_mip (str | list[str] | None): The target MIP.
    """

    institution_id: StrOrList = None
    variable_id: StrOrList = None
    target_mip: StrOrList = None


@dataclass(frozen=True)
class CMIP6Constraints(BaseSearchConstraints):
    """
    Constraints specific to CMIP6 searches.

    Attributes:
        experiment_id (str | list[str] | None): The experiment ID.
        source_id (str | list[str] | None): The source model ID.
        variant_label (str | list[str] | None): The variant label (ensemble member).
    """

    experiment_id: StrOrList = None
    source_id: StrOrList = None
    variant_label: StrOrList = None
