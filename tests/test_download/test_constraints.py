from climateset.download.constraints import (
    BaseSearchConstraints,
    CMIP6Constraints,
    Input4MIPsConstraints,
)


def test_base_constraints_to_esgf_params():
    """Test that BaseSearchConstraints correctly converts to ESGF parameters."""
    constraints = BaseSearchConstraints(project="CMIP6", variable="tas", frequency="mon")

    params = constraints.to_esgf_params()

    assert params == {"project": "CMIP6", "variable": "tas", "frequency": "mon"}

    # Verify strict None filtering
    assert "grid_label" not in params
    assert "version" not in params


def test_cmip6_constraints_inheritance():
    """Test that CMIP6Constraints includes fields from Base and its own."""
    constraints = CMIP6Constraints(
        project="CMIP6", experiment_id="ssp585", source_id="NorESM2-LM", variant_label="r1i1p1f1"
    )

    params = constraints.to_esgf_params()

    expected = {"project": "CMIP6", "experiment_id": "ssp585", "source_id": "NorESM2-LM", "variant_label": "r1i1p1f1"}
    assert params == expected


def test_input4mips_constraints_inheritance():
    """Test Input4MIPsConstraints serialization."""
    constraints = Input4MIPsConstraints(
        project="input4MIPs",
        target_mip="CMIP",
        institution_id="UoM",
    )

    params = constraints.to_esgf_params()

    expected = {
        "project": "input4MIPs",
        "target_mip": "CMIP",
        "institution_id": "UoM",
    }
    assert params == expected


def test_constraints_immutability():
    """Test that constraints are immutable (frozen)."""
    constraints = BaseSearchConstraints(project="CMIP6")

    try:
        constraints.project = "CMIP5"
        assert False, "Should have raised AttributeError"
    except AttributeError:
        pass  # Expected behavior


def test_base_constraints_multi_value_esgpull():
    """Test that constraints support list values and serialization to esgpull queries."""
    constraints = BaseSearchConstraints(project=["CMIP6", "input4MIPs"], variable=["tas", "pr"], frequency="mon")

    # The serialization should output the list values directly.
    esgf_params = constraints.to_esgf_params()
    esgpull_params = constraints.to_esgpull_query()

    expected = {"project": ["CMIP6", "input4MIPs"], "variable": ["tas", "pr"], "frequency": "mon"}

    assert esgf_params == expected
    assert esgpull_params == expected
