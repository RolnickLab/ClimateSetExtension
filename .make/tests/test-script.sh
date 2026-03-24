#!/usr/bin/env bash

# Set strict mode: exit on error, use unbound variables, fail on pipe errors
set -euo pipefail

# --- 1. Project Path Determination ---

SCRIPT_PATH="$(readlink -f "$0")"

# Get the directory of the script, its parent, and the project root
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
MAKE_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_DIR="$(dirname "$MAKE_DIR")"

# --- 2. Shared Configuration Variables ---

TEST_CONDA_ENV="lab-advanced-template-testing"
TEST_VENV_PATH="$PROJECT_DIR/.testvenv"
PIPX_TEST_ENV="$PROJECT_DIR/.testvenvpipx"

# Reusable Makefile argument strings
BASE_MAKEFILE_ARGS="-f $PROJECT_DIR/.make/base.make"
LINT_MAKEFILE_ARGS="-f $PROJECT_DIR/.make/lint.make"
CONDA_MAKEFILE_ARGS="-f $PROJECT_DIR/.make/conda.make"
POETRY_MAKEFILE_ARGS="-f $PROJECT_DIR/.make/poetry.make"
UV_MAKEFILE_ARGS="-f $PROJECT_DIR/.make/uv.make"

# Overrides for different configurations (now more focused)
MAKEFILE_CONDA_OVERRIDE="$BASE_MAKEFILE_ARGS $CONDA_MAKEFILE_ARGS $LINT_MAKEFILE_ARGS"
MAKEFILE_POETRY_OVERRIDE="$BASE_MAKEFILE_ARGS $CONDA_MAKEFILE_ARGS $POETRY_MAKEFILE_ARGS $LINT_MAKEFILE_ARGS"
MAKEFILE_UV_OVERRIDE="$BASE_MAKEFILE_ARGS $UV_MAKEFILE_ARGS $LINT_MAKEFILE_ARGS"

# Tool-specific arguments
POETRY_ARGS="DEFAULT_INSTALL_ENV=poetry DEFAULT_BUILD_TOOL=poetry VENV_PATH=$PROJECT_DIR/.venv"
UV_ARGS="DEFAULT_INSTALL_ENV=uv DEFAULT_BUILD_TOOL=uv"
POETRY_CONDA_ARGS="DEFAULT_INSTALL_ENV=conda DEFAULT_BUILT_TOOL=poetry CONDA_ENVIRONMENT=$TEST_CONDA_ENV"

# --- 3. Core Helper Functions ---

# Helper function to print section headers
print_header() {
  echo ""
  echo "###"
  echo "### $1"
  echo "###"
  echo ""
}

# Generic function to execute a make target with arguments
# Usage: make_test "Target Name" "MAKEFILE_ARGS" "CUSTOM_VARS" "make_target"
make_test() {
  local target_name="$1"
  local makefile_args="$2"
  local custom_vars="$3"
  local make_target="$4"

  echo "### Running 'make $make_target' ($target_name): ###"
  echo ""

  # Use command grouping (subshell) to ensure 'cd' doesn't affect the main script's environment
  (
    cd "$PROJECT_DIR" && \
    make $makefile_args $custom_vars $make_target
  )
}

# Generic function to run a sequence of tests
# Usage: run_test_suite "Suite Description" "MAKEFILE_ARGS" "CUSTOM_VARS" "TARGET_LIST"
run_test_suite() {
    local suite_description="$1"
    local makefile_args="$2"
    local custom_vars="$3"
    # Shift arguments to treat all subsequent args as the target list
    shift 3
    local targets=("$@")

    print_header "Test Suite: $suite_description"

    for target in "${targets[@]}"; do
        make_test "$suite_description" "$makefile_args" "$custom_vars" "$target"
    done
}


# --- 4. Setup and Cleanup Functions ---

test_cleanup() {
  print_header "Cleaning potential test environments"
  rm -rf "$PROJECT_DIR/.venv" "$TEST_VENV_PATH" "$PIPX_TEST_ENV" "$PROJECT_DIR/poetry.lock" "$PROJECT_DIR/uv.lock"
}

# --- 5. Specific Test Functions ---

base() {
  print_header "Test base targets"

  # Simple targets without special environment setup
  make_test "Base Targets" "" "" "targets"
  make_test "Base Targets" "" "" "version"

  # Test venv creation and removal with a temporary path
  run_test_suite \
      "Test temporary venv creation/removal" \
      "$MAKEFILE_POETRY_OVERRIDE" \
      "$POETRY_CONDA_ARGS VENV_PATH='$TEST_VENV_PATH'" \
      "venv-create" \
      "venv-remove AUTO_INSTALL=true"

  # Test versioning (bump targets) within a Conda environment managed by Poetry
  local version_targets=(
      "conda-create-env CONDA_ENVIRONMENT_FILE='$SCRIPT_DIR/test_environment.yml' CONDA_YES_OPTION='-y'"
      "conda-env-info"
      "conda-activate"
      "install AUTO_INSTALL=true"
      "bump-major dry"
      "bump-minor dry"
      "bump-patch dry"
      "conda-clean-env AUTO_INSTALL=true CONDA_YES_OPTION='-y'"
  )
  run_test_suite \
      "Test version bump targets (Poetry/Conda)" \
      "$MAKEFILE_POETRY_OVERRIDE" \
      "$POETRY_CONDA_ARGS" \
      "${version_targets[@]}"
}

conda() {
  local conda_targets=(
      "conda-create-env CONDA_ENVIRONMENT_FILE='$SCRIPT_DIR/test_environment.yml' CONDA_YES_OPTION='-y'"
      "conda-env-info"
      "conda-activate"
      "conda-clean-env AUTO_INSTALL=true CONDA_YES_OPTION='-y'"
  )
  run_test_suite \
      "Test core conda environment targets" \
      "$MAKEFILE_CONDA_OVERRIDE" \
      "$POETRY_CONDA_ARGS" \
      "${conda_targets[@]}"
}

# The 'lint' function is the most complex due to the matrix of combinations.
lint() {
    # Combination 1: Poetry (Default Venv Path)
    local poetry_venv_lint_targets=(
        "poetry-create-env" "install" "check-lint" "pylint" "check-complexity"
        "fix-lint" "precommit" "ruff" "ruff-fix" "ruff-format"
        "poetry-remove-env AUTO_INSTALL=true" "venv-remove AUTO_INSTALL=true"
    )
    run_test_suite \
        "Lint targets for Poetry (default venv)" \
        "$MAKEFILE_POETRY_OVERRIDE" \
        "$POETRY_ARGS" \
        "${poetry_venv_lint_targets[@]}"

    # Combination 2: Poetry (Conda-managed)
    local poetry_conda_lint_targets=(
        "conda-create-env CONDA_ENVIRONMENT_FILE='$SCRIPT_DIR/test_environment.yml' CONDA_YES_OPTION='-y'"
        "install" "check-lint" "pylint" "check-complexity"
        "fix-lint" "precommit" "ruff" "ruff-fix" "ruff-format"
        "conda-clean-env AUTO_INSTALL=true CONDA_YES_OPTION='-y'"
    )
    run_test_suite \
        "Lint targets for Poetry (Conda-managed)" \
        "$MAKEFILE_POETRY_OVERRIDE" \
        "$POETRY_CONDA_ARGS" \
        "${poetry_conda_lint_targets[@]}"

    # Combination 3: uv (Default Venv Path)
    local uv_lint_targets=(
        "uv-create-env" "install" "check-lint" "pylint" "check-complexity"
        "fix-lint" "precommit" "ruff" "ruff-fix" "ruff-format"
        "uv-remove-env AUTO_INSTALL=true"
    )
    run_test_suite \
        "Lint targets for uv" \
        "$MAKEFILE_UV_OVERRIDE" \
        "$UV_ARGS" \
        "${uv_lint_targets[@]}"
}

poetry() {
    # Test 1: Poetry managed environment
    local poetry_targets=(
        "poetry-create-env" "install" "poetry-remove-env AUTO_INSTALL=true" "venv-remove AUTO_INSTALL=true"
    )
    run_test_suite \
        "Poetry managed environment (poetry/venv)" \
        "$MAKEFILE_POETRY_OVERRIDE" \
        "$POETRY_ARGS" \
        "${poetry_targets[@]}"

    # Test 2: Conda managed environment
    local conda_targets=(
        "conda-create-env CONDA_ENVIRONMENT_FILE='$SCRIPT_DIR/test_environment.yml' CONDA_YES_OPTION='-y'"
        "install" "conda-clean-env AUTO_INSTALL=true"
    )
    run_test_suite \
        "Poetry managed environment (Conda)" \
        "$MAKEFILE_POETRY_OVERRIDE" \
        "$POETRY_CONDA_ARGS" \
        "${conda_targets[@]}"
}

poetry-pipx(){
  local pipx_poetry_targets=(
        "poetry-install-venv PIPX_VENV_PATH='$PIPX_TEST_ENV'"
        "poetry-uninstall-venv PIPX_VENV_PATH='$PIPX_TEST_ENV' AUTO_INSTALL=true"
        "poetry-install-venv" # Test default pipx path
  )
  run_test_suite \
        "Test pipx Poetry targets" \
        "$MAKEFILE_POETRY_OVERRIDE" \
        "$POETRY_ARGS" \
        "${pipx_poetry_targets[@]}"
}

uv(){
  local uv_targets=(
        "uv-create-env"
        "install AUTO_INSTALL=true"
        "uv-remove-env AUTO_INSTALL=true"
  )
  # Cleanup is done outside the suite for the specific cleanup step (rm -rf pyproject.toml.uv.backup)
  run_test_suite \
        "Test uv managed installs and migration" \
        "$MAKEFILE_UV_OVERRIDE" \
        "$UV_ARGS" \
        "${uv_targets[@]}"
}

uv-pipx(){
  local uv_pipx_targets=(
        "uv-install-venv PIPX_VENV_PATH='$PIPX_TEST_ENV'"
        "uv-uninstall-venv PIPX_VENV_PATH='$PIPX_TEST_ENV' AUTO_INSTALL=true"
  )
  run_test_suite \
        "Test pipx uv targets" \
        "$MAKEFILE_UV_OVERRIDE" \
        "$UV_ARGS" \
        "${uv_pipx_targets[@]}"
}

# --- 6. Execution Control Functions ---

all() {
  base
  conda
  lint
  poetry
  uv
}

list () {
    echo
    echo "!!! Do not run this script outside of the 'lab-advanced-template' repository."
    echo "This script exists only to test the makefiles for integrity when adding or"
    echo "modifying targets !!!"
    echo
    echo " List of available tests:"
    echo
    echo "    - base        : Test 'base' targets"
    echo "    - conda       : Test 'conda' targets"
    echo "    - lint        : Test 'linting' targets (Poetry/Conda/uv matrix)"
    echo "    - poetry      : Test 'poetry' targets (venv/Conda matrix)"
    echo "    - poetry-pipx : Test 'poetry' targets related to pipx"
    echo "    - uv          : Test 'uv' targets"
    echo "    - uv-pipx     : Test 'uv' targets related to pipx"

    echo
    echo " Full test suite"
    echo
    echo "    - all         : Run most tests, except '*-pipx'"
    echo

}

check_and_run_cleanup() {
    echo ""
    echo "🚨 WARNING: Executing cleanup will **permanently remove environments** (e.g., .venv, conda) and **lock files** (e.g., poetry.lock, requirements.lock)."
    read -r -p "Are you sure you want to continue with the cleanup? (y/n): " confirm_cleanup

    if [[ "$confirm_cleanup" =~ ^[Yy]$ ]]; then
        test_cleanup
        echo "Cleanup complete!"
    else
        echo "Cleanup aborted by user. Skipping environment removal."
    fi
}

# --- 7. Main Script Execution ---

if [[ "$#" -eq 0 ]]; then
    list
    exit 0
else
  check_and_run_cleanup
fi

for var in "$@"
do
    # Use a case statement to map arguments to function calls
    case "$var" in
        "list")
            list
            ;;
        "base")
            base
            ;;
        "conda")
            conda
            ;;
        "lint")
            lint
            ;;
        "poetry")
            poetry
            ;;
        "poetry-pipx")
            poetry-pipx
            ;;
        "uv")
            uv
            ;;
        "uv-pipx")
            uv-pipx
            ;;
        "all")
            all
            ;;
        *)
            echo "* * * * * * * * * * * * * * * * * * * * * * * * * "
            echo "* ""$var"" is not a valid input "
            echo "* Use the list command to see available inputs"
            echo "* * * * * * * * * * * * * * * * * * * * * * * * *"
            echo ""
            list
            exit 1
    esac
done
test_cleanup
