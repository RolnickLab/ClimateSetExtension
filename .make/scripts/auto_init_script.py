#!/usr/bin/env python3
"""
Initialization Script for Lab Advanced Template.

This script customizes the project based on user input and removes template-specific placeholders.
It handles:
- User input via CLI flags or interactive prompts.
- Automatic detection of repository URL.
- Variable replacement in configuration files.
- Directory renaming and import updates.
- Configuration updates in Makefile.variables.
- README customization.
- Link check configuration.
"""

import argparse
import json
import re
import shutil
import subprocess
from argparse import Namespace
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from typing import LiteralString
except ImportError:
    # Workaround for python<3.11
    LiteralString = str

# --- Constants ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Placeholders (These will be updated by the script itself after the first run)
PLACEHOLDER_PACKAGE_NAME = "core"
PLACEHOLDER_IMPORT_NAME = "my_awesome_project"
PLACEHOLDER_PROJECT_NAME = "My Awesome Project"
PLACEHOLDER_README_PROJECT_NAME = "\\<YOUR_PROJECT_NAME_HERE>"
PLACEHOLDER_DESCRIPTION_TOML = ""
PLACEHOLDER_DESCRIPTION_README = "[Provide a brief, one-sentence description of your project here.]"
PLACEHOLDER_PYTHON_VERSION = "3.12"
PLACEHOLDER_REPO_URL = "REPOSITORY_URL"
PLACEHOLDER_AUTHOR = "Author"
PLACEHOLDER_EMAIL = "author@example.com"

DEFAULT_PYTHON_VERSION = "3.12"
VALID_PYTHON_VERSIONS = ["3.11", "3.12", "3.13", "3.14"]
DEFAULT_INSTALL_ENV = "uv"
DEFAULT_BUILD_TOOL = "uv"
VALID_INSTALL_ENVS = ["uv", "poetry", "conda", "venv"]
VALID_BUILD_TOOLS = ["uv", "poetry"]

# Files to modify
PYPROJECT_TOML = PROJECT_ROOT / "pyproject.toml"
MAKEFILE_VARIABLES = PROJECT_ROOT / "Makefile.variables"
README_MD = PROJECT_ROOT / "README.md"
CHANGES_MD = PROJECT_ROOT / "CHANGES.md"
MARKDOWN_LINK_CHECK = PROJECT_ROOT / ".markdown-link-check.json"
INIT_MARKER_FILE = PROJECT_ROOT / ".make" / ".init_completed"

# Useful string variables
UV = "uv"
POETRY = "poetry"
CONDA = "conda"


# --- Helper Functions ---


def run_command(command: List[str], cwd: Path = PROJECT_ROOT, capture_output: bool = True) -> Optional[str]:
    """
    Runs a shell command and returns the output.

    Args:
        command: A list of strings representing the command to run.
        cwd: The working directory for the command. Defaults to PROJECT_ROOT.
        capture_output: Whether to capture stdout/stderr. Defaults to True.

    Returns:
        The stripped stdout string if successful, or None if the command fails.
    """
    try:
        result = subprocess.run(args=command, cwd=cwd, capture_output=capture_output, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        print(f"Error running command {' '.join(command)}")
        return None


def get_git_remote_url() -> Optional[str]:
    """
    Detects the git remote origin URL and converts it to HTTPS format.

    Returns:
        The HTTPS URL of the git remote, or None if not found.
    """
    url = run_command(command=["git", "remote", "get-url", "origin"])
    if not url:
        return None

    # Convert SSH to HTTPS
    # git@github.com:User/Repo.git -> https://github.com/User/Repo
    if url.startswith("git@"):
        url = url.replace(":", "/").replace("git@", "https://")

    if url.endswith(".git"):
        url = url[:-4]

    return url


def prompt_user(prompt: str, default: Optional[str] = None, choices: Optional[List[str]] = None) -> Optional[str]:
    """
    Prompts the user for input, with optional default and validation.

    Args:
        prompt: The text to display to the user.
        default: The default value to return if the user enters nothing.
        choices: A list of valid choices. If provided, input is validated against this list.

    Returns:
        The user's input or the default value.
    """
    while True:
        prompt_text = f"{prompt}"
        if default:
            prompt_text += f" - [{default}]"
        if choices:
            prompt_text += f" - Available choices: ({', '.join(choices)})"
        prompt_text += ": "

        value = input(prompt_text).strip()

        if not value and default:
            return default

        if not value and not default:
            print("Value is required.")
            continue

        if choices and value not in choices:
            print(f"Invalid choice. Must be one of: {', '.join(choices)}")
            continue

        return value


def replace_in_file(filepath: Path, replacements: Dict[str, str], dry_run: bool = False) -> None:
    """
    Replaces text in a file based on a dictionary of replacements.

    Args:
        filepath: Path to the file to modify.
        replacements: A dictionary where keys are strings to find and values are replacements.
        dry_run: If True, prints what would happen without modifying the file.
    """
    if not filepath.exists():
        print(f"Warning: File not found: {filepath}")
        return

    content = filepath.read_text(encoding="utf-8")

    new_content = content
    for search, replace in replacements.items():
        new_content = new_content.replace(search, replace)

    if content != new_content:
        if dry_run:
            print(f"[Dry Run] Would update {filepath}")
        else:
            filepath.write_text(data=new_content)
            print(f"Updated {filepath}")


def update_makefile_variables(
    filepath: Path,
    install_env: str,
    build_tool: str,
    python_version: str,
    app_name: str,
    conda_env: Optional[str],
    conda_tool: Optional[str],
    dry_run: bool,
) -> None:
    """
    Updates specific variables in Makefile.variables.

    Args:
        filepath: Path to the Makefile.variables file.
        install_env: The installation environment (e.g., 'uv', 'poetry', 'conda').
        build_tool: The build tool (e.g., 'uv', 'poetry').
        python_version: The Python version string.
        app_name: The application name.
        conda_env: The conda environment name (optional).
        conda_tool: The conda tool name (optional).
        dry_run: If True, prints changes without writing to file.
    """
    if not filepath.exists():
        print(f"Warning: File not found: {filepath}")
        return

    updates = {
        "DEFAULT_INSTALL_ENV :=": install_env,
        "DEFAULT_BUILD_TOOL :=": build_tool,
        "PYTHON_VERSION :=": python_version,
        "APPLICATION_NAME :=": app_name,
    }

    if not conda_env:
        conda_env = app_name
    updates["CONDA_ENVIRONMENT :="] = conda_env
    if conda_tool:
        updates["CONDA_TOOL :="] = conda_tool

    lines = filepath.read_text(encoding="utf-8").splitlines(keepends=True)

    new_lines = []
    for line in lines:
        matched_key = next((key for key in updates if line.startswith(key)), None)
        if matched_key:
            new_lines.append(f"{matched_key} {updates[matched_key]}\n")
        else:
            new_lines.append(line)

    if dry_run:
        print(f"[Dry Run] Would update {filepath} with configuration settings.")
    else:
        filepath.write_text(data="".join(new_lines), encoding="utf-8")
        print(f"Updated {filepath}")


def comment_block(match: re.Match[str]) -> LiteralString:
    """Utility function to help comment out sections of the pyproject.toml file, depending on the build tool
    selected."""
    block = match.group(0)
    lines = block.splitlines()
    new_lines = []
    for line in lines:
        if line.strip() == "":
            new_lines.append(line)
        elif line.startswith("# "):
            new_lines.append(line)
        else:
            new_lines.append(f"# {line}")
    return "\n".join(new_lines)


def uncomment_block(match: re.Match[str]) -> LiteralString:
    """Utility function to help un-comment sections of the pyproject.toml file, depending on the build tool selected."""
    block = match.group(0)
    lines = block.splitlines()
    new_lines = []
    for line in lines:
        if line.startswith("# "):
            new_lines.append(line[2:])
        elif line.startswith("#"):
            new_lines.append(line[1:])
        else:
            new_lines.append(line)
    return "\n".join(new_lines)


def update_pyproject_toml(
    filepath: Path,
    package_name: str,
    description: str,
    author: str,
    email: str,
    python_version: str,
    repo_url: Optional[str],
    build_tool: str,
    dry_run: bool,
) -> None:
    """
    Updates pyproject.toml with project metadata.

    Args:
        filepath: Path to pyproject.toml.
        package_name: The new package name.
        description: Project description.
        author: Author name.
        email: Author email.
        python_version: Python version string.
        repo_url: Repository URL.
        build_tool: The build tool to use ('uv' or 'poetry').
        dry_run: If True, prints changes without writing to file.
    """
    if not filepath.exists():
        print(f"Warning: File not found: {filepath}")
        return

    content = filepath.read_text(encoding="utf-8")

    # Basic replacements
    # Note: Using regex for more robust replacement where simple string replace might be ambiguous

    # Update name = "core" -> name = "package_name"
    # We look for name = "core" specifically in the [project] section usually at the top
    content = re.sub(
        pattern=f'name = "{PLACEHOLDER_PACKAGE_NAME}"', repl=f'name = "{package_name}"', string=content, count=1
    )

    # Update description
    content = re.sub(
        pattern=f'description = "{PLACEHOLDER_DESCRIPTION_TOML}"',
        repl=f'description = "{description}"',
        string=content,
        count=1,
    )

    # Update authors
    # authors = [{ name = "Francis Pelletier", email = "fplt.softwaredeveloper@gmail.com" }]
    new_authors = f'authors = [{{ name = "{author}", email = "{email}" }}]'
    content = re.sub(pattern=r"authors = \[.*]", repl=new_authors, string=content, count=1)

    # Update python version
    # requires-python = ">=3.12,<3.13"
    # We assume the user provides "3.12", we want ">=3.12,<3.13" logic or just update the base
    # For simplicity, let's try to construct the range if it matches X.Y format
    match = re.match(pattern=r"(\d+)\.(\d+)", string=python_version)
    if match:
        major, minor = map(int, match.groups())
        next_minor = minor + 1
        new_requires = f'requires-python = ">={major}.{minor},<{major}.{next_minor}"'
        content = re.sub(pattern=r'requires-python = ".*"', repl=new_requires, string=content, count=1)

    # Update tool.hatch.build.targets.wheel packages
    # packages = ["src/core"]
    content = content.replace(f'packages = ["src/{PLACEHOLDER_PACKAGE_NAME}"]', f'packages = ["src/{package_name}"]')

    # Update tool.poetry packages
    # packages = [
    #     { include = "core", from = "src" }
    # ]
    content = content.replace(f'include = "{PLACEHOLDER_PACKAGE_NAME}"', f'include = "{package_name}"')

    # Toggle Build Systems
    # Hatchling Block Pattern
    hatchling_pattern = (
        r"(?:# )?\[build-system\]\n"
        r"(?:# )?requires = \[\"hatchling\"\]\n"
        r"(?:# )?build-backend = \"hatchling\.build\"\n"
        r"\n"
        r"(?:# )?\[tool\.hatch\.build\.targets\.wheel\]\n"
        r"(?:# )?packages = \[\"src/[^\"]+\"\]"
    )

    # Poetry Block Pattern
    poetry_pattern = (
        r"(?:# )?\[build-system\]\n"
        r"(?:# )?requires = \[\"poetry-core\"\]\n"
        r"(?:# )?build-backend = \"poetry\.core\.masonry\.api\"\n"
        r"\n"
        r"(?:# )?\[tool\.poetry\]\n"
        r"(?:# )?packages = \[\n"
        r"(?:# )?    \{ include = \"[^\"]+\", from = \"src\" \}\n"
        r"(?:# )?\]"
    )

    if build_tool == "uv":
        # Enable Hatchling, Disable Poetry
        content = re.sub(pattern=hatchling_pattern, repl=uncomment_block, string=content, flags=re.MULTILINE)
        content = re.sub(pattern=poetry_pattern, repl=comment_block, string=content, flags=re.MULTILINE)
    elif build_tool == "poetry":
        # Disable Hatchling, Enable Poetry
        content = re.sub(pattern=hatchling_pattern, repl=comment_block, string=content, flags=re.MULTILINE)
        content = re.sub(pattern=poetry_pattern, repl=uncomment_block, string=content, flags=re.MULTILINE)

    # Update bumpversion files
    if repo_url:
        content = content.replace(PLACEHOLDER_REPO_URL, repo_url)

    # Update black target-version
    new_target_version = f'target-version = ["py{python_version.replace(".", "")}"]'
    content = re.sub(pattern=r"target-version = .*", repl=new_target_version, string=content, count=1)

    # Update mypy target-version
    new_target_version = f'python_version = "{python_version}"'
    content = re.sub(pattern=r"python_version = .+", repl=new_target_version, string=content, count=1)

    if dry_run:
        print(f"[Dry Run] Would update {filepath} with project metadata.")
    else:
        filepath.write_text(data=content, encoding="utf-8")
        print(f"Updated {filepath}")


def rename_package_directory(package_name: str, dry_run: bool) -> None:
    """
    Renames src/core to src/<package_name> and updates imports.

    Args:
        package_name: The new name for the package directory.
        dry_run: If True, prints changes without moving files or updating imports.
    """
    src_previous = PROJECT_ROOT / "src" / PLACEHOLDER_PACKAGE_NAME
    src_new = PROJECT_ROOT / "src" / package_name

    if not src_previous.exists():
        print(f"Warning: {src_previous} does not exist. Skipping rename.")
        return

    if dry_run:
        print(f"[Dry Run] Would rename {src_previous} to {src_new}")
    else:
        shutil.move(src=src_previous, dst=src_new)
        print(f"Renamed {src_previous} to {src_new}")

    # Update imports in all .py files
    # We need to walk through the project and replace "from my_awesome_project" with "from package_name"
    # and "import my_awesome_project" with "import package_name"

    # Directories to skip
    skip_dirs = {".git", ".venv", "__pycache__", ".nox", ".idea"}

    for file_path in PROJECT_ROOT.rglob("*.py"):
        # Check if any part of the path is in skip_dirs
        if any(part in skip_dirs for part in file_path.parts):
            continue

        # Skip this script itself to avoid self-modification during this step
        if file_path.resolve() == Path(__file__).resolve():
            continue

        content = file_path.read_text(encoding="utf-8")

        # Simple replacements for imports
        # This is a basic heuristic and might need refinement for complex cases
        new_content = content.replace(f"from {PLACEHOLDER_IMPORT_NAME}", f"from {package_name}")
        new_content = new_content.replace(f"import {PLACEHOLDER_IMPORT_NAME}", f"import {package_name}")

        if content != new_content:
            if dry_run:
                print(f"[Dry Run] Would update imports in {file_path}")
            else:
                file_path.write_text(data=new_content, encoding="utf-8")
                print(f"Updated imports in {file_path}")


def search_string(content: str, pattern: str) -> re.Match[str] | None:
    result = re.search(pattern=pattern, string=content)
    return result


def open_close_detail_sections(content: str, install_env: str) -> str:
    """
    Helper utility to handle changing the opened and closed detail sections.

    Args:
        content: The content to modify.
        install_env: The installation environment used.

    Returns:
        The modified content.
    """
    # Patterns to capture the details block.
    # Group 1: The opening tag (<details> or <details open>)
    # Group 2: The content (summary ... /details) including trailing newlines

    uv_pattern = r"(<details(?: open)?>)(\n<summary><strong>Stack: uv </strong></summary>[\s\S]*?</details>\n*)"
    poetry_pattern = r"(<details(?: open)?>)( <summary><strong>Stack: Poetry</strong></summary>[\s\S]*?</details>\n*)"
    conda_pattern = (
        r"(<details(?: open)?>)( <summary><strong>Stack: Poetry \+ Conda</strong></summary>[\s\S]*?</details>\n*)"
    )

    patterns = {
        UV: uv_pattern,
        POETRY: poetry_pattern,
        CONDA: conda_pattern,
    }

    working_content = content

    for env, pattern in patterns.items():
        if env == install_env:
            # Ensure Open
            replacement = r"<details open>\2"
        else:
            # Ensure Closed
            replacement = r"<details>\2"

        working_content = re.sub(pattern=pattern, repl=replacement, string=working_content)

    return working_content


def update_readme(
    readme_path: Path, project_name: str, description: str, install_env: str, python_version: str, dry_run: bool
) -> None:
    """
    Updates README.md content.

    Args:
        readme_path: Path to README.md.
        project_name: The project name.
        description: The project description.
        install_env: The installation environment used.
        python_version: The Python version used.
        dry_run: If True, prints changes without writing to file.
    """
    if not readme_path.exists():
        return

    content = readme_path.read_text(encoding="utf-8")

    # Replace Title and Description
    content = re.sub(
        pattern=f"# {re.escape(PLACEHOLDER_README_PROJECT_NAME)}", repl=f"# {project_name}", string=content
    )
    content = re.sub(pattern=re.escape(PLACEHOLDER_DESCRIPTION_README), repl=description, string=content)

    # Remove Template Initialization Section
    # We look for the section start and end
    start_marker = "## 🚀 Template Initialization"
    end_marker = "## 🐍 Python Version"

    start_idx = content.find(start_marker)
    end_idx = content.find(end_marker)

    if start_idx != -1 and end_idx != -1:
        # Keep the end marker section
        content = content[:start_idx] + content[end_idx:]

    # Update Python Version section text
    content = re.sub(
        pattern=f"This project uses \\*\\*Python {re.escape(PLACEHOLDER_PYTHON_VERSION)}\\*\\*",
        repl=f"This project uses **Python {python_version}**",
        string=content,
    )

    # Dynamic Content Removal based on stack
    content = open_close_detail_sections(content=content, install_env=install_env)

    if dry_run:
        print(f"[Dry Run] Would update {readme_path}")
    else:
        readme_path.write_text(data=content, encoding="utf-8")
        print(f"Updated {readme_path}")


def update_link_check(filepath: Path, repo_url: Optional[str], dry_run: bool) -> None:
    """
    Adds repository URL to .markdown-link-check.json ignore patterns.

    Args:
        filepath: Path to .markdown-link-check.json.
        repo_url: The repository URL to add to ignore patterns.
        dry_run: If True, prints changes without writing to file.
    """
    if not repo_url or not filepath.exists():
        return

    try:
        data = json.loads(filepath.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        print(f"Error reading {filepath}")
        return

    # Check if already exists
    patterns = data.get("ignorePatterns", [])
    if not any(p.get("pattern") == repo_url for p in patterns):
        patterns.append({"pattern": repo_url})
        data["ignorePatterns"] = patterns

        # Also replace REPOSITORY_URL placeholder if present
        new_patterns = []
        for p in patterns:
            if p.get("pattern") == f"{PLACEHOLDER_REPO_URL}/tree/main":
                new_patterns.append({"pattern": f"{repo_url}/tree/main"})
            else:
                new_patterns.append(p)
        data["ignorePatterns"] = new_patterns

        if dry_run:
            print(f"[Dry Run] Would update {filepath} with repo URL.")
        else:
            filepath.write_text(data=json.dumps(obj=data, indent=2), encoding="utf-8")
            print(f"Updated {filepath}")


def update_self(script_path: Path, replacements: Dict[str, str], dry_run: bool) -> None:
    """Updates the script's own constants to match the new project state."""
    if not script_path.exists():
        return

    content = script_path.read_text(encoding="utf-8")

    for key, value in replacements.items():
        # Look for KEY = "..." or KEY = '...'
        # We use json.dumps to generate a safe string representation (e.g. "value")
        # and replace the existing assignment.
        # We assume the constants are defined at the top level.

        # Regex explanation:
        # ^: Start of line
        # {key}: The constant name
        # \s*=\s*: Assignment operator with optional whitespace
        # .*$: Match the rest of the line (the value)

        new_line = f"{key} = {json.dumps(value)}"
        content = re.sub(pattern=rf"^{key}\s*=\s*.*$", repl=new_line, string=content, flags=re.MULTILINE)

    if dry_run:
        print(f"[Dry Run] Would update {script_path} constants.")
    else:
        script_path.write_text(data=content, encoding="utf-8")
        print(f"Updated {script_path} constants for future runs.")


def self_update_for_next_run_of_script(
    args: Namespace,
    author: str | Any,
    description: str | Any,
    email: str | Any,
    package_name: str | Any,
    project_name: str | Any,
    python_version: str | Any,
    repo_url: str | None,
) -> None:
    self_replacements = {
        "PLACEHOLDER_PACKAGE_NAME": package_name,
        "PLACEHOLDER_IMPORT_NAME": package_name,
        "PLACEHOLDER_PROJECT_NAME": project_name,
        "PLACEHOLDER_README_PROJECT_NAME": project_name,
        "PLACEHOLDER_DESCRIPTION_TOML": description,
        "PLACEHOLDER_PYTHON_VERSION": python_version,
        "DEFAULT_PYTHON_VERSION": python_version,
        "PLACEHOLDER_AUTHOR": author,
        "PLACEHOLDER_EMAIL": email,
    }

    if repo_url:
        self_replacements["PLACEHOLDER_REPO_URL"] = repo_url

    update_self(script_path=Path(__file__), replacements=self_replacements, dry_run=args.dry_run)


def gather_build_and_env_fields(args: Namespace, package_name: str | Any) -> tuple[Any, Any, Any, Any]:
    install_env = args.install_env
    if not install_env:
        install_env = prompt_user(prompt="Install Environment", default=DEFAULT_INSTALL_ENV, choices=VALID_INSTALL_ENVS)

    # Adjust build tool choices based on install env
    available_build_tools = VALID_BUILD_TOOLS
    if install_env == "conda" and "uv" in available_build_tools:
        # Logic from plan: If conda is selected remove uv from the choices of build-tools
        # and add explanation that only poetry is available when using conda
        print(
            "Note: When using Conda for environment management, 'poetry' is the only supported build tool in this "
            "template."
        )
        available_build_tools = ["poetry"]
    if install_env == "uv" and "poetry" in available_build_tools:
        # Logic from plan: If conda is selected remove uv from the choices of build-tools
        # and add explanation that only poetry is available when using conda
        print("Note: When using UV for environment management, 'uv' is the only supported build tool in this template.")
        available_build_tools = ["uv"]
    if install_env == "poetry" and "uv" in available_build_tools:
        # Logic from plan: If conda is selected remove uv from the choices of build-tools
        # and add explanation that only poetry is available when using conda
        print(
            "Note: When using Poetry for environment management, 'poetry' is the only supported build tool in this "
            "template."
        )
        available_build_tools = ["poetry"]

    build_tool = args.build_tool
    if not build_tool:
        default_bt = "poetry" if install_env in ["conda", "poetry"] else DEFAULT_BUILD_TOOL
        build_tool = prompt_user(prompt="Build Tool", default=default_bt, choices=available_build_tools)

    # Validate build_tool with install_env
    if install_env == "conda" and build_tool == "uv":
        print("Warning: 'uv' build tool is not supported with 'conda' environment in this template.")
        print("Switching build tool to 'poetry'.")
        build_tool = "poetry"

    conda_env_name = None
    conda_tool = None
    if install_env == "conda":
        conda_env_name = args.conda_env_name or prompt_user(prompt="Conda Environment Name", default=f"{package_name}")
        conda_tool = args.conda_tool or prompt_user(prompt="Conda Tool", default="mamba", choices=["mamba", "conda"])
    return build_tool, conda_env_name, conda_tool, install_env


def gather_metadata_fields(args: Namespace) -> tuple[str | Any, str | Any, str | Any, str | Any, str | Any, str | Any]:
    project_name = args.project_name or prompt_user(prompt="Project Name", default=PLACEHOLDER_PROJECT_NAME)
    package_name_args = [
        "Package Name (snake_case)",
        PLACEHOLDER_IMPORT_NAME.lower().replace(" ", "_").replace("-", "_"),
    ]
    if project_name != PLACEHOLDER_PROJECT_NAME and project_name:
        package_name_args = ["Package Name (snake_case)", project_name.lower().replace(" ", "_").replace("-", "_")]
    package_name = args.package_name or prompt_user(prompt=package_name_args[0], default=package_name_args[1])
    description = args.description or prompt_user(prompt="Project Description", default=PLACEHOLDER_DESCRIPTION_TOML)
    author = args.author or prompt_user(prompt="Author Name", default=PLACEHOLDER_AUTHOR)
    email = args.email or prompt_user(prompt="Author Email", default=PLACEHOLDER_EMAIL)
    python_version = args.python_version or prompt_user(
        prompt="Python Version", default=DEFAULT_PYTHON_VERSION, choices=VALID_PYTHON_VERSIONS
    )
    return author, description, email, package_name, project_name, python_version


def generate_args_and_parser() -> Namespace:
    parser = argparse.ArgumentParser(description="Initialize the project from the template.")

    # Project Metadata
    parser.add_argument("--project-name", help="Name of the project")
    parser.add_argument("--package-name", help="Python package name (snake_case)")
    parser.add_argument("--description", help="Brief project description")
    parser.add_argument("--author", help="Author's name")
    parser.add_argument("--email", help="Author's email")

    # Technical Configuration
    parser.add_argument("--python-version", choices=VALID_PYTHON_VERSIONS, help="Target Python version")
    parser.add_argument("--install-env", choices=VALID_INSTALL_ENVS, help="Tool for virtual environment management")
    parser.add_argument("--build-tool", choices=VALID_BUILD_TOOLS, help="Tool for dependency/build management")

    # Conda Specifics
    parser.add_argument("--conda-env-name", help="Name of the conda environment")
    parser.add_argument("--conda-tool", choices=["mamba", "conda"], help="Tool to use (mamba or conda)")

    # Flags
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without writing to disk")

    args = parser.parse_args()
    return args


# --- Main Execution ---


def main() -> None:
    """
    Main execution function for the initialization script.

    Parses arguments, gathers user input, and triggers file updates.
    """
    args = generate_args_and_parser()

    if INIT_MARKER_FILE.exists() and not args.dry_run:
        print("\n⚠️  WARNING: It looks like this project has already been initialized.")
        print(f"Marker file exists: {INIT_MARKER_FILE}")
        print("Re-running this script might overwrite your changes or cause unexpected behavior.")

        should_continue = prompt_user(prompt="Do you want to continue anyway?", default="no", choices=["yes", "no"])
        if should_continue != "yes":
            print("Aborting initialization.")
            return

    print("🚀 Starting Project Initialization...")

    # 1. Gather Inputs
    author, description, email, package_name, project_name, python_version = gather_metadata_fields(args=args)

    build_tool, conda_env_name, conda_tool, install_env = gather_build_and_env_fields(
        args=args, package_name=package_name
    )

    # 2. Automatic Repo Detection
    repo_url = get_git_remote_url()
    if repo_url:
        print(f"Detected Repository URL: {repo_url}")
    else:
        print("Could not detect repository URL. Placeholders will remain.")

    # 3. Execution
    print("\nApplying changes...")

    # Update Makefile.variables
    update_makefile_variables(
        filepath=MAKEFILE_VARIABLES,
        install_env=install_env,
        build_tool=build_tool,
        python_version=python_version,
        app_name=package_name,
        conda_env=conda_env_name,
        conda_tool=conda_tool,
        dry_run=args.dry_run,
    )

    # Update pyproject.toml
    update_pyproject_toml(
        filepath=PYPROJECT_TOML,
        package_name=package_name,
        description=description,
        author=author,
        email=email,
        python_version=python_version,
        repo_url=repo_url,
        build_tool=build_tool,
        dry_run=args.dry_run,
    )

    # Rename Directory and Update Imports
    if package_name != PLACEHOLDER_PACKAGE_NAME:
        rename_package_directory(package_name=package_name, dry_run=args.dry_run)

    # Update README.md
    update_readme(
        readme_path=README_MD,
        project_name=project_name,
        description=description,
        install_env=install_env,
        python_version=python_version,
        dry_run=args.dry_run,
    )

    # Update CHANGES.md
    replace_in_file(
        filepath=CHANGES_MD,
        replacements={
            PLACEHOLDER_PROJECT_NAME: project_name,
            PLACEHOLDER_README_PROJECT_NAME: project_name,
            PLACEHOLDER_AUTHOR: author,
            PLACEHOLDER_REPO_URL: repo_url if repo_url else PLACEHOLDER_REPO_URL,
        },
        dry_run=args.dry_run,
    )

    # Update .markdown-link-check.json
    update_link_check(filepath=MARKDOWN_LINK_CHECK, repo_url=repo_url, dry_run=args.dry_run)

    # 4. Update self (the script itself) to prepare for next run
    self_update_for_next_run_of_script(
        args=args,
        author=author,
        description=description,
        email=email,
        package_name=package_name,
        project_name=project_name,
        python_version=python_version,
        repo_url=repo_url,
    )

    print("\n✅ Initialization Complete!")
    if args.dry_run:
        print("(This was a dry run. No files were modified.)")
    else:
        INIT_MARKER_FILE.touch()
        print("\n🔍  Please review the changes and commit them.")
        print("\n📦  You can now process to installing the package using the 'make install' command.")


if __name__ == "__main__":
    main()
