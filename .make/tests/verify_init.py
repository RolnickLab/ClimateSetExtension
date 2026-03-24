#!/usr/bin/env python3
import importlib.util
import io
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

PROJECT_ROOT_ACTUAL = Path(__file__).resolve().parent.parent.parent
SCRIPT_REL_PATH = Path(".make/scripts/auto_init_script.py")
MAKE_FOLDER = ".make"


class TestAutoInit(unittest.TestCase):
    @staticmethod
    def _setup_test_env(temp_dir_path):
        """Copies essential project files to the temp directory."""
        temp_root = Path(temp_dir_path)

        items_to_copy = [
            "pyproject.toml",
            "Makefile.variables",
            "README.md",
            ".markdown-link-check.json",
            MAKE_FOLDER,
            "src",
            "CHANGES.md",
        ]

        for item in items_to_copy:
            src = PROJECT_ROOT_ACTUAL / item
            dst = temp_root / item
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True)
            elif src.exists():
                shutil.copy2(src, dst)

        # Create a dummy python file to test import replacement
        # We use 'my_awesome_project' as the import because that's what PLACEHOLDER_IMPORT_NAME is in the script
        dummy_file = temp_root / "src" / "dummy_script.py"
        dummy_file.write_text(
            "import my_awesome_project\nfrom my_awesome_project.utils import create_logger\n", encoding="utf-8"
        )

        # Initialize dummy git repo to test git detection
        try:
            subprocess.run(["git", "init"], cwd=temp_root, check=True, capture_output=True)
            subprocess.run(
                [
                    "git",
                    "remote",
                    "add",
                    "origin",
                    "https://githubbogusexample.com/testuser/testrepo.git",
                ],
                cwd=temp_root,
                check=True,
                capture_output=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Warning: Failed to initialize git repo in temp dir. Git tests might fail.")

        return temp_root

    @staticmethod
    def _load_module_from_path(name, path):
        """Loads a module dynamically from a file path."""
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load script from {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    def _verify_markdown_link_check(self, temp_root: Path):
        link_check = (temp_root / ".markdown-link-check.json").read_text()
        self.assertIn("https://githubbogusexample.com/testuser/testrepo", link_check)

    def _verify_self_update(self, python_version, script_path: Path, test_package_name: str):
        script_content = script_path.read_text()
        self.assertIn(f'PLACEHOLDER_PACKAGE_NAME = "{test_package_name}"', script_content)
        self.assertIn(
            'PLACEHOLDER_REPO_URL = "https://githubbogusexample.com/testuser/testrepo"',
            script_content,
        )
        self.assertIn(f'PLACEHOLDER_PYTHON_VERSION = "{python_version}"', script_content)

    def _verify_package_imports(self, temp_root: Path, test_package_name: str):
        dummy_file_after = temp_root / "src" / "dummy_script.py"
        dummy_content = dummy_file_after.read_text()
        self.assertIn(f"import {test_package_name}", dummy_content)
        self.assertIn(f"from {test_package_name}.utils import create_logger", dummy_content)

    def _verify_changes(self, temp_root: Path):
        changes = (temp_root / "CHANGES.md").read_text()
        self.assertIn("https://githubbogusexample.com/testuser/testrepo", changes)

    def _verify_readme(self, python_version, temp_root: Path, test_project_name: str):
        readme = (temp_root / "README.md").read_text()
        self.assertIn(f"# {test_project_name}", readme)
        self.assertNotIn("## 🚀 Template Initialization", readme)
        self.assertIn(f"This project uses **Python {python_version}**", readme)

    def _verify_makefile_variables(self, build_tool, install_env: str, python_version, temp_root: Path):
        makefile_vars = (temp_root / "Makefile.variables").read_text()
        self.assertIn(f"DEFAULT_INSTALL_ENV := {install_env}", makefile_vars)

        expected_build_tool = build_tool
        if install_env == "conda" and build_tool == "uv":
            expected_build_tool = "poetry"
        self.assertIn(f"DEFAULT_BUILD_TOOL := {expected_build_tool}", makefile_vars)
        self.assertIn(f"PYTHON_VERSION := {python_version}", makefile_vars)

    def _verify_pyprojecttoml(
        self, build_tool, install_env: str, python_version, temp_root: Path, test_package_name: str
    ):
        pyproject = (temp_root / "pyproject.toml").read_text()
        self.assertIn(f'name = "{test_package_name}"', pyproject)
        self.assertIn("https://githubbogusexample.com/testuser/testrepo", pyproject)

        # Build System Assertions
        if build_tool == "uv" and install_env != "conda":
            self.assertRegex(pyproject, r'(?m)^build-backend = "hatchling\.build"')
            self.assertNotRegex(pyproject, r'(?m)^build-backend = "poetry\.core\.masonry\.api"')
        elif build_tool == "poetry" or install_env == "conda":
            self.assertRegex(pyproject, r'(?m)^build-backend = "poetry\.core\.masonry\.api"')
            self.assertNotRegex(pyproject, r'(?m)^build-backend = "hatchling\.build"')

        # Python Version Assertion
        major, minor = python_version.split(".")
        next_minor = int(minor) + 1
        expected_requires = f'requires-python = ">={major}.{minor},<{major}.{next_minor}"'
        self.assertIn(expected_requires, pyproject)

        expected_python_target = f'target-version = ["py{major}{minor}"]'
        self.assertIn(expected_python_target, pyproject)

        expected_python_version = f'python_version = "{major}.{minor}"'
        self.assertIn(expected_python_version, pyproject)

    def _verify_directory_structure(self, temp_root: Path, test_package_name: str):
        expected_src = temp_root / "src" / test_package_name
        self.assertTrue(expected_src.exists(), f"Source directory {expected_src} should exist")
        self.assertFalse(
            (temp_root / "src" / "core").exists(),
            "Old source directory 'core' should be gone",
        )

    def run_init_scenario(self, install_env, build_tool, python_version, dry_run=False):
        """Runs a single initialization scenario in a fresh temporary directory."""

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = self._setup_test_env(temp_dir)
            script_path = temp_root / SCRIPT_REL_PATH

            module_name = (
                f"auto_init_{install_env}_{build_tool}_{python_version.replace('.', '')}_{'dry' if dry_run else 'full'}"
            )

            script_module = self._load_module_from_path(module_name, script_path)

            # Mock sys.argv
            test_package_name = f"test_pkg_{install_env}_{build_tool}"
            test_project_name = "Test Project"
            args = [
                str(script_path),
                "--project-name",
                test_project_name,
                "--package-name",
                test_package_name,
                "--description",
                "Test Description",
                "--author",
                "Test Author",
                "--email",
                "test@example.com",
                "--python-version",
                python_version,
                "--install-env",
                install_env,
            ]

            if install_env == "conda":
                args.extend(["--conda-env-name", "test-conda-env", "--conda-tool", "mamba"])

            args.extend(["--build-tool", build_tool])

            if dry_run:
                args.append("--dry-run")

            original_argv = sys.argv
            sys.argv = args

            f = io.StringIO()
            try:
                with redirect_stdout(f):
                    script_module.main()
            except SystemExit as e:
                if e.code != 0:
                    self.fail(f"Script exited with code {e.code}. Output:\n{f.getvalue()}")
            except Exception as e:
                self.fail(f"Script crashed: {e}. Output:\n{f.getvalue()}")
            finally:
                sys.argv = original_argv

            # --- Assertions ---
            if dry_run:
                # In dry run, 'core' should still exist
                self.assertTrue((temp_root / "src" / "core").exists())
                self.assertFalse((temp_root / "src" / test_package_name).exists())
                self.assertFalse((temp_root / MAKE_FOLDER / ".init_completed").exists())
                return

            # 1. Check Directory Structure
            self._verify_directory_structure(temp_root, test_package_name)

            # 2. Check pyproject.toml
            self._verify_pyprojecttoml(build_tool, install_env, python_version, temp_root, test_package_name)

            # 3. Check Makefile.variables
            self._verify_makefile_variables(build_tool, install_env, python_version, temp_root)

            # 4. Check README.md
            self._verify_readme(python_version, temp_root, test_project_name)

            # 5. Check CHANGES.md - Relaxed assertion
            self._verify_changes(temp_root)

            # 6. Check Import Updates
            self._verify_package_imports(temp_root, test_package_name)

            # 7. Check Marker File
            self.assertTrue(
                (temp_root / MAKE_FOLDER / ".init_completed").exists(),
                "Init marker should exist",
            )

            # 8. Check Script Self-Update
            self._verify_self_update(python_version, script_path, test_package_name)

            # 9. Check .markdown-link-check.json
            self._verify_markdown_link_check(temp_root)

    def test_matrix(self):
        """Iterate through all valid combinations."""
        combinations = [
            ("uv", "uv", "3.12"),
            ("poetry", "poetry", "3.13"),
            ("conda", "poetry", "3.12"),
            ("venv", "uv", "3.12"),
        ]
        print("")  # Spacing
        for install_env, build_tool, py_ver in combinations:
            with self.subTest(install_env=install_env, build_tool=build_tool, python_version=py_ver):
                print(
                    f"Testing combination: Env={install_env}, Build={build_tool}, Py={py_ver} ...",
                    end="",
                    flush=True,
                )
                self.run_init_scenario(install_env, build_tool, py_ver)
                print(" OK")

    def test_dry_run(self):
        """Tests that dry-run does not modify files."""
        print("Testing Dry Run ...", end="", flush=True)
        self.run_init_scenario("uv", "uv", "3.12", dry_run=True)
        print(" OK")


if __name__ == "__main__":
    unittest.main()
