# TASK-002: Lift `isolated_esgpull_context` into V2 downloader `download()` methods

## 1. Goal
Move the `isolated_esgpull_context` lifecycle from the utility functions (removed in TASK-001) into the V2 downloader `download()` methods, so each full download run opens exactly one isolated context and threads `esg` through all iterations.

## 2. Context & References
- **Plan section:** Steps 2 and 3 from `plan-refactor-esgpull-isolation.md`:
  > 2. Lift context creation into V2 downloader `download()` methods — wrap the entire iteration loop in a single `with isolated_esgpull_context(self.config.data_dir) as esg:` block.
  > 3. Update intermediate methods to thread `esg` through — Modify ... to accept and forward the `esg` parameter.
- **Upstream tasks:** TASK-001 — after which the three `esgpull_search_and_download_*` functions have this signature pattern:
  ```python
  def esgpull_search_and_download_esgf_model_single_var(
      esg: Esgpull,       # <-- NEW: injected instance
      model: str,
      variable: str,
      experiment: str,
      project: str,
      default_grid_label: str,
      default_frequency: str,
      preferred_version: str,
      max_ensemble_members: int,
      ensemble_members: list[str],
      data_dir: Path,     # <-- now Path only, used for dest_dir
      distrib: bool = False,
      logger: logging.Logger = LOGGER,
  ): ...
  ```
- **Key files:**
  - `climateset/download/cmip6_downloader.py` (lines 86–156) — `CMIP6DownloaderV2`
  - `climateset/download/input4mips_downloader.py` (lines 122–224) — `Input4MipsDownloaderV2`
- **Relevant skills:** `python` (DIP, composition), `systemdesign` (lifecycle management, SRP)

## 3. Subtasks
- [x] 1. Add `from climateset.download.esgpull_utils import isolated_esgpull_context` to `cmip6_downloader.py` (if not already imported)
- [x] 2. Wrap the triple-nested loop in `CMIP6DownloaderV2.download()` with `with isolated_esgpull_context(self.config.data_dir) as esg:`
- [x] 3. Add `esg: Esgpull` parameter to `CMIP6DownloaderV2.download_from_model_single_var()` and pass it through to `esgpull_search_and_download_esgf_model_single_var`
- [x] 4. Update the call site in `CMIP6DownloaderV2.download()` to pass `esg` to `self.download_from_model_single_var()`
- [x] 5. Add `from climateset.download.esgpull_utils import isolated_esgpull_context` to `input4mips_downloader.py` (if not already imported)
- [x] 6. Wrap all iteration loops in `Input4MipsDownloaderV2.download()` with a single `with isolated_esgpull_context(self.config.data_dir) as esg:` block
- [x] 7. Add `esg: Esgpull` parameter to `Input4MipsDownloaderV2.download_raw_input_single_var()` and `Input4MipsDownloaderV2.download_meta_historic_biomassburning_single_var()` — forward to the corresponding `esgpull_utils` function
- [x] 8. Update all call sites in `Input4MipsDownloaderV2.download()` to pass `esg`
- [x] 9. Verify syntax validity of both modified files

## 4. Requirements & Constraints
- **Technical:**
  - The `with` block must encompass **all** iteration loops in `download()`, not just one sub-loop. For `Input4MipsDownloaderV2`, this includes the raw variables loop, biomass loop, and metafiles loops.
  - `self.config.data_dir` may be `str` or `Path` — `isolated_esgpull_context` already handles coercion, so pass it directly.
  - When passing `data_dir` to the refactored `esgpull_utils` functions, coerce to `Path` at the call site: `Path(self.config.data_dir)`.
- **Business:** Download behaviour must be identical — same files, same destination paths, same error handling.
- **Out of scope:** Modifying V1 downloaders (`CMIP6Downloader`, `Input4MipsDownloader`), tests (TASK-004), or the `isolated_esgpull_context` function itself.

## 5. Acceptance Criteria
- [x] AC-1: `CMIP6DownloaderV2.download()` opens exactly one `isolated_esgpull_context` wrapping all iterations.
- [x] AC-2: `Input4MipsDownloaderV2.download()` opens exactly one `isolated_esgpull_context` wrapping all iterations.
- [x] AC-3: `download_from_model_single_var`, `download_raw_input_single_var`, and `download_meta_historic_biomassburning_single_var` all accept and forward `esg: Esgpull`.
- [x] AC-4: No direct calls to `isolated_esgpull_context` exist outside the two `download()` methods.
- [x] AC-5: `python3 -c "import ast; ast.parse(open('climateset/download/cmip6_downloader.py').read())"` exits 0.
- [x] AC-6: `python3 -c "import ast; ast.parse(open('climateset/download/input4mips_downloader.py').read())"` exits 0.
- [x] AC-7: `make precommit` and `make pylint` exit 0 for modified files.

## 6. Testing & Validation
```bash
# Syntax check both files
python3 -c "import ast; ast.parse(open('climateset/download/cmip6_downloader.py').read())"
python3 -c "import ast; ast.parse(open('climateset/download/input4mips_downloader.py').read())"
# Expected: exit 0, no output

# Lint
make precommit
make pylint
# Expected: exit 0

# Verify context is opened exactly once per download() method
grep -n "isolated_esgpull_context" climateset/download/cmip6_downloader.py
# Expected: import line + one call inside download()

grep -n "isolated_esgpull_context" climateset/download/input4mips_downloader.py
# Expected: import line + one call inside download()

# Verify esg is threaded through intermediate methods
grep -n "esg:" climateset/download/cmip6_downloader.py
# Expected: parameter in download_from_model_single_var signature

grep -n "esg:" climateset/download/input4mips_downloader.py
# Expected: parameter in download_raw_input_single_var and download_meta_historic_biomassburning_single_var
```

## 7. Completion Protocol
1. Verify every AC is checked off in Section 5.
2. Run all commands in Section 6 and confirm expected output.
3. Stage and commit with a scoped message:
   ```bash
   git add climateset/download/cmip6_downloader.py climateset/download/input4mips_downloader.py
   git commit -m "refactor(downloaders): lift esgpull context to batch scope in V2 downloaders — closes TASK-002"
   ```
4. Update this file: check off completed subtasks and ACs, note any deviations.
5. Notify the user with a concise summary and request approval before proceeding to TASK-003.
