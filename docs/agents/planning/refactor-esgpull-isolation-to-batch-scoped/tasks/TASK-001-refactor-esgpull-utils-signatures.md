# TASK-001: Refactor search-and-download functions to accept an injected `Esgpull` instance

## 1. Goal
Replace the per-call `isolated_esgpull_context` usage inside the three `esgpull_search_and_download_*` functions with an externally provided `esg: Esgpull` parameter, so callers control the context lifecycle.

## 2. Context & References
- **Plan section:** Steps 1 and 4 from `plan-refactor-esgpull-isolation.md`:
  > 1. Add `esg: Esgpull` parameter to search-and-download functions — Modify ... to accept an `esg: Esgpull` parameter instead of `data_dir: Path | str`. Remove the internal `with isolated_esgpull_context(...)` block from each function body.
  > 4. Re-derive `dest_dir` from `data_dir` inside search-and-download functions — ... add `data_dir` as a separate parameter alongside `esg` ... Option (a) is simpler and more explicit — prefer it.
- **Upstream tasks:** None — this is the first task.
- **Key files:**
  - `climateset/download/esgpull_utils.py` (lines 56–167) — the three functions to modify
- **Relevant skills:** `python` (DIP, strict typing), `systemdesign` (decoupling context lifecycle from business logic)

### Current signatures (to be changed)
```python
# esgpull_utils.py:56
def esgpull_search_and_download_esgf_raw_single_var(
    variable: str, institution_id: str, project: str,
    default_grid_label: str, default_frequency: str, preferred_version: str,
    data_dir: Path | str,  # <-- drives both isolation AND dest_dir
    distrib: bool = False, logger: logging.Logger = LOGGER,
): ...

# esgpull_utils.py:94
def esgpull_search_and_download_esgf_biomass_single_var(
    variable: str, variable_id: str, institution_id: str, project: str,
    default_grid_label: str, default_frequency: str, preferred_version: str,
    data_dir: Path | str,
    distrib: bool = False, logger: logging.Logger = LOGGER,
): ...

# esgpull_utils.py:129
def esgpull_search_and_download_esgf_model_single_var(
    model: str, variable: str, experiment: str, project: str,
    default_grid_label: str, default_frequency: str, preferred_version: str,
    max_ensemble_members: int, ensemble_members: list[str],
    data_dir: Path | str,
    distrib: bool = False, logger: logging.Logger = LOGGER,
): ...
```

## 3. Subtasks
- [x] 1. Replace `data_dir: Path | str` with `esg: Esgpull` and `data_dir: Path` (keeping `data_dir` for `dest_dir` derivation) in `esgpull_search_and_download_esgf_raw_single_var`
- [x] 2. Remove the `with isolated_esgpull_context(data_dir) as esg:` block — dedent the function body one level so it uses the injected `esg` directly
- [x] 3. Repeat subtasks 1–2 for `esgpull_search_and_download_esgf_biomass_single_var`
- [x] 4. Repeat subtasks 1–2 for `esgpull_search_and_download_esgf_model_single_var`
- [x] 5. Verify `isolated_esgpull_context` is no longer called inside any of the three functions (but remains defined in the module for external use)
- [x] 6. Run `python3 -c "import ast; ast.parse(open('climateset/download/esgpull_utils.py').read())"` to confirm syntax validity

## 4. Requirements & Constraints
- **Technical:**
  - `esg` parameter must be typed as `Esgpull` (imported from `esgpull`).
  - `data_dir` parameter changes type from `Path | str` to `Path` — the `str` coercion responsibility moves to the caller (the V2 downloader).
  - `isolated_esgpull_context()` definition must remain in the module — it is still the public API for creating contexts. Only its **call sites** inside the three functions are removed.
- **Business:** The function bodies must remain functionally identical — only the context acquisition changes.
- **Out of scope:** Modifying V2 downloader classes or tests — those are TASK-002 and TASK-004.

## 5. Acceptance Criteria
- [x] AC-1: All three `esgpull_search_and_download_*` functions accept `esg: Esgpull` as their first parameter.
- [x] AC-2: All three functions retain `data_dir: Path` as a parameter (used only for `dest_dir` derivation).
- [x] AC-3: No call to `isolated_esgpull_context` exists inside any of the three function bodies.
- [x] AC-4: `isolated_esgpull_context` remains defined and importable from `esgpull_utils`.
- [x] AC-5: `python3 -c "import ast; ast.parse(open('climateset/download/esgpull_utils.py').read())"` exits 0.
- [x] AC-6: `make precommit` and `make pylint` exit 0 for `climateset/download/esgpull_utils.py`.

## 6. Testing & Validation
```bash
# Syntax check
python3 -c "import ast; ast.parse(open('climateset/download/esgpull_utils.py').read())"
# Expected: exit 0, no output

# Lint
make precommit
make pylint
# Expected: exit 0 (or only pre-existing warnings unrelated to this change)

# Verify isolated_esgpull_context is still defined
grep -n "def isolated_esgpull_context" climateset/download/esgpull_utils.py
# Expected: one match (the definition)

# Verify no internal calls remain
grep -n "isolated_esgpull_context" climateset/download/esgpull_utils.py
# Expected: only the definition line and any import/export — no calls inside the three functions
```

## 7. Completion Protocol
1. Verify every AC is checked off in Section 5.
2. Run all commands in Section 6 and confirm expected output.
3. Stage and commit with a scoped message:
   ```bash
   git add climateset/download/esgpull_utils.py
   git commit -m "refactor(esgpull): inject Esgpull instance into search-and-download functions — closes TASK-001"
   ```
4. Update this file: check off completed subtasks and ACs, note any deviations.
5. Notify the user with a concise summary and request approval before proceeding to TASK-002.
