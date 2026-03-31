# TASK-003: Verify `esg.context` statefulness across sequential queries

## 1. Goal
Prove that reusing a single `Esgpull` instance across multiple sequential search-and-download calls does not cause facet bleed — i.e., constraints from call N do not leak into call N+1.

## 2. Context & References
- **Plan section:** Step 5 from `plan-refactor-esgpull-isolation.md`:
  > 5. Verify `esg.context` statefulness — Write a focused integration test that calls `esgpull_search_and_download_esgf_model_single_var` twice with different `(variable, experiment)` pairs using the **same** `esg` instance, and asserts that each call returns files matching only its own constraints.
- **FMEA risk:** The plan identifies "`esg.context` carries stale facet state between variables" as **medium likelihood**.
- **Upstream tasks:** TASK-001 (signatures changed), TASK-002 (V2 downloaders updated). The refactored function signature from TASK-001:
  ```python
  def esgpull_search_and_download_esgf_model_single_var(
      esg: Esgpull, model: str, variable: str, experiment: str,
      project: str, default_grid_label: str, default_frequency: str,
      preferred_version: str, max_ensemble_members: int,
      ensemble_members: list[str], data_dir: Path,
      distrib: bool = False, logger: logging.Logger = LOGGER,
  ) -> list[Path] | None: ...
  ```
- **Key files:**
  - `climateset/download/esgpull_utils.py` — functions under test
  - `tests/test_download/test_utils.py` — existing test file (extend or add sibling)
- **Relevant skills:** `tdd` (Red-Green-Refactor), `python` (pytest fixtures)

## 3. Subtasks
- [x] 1. Create `tests/test_download/test_esgpull_context_reuse.py` with a focused integration test class
- [x] 2. Write a test that opens one `isolated_esgpull_context`, calls `esgpull_search_and_download_esgf_model_single_var` twice with different `(variable, experiment)` pairs, and asserts each call's `Query` was constructed with the correct constraints (TDD Red — expect this to pass if esgpull creates fresh `Query` objects per call, or fail if state bleeds)
- [x] 3. If facet bleed is detected: add a `Query` reset or fresh `Query()` construction guard at the top of each function — then re-run to green
- [x] 4. Add a lifecycle test asserting that `isolated_esgpull_context` creates exactly one `.esgpull_jobs/<UUID>` directory and cleans it up on exit
- [x] 5. Run the new tests in isolation to confirm they pass

## 4. Requirements & Constraints
- **Technical:**
  - Integration tests that hit real ESGF nodes should be marked with `@pytest.mark.integration` (or `@pytest.mark.slow`) so they can be skipped in CI fast-path runs.
  - Mock `esg.download` (the async download) to avoid bandwidth usage — only the search/hints layer should hit the network.
  - Use `tmp_path` fixture for the `data_dir` / `raw_data_path` to avoid polluting the workspace.
- **Business:** If facet bleed is confirmed, the fix must be minimal (e.g., ensuring `Query` is freshly constructed per call — which the current code already does). Do not refactor `esgpull` internals.
- **Out of scope:** Testing V1 downloaders. Performance benchmarking of shared vs. per-variable contexts.

## 5. Acceptance Criteria
- [x] AC-1: A test exists that reuses one `Esgpull` instance across 2+ calls to `esgpull_search_and_download_esgf_model_single_var` with different parameters.
- [x] AC-2: The test asserts that each call's search results contain only files matching its own constraints (no cross-contamination).
- [x] AC-3: A lifecycle test asserts exactly one UUID directory is created and cleaned up per `isolated_esgpull_context` usage.
- [x] AC-4: All new tests pass: `make test-custom TEST_ARGS="tests/test_download/test_esgpull_context_reuse.py -v"` exits 0.
- [x] AC-5: `make precommit` exits 0.

## 6. Testing & Validation
```bash
# Run the new test file
make test-custom TEST_ARGS="tests/test_download/test_esgpull_context_reuse.py -v"
# Expected: all tests pass

# Run with integration marker (if tests hit real nodes)
make test-marker TEST_ARGS="integration"
# Expected: integration tests pass (may be slow)

# Lint
make precommit
# Expected: exit 0

# Verify no regressions in existing tests
make test-custom TEST_ARGS="tests/test_download/ -v"
# Expected: all tests pass
```

## 7. Completion Protocol
1. Verify every AC is checked off in Section 5.
2. Run all commands in Section 6 and confirm expected output.
3. Stage and commit with a scoped message:
   ```bash
   git add tests/test_download/test_esgpull_context_reuse.py
   git commit -m "test(esgpull): add context reuse and lifecycle tests — closes TASK-003"
   ```
4. Update this file: check off completed subtasks and ACs, note any deviations.
5. Notify the user with a concise summary and request approval before proceeding to TASK-004.
