# TASK-004: Update existing tests and run full validation suite

## 1. Goal
Align all existing test files with the refactored function signatures from TASK-001/TASK-002, add a lifecycle assertion for single-context-per-batch, and confirm zero regressions across the full test suite and linting pipeline.

## 2. Context & References
- **Plan section:** Steps 6 and 7 from `plan-refactor-esgpull-isolation.md`:
  > 6. Update existing unit tests — Adjust `test_esgpull_utils.py` and any V2 downloader tests to reflect the new signatures. Add a lifecycle assertion that `isolated_esgpull_context` is entered exactly once per `download()` invocation.
  > 7. Run full test suite and linting — Execute `make test`, `make precommit`, `make pylint`, and `make mypy` to confirm no regressions.
- **Upstream tasks:**
  - TASK-001: `esgpull_search_and_download_*` now accept `esg: Esgpull` + `data_dir: Path` instead of `data_dir: Path | str` (no internal context).
  - TASK-002: V2 downloader `download()` methods wrap loops in `isolated_esgpull_context`; intermediate methods accept `esg: Esgpull`.
  - TASK-003: New test file `test_esgpull_context_reuse.py` exists and passes.
- **Key files:**
  - `tests/test_download/test_utils.py` — likely contains tests for esgpull utility functions
  - `tests/test_download/test_downloader.py` — likely contains tests for V2 downloaders
  - `tests/test_download/test_constraints.py` — should be unaffected (no signature changes)
  - `tests/test_download/test_search_client.py` — should be unaffected
- **Relevant skills:** `tdd` (green/refactor), `qa` (full validation pass), `python` (testing, type checking)

## 3. Subtasks
- [ ] 1. Read `tests/test_download/test_utils.py` and identify all tests that call `esgpull_search_and_download_*` functions directly — update call sites to pass a mock `Esgpull` instance as `esg` and `Path` as `data_dir`
- [ ] 2. Read `tests/test_download/test_downloader.py` and identify all tests that instantiate `CMIP6DownloaderV2` or `Input4MipsDownloaderV2` — update mocks/patches for the new `isolated_esgpull_context` call site (now in `download()`, not in `esgpull_utils`)
- [ ] 3. Add a lifecycle assertion: patch `isolated_esgpull_context` and assert it is called exactly **once** during a `download()` invocation with multiple variables
- [ ] 4. Run `make test-custom TEST_ARGS="tests/test_download/ -v"` and fix any failures
- [ ] 5. Run `make test` (full suite) and confirm no regressions
- [ ] 6. Run `make precommit`, `make pylint`, and `make mypy` across modified files

## 4. Requirements & Constraints
- **Technical:**
  - When mocking `Esgpull`, mock only the instance — not the `isolated_esgpull_context` context manager — in tests that exercise the utility functions directly. For downloader-level tests, mock the context manager to avoid filesystem side effects.
  - Do not modify V1 downloader tests — V1 code was not touched by this refactor.
- **Business:** All tests that passed before TASK-001 must still pass. Zero regressions.
- **Out of scope:** Writing new integration tests against live ESGF nodes (covered in TASK-003). Performance benchmarking.

## 5. Acceptance Criteria
- [ ] AC-1: All tests in `tests/test_download/test_utils.py` pass with the updated signatures.
- [ ] AC-2: All tests in `tests/test_download/test_downloader.py` pass with the updated V2 downloader code.
- [ ] AC-3: A test asserts `isolated_esgpull_context` is entered exactly once per `download()` call.
- [ ] AC-4: `make test` (full suite) exits 0 with no failures or errors.
- [ ] AC-5: `make precommit` exits 0.
- [ ] AC-6: `make pylint` exits 0 for modified files.
- [ ] AC-7: `make mypy` exits 0 for `climateset/download/esgpull_utils.py climateset/download/cmip6_downloader.py climateset/download/input4mips_downloader.py` (or only pre-existing errors unrelated to this change).

## 6. Testing & Validation
```bash
# Targeted test run
make test-custom TEST_ARGS="tests/test_download/ -v"
# Expected: all tests pass

# Full test suite
make test
# Expected: exit 0, no failures

# Lint
make precommit
make pylint
# Expected: exit 0

# Type check
make mypy
# Expected: exit 0 or only pre-existing errors
```

## 7. Completion Protocol
1. Verify every AC is checked off in Section 5.
2. Run all commands in Section 6 and confirm expected output.
3. Stage and commit with a scoped message:
   ```bash
   git add tests/test_download/
   git commit -m "test(esgpull): update tests for batch-scoped context refactor — closes TASK-004"
   ```
4. Update this file: check off completed subtasks and ACs, note any deviations.
5. Notify the user with a concise summary. This is the final task — the refactor is complete.
