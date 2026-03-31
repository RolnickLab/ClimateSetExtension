# TASK-05: Update Tests

## Goal
Refactor the test suite in `tests/test_download/test_downloader.py` (and any related test files) to instantiate the base downloaders with the appropriate `ClientType` config rather than importing and testing `*V2` classes separately.

## Context & References
- **Source Plan**: docs/agents/planning/refactor-download-compose.md
- **Relevant Specs**: N/A
- **Existing Code**: 
  - `tests/test_download/test_downloader.py`
  - Any mocks inside `tests/test_download/`

## Subtasks
1. [ ] Find tests covering `CMIP6DownloaderV2` and `Input4MipsDownloaderV2`.
2. [ ] Change imports to the base classes.
3. [ ] Modify the fixture or config initialization to explicitly set `client_type=ClientType.ESGPULL` for testing the `esgpull` branches.
4. [ ] Use `pytest.mark.parametrize` to run shared tests for both `PYESGF` and `ESGPULL` where appropriate to ensure identical behavior.

## Requirements & Constraints
- Coverage for both `pyesgf` and `esgpull` code paths must be maintained.
- Tests should not hang or fail due to improper mocking of the selected backend.

## Acceptance Criteria (AC)
- [ ] AC 1: All V2 tests are ported to use the base classes with config injection.
- [ ] AC 2: `make test` passes with 100% success rate across both `pyesgf` and `esgpull` configurations.

## Testing & Validation
- **Command**: `make test-custom TEST_ARGS="tests/test_download/test_downloader.py -v"`
- **Success State**: All tests pass.
- **Manual Verification**: Verify test logs run both branches of the conditional dispatch.

## Completion Protocol
1. [ ] All ACs are met.
2. [ ] Tests pass without regressions.
3. [ ] Code is linted via `make precommit` and `make pylint`.
4. [ ] Documentation updated (if applicable).
5. [ ] Commit work: `git commit -m "test(download): task 05 - update test suite for harmonized downloaders"`
6. [ ] Update this document: Mark as COMPLETE.