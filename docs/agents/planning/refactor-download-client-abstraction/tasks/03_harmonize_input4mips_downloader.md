# TASK-03: Harmonize Input4MipsDownloader

## Goal
Harmonize the Input4Mips downloader so that it dispatches the actual download steps to the correct utility backends based on the `client_type` configuration, effectively removing the need for `Input4MipsDownloaderV2`.

## Context & References
- **Source Plan**: docs/agents/planning/refactor-download-compose.md
- **Relevant Specs**: N/A
- **Existing Code**: 
  - `climateset/download/input4mips_downloader.py`

## Subtasks
1. [ ] Update `Input4MipsDownloader.download()` to conditionally manage the `isolated_esgpull_context` when `ClientType.ESGPULL` is selected.
2. [ ] Modify `Input4MipsDownloader.download_raw_input_single_var` to accept an optional `esg: Esgpull = None`.
3. [ ] Modify `Input4MipsDownloader.download_meta_historic_biomassburning_single_var` to accept an optional `esg: Esgpull = None`.
4. [ ] In `ESGPULL` mode, forward `esg` to the corresponding `esgpull_utils` search-and-download functions.
5. [ ] Delete `Input4MipsDownloaderV2` class from `input4mips_downloader.py`.

## Requirements & Constraints
- Both utility functions for both standard inputs and biomass must accept identical parameters from the caller.

## Acceptance Criteria (AC)
- [ ] AC 1: `Input4MipsDownloader` acts as a unified facade for both clients.
- [ ] AC 2: `Input4MipsDownloaderV2` no longer exists in the codebase.
- [ ] AC 3: Dispatches to standard and biomass endpoints work for both client types.

## Testing & Validation
- **Command**: `make test-custom TEST_ARGS="tests/test_download/test_downloader.py"` (specifically tests targeting `Input4MipsDownloader`)
- **Success State**: The harmonized class works for both client types.
- **Manual Verification**: Verify `download_raw_input_single_var` branches correctly using `make test`.

## Completion Protocol
1. [ ] All ACs are met.
2. [ ] Tests pass without regressions.
3. [ ] Code is linted via `make precommit` and `make pylint`.
4. [ ] Documentation updated (if applicable).
5. [ ] Commit work: `git commit -m "refactor(download): task 03 - harmonize Input4MipsDownloader and remove V2"`
6. [ ] Update this document: Mark as COMPLETE.