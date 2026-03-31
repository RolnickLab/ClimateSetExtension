# TASK-02: Harmonize CMIP6Downloader

## Goal
Harmonize the CMIP6 downloader so that it dispatches the actual download step to the correct utility backend based on the `client_type` configuration, effectively removing the need for `CMIP6DownloaderV2`.

## Context & References
- **Source Plan**: docs/agents/planning/refactor-download-compose.md
- **Relevant Specs**: N/A
- **Existing Code**: 
  - `climateset/download/cmip6_downloader.py`

## Subtasks
1. [ ] Update `CMIP6Downloader.download()` to manage the `isolated_esgpull_context` if `client_type == ClientType.ESGPULL`.
2. [ ] Refactor the inner loops of `CMIP6Downloader.download()` into a helper (or thread `esg` through) so that `download_from_model_single_var` accepts an optional `esg: Esgpull = None`.
3. [ ] If `ESGPULL` mode, call `esgpull_search_and_download_esgf_model_single_var` passing the injected `esg`.
4. [ ] Delete `CMIP6DownloaderV2` class from `cmip6_downloader.py`.

## Requirements & Constraints
- Both utility functions (`search_and_download_esgf_model_single_var` and its esgpull counterpart) must have consistent argument signatures as called from the class.

## Acceptance Criteria (AC)
- [ ] AC 1: `CMIP6Downloader` acts as a unified facade for both clients.
- [ ] AC 2: `CMIP6DownloaderV2` no longer exists in the codebase.
- [ ] AC 3: `CMIP6Downloader` successfully passes config parameters to both backends without TypeErrors.

## Testing & Validation
- **Command**: `make test-custom TEST_ARGS="tests/test_download/test_downloader.py"` (specifically tests targeting `CMIP6Downloader`)
- **Success State**: The harmonized class works for both client types.
- **Manual Verification**: Check if `make precommit` and `make pylint` pass after deleting `CMIP6DownloaderV2`.

## Completion Protocol
1. [ ] All ACs are met.
2. [ ] Tests pass without regressions.
3. [ ] Code is linted via `make precommit` and `make pylint`.
4. [ ] Documentation updated (if applicable).
5. [ ] Commit work: `git commit -m "refactor(download): task 02 - harmonize CMIP6Downloader and remove V2"`
6. [ ] Update this document: Mark as COMPLETE.