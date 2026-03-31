# TASK-01: Client Configuration

## Goal
Introduce a `ClientType` Enum to configure which ESGF client (`pyesgf` or `esgpull`) the downlaoder should use, and add it to `BaseDownloaderConfig`.

## Context & References
- **Source Plan**: docs/agents/planning/refactor-download-compose.md
- **Relevant Specs**: N/A
- **Existing Code**: 
  - `climateset/download/downloader_config.py`

## Subtasks
1. [ ] Define a `ClientType` Enum (e.g. `PYESGF` and `ESGPULL`).
2. [ ] Add `client_type: ClientType = ClientType.PYESGF` property to `BaseDownloaderConfig`'s `__init__`.
3. [ ] Update YAML configuration parsing or ensure that standard parsing can handle mapping strings like "pyesgf" or "esgpull" to the `ClientType` Enum if provided in the kwargs.

## Requirements & Constraints
- Must default to `ClientType.PYESGF` to avoid breaking changes in behavior if `client_type` is missing from the config files.
- Configuration must fail fast if an invalid client string is provided in the YAML.

## Acceptance Criteria (AC)
- [ ] AC 1: `BaseDownloaderConfig` successfully initializes with either `PYESGF` or `ESGPULL` string/enum.
- [ ] AC 2: ValueError or similar validation error is raised when an unsupported client type string is loaded from config.

## Testing & Validation
- **Command**: `make test-custom TEST_ARGS="tests/test_download/test_downloader.py"`
- **Success State**: Configuration parsing tests pass or downloader tests load default successfully.
- **Manual Verification**: Run a quick script to load `BaseDownloaderConfig` with a test yaml setting `client_type` explicitly.

## Completion Protocol
1. [ ] All ACs are met.
2. [ ] Tests pass without regressions.
3. [ ] Code is linted via `make precommit` and `make pylint`.
4. [ ] Documentation updated (if applicable).
5. [ ] Commit work: `git commit -m "feat(download): task 01 - add ClientType configuration"`
6. [ ] Update this document: Mark as COMPLETE.