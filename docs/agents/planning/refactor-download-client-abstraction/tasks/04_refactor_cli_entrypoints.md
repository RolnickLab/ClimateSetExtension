# TASK-04: Refactor CLI / Entrypoints

## Goal
Update the CLI and example scripts to use the unified classes (`Input4MipsDownloader` and `CMIP6Downloader`) and remove references to the `*V2` classes.

## Context & References
- **Source Plan**: docs/agents/planning/refactor-download-compose.md
- **Relevant Specs**: N/A
- **Existing Code**: 
  - `climateset/cli.py`
  - `scripts/download_example.py`
  - Any config yaml references if applicable.

## Subtasks
1. [ ] Replace imports of `CMIP6DownloaderV2` with `CMIP6Downloader` in `climateset/cli.py` and `scripts/download_example.py`.
2. [ ] Replace imports of `Input4MipsDownloaderV2` with `Input4MipsDownloader` in `climateset/cli.py` and `scripts/download_example.py`.
3. [ ] Update initialization logic in these entrypoints if `client_type` needs to be explicitly passed or parsed from CLI arguments (if desired), or rely on config parsing doing it.

## Requirements & Constraints
- The user-facing CLI behavior should stay identical, but internally it routes through the unified class rather than picking between V1/V2 classes at the CLI layer.

## Acceptance Criteria (AC)
- [ ] AC 1: No references to `V2` downloaders exist in the entrypoint scripts.
- [ ] AC 2: CLI scripts run successfully using the new consolidated class.

## Testing & Validation
- **Command**: `python scripts/download_example.py --help` (or equivalent test execution)
- **Success State**: Script parses without import errors.
- **Manual Verification**: Run `grep -r "DownloaderV2"` to ensure all usages are purged.

## Completion Protocol
1. [ ] All ACs are met.
2. [ ] Tests pass without regressions.
3. [ ] Code is linted via `make precommit` and `make pylint`.
4. [ ] Documentation updated (if applicable).
5. [ ] Commit work: `git commit -m "refactor(cli): task 04 - replace V2 downloader references in CLI"`
6. [ ] Update this document: Mark as COMPLETE.