# Task 2: Isolated Context Engine Foundation

**Status:** Completed
**Sequence:** 2

## Goal
Build a safe, isolated execution environment for `esgpull` to avoid file lock SQLite collisions during parallel SLURM batch jobs and prevent pollution of the user's `$HOME` directory.

## Sub-tasks
- [x] Construct a context manager/helper named `isolated_esgpull_context(raw_data_path)` inside `climateset/download/utils.py` (or a newly created `esgpull_utils.py`).
- [x] Implement path logic to create a uniquely hashed path mapping to `RAW_DATA/.esgpull_jobs/<UUID>`.
- [x] Ensure that `esgpull` initializes in this unique path (`esg = Esgpull(path=hash_path)`), which will generate its local SQLite directory and configuration, bypassing `~/.esgpull`.
- [x] Implement a strict `try/finally` block that wraps `shutil.rmtree` to tear down and safely purge the isolation folder and its SQLite DB, regardless of whether the download succeeds or fails.

## Constraints & Assumptions
- **Isolated Contexts**: `esgpull` must not initialize in `$HOME`. It must initialize in a unique, isolated path natively inside `RAW_DATA` (e.g. `RAW_DATA/.esgpull_jobs/<hash>`) to avoid file lock SQLite collisions between parallel SLURM batch jobs.
- `esgpull` strictly requires an installation path where it generates a local SQLite directory. The context must manage this lifecycle safely.

## Acceptance Criteria
- **AC4:** The `esgpull` execution context initializes in an isolated cluster hash folder within `RAW_DATA`, completely bypassing `~/.esgpull`.

## Notes
- *Delegation:* systemdesign
- Requires Task 1 to be completed.
