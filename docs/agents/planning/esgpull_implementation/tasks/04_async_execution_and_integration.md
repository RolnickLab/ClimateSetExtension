# Task 4: Asynchronous Execution & Integration

**Status:** Completed
**Sequence:** 4

## Goal
Execute the asynchronous download of files tracked in the isolated `esgpull` database and move the final artifacts into the project's standard Data Reference Syntax schema. This replaces the `subprocess.run` bash scripts with native Python `asyncio`.

## Sub-tasks
- [x] Ensure the new `esgpull` client does NOT rely on bash-generation scripts (like `_download_result` and `_download_process`).
- [x] Implement the two-step asynchronous tracker model: First, add the tracked search result files to the local database: `esg.db.add(*files)`.
- [x] Wrap the retrieval execution using `asyncio`: `downloaded, errors = asyncio.run(esg.download(files, show_progress=False))`.
- [x] Implement the finalization pipeline: `esgpull` naturally handles standard DRS. We must physically execute `shutil.move()` to safely transfer `.nc` chunk artifacts from the isolated DB cache DRS tree to the final target `RAW_DATA` path matching the project's specific directory schema.
- [x] Ensure extraction is fully native Python (`asyncio`) without invoking `subprocess.run(["bash", ...])`.

## Constraints & Assumptions
- **Storage Cleanup**: `esgpull` downloads files natively to its internal cache. We must `shutil.move()` them to our strict local directory schema and cleanly delete the `.esgpull_jobs/<hash>` context immediately afterward to avoid cluttering disk space.

## Acceptance Criteria
- **AC5:** A finalization block safely transfers `.nc` files from the isolation folder to the target Data Reference Syntax and purges the isolation folder afterward.
- **AC6:** End-to-end extraction in the new client is native `asyncio.run(esg.download())` without invoking `subprocess.run(["bash", ...])`.

## Notes
- *Delegation:* python
- Requires Task 3 to be completed.
