# Task 5: Verification & Testing

**Status:** Completed
**Sequence:** 5

## Goal
Verify the robustness, correctness, and isolation of the new `esgpull` client while ensuring backward compatibility with the existing `esgf-pyclient` downloader.

## Sub-tasks
- [x] Add new automated tests specifically for the `esgpull` client.
- [x] Implement new mocks for unit tests, accounting for the fact that `esgpull` `search()` returns `File` classes rather than `PyESGF` `ResultSet` objects.
- [x] Ensure existing tests for the `esgf-pyclient` downloader continue to pass unmodified.
- [x] **Critical:** Ensure the *search querying phase* is NEVER mocked, so queries touch real node indices and actual `File` arrays are returned.
- [x] Implement `unittest.mock` intercept **only** for the `Esgpull.download` asynchronous function to prevent massive data bandwidth usage in CI.
- [x] Perform manual end-to-end verification of storage independence: run a script with the new client and verify via terminal `ls` that `RAW_DATA/.esgpull_jobs/` spawns correctly and cleans up cleanly. Verify `~/.esgpull` is untouched.
- [x] Perform manual verification of subprocess elimination: inspect output logs to ensure Python asynchronously tracks chunks rather than arbitrary `wget` process logs. Confirm `.nc` files are correctly structured in the output directory.
- [x] Perform manual regression testing using the existing `esgf-pyclient` implementation.

## Constraints & Assumptions
- **Deterministic Search Testing**: Standard unit testing uses heavy mocking. **CRITICAL:** the *search querying phase* mechanism of this application must NEVER be mocked to ensure constraints map accurately against real ESGF nodes.

## Acceptance Criteria
- All previous Acceptance Criteria (AC1-AC6) are proven correct by fulfilling this verification plan.

## Notes
- Some tests can be written in parallel with Tasks 3 and 4. Manual verification requires Tasks 1-4 to be fully implemented.
