# Formal Design Document: Refactor Per-Variable Esgpull Isolation to Batch-Scoped Contexts

## 1. 🎯 Scope & Context

`isolated_esgpull_context()` in `esgpull_utils.py` currently creates a new UUID directory, SQLite database, and full `Esgpull(install=True)` instance for **every single variable download call**. In the V2 downloaders, the triple-nested loop `models x variables x experiments` means a configuration with 5 models, 10 variables, and 4 experiments spawns **200 isolated contexts** — each paying the full cost of filesystem setup, SQLite initialization, and teardown via `shutil.rmtree`. This document proposes lifting the isolation boundary from per-variable to per-batch (i.e., per `download()` invocation), so a single `Esgpull` instance is reused across all iterations within one downloader run. The constraint from `ESGPULL_CLIENT_IMPLEMENTATION.md` — isolation from `$HOME` and between parallel SLURM jobs — must be preserved.

## 2. 🧠 Architectural Approach (Trade-offs & Strategy)

- **Chosen pattern: Dependency Injection of a shared context.** Instead of each `esgpull_search_and_download_*` function creating its own context internally, the V2 downloader's `download()` method opens **one** `isolated_esgpull_context()` and passes the `Esgpull` instance into each search-and-download call. This applies the **Dependency Inversion Principle (DIP)** — the download functions depend on an abstract `Esgpull` handle, not on how it was constructed.

- **Why not a module-level singleton?** A singleton would violate the SLURM parallel-job isolation requirement documented in `ESGPULL_CLIENT_IMPLEMENTATION.md` (Constraint 1). Two concurrent SLURM jobs sharing a single SQLite file would deadlock. The batch-scoped context preserves one-context-per-OS-process semantics while eliminating the per-variable churn.

- **Why not a connection pool / context cache keyed by `data_dir`?** This adds complexity (LRU eviction, thread-safety, lifecycle management) that is not justified. The current execution model is single-threaded and sequential within a downloader. A pool solves a problem that does not yet exist.

- **Accepted trade-off: SQLite DB grows within a batch.** The `esg.db.add(*files)` calls accumulate rows across the entire batch instead of starting fresh each time. This is acceptable because (a) the DB is still ephemeral and deleted on context exit, (b) file deduplication already happens via `unique_files` in `_download_and_move_files`, and (c) the total row count for a realistic batch (hundreds of files) is trivial for SQLite.

- **Accepted trade-off: blast radius of a mid-batch failure increases.** If the process crashes at variable 150 of 200, the single context's `finally` block cleans up the isolation directory — including any already-moved files' cache remnants. This is identical to the current behavior per-variable, but now the cleanup is deferred. Already-moved `.nc` files in their final destination are unaffected by the teardown, so no data loss occurs.

- **Principle: Easier to Change (ETC).** The refactored functions accept an `Esgpull` instance as a parameter, making them testable with a mock or stub `Esgpull` without needing to patch the context manager. This also unblocks future parallelization — a caller could open N contexts and distribute work across them.

## 3. 🛡️ Verification & Failure Modes (FMEA)

### Test Strategy

- **Unit tests (`test_esgpull_utils.py`):** Verify that refactored `esgpull_search_and_download_*` functions accept an `Esgpull` instance parameter and invoke `esg.context.search`, `esg.db.add`, and `esg.download` on it. Mock only the `Esgpull` object — never the context manager boundary.
- **Integration tests (existing real-node tests):** Run the V2 downloaders against a live ESGF node for a minimal config (1 model, 1 variable, 1 experiment) and confirm files land in the correct directory tree. This validates that a shared context does not corrupt search state across iterations.
- **Context lifecycle test:** Assert that exactly **one** UUID directory is created under `.esgpull_jobs/` during a multi-variable `download()` call, and that it is removed after the call completes (both on success and on exception).
- **Regression:** Existing V1 downloader tests must pass unchanged — V1 code is not touched by this refactor.

### Known Risks

| Failure Mode | Likelihood | Impact | Mitigation |
|---|---|---|---|
| **SQLite row accumulation causes query slowdown** | Low — realistic batches are < 1 000 files | Minor latency increase | Monitor; if observed, add `esg.db` pruning between iterations |
| **`esg.context` carries stale facet state between variables** | Medium — depends on esgpull internals | Wrong files downloaded | Integration test with 2+ variables asserts correct file counts per variable |
| **`asyncio.run()` called multiple times on same event loop** | Low — each `_download_and_move_files` call uses `asyncio.run()` which creates a fresh loop | RuntimeError if loop is already running | Wrap in `asyncio.new_event_loop()` + `loop.run_until_complete()` if needed; verify in unit test |
| **Interrupted batch leaves orphan `.esgpull_jobs/<UUID>` dir** | Same as current | Disk space leak | No change — `finally` block + `shutil.rmtree` already handles this. Document manual cleanup in ops runbook |

## 4. 📋 Granular Implementation Steps

1. **Add `esg: Esgpull` parameter to search-and-download functions** — Modify `esgpull_search_and_download_esgf_raw_single_var`, `esgpull_search_and_download_esgf_biomass_single_var`, and `esgpull_search_and_download_esgf_model_single_var` in `esgpull_utils.py` to accept an `esg: Esgpull` parameter instead of `data_dir: Path | str`. Remove the internal `with isolated_esgpull_context(...)` block from each function body.

2. **Lift context creation into V2 downloader `download()` methods** — In `CMIP6DownloaderV2.download()` and `Input4MipsDownloaderV2.download()`, wrap the entire iteration loop in a single `with isolated_esgpull_context(self.config.data_dir) as esg:` block. Pass `esg` through to each `download_from_model_single_var` / `download_raw_input_single_var` / `download_meta_historic_biomassburning_single_var` call.

3. **Update intermediate methods to thread `esg` through** — Modify `CMIP6DownloaderV2.download_from_model_single_var`, `Input4MipsDownloaderV2.download_raw_input_single_var`, and `Input4MipsDownloaderV2.download_meta_historic_biomassburning_single_var` to accept and forward the `esg` parameter to the corresponding `esgpull_utils` function.

4. **Re-derive `dest_dir` from `data_dir` inside search-and-download functions** — Since `data_dir` is no longer passed directly, either (a) add `data_dir` as a separate parameter alongside `esg`, or (b) derive it from `esg.config.paths` (the esgpull instance knows its root path's parent). Option (a) is simpler and more explicit — prefer it.

5. **Verify `esg.context` statefulness** — Write a focused integration test that calls `esgpull_search_and_download_esgf_model_single_var` twice with different `(variable, experiment)` pairs using the **same** `esg` instance, and asserts that each call returns files matching only its own constraints. This guards against facet bleed between iterations.

6. **Update existing unit tests** — Adjust `test_esgpull_utils.py` and any V2 downloader tests to reflect the new signatures. Add a lifecycle assertion that `isolated_esgpull_context` is entered exactly once per `download()` invocation.

7. **Run full test suite and linting** — Execute `make test`, `make precommit`, `make pylint`, and `make mypy` to confirm no regressions.

## 5. ⏭️ Next Step

> Shall I proceed with Step 1 — adding the `esg: Esgpull` parameter to the three search-and-download functions in `esgpull_utils.py` and removing their internal `isolated_esgpull_context` blocks?
