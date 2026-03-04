# SPEC: Add ESGF Download Client using `esgpull`

## 1. Goal Description
The objective is to implement a new asynchronous, stateful `esgpull` downloader client alongside the existing `esgf-pyclient` implementation for ESGF download operations within `climateset`. 

This addition introduces an architecture paradigm shift: providing a database-backed approach driven by `asyncio` as an alternative to the existing dynamic, stateless HTTP search scripts operating over `subprocess`. As the **Orchestrator**, the goal is a contract-first implementation integrating asynchronous Python fetching with strict SLURM cluster isolation requirements, while maintaining backwards compatibility with the existing client.

## 2. Constraints & Assumptions
### Non-Functional Requirements
1. **Isolated Contexts**: `esgpull` must not initialize in `$HOME`. It must initialize in a unique, isolated path natively inside `RAW_DATA` (e.g. `RAW_DATA/.esgpull_jobs/<hash>`) to avoid file lock SQLite collisions between parallel SLURM batch jobs.
2. **Storage Cleanup**: `esgpull` downloads files natively to its internal cache. We must `shutil.move()` them to our strict local directory schema and cleanly delete the `.esgpull_jobs/<hash>` context immediately afterward to avoid cluttering disk space.
3. **Deterministic Search Testing**: Standard unit testing uses heavy mocking. **CRITICAL:** the *search querying phase* mechanism of this application must NEVER be mocked to ensure constraints map accurately against real ESGF nodes. 

## 3. Acceptance Criteria (AC)
- **AC1:** A new `EsgpullDownloader` (or similar) is added as an alternative search and download engine, co-existing with the current `esgf-pyclient` based downloader.
- **AC2:** Existing download clients and their related helper scripts/wrappers are preserved without breaking changes.
- **AC3:** Complex combinations of multi-value constraints and wildcards route through `esgpull.models.Query` successfully in the new client.
- **AC4:** The `esgpull` execution context initializes in an isolated cluster hash folder within `RAW_DATA`, completely bypassing `~/.esgpull`.
- **AC5:** A finalization block safely transfers `.nc` files from the isolation folder to the target Data Reference Syntax and purges the isolation folder afterward.
- **AC6:** End-to-end extraction in the new client is native `asyncio.run(esg.download())` without invoking `subprocess.run(["bash", ...])`.

---

## 4. Orchestrator Implementation Phases

### Blueprint Phase: Environment & Interfaces
* **Delegation**: `python`, `systemdesign`
* **Implement**: 
  1. Modify `pyproject.toml` or `environment.yml`: Add `esgpull` to dependencies while retaining `esgf-pyclient`.
  2. Map out `climateset/download/constraints.py`. Ensure serialization outputs can be seamlessly adapted into dictionaries compatible with `esgpull.models.Query(selection=...)` for the new client, without breaking existing serialization.

### Foundation Phase: The Isolated Context Engine
* **Delegation**: `systemdesign`
* **Implement**: 
  1. Construct a context manager/helper inside `climateset/download/utils.py` (or a new `esgpull_utils.py`) named `isolated_esgpull_context(raw_data_path)`.
  2. Implement uniquely hashed path logic mapping to `RAW_DATA/.esgpull_jobs/<UUID>`.
  3. Ensure a strict `try/finally` block wraps `shutil.rmtree` to tear down the environment safely regardless of download success/failure.

### Implementation Phase A: The Search Contract
* **Delegation**: `python`
* **Implement**: 
  1. Create a new downloader module (e.g., `climateset/download/esgpull_downloader.py`) that inherits from `AbstractDownloader` (if applicable) or implements the necessary download interface.
  2. Implement dynamic facet lookups (e.g., `get_grid_label`) in the new client using `hints = esg.context.hints(query)`.
  3. Implement `search_and_download_*` equivalents that instantiate `esg = Esgpull(path=hash_path)`.
  4. Trigger searches to return native tracking files: `files = esg.context.search(query, file=True)`.

### Implementation Phase B: Asynchronous Execution & Integration
* **Delegation**: `python`
* **Implement**: 
  1. Ensure the new client does not rely on bash-generation scripts like `_download_result` and `_download_process`.
  2. Inside the new download methods, add tracked files: `esg.db.add(*files)`.
  3. Wrap the retrieval execution: `downloaded, errors = asyncio.run(esg.download(files))`.
  4. Implement the file-moving pipeline transferring `.nc` chunk artifacts from the isolated DB cache to the final `RAW_DATA` path.

---

## 5. Verification Plan

Every Acceptance Criterion must be verified natively.

### Automated Verification
- **Test Alignment (AC1, AC2, AC3)**: Add new tests specifically for the `esgpull` client. Ensure existing tests for the `esgf-pyclient` downloader continue to pass unmodified.
- **Search Isolation Guarantee (AC3)**: Review test logs to confirm queries touch real node indices and that actual `File` arrays are returned for the new client.
- **Download Mocking**: Ensure the `Esgpull.download` asynchronous function is the **only** layer intercepted by `unittest.mock` to prevent massive data bandwidth usage in CI.

### Manual End-To-End Verification
- **Storage Independence (AC4, AC5)**: Execute a download script configured to use the new `esgpull` client. During runtime, verify (via terminal `ls`) that `RAW_DATA/.esgpull_jobs/` spawns the correct tracked directory, and verify it deletes cleanly when the script finishes. Confirm `~/.esgpull` is untouched.
- **Subprocess Elimination (AC6)**: Inspect the output of the `esgpull` download script to ensure Python logs asynchronously track chunks rather than arbitrary `wget` process logs. Confirm the `.nc` files correctly structure themselves in the output directory.
- **Regression Testing**: Execute a download using the existing `esgf-pyclient` implementation to verify it still functions as expected.
