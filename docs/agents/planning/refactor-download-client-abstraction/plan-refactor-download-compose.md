# 🎯 Scope & Context
The project currently has duplicated classes for downloading data (`Input4MipsDownloader` vs `Input4MipsDownloaderV2` and `CMIP6Downloader` vs `CMIP6DownloaderV2`), where the base versions use `pyesgf` and the `V2` versions use `esgpull`. We need to harmonize these into single `Input4MipsDownloader` and `CMIP6Downloader` classes. The user should instead be able to choose the underlying client implementation (`pyesgf` vs `esgpull`) via configuration, consolidating the shared business logic (e.g., variable iteration, metadata downloading) into a single location.

# 🏛️ Architectural Approach
We will utilize the **Strategy Pattern** and **Dependency Injection** principles. Instead of duplicating the iterator loop logic via inheritance, the core downloaders will switch between backend functions based on a new configuration field. 
- A `ClientType` enumerator (`PYESGF` / `ESGPULL`) will be added to `BaseDownloaderConfig` to formalize the choice.
- **Backend Lifecycle Management:** For `ESGPULL`, the downloader's `download()` method will manage a single `isolated_esgpull_context`. This shared `Esgpull` instance will be injected into all search-and-download calls within the batch to ensure isolation from `$HOME` and performance efficiency (lifting the isolation boundary from per-variable to per-batch).
- The downloader classes will act as orchestrators, reading the configuration and routing the payload parameters to the respective utility functions.

# 🧪 Verification & FMEA
- **Verification Strategy:** Run the existing unit test suites for `test_downloader.py` using `make test`. We will parametrize the tests to run against both `ClientType.PYESGF` and `ClientType.ESGPULL`.
- **FMEA (Failure Mode and Effects Analysis):**
  - **Failure:** User passes an invalid client string in their YAML config. 
    **Mitigation:** `BaseDownloaderConfig` uses Pydantic or strict Enum validation to fail fast.
  - **Failure:** `esg.context` carries stale facet state between variables in `ESGPULL` mode.
    **Mitigation:** Verify search results in TASK-05 to ensure fresh `Query` objects are used per call.
  - **Failure:** Inconsistent function signatures between `utils.py` and `esgpull_utils.py`.
    **Mitigation:** Standardize the orchestrator-to-utility interface to always accept an optional `esg: Esgpull = None` parameter.

# 🪜 Implementation Steps
1. **Add Client Configuration:** Introduce a `ClientType` Enum and add a `client_type: ClientType = ClientType.PYESGF` property to `BaseDownloaderConfig` in `climateset/download/downloader_config.py`. Update YAML parsing to handle this field.
2. **Harmonize `CMIP6Downloader`:** Update `CMIP6Downloader.download_from_model_single_var` to conditionally call `search_and_download_esgf_model_single_var` (for pyesgf) or `esgpull_search_and_download_esgf_model_single_var` (for esgpull) based on `self.config.client_type`. Delete `CMIP6DownloaderV2`.
3. **Harmonize `Input4MipsDownloader`:** Update `Input4MipsDownloader` functions to conditionally dispatch to the `pyesgf` or `esgpull` utility equivalents. Delete `Input4MipsDownloaderV2`.
4. **Refactor CLI / Entrypoints:** Update `cli.py` and any `scripts/download_example.py` to utilize the unified classes, passing the chosen client type via configuration.
5. **Update Tests:** Refactor the test suite in `tests/test_download/test_downloader.py` to instantiate the base downloaders with the appropriate `ClientType` config rather than importing `*V2` classes.

# 🤝 Next Step
Are you ready to approve Step 1 of the implementation to add the `ClientType` configuration?