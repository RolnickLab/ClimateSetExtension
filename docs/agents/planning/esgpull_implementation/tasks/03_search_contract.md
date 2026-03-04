# Task 3: Search Contract Implementation

**Status:** Completed
**Sequence:** 3

## Goal
Implement the new `esgpull` downloader class to handle search querying against ESGF nodes, producing `File` arrays. This replaces the stateless, iterative facet exploration of `esgf-pyclient` with `esgpull`'s bulk, stateful query system.

## Sub-tasks
- [x] Create a new downloader module (e.g., `climateset/download/esgpull_downloader.py`).
- [x] Implement the `EsgpullDownloader` class, ensuring it adheres to the existing downloader interface (e.g., inheriting from `AbstractDownloader`).
- [x] Implement `search_and_download_*` equivalents that instantiate `esg = Esgpull(path=hash_path)`.
- [x] Replace the iterative `get_available_facets`-based fallback logic used by `esgf-pyclient` with bulk queries, utilizing `esgpull`'s wildcard and multi-value support where appropriate. If dynamic lookup is still required, use `hints = esg.context.hints(query)`.
- [x] Implement distributed search handling (`distrib=True/False`) natively through `esgpull` options (`query(options=Options(distrib=True))`).
- [x] Trigger search queries using the new client that return native tracking files: `files = esg.context.search(query, file=True)`.
- [x] Ensure complex combinations of multi-value constraints and wildcards route through `esgpull.models.Query` successfully.

## Constraints & Assumptions
- Constraints must seamlessly translate to `esgpull.models.Query` ensuring real ESGF nodes are queried correctly.
- `esgpull` handles failover and HTTP management internally, bypassing the manual `_rotate_node()` logic used by `esgf-pyclient`.

## Acceptance Criteria
- **AC1:** A new `EsgpullDownloader` (or similar) is added as an alternative search and download engine, co-existing with the current `esgf-pyclient` based downloader.
- **AC3:** Complex combinations of multi-value constraints and wildcards route through `esgpull.models.Query` successfully in the new client.

## Notes
- *Delegation:* python
- Requires Task 2 to be completed.
