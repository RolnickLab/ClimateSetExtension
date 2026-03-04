# Comprehensive Mapping: `esgf-pyclient` vs. `esgpull`

This document serves as a direct technical translation guide for refactoring `climateset/download/utils.py` and `climateset/download/client.py` to use `esgpull`. 

## 1. High-Level Architecture Mapping

| Concept                | `esgf-pyclient` (Current)                 | `esgpull` (New)                                                                                         |
| :--------------------- | :---------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| **Client Instance**    | `pyesgf.search.SearchConnection(url=...)` | `esgpull.Esgpull(path=data_dir)`                                                                        |
| **Session / State**    | Custom `SearchClient` + `SearchSession`   | `esg.context` (`esgpull.context.Context`)                                                               |
| **Nodes & Failover**   | Manually looping `NODE_LINK_URLS`         | Handled natively by configuring `esg.config.api.index_node` or passing `index_node=` to search methods. |
| **Distributed Search** | `SearchConnection(distrib=True)`          | `esg.context.search(..., query(options=Options(distrib=True)))`                                         |
| **Constraints/Facets** | `ctx.constrain(project="CMIP6", ...)`     | `query = Query(selection=dict(project=["CMIP6"], ...))`                                                 |

## 2. Code-Level Translation Guide

### A. Initialization
**Current `esgf-pyclient`:**
```python
from climateset.download.client import SearchClient

with SearchClient(distrib=False) as client:
    session = client.new_session()
    # ...
```

**New `esgpull` Equivalent:**
```python
from esgpull import Esgpull

# Esgpull requires a path to initialize its SQLite DB and configuration
esg = Esgpull(path=data_dir) 
```

### B. Applying Constraints (Search Definition)
**Current `esgf-pyclient`:**
```python
from climateset.download.constraints import Input4MIPsConstraints

initial_constraints = Input4MIPsConstraints(
    project="input4MIPs", 
    institution_id="VUA", 
    variable="CO2"
)
session.constrain(initial_constraints)
```

**New `esgpull` Equivalent:**
```python
from esgpull.models import Query

# Constraints are passed as lists to allow multiple arguments natively
query = Query(
    selection=dict(
        project=["input4MIPs"],
        institution_id=["VUA"],
        variable=["CO2"] # Maps to 'variable_id' in newer schemas, verify ESGF node schema mapping
    )
)
```

### C. Discovering Facets dynamically
**Current `esgf-pyclient`:**
```python
# To check what grid labels exist before applying one
grid_labels = session.get_available_facets("grid_label") 
if default_grid_label in grid_labels:
    session.constrain(grid_label=default_grid_label)
```

**New `esgpull` Equivalent:**
```python
# Use the context's hints method to fetch available facets
hints = esg.context.hints(query, file=False, facets=["grid_label"])

# hints usually returns a list of dictionaries mapping node -> facets -> counts
if hints and "grid_label" in hints[0]:
    available_grids = list(hints[0]["grid_label"].keys())
    # Apply logic...
```
*(Note: Because `esgpull` supports wildcards and native multi-value arguments like `grid_label=["gn", "gr"]`, dynamically checking facets before querying is often no longer strictly necessary, greatly simplifying `utils.py`)*

### D. Executing the Search
**Current `esgf-pyclient`:**
```python
# Executes against the server and returns a ResultSet of THREDDS pointers
results = session.search() 
```

**New `esgpull` Equivalent:**
```python
# Executes the Solr query and returns list of `File` objects directly mapped to SQLite
files = esg.context.search(query, file=True) 

print(f"Found {len(files)} files.")
```

### E. Executing the Download
**Current `esgf-pyclient`:**
```python
# Iterates over ResultSet, extracts wget bash script, and runs subprocess
for result in search_results:
    wget_script_content = result.file_context().get_download_script()
    subprocess.run(["bash", "-c", wget_script_content, "download", "-s"], cwd=temp_download_path)
```

**New `esgpull` Equivalent:**
```python
import asyncio

# esgpull completely replaces THREDDS bash scripts with native asynchronous fetching.
# Note: files must be added to the internal queue/database before downloading.
esg.db.add(*files)

# Run the async download concurrently
async def run_downloads():
    # Will download tracked files concurrently into `esg.path / data` 
    downloaded, errors = await esg.download(files, show_progress=False)
    return downloaded

# Since utils.py is synchronous right now, we must wrap it:
downloaded_files = asyncio.run(run_downloads())
```

## 3. Notable Edge Cases

1. `esgpull` strictly requires an installation path where it generates a local SQLite directory (`.esgpull/`). The path in `climateset/download/utils.py` operations (`data_dir`) must be initialized properly so `esg.db` does not throw an error.
2. The current implementation creates folders natively like `RAW_DATA / "{project}/{model_id}/{variable}"`. `esgpull` naturally handles standard DRS (Data Reference Syntax) directory generation inside its configuration, but if we need custom `temp_download_path` mapping, we will need to intercept the downloaded files and physically execute `shutil.move()` or configure `esgpull`'s internal data tree path structures.
3. Because `esgpull` `search()` returns `File` classes rather than PyESGF `ResultSet` objects, existing unit tests evaluating `climateset` modules will require new mocks.
