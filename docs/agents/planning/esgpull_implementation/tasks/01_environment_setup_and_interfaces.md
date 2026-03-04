# Task 1: Environment Setup & Interfaces Blueprint

**Status:** Completed
**Sequence:** 1

## Goal
Update project dependencies to include `esgpull` and prepare the existing constraint definitions for compatibility with `esgpull`'s querying system, without breaking the existing `esgf-pyclient` setup. This is part of the overarching goal to implement a new asynchronous, stateful `esgpull` downloader client alongside the existing implementation.

## Sub-tasks
- [x] Update `pyproject.toml` or `environment.yml`: Add `esgpull` to the project's dependencies while retaining `esgf-pyclient`.
- [x] Review and map out `climateset/download/constraints.py`.
- [x] Ensure serialization outputs from constraints can be seamlessly adapted into dictionaries compatible with `esgpull.models.Query(selection=...)`. Handle the transition from strict Pydantic models to `esgpull`'s native multi-value lists (e.g., `project=["CMIP6"]`).
- [x] Verify that these modifications do not break the existing constraint serialization logic used by `esgf-pyclient`.

## Constraints & Assumptions
- Existing code and configurations must be preserved for backward compatibility.
- Transition from iterative constraint building to bulk defining searches where possible, leveraging `esgpull`'s native multi-value, exclusion (`!`), and wildcard (`*`) capabilities.

## Acceptance Criteria
- **AC2:** Existing download clients and their related helper scripts/wrappers are preserved without breaking changes.

## Notes
- *Delegation:* python, systemdesign
