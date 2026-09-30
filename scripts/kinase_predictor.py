# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["kinase-library==1.8.0"]
# ///
"""Kinase Library predictions, run in their own environment.

Every kinase-library release pins numpy ~1.26 / pandas ~2.2, which can't share
an environment with this project's numpy >= 2 / pandas >= 3 -- installing it
into the project environment downgrades them, and the next `uv sync` removes
it again (that silently broke kinase predictions from 2026-07 until this
helper). So step 4 runs this script with `uv run --script`, which builds and
caches a separate environment from the dependency block above (versions
pinned in kinase_predictor.py.lock) without touching the project's.

Protocol: phosphosite windows (15-mers, center residue lowercase) on stdin,
one per line; one JSON object per window on stdout, in order:
    {"window": "...", "prediction": "CDK1(3.50,99.0%); ..."}
    {"window": "...", "error": "ValueError: ..."}
`--check` imports the library and exits, to build/verify the environment.

The Kinase Library (Johnson et al., Nature 2023) is licensed CC BY-NC-SA 3.0.
"""
from __future__ import annotations

import json
import sys
import warnings
from functools import lru_cache

TOP_K = 5


def speed_up_kinase_library() -> None:
    """Memoize kinase_library's reference-data loaders for this process.

    Substrate.predict() re-reads the same kinome/matrix reference files from
    disk on every call instead of caching them; wrapping the two most-called
    loaders in an lru_cache eliminates that redundant I/O. Idempotent.
    """
    import kinase_library.modules.data as kl_data
    if getattr(kl_data, "_cluster_scout_cached", False):
        return
    kl_data.get_kinase_list = lru_cache(maxsize=None)(kl_data.get_kinase_list)
    kl_data.get_kinome_info = lru_cache(maxsize=None)(kl_data.get_kinome_info)
    kl_data._cluster_scout_cached = True


def predict_kinases(window: str) -> str:
    """Run the Kinase Library on a 15-mer window and return a formatted top-5 string."""
    import kinase_library as kl
    speed_up_kinase_library()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sub = kl.Substrate(window)
        result = sub.predict()
    top = result.head(TOP_K)
    parts = []
    for kinase, row in top.iterrows():
        parts.append(f"{kinase}({row['Score']:.2f},{row['Percentile']:.1f}%)")
    return "; ".join(parts)


def main() -> None:
    if "--check" in sys.argv:
        import kinase_library  # noqa: F401 -- building/importing the environment is the check
        speed_up_kinase_library()
        print(json.dumps({"ok": True}), flush=True)
        return
    for line in sys.stdin:
        window = line.strip()
        if not window:
            continue
        try:
            out = {"window": window, "prediction": predict_kinases(window)}
        except Exception as exc:
            out = {"window": window, "error": f"{type(exc).__name__}: {exc}"}
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
