# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["aiupred @ https://github.com/doszilab/AIUPred/archive/c5cd78d6d485c09001964e6a98c08b993907eca1.zip"]
# ///
"""AIUPred disorder / binding predictions, run locally in their own environment.

AIUPred needs PyTorch, a large download that the rest of the project doesn't use,
so step 4 runs this script with `uv run --script`, which builds and caches a
separate environment from the dependency block above (versions pinned in
aiupred_predictor.py.lock) without touching the project's. The model weights
ship inside the aiupred package. It is pinned to a commit and fetched as a zip
so no `git` is needed.

Protocol: one `UNIPROT_ID<TAB>SEQUENCE` per line on stdin; one JSON object per
protein on stdout, in order:
    {"id": "P04637", "general": [0.97, ...], "binding": [0.87, ...]}
    {"id": "P04637", "error": "ValueError: ..."}
Scores are per residue, in sequence order (position 1 first). `--check` loads the
model and exits, to build/verify the environment.

AIUPred (Erdős & Dosztányi, 2024) is MIT licensed:
https://github.com/doszilab/AIUPred
"""
from __future__ import annotations

import json
import sys


def load_model():
    import torch
    from aiupred import AIUPred
    return AIUPred(force_cpu=not torch.cuda.is_available())


def main() -> None:
    model = load_model()
    if "--check" in sys.argv:
        print(json.dumps({"ok": True}), flush=True)
        return
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        uniprot_id, _, sequence = line.partition("\t")
        try:
            out = {
                "id": uniprot_id,
                "general": [float(x) for x in model.predict_disorder(sequence)],
                "binding": [float(x) for x in model.predict_binding(sequence)],
            }
        except Exception as exc:
            out = {"id": uniprot_id, "error": f"{type(exc).__name__}: {exc}"}
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
