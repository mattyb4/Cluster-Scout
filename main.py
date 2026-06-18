from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
from pipeline_utils import PTM_PROXIMITY_STEPS, MUTATION_CLUSTERING_STEPS  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
INPUT_TSV = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
MODELS_DIR = PROJECT_ROOT / "cif_models"

RUN_ONLY_UNIPROT: str | None = None


def _bar(char: str = "─", width: int = 60) -> str:
    return char * width


def run_step(label: str, step_num: int, total: int, cmd: list[str]) -> float:
    print()
    print(_bar("─"))
    print(f"  Step {step_num}/{total}: {label}")
    print(f"  Command: {' '.join(cmd)}")
    print(_bar("─"))
    t0 = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - t0
    print(_bar("─"))
    if result.returncode == 0:
        print(f"  ✓  Step {step_num} finished in {elapsed:.1f}s")
    else:
        print(f"  ✗  Step {step_num} FAILED (exit code {result.returncode}) after {elapsed:.1f}s")
        print(_bar("═"))
        sys.exit(result.returncode)
    print(_bar("─"))
    return elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the mutation proximity pipeline.")
    parser.add_argument(
        "--mode",
        choices=["ptm-proximity", "mutation-clustering"],
        default="ptm-proximity",
        help=(
            "Pipeline mode. 'ptm-proximity' (default) finds recurrent cancer mutations "
            "that cluster in 3D space near disease-associated PTM sites. "
            "'mutation-clustering' finds recurrent mutations that cluster together in 3D "
            "space without any PTM requirement."
        ),
    )
    args = parser.parse_args()
    mode = args.mode

    STEPS = PTM_PROXIMITY_STEPS if mode == "ptm-proximity" else MUTATION_CLUSTERING_STEPS

    print()
    print(_bar("═"))
    print("       Bio465 Capstone Pipeline")
    print(_bar("═"))
    print(f"  Mode         : {mode}")
    print(f"  Project root : {PROJECT_ROOT}")
    print(f"  Input TSV    : {INPUT_TSV}")
    print(f"  Models dir   : {MODELS_DIR}")
    if RUN_ONLY_UNIPROT:
        print(f"  Limiting step 3 to UniProt: {RUN_ONLY_UNIPROT}")
    print()
    print("  Steps to run:")
    for i, label in enumerate(STEPS, 1):
        print(f"    {i}. {label}")
    print(_bar("═"))

    python_exe = sys.executable

    step1_cmd = [python_exe, str(SCRIPTS_DIR / "1_filter.py"), "--mode", mode]

    step2_cmd = [
        python_exe,
        str(SCRIPTS_DIR / "2_download_structures.py"),
        str(INPUT_TSV),
        "--id_column", "uniprot_id",
        "--out_dir", str(MODELS_DIR),
        "--prefer", "cif",
        "--delay", "0.1",
        "--also_pae",
        "--logs_dir", str(PROJECT_ROOT / "Output" / "logs"),
    ]

    step3_cmd = [python_exe, str(SCRIPTS_DIR / "3_find_nearby_mutations.py"), "--mode", mode]
    if RUN_ONLY_UNIPROT:
        step3_cmd.extend(["--uniprot", RUN_ONLY_UNIPROT])

    pipeline_start = time.time()

    t1 = run_step(STEPS[0], 1, len(STEPS), step1_cmd)
    t2 = run_step(STEPS[1], 2, len(STEPS), step2_cmd)
    t3 = run_step(STEPS[2], 3, len(STEPS), step3_cmd)

    step4_cmd = [python_exe, str(SCRIPTS_DIR / "4_annotate_1433pred.py")]
    step5_cmd = [python_exe, str(SCRIPTS_DIR / "5_annotate_polyphen.py")]

    if mode == "ptm-proximity":
        t4 = run_step(STEPS[3], 4, len(STEPS), step4_cmd)
        t5 = run_step(STEPS[4], 5, len(STEPS), step5_cmd)

    total = time.time() - pipeline_start
    print()
    print(_bar("═"))
    print("  Pipeline complete!")
    print()
    print(f"  Step 1 ({STEPS[0]}): {t1:.1f}s")
    print(f"  Step 2 ({STEPS[1]}): {t2:.1f}s")
    print(f"  Step 3 ({STEPS[2]}): {t3:.1f}s")
    if mode == "ptm-proximity":
        print(f"  Step 4 ({STEPS[3]}): {t4:.1f}s")
        print(f"  Step 5 ({STEPS[4]}): {t5:.1f}s")
    print(f"  Total elapsed: {total:.1f}s")
    print(_bar("═"))
    print()


if __name__ == "__main__":
    main()
