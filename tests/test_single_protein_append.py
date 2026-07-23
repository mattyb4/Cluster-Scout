"""Tests for single-protein append/replace logic in the proximity database."""
import csv
from pathlib import Path

import pandas as pd
import pytest

from conftest import import_script

mod = import_script("analyze_single_cif_nearby_mutations.py")

OUTPUT_COLUMNS = mod.OUTPUT_COLUMNS_STEP3


def _write_db(path: Path, rows: list[dict]) -> None:
    """Write a minimal proximity DB in UTF-16 TSV format."""
    with path.open("w", encoding="utf-16", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_db(path: Path) -> list[dict]:
    """Read back all rows from a UTF-16 proximity DB."""
    df = pd.read_csv(path, sep="\t", encoding="utf-16", dtype=str, keep_default_na=False)
    return df.to_dict("records")


def _make_row(uniprot: str, ptm_site: str, **overrides) -> dict:
    """Build a minimal proximity DB row with defaults for all required columns."""
    row = {col: "" for col in OUTPUT_COLUMNS}
    row["UniProt"] = uniprot
    row["ptm_site"] = ptm_site
    row.update(overrides)
    return row


class TestAppendRowsNoDuplicates:
    def test_appends_to_empty_file(self, tmp_path):
        db = tmp_path / "db.tsv"  # doesn't exist yet -- exercises the "no prior file" path
        rows = [_make_row("P12345", "S100"), _make_row("P12345", "T200")]
        added = mod.append_rows_no_duplicates(db, rows)
        assert added == 2, f"both new rows should be reported as added when the DB starts empty, got {added}"
        result = _read_db(db)
        assert len(result) == 2, f"the DB file should now contain both rows, got {len(result)}"

    def test_skips_exact_duplicates(self, tmp_path):
        db = tmp_path / "db.tsv"
        row = _make_row("P12345", "S100")
        _write_db(db, [row])
        added = mod.append_rows_no_duplicates(db, [row])
        assert added == 0, (
            f"appending a row identical to one already in the DB (same UniProt+ptm_site "
            f"key) must be reported as 0 added, not silently duplicated, got {added}"
        )
        assert len(_read_db(db)) == 1, "the DB should still contain exactly one row after the no-op append"

    def test_appends_new_ptm_for_same_uniprot(self, tmp_path):
        db = tmp_path / "db.tsv"
        existing = _make_row("P12345", "S100")
        _write_db(db, [existing])
        new = _make_row("P12345", "T200")
        added = mod.append_rows_no_duplicates(db, [new])
        assert added == 1, (
            f"a new PTM site for an already-present protein is not a duplicate (dedup key "
            f"is UniProt+ptm_site, not just UniProt) and should be added, got {added}"
        )
        result = _read_db(db)
        assert len(result) == 2, f"the DB should now have both PTM sites for P12345, got {len(result)} row(s)"
        sites = {r["ptm_site"] for r in result}
        assert sites == {"S100", "T200"}, f"both PTM sites should be present, got {sites}"


class TestRemoveAndReplace:
    """Test the replace flow: remove existing rows for a UniProt, then append new ones."""

    def _remove_uniprot_rows(self, db_path: Path, uniprot: str) -> int:
        """Replicate the app's _remove_uniprot_rows logic for testing."""
        df = pd.read_csv(db_path, sep="\t", encoding="utf-16", dtype=str,
                         keep_default_na=False)
        before = len(df)
        df = df[df["UniProt"] != uniprot]
        df.to_csv(db_path, sep="\t", index=False, encoding="utf-16")
        return before - len(df)

    def test_replace_removes_old_and_appends_new(self, tmp_path):
        db = tmp_path / "db.tsv"
        old_rows = [
            _make_row("P12345", "S100", gene="TP53"),
            _make_row("P12345", "T200", gene="TP53"),
            _make_row("Q99999", "S50", gene="OTHER"),
        ]
        _write_db(db, old_rows)
        assert len(_read_db(db)) == 3, "sanity check: all 3 seed rows should be present before the replace"

        removed = self._remove_uniprot_rows(db, "P12345")
        assert removed == 2, f"both of P12345's rows should be removed, got {removed}"
        remaining = _read_db(db)
        assert len(remaining) == 1, f"only Q99999's row should remain after removing P12345, got {len(remaining)}"
        assert remaining[0]["UniProt"] == "Q99999", (
            f"the surviving row must belong to the untouched protein, got {remaining[0]['UniProt']}"
        )

        new_rows = [
            _make_row("P12345", "S100", gene="TP53", mutation_at_ptm_site="yes"),
            _make_row("P12345", "T200", gene="TP53"),
            _make_row("P12345", "Y300", gene="TP53"),
        ]
        added = mod.append_rows_no_duplicates(db, new_rows)
        assert added == 3, f"all 3 fresh rows should be added post-removal (no stale dupes left), got {added}"
        result = _read_db(db)
        assert len(result) == 4, f"expected Q99999's 1 row + P12345's 3 new rows = 4 total, got {len(result)}"
        p12345_rows = [r for r in result if r["UniProt"] == "P12345"]
        assert len(p12345_rows) == 3, f"P12345 should have exactly its 3 re-analyzed PTM sites, got {len(p12345_rows)}"
        sites = {r["ptm_site"] for r in p12345_rows}
        assert sites == {"S100", "T200", "Y300"}, f"expected all 3 new PTM sites, got {sites}"

    def test_replace_with_no_existing_rows_is_harmless(self, tmp_path):
        db = tmp_path / "db.tsv"
        _write_db(db, [_make_row("Q99999", "S50")])
        removed = self._remove_uniprot_rows(db, "P12345")
        assert removed == 0, (
            f"removing a protein that was never in the DB should be a harmless no-op "
            f"(0 removed), not raise a KeyError, got {removed}"
        )
        assert len(_read_db(db)) == 1, "the unrelated existing row must be untouched"

    def test_other_proteins_are_untouched(self, tmp_path):
        db = tmp_path / "db.tsv"
        rows = [
            _make_row("P12345", "S100"),
            _make_row("Q11111", "T50"),
            _make_row("Q22222", "S75"),
        ]
        _write_db(db, rows)
        self._remove_uniprot_rows(db, "P12345")
        result = _read_db(db)
        uniprots = {r["UniProt"] for r in result}
        assert uniprots == {"Q11111", "Q22222"}, (
            f"removing P12345 must leave every other protein's rows exactly as they "
            f"were, got remaining UniProts {uniprots}"
        )
