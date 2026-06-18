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
        db = tmp_path / "db.tsv"
        rows = [_make_row("P12345", "S100"), _make_row("P12345", "T200")]
        added = mod.append_rows_no_duplicates(db, rows)
        assert added == 2
        result = _read_db(db)
        assert len(result) == 2

    def test_skips_exact_duplicates(self, tmp_path):
        db = tmp_path / "db.tsv"
        row = _make_row("P12345", "S100")
        _write_db(db, [row])
        added = mod.append_rows_no_duplicates(db, [row])
        assert added == 0
        assert len(_read_db(db)) == 1

    def test_appends_new_ptm_for_same_uniprot(self, tmp_path):
        db = tmp_path / "db.tsv"
        existing = _make_row("P12345", "S100")
        _write_db(db, [existing])
        new = _make_row("P12345", "T200")
        added = mod.append_rows_no_duplicates(db, [new])
        assert added == 1
        result = _read_db(db)
        assert len(result) == 2
        sites = {r["ptm_site"] for r in result}
        assert sites == {"S100", "T200"}


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
        assert len(_read_db(db)) == 3

        removed = self._remove_uniprot_rows(db, "P12345")
        assert removed == 2
        remaining = _read_db(db)
        assert len(remaining) == 1
        assert remaining[0]["UniProt"] == "Q99999"

        new_rows = [
            _make_row("P12345", "S100", gene="TP53", mutation_at_ptm_site="yes"),
            _make_row("P12345", "T200", gene="TP53"),
            _make_row("P12345", "Y300", gene="TP53"),
        ]
        added = mod.append_rows_no_duplicates(db, new_rows)
        assert added == 3
        result = _read_db(db)
        assert len(result) == 4
        p12345_rows = [r for r in result if r["UniProt"] == "P12345"]
        assert len(p12345_rows) == 3
        sites = {r["ptm_site"] for r in p12345_rows}
        assert sites == {"S100", "T200", "Y300"}

    def test_replace_with_no_existing_rows_is_harmless(self, tmp_path):
        db = tmp_path / "db.tsv"
        _write_db(db, [_make_row("Q99999", "S50")])
        removed = self._remove_uniprot_rows(db, "P12345")
        assert removed == 0
        assert len(_read_db(db)) == 1

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
        assert uniprots == {"Q11111", "Q22222"}
