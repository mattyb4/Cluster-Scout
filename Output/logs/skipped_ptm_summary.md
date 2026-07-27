# Skipped PTM Summary

Regenerated 2026-07-27 against the current `ptm_skipped.tsv` / `download_errors.tsv` (July 23 run). Supersedes the June 8 version — most protein-level findings carry over unchanged, with some new entries reflecting data/threshold drift since then.

68 proteins had PTMs excluded from analysis for non-biological reasons (26 missing/incomplete AlphaFold models, 29 residue mismatches, 5 positions beyond canonical length — some proteins appear in more than one category). Each entry lists the skip reason, PTM count, and what to do if you want to recover this data.

---

## No AFDB Entry (17 proteins)
These proteins have no AlphaFold DB record at all. Almost all are very large proteins likely excluded during AFDB database construction. To recover: generate structures via the AlphaFold Server (https://alphafoldserver.com, limit 5000 AA) or ESMFold, place CIF in `cif_models/{UniProt}/AF-{UniProt}-F1-model_v6.cif`.

| UniProt | Gene | Skipped PTMs | Length | Notes |
|---------|------|-------------|--------|-------|
| Q15149 | PLEC | 328 | 4684 AA | fits AlphaFold Server — new since the June summary |
| P49792 | RANBP2 | 306 | 3224 AA | fits AlphaFold Server |
| Q8WXI7 | MUC16 | 278 | 14507 AA | too large for AlphaFold Server |
| P25054 | APC | 186 | 2843 AA | major tumor suppressor, worth modeling |
| P51587 | BRCA2 | 183 | 3418 AA | major breast/ovarian cancer suppressor, worth modeling |
| Q13315 | ATM | 132 | 3056 AA | major DNA damage checkpoint kinase, worth modeling |
| O14686 | KMT2D | 127 | 5537 AA | too large for AlphaFold Server |
| Q96T58 | SPEN | 119 | 3664 AA | fits AlphaFold Server |
| Q03164 | KMT2A | 110 | 3969 AA | recurrent leukemia fusion gene |
| O95071 | UBR5 | 85 | 2799 AA | fits AlphaFold Server |
| Q14517 | FAT1 | 76 | 4588 AA | fits AlphaFold Server |
| Q8NEZ4 | KMT2C | 76 | 4911 AA | fits AlphaFold Server |
| Q9Y4A5 | TRRAP | 60 | 3859 AA | fits AlphaFold Server |
| Q9NR09 | BIRC6 | 54 | 4857 AA | fits AlphaFold Server |
| Q15911 | ZFHX3 | 38 | 3703 AA | fits AlphaFold Server |
| Q6V0I7 | FAT4 | 13 | 4981 AA | fits AlphaFold Server |
| Q9NZR2 | LRP1B | 8 | 4599 AA | fits AlphaFold Server |

---

## No Canonical AFDB Model (9 proteins)
AFDB has entries for these proteins but only modeled specific isoforms — never the canonical sequence. To recover: same approach as above, model the canonical sequence manually.

| UniProt | Gene | Skipped PTMs | Canonical Length | AFDB isoform models |
|---------|------|-------------|-------------------|---------------------|
| Q8IZT6 | ASPM | 193 | 3477 AA | Q8IZT6-2 (1892 AA) — new since June |
| Q63HN8 | RNF213 | 100 | 5207 AA — too large for AlphaFold Server | Q63HN8-6 (557 AA), Q63HN8-5 (1063 AA) |
| P21359 | NF1 | 73 | 2839 AA | P21359-5 (593), P21359-4 (1598), P21359-3 (551 AA) — major RAS pathway tumor suppressor, worth modeling |
| O15417 | TNRC18 | 50 | 2968 AA | O15417-4 (314), O15417-2 (2256 AA) — new since June |
| Q8NFP9 | NBEA | 20 | 2946 AA | Q8NFP9-3 (739), Q8NFP9-2 (524 AA) |
| Q7Z407 | CSMD3 | 4 | 3707 AA | Q7Z407-5 (2675 AA) |
| Q99102 | MUC4 | 4 | 5412 AA — too large for AlphaFold Server | Q99102-13 (1176), Q99102-12 (1125 AA) |
| Q99996 | AKAP9 | 2 | 3907 AA | Q99996-4 (1643 AA) |
| Q8TDW7 | FAT3 | 1 | 4557 AA | Q8TDW7-2 (1254 AA) |

---

## Residue Mismatch (29 proteins)
The amino acid at the PTM position in the canonical AlphaFold structure does not match the residue in the PTMD annotation. Causes: PTMD annotated against an outdated UniProt sequence version, an alternatively spliced exon in the source isoform, or the entire annotation being from a different isoform. PTMD does not record which sequence version or isoform it used. These PTMs cannot be mapped to the canonical structure without resolving the source sequence.

Ratios reflect (mismatched PTMs) / (total PTMs attempted for that protein) — a high ratio suggests the whole annotation is on the wrong isoform; a low, scattered ratio suggests an isoform-specific region or occasional sequence-version drift.

| UniProt | Gene | Mismatched / Total | Likely cause |
|---------|------|---------------------|--------------|
| Q8NFD5 | ARID1B | 33/35 | entire annotation on wrong isoform |
| Q00987 | MDM2 | 31/35 | mostly wrong isoform, scattered |
| P52948 | NUP98 | 20/26 | likely isoform-specific region |
| P15941 | MUC1 | 13/15 | mostly wrong isoform — new since June |
| Q96RK0 | CIC | 11/68 | scattered, likely isoform-specific exons |
| P78549 | NTHL1 | 10/10 | entire annotation on wrong isoform |
| Q9HB09 | BCL2L12 | 9/9 | entire annotation on wrong isoform |
| Q99836 | MYD88 | 8/13 | moderate mismatch, scattered — new since June (also see length-beyond-canonical below) |
| P38936 | CDKN1A | 7/24 | scattered — also has position_not_in_structure, see below |
| P16455 | MGMT | 5/9 | moderate mismatch — new since June (also see length-beyond-canonical below) |
| Q6NWY9 | PRPF40B | 5/6 | mostly wrong isoform |
| O75030 | MITF | 5/13 | scattered mismatches — isoform-specific |
| O00255 | MEN1 | 5/17 | mixed — likely partially isoform-specific |
| P55197 | MLLT10 | 4/10 | positions 696–709, clustered — alternatively spliced exon |
| Q9UPS6 | SETD1B | 3/14 | scattered — new since June |
| P19544 | WT1 | 3/9 | also has position_not_in_structure; isoform-specific |
| Q9HBE5 | IL21R | 3/3 | entire annotation on wrong isoform |
| O15013 | ARHGEF10 | 3/9 | clustered near C-terminus — alternative C-terminal exon |
| P16220 | CREB1 | 2/2 | entire annotation on wrong isoform |
| Q9NVI1 | FANCI | 2/8 | scattered — new since June |
| P30622 | CLIP1 | 2/23 | spread out — sequence version drift |
| P10275 | AR | 1/41 | K630→R630 — conservative substitution, sequence version update |
| P08575 | PTPRC | 1/13 | S973→R973 — sequence version update |
| P01106 | MYC | 1/13 | T58→P58 — likely sequence version update |
| Q03112 | MECOM | 1/2 | S538→T538 — conservative substitution, sequence version update |
| Q5H9F3 | BCORL1 | 1/25 | near C-terminus — alternative C-terminal exon |
| Q16620 | NTRK2 | 1/23 | Y722→P722 — new since June |
| Q15910 | EZH2 | 1/46 | K27→R27 — conservative substitution, sequence version update |
| Q86Y26 | NUTM1 | 1/2 | near C-terminus — alternative C-terminal exon |

---

## Position Beyond Canonical Length (5 proteins)
PTM positions in PTMD exceed the canonical protein length entirely. These positions do not exist in any known sequence of the protein. Likely PTMD annotation errors or annotations from an undocumented longer isoform.

| UniProt | Gene | Skipped PTMs | Detail |
|---------|------|-------------|--------|
| P38936 | CDKN1A | 10 | Canonical length 164 AA; PTMs annotated at positions 175–197 — positions do not exist |
| P16455 | MGMT | 3 | Canonical length 207 AA; PTMs at positions 209, 224, 232 — just beyond end, new since June |
| P19544 | WT1 | 1 | Canonical length 449 AA; S461 annotated beyond end |
| P31749 | AKT1 | 1 | Canonical length 480 AA; T554 annotated beyond end |
| Q99836 | MYD88 | 1 | Canonical length 296 AA; position 304 annotated beyond end, new since June |

Note: GNAS, which appeared in this category in June (28 skipped PTMs, likely long-isoform XLas), no longer appears in the current skip logs — its AlphaFold model or the underlying gene set has likely changed since then.

---

## Priority Candidates for Manual Recovery

If the group wants to manually obtain structures and recover data, these are the highest-value targets (all fit AlphaFold Server's 5000-residue limit):

1. **PLEC (Q15149)** — 328 PTMs, 4684 AA, newly missing since June
2. **RANBP2 (P49792)** — 306 PTMs, 3224 AA
3. **ASPM (Q8IZT6)** — 193 PTMs, 3477 AA, AFDB has only a 1892 AA isoform model
4. **APC (P25054)** — 186 PTMs, 2843 AA, major cancer gene
5. **BRCA2 (P51587)** — 183 PTMs, 3418 AA, major cancer gene
6. **ATM (Q13315)** — 132 PTMs, 3056 AA, major cancer gene
7. **SPEN (Q96T58)** — 119 PTMs, 3664 AA
8. **KMT2A (Q03164)** — 110 PTMs, 3969 AA, recurrent leukemia gene
9. **NF1 (P21359)** — 73 PTMs, 2839 AA, isoform-only in AFDB
