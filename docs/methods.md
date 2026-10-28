# Methods

This document records the scientific methods behind Cluster-Scout's analysis, in enough detail to describe them in a paper or reproduce them independently. See `README.md` for how to run the pipeline, and `docs/help.md` for the in-app user guide.

---

## PolyPhen-2 predictions

Implemented in `scripts/4_annotate.py` (step 4). Each mutation gets the PolyPhen-2 **HDIV** prediction from dbNSFP, served by the myvariant.info API.

- **Which records:** the dbNSFP records for the mutation's gene whose reference residue and alternate residue match the mutation and which list its position on *any* transcript. dbNSFP gives one position per Ensembl transcript.
- **Which prediction:** across those records and their transcripts, the most severe class (probably damaging > possibly damaging > benign) wins. Within that class, the highest score is reported. Before 2026-09-29, the first score encountered in that class was reported, which depended on the order records came back. Classes were unaffected; some cached scores from before then differ slightly from what a fresh lookup reports.
- **Batching:** myvariant.info is queried gene by gene, for up to 100 positions per request, and each mutation is matched to its records locally. This is equivalent to querying each variant on its own: tested against 1,184 single-variant answers, all agreed once the tie-break above was applied. It's about 15 times faster. A request returning more than 1,000 records (myvariant.info's page size) is split. Requests are sequential and paced, because concurrent requests triggered the service's rate limiting.
- **Caching:** results are cached in `data/cache/polyphen.tsv` by (gene, mutation). A failed request isn't cached, so it's retried on the next run.

---

## ClinVar as a mutation source

Source file: ClinVar's `variant_summary.txt` (one row per variant per genome assembly).

1. **Variants kept.** Rows on the GRCh38 assembly whose `Name` ends in a protein change `(p.Xaa123Yaa)` between two different standard amino acids. Synonymous (`=`), stop (`Ter`), frameshift, extension, indel and selenocysteine/pyrrolysine changes are dropped. The three-letter codes are converted to one-letter notation. For variants spanning several genes (`GeneSymbol` of the form `A;B`), the gene named in the variant's transcript is used.
2. **Classes.** A variant is included when any part of its `ClinicalSignificance` (split on `/`, so `Pathogenic/Likely pathogenic` counts as both) belongs to a selected class: Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign, or Conflicting classifications of pathogenicity. The default is Pathogenic and Likely pathogenic.
3. **Counts.** A mutation's count is the sum of `NumberSubmitters` over every included variant of the same gene and protein change (a variant with 0 submitters counts as 1). This stands in for COSMIC's number of affected samples; the hotspot threshold (default 3) is applied to it, and the per-gene total is the sum over all of the gene's included missense variants before that threshold.
4. **Numbering.** ClinVar numbers a mutation along its RefSeq transcript's protein. Positions are then put on the canonical AlphaFold-modeled sequence by the same procedure as COSMIC's (next section).

ClinVar submitter counts are not independent patients: several submitters can report the same patient, and ClinVar's clinical focus means well-studied genes have more submissions. Counts therefore compare variants within a ClinVar run, not against COSMIC sample counts.

## Mapping COSMIC mutations onto the canonical protein sequence

Implemented in `scripts/cosmic_numbering.py`. It runs in pipeline step 1 (`scripts/1_filter.py`, `apply_canonical_numbering`) for every pipeline mode and PTM source.

### The problem

Cluster-Scout measures distances on AlphaFold DB models, which are built from each protein's **canonical** UniProt sequence. PTM positions from PTMD and PhosphoSitePlus are also given in canonical numbering. COSMIC, however, numbers each mutation along the transcript it annotated the gene with, and for some genes that transcript encodes a **different isoform**. A mutation COSMIC calls R248Q may then be residue 249 of the canonical protein, or it may lie in an exon the canonical protein doesn't have at all.

Taking COSMIC positions at face value puts these mutations on the wrong residue of the structure. That makes their 3D distances, sequence distances, and PTM-proximity calls wrong.

### Why COSMIC and canonical numbering disagree

Most genes produce several protein isoforms through alternative splicing. COSMIC reports every mutation in a gene along one transcript chosen for that gene. UniProt designates one isoform per protein as canonical, by its own criteria. AlphaFold DB models that canonical sequence, and PTMD and PhosphoSitePlus number PTM sites on it. For most genes the two choices agree. When they don't, COSMIC's numbering follows a different isoform than everything else in the analysis.

In the PTMD-based results (COSMIC Mutant Census v104, minimum 3 samples), this applied to 103 of 712 proteins. COSMIC's isoform was the longer one for 87 of them. That is consistent with COSMIC having historically preferred the longest transcript, but it is an inference from the data, not a documented COSMIC rule. Aligning each COSMIC isoform to the canonical sequence showed these differences:

| Difference from the canonical isoform | Proteins | Examples |
|---|---|---|
| One inserted or skipped segment in the middle (an alternatively spliced exon) | 54 | MEN1, NFIB, SGK1 |
| Two or more internal differences | 21 | PTPRT, BRCA1, FGFR2, TPM4 |
| Longer N-terminus (alternative start or first exon) | 15 | HSP90AA1, TMPRSS2, MYCL |
| Different C-terminus (alternative last exon) | 3 | CDKN2A, RGS7 |
| Other combinations, including a shorter N-terminus | 10 | FGFR1, WT1, ASXL2 |

Some isoforms also differ by a few substituted residues. Positions are counted from the start of the protein, so any inserted or removed segment shifts the numbering of every residue after it. An extra 5-residue exon at position 100, for example, makes COSMIC's R248 the canonical R243. An N-terminal difference shifts the entire protein.

Some mutations can't be placed on the canonical sequence, for two reasons:

- **Isoform-only sequence.** They lie inside a segment only COSMIC's isoform contains, so neither the canonical protein nor its AlphaFold model has that residue. This was the case in 43 of the 103 renumbered proteins.
- **Mismatched reference residues on canonical-numbered proteins.** In 19 proteins, COSMIC's numbering matches the canonical sequence overall, but a few reference residues don't fit. Possible causes are sequence revisions between the releases COSMIC and UniProt used, a few samples annotated on another transcript, or data errors. This group hasn't been investigated individually.

### Why transcript cross-references aren't used

The obvious approach is to look up COSMIC's transcript ID and find the matching UniProt isoform. On real data this failed in two ways:

- **Missing cross-references.** UniProt has no cross-reference to many of COSMIC's Ensembl transcripts. In a 2026-09 test, 209 of the 1,269 affected mutation–PTM pairs were on transcripts UniProt didn't link to any entry.
- **Retired or changed transcripts.** Fetching the transcript's protein directly from Ensembl returned nothing for 32 of 118 transcripts, because Ensembl has retired them. For another 79 mutations, COSMIC's reference residue didn't match Ensembl's current version of the transcript, so the transcript had changed since COSMIC used it.

The method below therefore uses the mutations themselves as evidence, rather than transcript identifiers.

### Procedure

For each protein, the pipeline runs these steps.

**1. Collect candidate sequences.** It fetches every isoform sequence UniProt holds for the protein's canonical accession, including the canonical sequence itself. The query is the UniProt REST `uniprotkb/stream` endpoint with `includeIsoform=true`, 25 accessions per request. Results are cached in `data/cache/uniprot_isoform_sequences.tsv`.

**2. Identify COSMIC's numbering.** Each candidate is scored by how many of the gene's COSMIC hotspot mutations have their reference residue at their stated position in that sequence. For example, R248Q scores a point if residue 248 is R. The highest-scoring isoform is taken as the sequence COSMIC numbered by. The canonical sequence wins ties, so another isoform must explain *strictly more* mutations to be chosen. If the canonical sequence wins, no mutation is moved.

This assumes COSMIC numbers all of a gene's mutations along a single transcript, which matches how the pipeline reads COSMIC: one transcript per gene.

**3. Align the chosen isoform to the canonical sequence.** A global pairwise alignment (Needleman–Wunsch with affine gap penalties, via Biopython's `PairwiseAligner`) is run with these scores:

| Parameter | Value |
|---|---|
| Match | +2 |
| Mismatch | −3 |
| Gap open | −5 |
| Gap extend | −0.1 |

Isoforms of one gene are nearly identical sequences broken up by whole inserted or skipped exons, and sometimes by alternative exons of similar length. Making gaps costly to open but cheap to extend lets an exon-sized difference become one long gap rather than many short ones. Otherwise the alignment would scatter small gaps and mismatches through otherwise identical sequence. The mismatch penalty is higher than the match reward, so substitutions are accepted only where the sequences genuinely differ.

**4. Build a position map.** Every aligned column, whether identical or substituted, gives an isoform position → canonical position pair (both 1-based). Isoform positions that fall in a gap, i.e. sequence only the isoform has, get no canonical position.

**5. Place each mutation.** For a mutation with reference residue *X* at COSMIC position *p*:

- **Mapped position checks out:** if *p* maps to canonical position *q* and the canonical residue at *q* is *X*, the mutation is **moved** to *q*.
- **Original position already fits:** otherwise, if the canonical residue at *p* itself is *X*, the mutation **keeps position *p***. The residue check is treated as stronger evidence than the alignment. In testing this applied to 2 of 413 mutations that already matched the canonical sequence.
- **Neither:** the mutation is **unmapped**. It lies in isoform-only sequence, or no isoform explains its reference residue.
- **Poor-fit protein:** checked before everything else, per protein. If even the best-fitting isoform leaves more than 1% of the protein's COSMIC reference residues unmatched, COSMIC numbered on a transcript UniProt doesn't hold. Every mutation of that protein is **unmapped** and keeps COSMIC's label, whether or not its own residue happens to match (see "Validation of the isoform choice" below).
- **Ambiguous placement:** checked next. If equally scoring optimal alignments (up to 256 are compared) put the mutation's isoform position on different canonical positions, or in a gap in some of them, the mutation is **unmapped** and also recorded as ambiguous. This happens where a gap can slide because the residues at its edges match, or where a repeat lets the gap fall in a different unit at the same score. The residue check cannot catch these: the wrong position holds the same residue as the right one.

The residue check in the first two cases makes a wrong move unlikely. A mutation only moves if the residue it claims to mutate is actually at its new position.

### How results are recorded

- **Labels.** Moved mutations are relabelled in canonical numbering (R248Q → R249Q) everywhere in the outputs. That's the numbering the PTM sites, structures, sequence distances and plots all use.
- **COSMIC's own label** is kept in the step-1 file's `cosmic_mutation_labels` column (as `R249Q=R248Q`). The long-format outputs carry it in `cosmic_mutation` (and `cosmic_anchor_mutation` for Mutation Clustering anchors), and the app shows it as "COSMIC label". Use it to look a mutation up in COSMIC. The same isoform check runs on ClinVar mutations (ClinVar also numbers a variant along its own RefSeq transcript's protein, which can be a different isoform), and for ClinVar results the app shows this column as "ClinVar label"; the column names in the data files stay `cosmic_mutation` and `cosmic_anchor_mutation`.
- **Which isoform COSMIC used** for each renumbered protein is recorded in the step-1 file's `cosmic_numbering_isoform` column.
- **Unmapped mutations** keep COSMIC's label and are listed in the step-1 file's `unmapped_mutations` column. Step 3 still measures them at COSMIC's position number but tags them `(isoform?)`. Step 3 also tags any mutation whose reference residue doesn't match the AlphaFold model's residue at its position, which catches differences between the model and the current UniProt sequence.

### Validation

These figures come from the PTMD-based PTM Proximity results (COSMIC Mutant Census v104, minimum 3 samples), compared before and after the change:

| | Before | After |
|---|---|---|
| Mutation–PTM pairs tagged `(isoform?)` | 1,269 | 176 |
| Proteins found to be numbered by a non-canonical isoform | — | 103 |

- **Wrong pairs removed:** 1,074 pairs present before are gone. They had been measured from the wrong residue.
- **Correct pairs added:** 1,031 pairs appear at the correct canonical positions.
- **Unchanged:** the remaining 12,955 pairs are identical.
- **Fit of the chosen isoform:** among the proteins with tagged mutations, the chosen isoform matched *every* COSMIC reference residue for 72 of 88.

Example: SGK1, where COSMIC's numbering follows a longer isoform. COSMIC's C144W is canonical C49W, 3.9 Å from the K50 PTM site. Before this method it was measured as if it were residue 144.

### Validation of the alignment against UniProt's isoform annotation

`scripts/validate_numbering.py` builds an alignment-free ground truth. UniProt records each isoform as a list of "Alternative sequence" (VAR_SEQ) edits with exact canonical coordinates. Applying them to the canonical sequence gives an isoform → canonical position map. Positions inside a replaced block have no unambiguous truth and are excluded. Run on the 1,504 isoforms of the 527 proteins in the isoform cache (every isoform rebuilt exactly from its edits):

| | Positions |
|---|---|
| With unambiguous truth | 1,183,161 |
| Alignment agrees | 1,182,834 (99.9724%) |
| Alignment disagrees | 327, in 189 isoforms |
| ...placed on a wrong canonical position | 280 |
| ...left unmapped (safe) | 47 |

All 280 wrong placements land on a canonical residue identical to the true one, so the reference-residue check would not catch them, and they sit in repeat-like regions (the largest cluster is Q9NNX6). The ambiguity check above flags all 280, and 586 correctly placed positions as well. In every one of the 189 isoforms the aligner's score exceeds that of UniProt's layout, so these are consequences of the scoring model, not implementation errors, and tuning the scores would not remove them.

Applied to the COSMIC mutations in proteins renumbered by a non-canonical isoform, the check withholds 38 of 39,473 mutations (0.1%) in the full COSMIC table, and 2 of 4,432 in the PTMD-based table.

What this does not test: whether the *chosen isoform* is the one COSMIC used. That is supported by the reference-residue fit above and needs a spot-check against COSMIC's own transcript annotation.

### Validation of the isoform choice against COSMIC's own transcript

Every gene in the Mutant Census v104 file is numbered on a single Ensembl transcript (756 of 756 genes), which confirms the one-numbering-per-gene assumption. `scripts/validate_isoform_choice.py` fetches each gene's Ensembl protein (current version only) and renumbers the same mutations with it as the sole candidate, then compares the labels with the pipeline's. Of 745 proteins, 681 could be compared; 51 had no Ensembl or UniProt sequence and 13 had an Ensembl sequence that no longer fits COSMIC's residues (version drift), so 64 (8.6%) are unvalidated.

| | |
|---|---|
| Proteins where the chosen sequence is identical to COSMIC's protein | 660 of 681, with **no** label disagreements |
| Proteins where COSMIC's protein is not any UniProt isoform | 21, of which 15 have mutations placed differently |
| Mutations compared | 257,968 |
| Same label | 250,810 (97.23%) |
| Both unmapped | 2,269 (0.88%) |
| Placed by Ensembl's numbering but unmapped by the pipeline | 4,889 (1.90%) |
| Placed at a different position than Ensembl's numbering | **0** |

When a UniProt isoform matches COSMIC's protein, the method is right every time. The failures are all in the 15 proteins whose COSMIC protein has no UniProt counterpart (for example AFDN, PTPRT, RUNX1T1, GNAS). Before the poor-fit rule, 316 mutations in 12 of these proteins were placed at a different position than Ensembl's numbering. They passed the reference-residue check at both positions, so that check could not catch them. The proteins are identifiable without Ensembl: the chosen numbering explains only 2-92% of their mutations' reference residues (median 18%), against at least 99% for all but 2 of the other 666 proteins (a threshold of 0.99 flags 17 proteins and catches all 15). That is the basis of the poor-fit rule: with it, no mutation is placed differently from Ensembl's numbering. The cost is that mutations Ensembl's sequence could place stay unmapped (up from 2,806 to 4,889 of the 257,968 compared), tagged `(isoform?)` in step 3 rather than excluded.

In the full COSMIC table the rule applies to 27 proteins and 11,743 of 280,845 mutations (4.2%); in the PTMD-based table, 26 proteins and 1,388 of 34,167 mutations (4.1%); in the PSP-based table, 32 proteins and 1,888 of 35,067 (5.4%). Every one of these has a best-fitting isoform explaining at most 93% of its residues, including the 10 not compared against Ensembl (best fit at most 59%), so none looks like a borderline case.

### Limitations

- **One numbering per gene.** The method assumes COSMIC used one isoform for all of a gene's mutations. Mutations from a second transcript would score poorly and be left unmapped rather than moved wrongly.
- **Reference-residue evidence can coincide.** A residue can match by chance, which could favor the wrong isoform when a gene has only one or two hotspot mutations. The tie-break toward canonical numbering and the per-mutation residue check limit the effect, but can't remove it entirely.
- **Sequence versions.** UniProt sequences, including canonical ones, occasionally change between releases, and the isoform cache is not refreshed automatically. AlphaFold DB models may also predate the current UniProt sequence. Step 3's per-residue check against the model flags these cases as `(isoform?)`. Deleting the cache file forces a refetch.
- **Isoform-only sequence has no structure.** Mutations in exons absent from the canonical protein have no position on the AlphaFold DB model. Placing them would need a structure of that isoform, for example from AlphaFold Server, analyzed with the Single Protein mode.
- **Unavailable sequences.** If UniProt fails to return a protein's isoforms after retries, that protein's mutations keep COSMIC's numbering for that run, with a warning in the step-1 log. They're retried on the next run.
