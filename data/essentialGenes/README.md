# Stanford yeast deletion-project reference gene lists

Two newline-separated SGD ORF lists, extracted from
`code/modelTests/essentialGenes.m`'s hardcoded `inviableORFs` /
`verifiedORFs` local functions so that both the MATLAB and Python
implementations of the essential-gene benchmark can read from a single
source.

## Files

- **`inviable_orfs.txt`** — 1191 entries (1122 unique). Essential ORFs
  from the Stanford yeast deletion project (14 aug 2011 snapshot).
  Original source:
  http://www-sequence.stanford.edu/group/yeast_deletion_project/downloads.html
- **`verified_orfs.txt`** — 5061 unique ORFs. The verified set from SGD
  (27 August 2013).
  Original source:
  http://www.yeastgenome.org/cgi-bin/search/featureSearch?featuretype=ORF&qualifier=Verified

Duplicates in `inviable_orfs.txt` come from the original list and are
preserved verbatim; both languages deduplicate at use time so the
benchmark counts match.

## Why curated, not regenerated each run

The Stanford collection was screened in complex media supplemented
with the auxotrophic markers required by the deletion strains; the
list is an imperfect-but-stable reference used here for comparative
benchmarking of model versions. Keep the file frozen unless deliberately
updating against a newer Stanford snapshot.
