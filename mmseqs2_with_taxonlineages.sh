#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# mmseqs_sting_taxonomy.sh
#
# Run mmseqs search for a set of query protein sequences against a target MMseqs
# protein DB, export alignments, then add NCBI TaxID + lineage using an
# accession->taxid mapping file and taxonkit.
#
# Requirements (must be in PATH):
#   - mmseqs
#   - taxonkit
#   - awk, cut, sort, join, mkdir, rm, etc.
#
# Defaults:
#   - Intermediate taxonomy files are REMOVED by default.
#   - Use --keep-tax to keep intermediate taxonomy files.
#   - Use --keep-tmp to keep MMseqs tmp directory.
###############################################################################

usage() {
  cat <<'EOF'
Usage:
  mmseqs_sting_taxonomy.sh \
    -p PREFIX \
    -q QUERY_FASTA \
    -T TARGET_MMSEQS_DB \
    -t THREADS \
    -n NUM_ITERATIONS \
    -e EVALUE \
    -m PROT_ACCESSION2TAXID \
    -d TAXONKIT_TAXDUMP_DIR \
    [-s SENSITIVITY] \
    [-M MAX_SEQS] \
    [-o OUTDIR] \
    [--keep-tmp] \
    [--keep-tax]

Required arguments:
  -p  PREFIX                  Prefix for naming outputs (e.g. sting_MM)
  -q  QUERY_FASTA             FASTA file with one or more protein queries
  -T  TARGET_MMSEQS_DB        Target MMseqs DB (e.g. MM_refseq_protein)
  -t  THREADS                 Number of threads for mmseqs
  -n  NUM_ITERATIONS          mmseqs --num-iterations (e.g. 3)
  -e  EVALUE                  mmseqs -e e-value cutoff (e.g. 1e-5)
  -m  PROT_ACCESSION2TAXID    Path to NCBI protein accession->taxid mapping file
                              (e.g. prot.accession2taxid or .gz)
  -d  TAXONKIT_TAXDUMP_DIR    Path to NCBI taxonomy dump directory for taxonkit
                              (contains nodes.dmp, names.dmp, etc.)

Optional arguments:
  -s  SENSITIVITY             mmseqs -s sensitivity (default: 7.5)
  -M  MAX_SEQS                mmseqs --max-seqs (default: 2000)
  -o  OUTDIR                  Output directory (default: .)
  --keep-tmp                  Do not delete MMseqs temporary directory
  --keep-tax                  Keep intermediate taxonomy files (default: delete)

Outputs (in OUTDIR):
  PREFIX_Res.tsv                   mmseqs convertalis output
  PREFIX_Res.with_taxonomy.tsv     final table with added taxid + lineage columns

If --keep-tax is set, these are also kept:
  PREFIX.accessions
  PREFIX.taxid
  PREFIX.accver_taxid.tsv
  PREFIX.taxids
  PREFIX.taxid_lineage.tsv
  PREFIX.accver_taxid_lineage.tsv

Example:
  mmseqs_sting_taxonomy.sh \
    -p sting_MM \
    -q sting.faa \
    -T MM_refseq_protein \
    -t 60 \
    -n 3 \
    -e 1e-5 \
    -m /thirdisk/public_databases/TaxonomyDB/prot.accession2taxid.gz \
    -d /thirdisk/public_databases/TaxonomyDB/ \
    -o results_sting
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

# Defaults
SENSITIVITY="7.5"
MAX_SEQS="2000"
OUTDIR="."
KEEP_TMP=0
KEEP_TAX=0

# Parse args
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p) PREFIX="$2"; shift 2 ;;
    -q) QUERY_FASTA="$2"; shift 2 ;;
    -T) TARGET_DB="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -n) NUM_ITER="$2"; shift 2 ;;
    -e) EVALUE="$2"; shift 2 ;;
    -m) MAPFILE="$2"; shift 2 ;;
    -d) TAXDUMP="$2"; shift 2 ;;
    -s) SENSITIVITY="$2"; shift 2 ;;
    -M) MAX_SEQS="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    --keep-tmp) KEEP_TMP=1; shift 1 ;;
    --keep-tax) KEEP_TAX=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (use -h for help)" ;;
  esac
done

# Validate required args
: "${PREFIX:?Missing -p PREFIX}"
: "${QUERY_FASTA:?Missing -q QUERY_FASTA}"
: "${TARGET_DB:?Missing -T TARGET_MMSEQS_DB}"
: "${THREADS:?Missing -t THREADS}"
: "${NUM_ITER:?Missing -n NUM_ITERATIONS}"
: "${EVALUE:?Missing -e EVALUE}"
: "${MAPFILE:?Missing -m PROT_ACCESSION2TAXID}"
: "${TAXDUMP:?Missing -d TAXONKIT_TAXDUMP_DIR}"

[[ -s "$QUERY_FASTA" ]] || die "Query FASTA not found or empty: $QUERY_FASTA"
[[ -e "$MAPFILE" ]]     || die "Mapping file not found: $MAPFILE"
[[ -d "$TAXDUMP" ]]     || die "Taxdump dir not found: $TAXDUMP"

command -v mmseqs >/dev/null 2>&1 || die "mmseqs not found in PATH"
command -v taxonkit >/dev/null 2>&1 || die "taxonkit not found in PATH"
command -v awk >/dev/null 2>&1 || die "awk not found"
command -v join >/dev/null 2>&1 || die "join not found"
command -v sort >/dev/null 2>&1 || die "sort not found"
command -v cut >/dev/null 2>&1 || die "cut not found"

mkdir -p "$OUTDIR"

# Standardized names
QUERY_DB="${OUTDIR}/${PREFIX}.queryDB"
RES_DB="${OUTDIR}/${PREFIX}.searchRes"
TMPDIR="${OUTDIR}/${PREFIX}.mmseqs_tmp"

RES_TSV="${OUTDIR}/${PREFIX}_Res.tsv"
RES_WITH_TAX="${OUTDIR}/${PREFIX}_Res.with_taxonomy.tsv"

# Intermediate taxonomy files (kept only with --keep-tax)
ACCESSIONS="${OUTDIR}/${PREFIX}.accessions"
TAXID_RAW="${OUTDIR}/${PREFIX}.taxid"
ACCVER_TAXID="${OUTDIR}/${PREFIX}.accver_taxid.tsv"
TAXIDS="${OUTDIR}/${PREFIX}.taxids"
TAXID_LINEAGE="${OUTDIR}/${PREFIX}.taxid_lineage.tsv"
ACCVER_TAXID_LINEAGE="${OUTDIR}/${PREFIX}.accver_taxid_lineage.tsv"

# Decide how to read mapping file (plain or gz)
MAPCAT="cat"
if [[ "$MAPFILE" =~ \.gz$ ]]; then
  command -v zcat >/dev/null 2>&1 || die "zcat not found, but mapping file is .gz"
  MAPCAT="zcat"
fi

echo "[INFO] Prefix:              $PREFIX"
echo "[INFO] Query FASTA:          $QUERY_FASTA"
echo "[INFO] Target MMseqs DB:      $TARGET_DB"
echo "[INFO] Threads:              $THREADS"
echo "[INFO] Iterations:           $NUM_ITER"
echo "[INFO] E-value:              $EVALUE"
echo "[INFO] Sensitivity (-s):     $SENSITIVITY"
echo "[INFO] Max seqs:             $MAX_SEQS"
echo "[INFO] Mapping file:         $MAPFILE"
echo "[INFO] Taxdump dir:          $TAXDUMP"
echo "[INFO] Output dir:           $OUTDIR"
echo "[INFO] Keep tmp:             $KEEP_TMP"
echo "[INFO] Keep tax intermediates:$KEEP_TAX"
echo

###############################################################################
# 1) Prepare query database
###############################################################################
echo "[INFO] Creating MMseqs query DB..."
mmseqs createdb "$QUERY_FASTA" "$QUERY_DB" --dbtype 1 --threads "$THREADS"

###############################################################################
# 2) Run mmseqs search
###############################################################################
echo "[INFO] Running mmseqs search..."
mmseqs search "$QUERY_DB" "$TARGET_DB" "$RES_DB" "$TMPDIR" \
  --threads "$THREADS" \
  --num-iterations "$NUM_ITER" \
  -s "$SENSITIVITY" \
  -e "$EVALUE" \
  --max-seqs "$MAX_SEQS" \
  --rescore-mode 0

###############################################################################
# 3) Export results as TSV
###############################################################################
echo "[INFO] Converting alignments to TSV..."
mmseqs convertalis "$QUERY_DB" "$TARGET_DB" "$RES_DB" "$RES_TSV" \
  --format-output "query,qlen,qstart,qend,alnlen,qcov,fident,target,tlen,tstart,tend,tcov,theader,evalue,bits" \
  --threads "$THREADS"

###############################################################################
# 4) Extract unique target accessions from column 8
###############################################################################
echo "[INFO] Extracting unique target accessions (column 8)..."
cut -f8 "$RES_TSV" | sort -u > "$ACCESSIONS"

###############################################################################
# 5) Filter mapping file to only these accessions (match accession.version in col2)
#    Mapping format (observed):
#      accession <tab> accession.version <tab> taxid <tab> gi
###############################################################################
echo "[INFO] Filtering mapping file to your accessions (streaming; RAM-safe)..."
awk -F'\t' 'NR==FNR{a[$1]=1; next} ($2 in a){print}' \
  "$ACCESSIONS" <($MAPCAT "$MAPFILE") > "$TAXID_RAW"

###############################################################################
# 6) Build accession.version -> taxid 2-column file
###############################################################################
echo "[INFO] Building accession.version -> taxid mapping..."
awk 'BEGIN{FS=OFS="\t"} {print $2,$3}' "$TAXID_RAW" > "$ACCVER_TAXID"

###############################################################################
# 7) Get lineage for unique taxids using taxonkit
###############################################################################
echo "[INFO] Extracting unique taxids..."
cut -f2 "$ACCVER_TAXID" | sort -u > "$TAXIDS"

echo "[INFO] Running taxonkit lineage..."
taxonkit lineage --data-dir "$TAXDUMP" "$TAXIDS" > "$TAXID_LINEAGE"
# Output: taxid <tab> lineage

###############################################################################
# 8) Join accession.version + taxid with lineage -> accession.version taxid lineage
###############################################################################
echo "[INFO] Joining accession->taxid with taxid->lineage..."
sort -k2,2 "$ACCVER_TAXID" > "${ACCVER_TAXID}.sorted"
sort -k1,1 "$TAXID_LINEAGE" > "${TAXID_LINEAGE}.sorted"

join -t $'\t' -1 2 -2 1 \
  "${ACCVER_TAXID}.sorted" \
  "${TAXID_LINEAGE}.sorted" \
| awk 'BEGIN{FS=OFS="\t"} {print $2,$1,$3}' \
> "$ACCVER_TAXID_LINEAGE"

rm -f "${ACCVER_TAXID}.sorted" "${TAXID_LINEAGE}.sorted"

###############################################################################
# 9) Append taxid + lineage to the original result table using column 8 accession
###############################################################################
echo "[INFO] Appending taxonomy columns to results..."
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR { tax[$1]=$2"\t"$3; next }
  {
    key=$8
    if (key in tax) print $0, tax[key]
    else            print $0, "NA", "NA"
  }' "$ACCVER_TAXID_LINEAGE" "$RES_TSV" \
> "$RES_WITH_TAX"

###############################################################################
# Cleanup intermediates (taxonomy)
###############################################################################
if [[ "$KEEP_TAX" -eq 0 ]]; then
  echo "[INFO] Removing intermediate taxonomy files (use --keep-tax to keep them)..."
  rm -f "$ACCESSIONS" "$TAXID_RAW" "$ACCVER_TAXID" "$TAXIDS" "$TAXID_LINEAGE" "$ACCVER_TAXID_LINEAGE"
else
  echo "[INFO] Keeping intermediate taxonomy files (--keep-tax set)."
fi

###############################################################################
# Cleanup tmp
###############################################################################
if [[ "$KEEP_TMP" -eq 0 ]]; then
  echo "[INFO] Cleaning MMseqs tmp directory: $TMPDIR"
  rm -rf "$TMPDIR"
else
  echo "[INFO] Keeping MMseqs tmp directory: $TMPDIR"
fi

echo
echo "[DONE] Results:"
echo "  - Raw mmseqs TSV:          $RES_TSV"
echo "  - Final TSV with taxonomy: $RES_WITH_TAX"
