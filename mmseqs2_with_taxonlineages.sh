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
# Logging:
#   - A combined log is written to: OUTDIR/PREFIX.run.log
#   - A summary report is written to: OUTDIR/PREFIX.summary.txt
#
# Defaults:
#   - Intermediate taxonomy files are REMOVED by default.
#   - MMseqs query/result DB artifacts are REMOVED by default.
#   - The intermediate MMseqs TSV (PREFIX_Res.tsv) is REMOVED by default.
#
# Flags:
#   - --keep-tax        keep intermediate taxonomy files
#   - --keep-mmseqs-db  keep MMseqs query/result DB artifacts
#   - --keep-raw        keep raw MMseqs TSV (PREFIX_Res.tsv)
#   - --keep-tmp        keep MMseqs tmp directory
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
    [--keep-tax] \
    [--keep-mmseqs-db] \
    [--keep-raw]

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

Optional flags:
  --keep-tmp                  Do not delete MMseqs temporary directory
  --keep-tax                  Keep intermediate taxonomy files (default: delete)
  --keep-mmseqs-db            Keep MMseqs query/result DB artifacts (default: delete)
  --keep-raw                  Keep raw MMseqs TSV (PREFIX_Res.tsv) (default: delete)

Outputs (in OUTDIR):
  PREFIX_Res.with_taxonomy.tsv     final table with added taxid + lineage columns
  PREFIX.run.log                   combined log (mmseqs + taxonkit + pipeline)
  PREFIX.summary.txt               summary report (hits, taxonomy coverage, species count)

If --keep-raw is set:
  PREFIX_Res.tsv                   raw mmseqs convertalis output

If --keep-tax is set, these are also kept:
  PREFIX.accessions
  PREFIX.taxid
  PREFIX.accver_taxid.tsv
  PREFIX.taxids
  PREFIX.taxid_lineage.tsv
  PREFIX.accver_taxid_lineage.tsv

If --keep-mmseqs-db is set, these are also kept:
  PREFIX.queryDB*                  MMseqs query DB artifacts
  PREFIX.searchRes*                MMseqs result DB artifacts

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
KEEP_MMSEQS_DB=0
KEEP_RAW=0

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
    --keep-mmseqs-db) KEEP_MMSEQS_DB=1; shift 1 ;;
    --keep-raw) KEEP_RAW=1; shift 1 ;;
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

mkdir -p "$OUTDIR"

# Logging files
LOGFILE="${OUTDIR}/${PREFIX}.run.log"
SUMMARY="${OUTDIR}/${PREFIX}.summary.txt"

# Start logging everything (stdout+stderr) to file AND screen
# NOTE: this must happen after OUTDIR exists
exec > >(tee -a "$LOGFILE") 2>&1

echo "============================================================"
echo "[RUN] mmseqs_sting_taxonomy.sh"
echo "[PREFIX] $PREFIX"
echo "[DATE]   $(date -Is)"
echo "============================================================"
echo "[ARGS]  QUERY_FASTA=$QUERY_FASTA"
echo "[ARGS]  TARGET_DB=$TARGET_DB"
echo "[ARGS]  THREADS=$THREADS"
echo "[ARGS]  NUM_ITER=$NUM_ITER"
echo "[ARGS]  EVALUE=$EVALUE"
echo "[ARGS]  SENSITIVITY=$SENSITIVITY"
echo "[ARGS]  MAX_SEQS=$MAX_SEQS"
echo "[ARGS]  MAPFILE=$MAPFILE"
echo "[ARGS]  TAXDUMP=$TAXDUMP"
echo "[ARGS]  OUTDIR=$OUTDIR"
echo "[ARGS]  KEEP_TMP=$KEEP_TMP"
echo "[ARGS]  KEEP_TAX=$KEEP_TAX"
echo "[ARGS]  KEEP_MMSEQS_DB=$KEEP_MMSEQS_DB"
echo "[ARGS]  KEEP_RAW=$KEEP_RAW"
echo "============================================================"
echo

# Standardized names
QUERY_DB="${OUTDIR}/${PREFIX}.queryDB"
RES_DB="${OUTDIR}/${PREFIX}.searchRes"
TMPDIR="${OUTDIR}/${PREFIX}.mmseqs_tmp"

# Raw TSV (intermediate) and final TSV
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

###############################################################################
# PRE-FLIGHT VALIDATION CHECKS
###############################################################################
echo "------------------------------------------------------------"
echo "[CHECKPOINT] Validating inputs and environment..."
echo "------------------------------------------------------------"

# 1) Check required executables
for tool in mmseqs taxonkit awk sort join cut tee; do
  command -v "$tool" >/dev/null 2>&1 || die "Required executable '$tool' not found in PATH."
done
echo "[OK] Required executables found."

# 2) Check query FASTA
[[ -f "$QUERY_FASTA" ]] || die "Query FASTA file not found: $QUERY_FASTA"
[[ -s "$QUERY_FASTA" ]] || die "Query FASTA file is empty: $QUERY_FASTA"
echo "[OK] Query FASTA exists and is non-empty."

# 3) Check target MMseqs DB (basic sanity: .dbtype file must exist)
[[ -f "${TARGET_DB}.dbtype" ]] || die "Target MMseqs DB not found or not formatted: $TARGET_DB (missing ${TARGET_DB}.dbtype)"
echo "[OK] Target MMseqs DB detected."

# 4) Check mapping file exists & readability (plain or gz)
[[ -f "$MAPFILE" ]] || die "Mapping file not found: $MAPFILE"
if [[ "$MAPFILE" =~ \.gz$ ]]; then
  zcat "$MAPFILE" | head -n 1 >/dev/null 2>&1 || die "Mapping file (.gz) appears unreadable: $MAPFILE"
else
  head -n 1 "$MAPFILE" >/dev/null 2>&1 || die "Mapping file appears unreadable: $MAPFILE"
fi
echo "[OK] Mapping file exists and is readable."

# 5) Check taxdump directory structure
[[ -f "$TAXDUMP/nodes.dmp" ]] || die "Taxdump invalid: missing $TAXDUMP/nodes.dmp"
[[ -f "$TAXDUMP/names.dmp" ]] || die "Taxdump invalid: missing $TAXDUMP/names.dmp"
echo "[OK] Taxdump directory validated."

# 6) Validate numeric parameters
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "THREADS must be an integer: $THREADS"
[[ "$NUM_ITER" =~ ^[0-9]+$ ]] || die "NUM_ITERATIONS must be an integer: $NUM_ITER"
echo "[OK] Numeric parameters validated."

echo "------------------------------------------------------------"
echo "[CHECKPOINT PASSED] All inputs validated."
echo "------------------------------------------------------------"
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
# 3) Export results as TSV (intermediate; deleted by default)
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
# SUMMARY REPORT
###############################################################################
echo "[INFO] Building summary report: $SUMMARY"

total_hits=$(wc -l < "$RES_WITH_TAX" | awk '{print $1}')
tax_hits=$(awk -F'\t' '$NF!="NA"{c++} END{print c+0}' "$RES_WITH_TAX")
unique_taxids=$(awk -F'\t' '$NF!="NA"{print $(NF-1)}' "$RES_WITH_TAX" | sort -u | wc -l | awk '{print $1}')
unique_species=$(awk -F'\t' '$NF!="NA"{print $NF}' "$RES_WITH_TAX" | sort -u | wc -l | awk '{print $1}')

{
  echo "============================================================"
  echo "MMseqs taxonomy run summary"
  echo "Prefix:                $PREFIX"
  echo "Date:                  $(date -Is)"
  echo "Query FASTA:            $QUERY_FASTA"
  echo "Target MMseqs DB:       $TARGET_DB"
  echo "Threads:               $THREADS"
  echo "Iterations:            $NUM_ITER"
  echo "E-value:               $EVALUE"
  echo "Sensitivity (-s):      $SENSITIVITY"
  echo "Max seqs:              $MAX_SEQS"
  echo "Mapping file:          $MAPFILE"
  echo "Taxdump dir:           $TAXDUMP"
  echo "Output directory:      $OUTDIR"
  echo "------------------------------------------------------------"
  echo "Total hits (rows):     $total_hits"
  echo "Hits with taxonomy:    $tax_hits"
  echo "Unique TaxIDs:         $unique_taxids"
  echo "Unique species(lineage): $unique_species"
  echo "------------------------------------------------------------"
  echo "Final output:          $RES_WITH_TAX"
  echo "Log file:              $LOGFILE"
  echo "============================================================"
} > "$SUMMARY"

echo "[INFO] Summary:"
cat "$SUMMARY"

###############################################################################
# Cleanup intermediates (taxonomy)
###############################################################################
if [[ "$KEEP_TAX" -eq 0 ]]; then
  echo "[INFO] Removing intermediate taxonomy files (use --keep-tax to keep them)..."
  rm -f "$ACCESSIONS" "$TAXID_RAW" "$ACCVER_TAXID" "$TAXIDS" "$TAXID_LINEAGE" "$ACCVER_TAXID_LINEAGE" 2>/dev/null || true
else
  echo "[INFO] Keeping intermediate taxonomy files (--keep-tax set)."
fi

###############################################################################
# Cleanup raw MMseqs TSV
###############################################################################
if [[ "$KEEP_RAW" -eq 0 ]]; then
  echo "[INFO] Removing raw MMseqs TSV (use --keep-raw to keep it): $RES_TSV"
  rm -f "$RES_TSV" 2>/dev/null || true
else
  echo "[INFO] Keeping raw MMseqs TSV (--keep-raw set): $RES_TSV"
fi

###############################################################################
# Cleanup MMseqs DB artifacts (query DB + search result DB)
###############################################################################
if [[ "$KEEP_MMSEQS_DB" -eq 0 ]]; then
  echo "[INFO] Removing MMseqs DB artifacts (use --keep-mmseqs-db to keep them)..."

  # Query DB artifacts
  rm -f \
    "${QUERY_DB}" \
    "${QUERY_DB}.dbtype" \
    "${QUERY_DB}.index" \
    "${QUERY_DB}.lookup" \
    "${QUERY_DB}.source" \
    "${QUERY_DB}_h" \
    "${QUERY_DB}_h.dbtype" \
    "${QUERY_DB}_h.index" \
    2>/dev/null || true

  # Search results DB artifacts
  rm -f \
    "${RES_DB}" \
    "${RES_DB}.dbtype" \
    "${RES_DB}.index" \
    2>/dev/null || true
else
  echo "[INFO] Keeping MMseqs DB artifacts (--keep-mmseqs-db set)."
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
echo "[DONE] Result:"
echo "  - Final TSV with taxonomy: $RES_WITH_TAX"
echo "  - Log file:                $LOGFILE"
echo "  - Summary report:          $SUMMARY"
