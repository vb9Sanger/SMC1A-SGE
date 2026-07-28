#!/usr/bin/env bash
# extract_clinvar_vcfs.sh
#
# For each targeton in the regions guide TSV, extracts variants from a
# ClinVar VCF using bcftools and writes a per-targeton bgzipped VCF + index.
#
# Usage:
#   bash extract_clinvar_vcfs.sh \
#       --regions  targeton_regions.tsv \
#       --clinvar  clinvar.vcf.gz \
#       --outdir   clinvar_vcfs/
#
# Requirements: bcftools (>=1.9), tabix (part of htslib)
#
# Input TSV must have headers: Targeton_ID, chrom, start, end
# (as produced by extract_targeton_regions.py)

set -euo pipefail

# ---------- defaults ----------
REGIONS_TSV=""
CLINVAR_VCF=""
OUTDIR="clinvar_vcfs"

# ---------- argument parsing ----------
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --regions)  REGIONS_TSV="$2"; shift 2 ;;
        --clinvar)  CLINVAR_VCF="$2"; shift 2 ;;
        --outdir)   OUTDIR="$2";      shift 2 ;;
        -h|--help)  usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ---------- validation ----------
if [[ -z "$REGIONS_TSV" || -z "$CLINVAR_VCF" ]]; then
    echo "ERROR: --regions and --clinvar are required."
    usage
fi

if [[ ! -f "$REGIONS_TSV" ]]; then
    echo "ERROR: regions file not found: $REGIONS_TSV"; exit 1
fi

if [[ ! -f "$CLINVAR_VCF" ]]; then
    echo "ERROR: ClinVar VCF not found: $CLINVAR_VCF"; exit 1
fi

if [[ ! -f "${CLINVAR_VCF}.tbi" ]]; then
    echo "ERROR: ClinVar VCF index not found: ${CLINVAR_VCF}.tbi"
    echo "       Run: tabix -p vcf ${CLINVAR_VCF}"
    exit 1
fi

for tool in bcftools tabix; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found on PATH"; exit 1
    fi
done

mkdir -p "$OUTDIR"

# ---------- main loop ----------
# Read header line to find column indices (handles any column order)
HEADER=$(grep -m1 "^Targeton_ID" "$REGIONS_TSV")
if [[ -z "$HEADER" ]]; then
    echo "ERROR: could not find header row starting with 'Targeton_ID' in $REGIONS_TSV"
    exit 1
fi

ID_COL=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^Targeton_ID$" | cut -d: -f1)
CHR_COL=$(echo "$HEADER"  | tr '\t' '\n' | grep -n "^chrom$"       | cut -d: -f1)
START_COL=$(echo "$HEADER"| tr '\t' '\n' | grep -n "^start$"       | cut -d: -f1)
END_COL=$(echo "$HEADER"  | tr '\t' '\n' | grep -n "^end$"         | cut -d: -f1)

echo "Column indices — Targeton_ID:$ID_COL  chrom:$CHR_COL  start:$START_COL  end:$END_COL"
echo ""

SUCCESS=0
EMPTY=0
FAIL=0

while IFS=$'\t' read -r -a FIELDS; do
    # Skip header
    [[ "${FIELDS[0]}" == "Targeton_ID" ]] && continue

    TARGETON="${FIELDS[$((ID_COL-1))]}"
    CHROM="${FIELDS[$((CHR_COL-1))]}"
    START="${FIELDS[$((START_COL-1))]}"
    END="${FIELDS[$((END_COL-1))]}"

    if [[ -z "$TARGETON" || -z "$CHROM" || -z "$START" || -z "$END" ]]; then
        echo "WARNING: skipping incomplete row: ${FIELDS[*]}"
        continue
    fi

    # Strip 'chr' prefix if present — ClinVar GRCh38 VCF uses plain contig names (e.g. X, not chrX)
    CHROM_BARE="${CHROM#chr}"
    REGION="${CHROM_BARE}:${START}-${END}"
    OUTFILE="${OUTDIR}/${TARGETON}_clinvar.vcf.gz"

    echo -n "  ${TARGETON}  ${REGION}  →  ${OUTFILE} ... "

    if bcftools view \
            --regions "$REGION" \
            --output-type z \
            --output "$OUTFILE" \
            "$CLINVAR_VCF" 2>/dev/null; then

        # Check if any variants were written (more than just the header)
        N_VARS=$(bcftools view -H "$OUTFILE" 2>/dev/null | wc -l)

        if [[ "$N_VARS" -eq 0 ]]; then
            echo "OK (0 variants)"
            EMPTY=$((EMPTY+1))
        else
            bcftools index --tbi "$OUTFILE"
            echo "OK (${N_VARS} variant(s))"
            SUCCESS=$((SUCCESS+1))
        fi
    else
        echo "FAILED"
        FAIL=$((FAIL+1))
    fi

done < "$REGIONS_TSV"

# ---------- summary ----------
echo ""
echo "Done."
echo "  With variants (indexed): $SUCCESS"
echo "  Empty (no ClinVar hits): $EMPTY"
echo "  Failed:                  $FAIL"

if [[ "$EMPTY" -gt 0 ]]; then
    echo ""
    echo "NOTE: Empty VCFs were not indexed. If you expected variants in those"
    echo "      regions, check that the contig naming matches (e.g. 'chrX' vs 'X')."
    echo "      ClinVar GRCh38 uses 'chr'-prefixed contigs; GRCh37 does not."
fi
