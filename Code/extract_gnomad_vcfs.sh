#!/usr/bin/env bash
# extract_gnomad_vcfs.sh
#
# For each targeton in the regions guide TSV, queries gnomAD v4.1 sites VCFs
# REMOTELY (no local download required) using bcftools + HTTPS range requests,
# and writes per-targeton bgzipped VCF + index, for BOTH exomes and genomes.
#
# Mirrors extract_clinvar_vcfs.sh in structure/output conventions, with two
# differences:
#   1. Two source files (exomes + genomes) instead of one, since assay
#      variants can fall in regions with uneven exome-capture coverage
#      (intronic/splice positions in particular).
#   2. gnomAD v4.1 GRCh38 VCFs use 'chr'-prefixed contigs, so — unlike the
#      ClinVar script — no contig stripping is needed; --regions is built
#      directly from the chrom column as given.
#
# Usage:
#   bash extract_gnomad_vcfs.sh \
#       --regions      targeton_regions.tsv \
#       --exomes-url   https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chrX.vcf.bgz \
#       --genomes-url  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrX.vcf.bgz \
#       --outdir       gnomad_vcfs/
#
# Requirements: bcftools (>=1.9, built with libcurl for remote HTTPS access)
#
# Input TSV must have headers: Targeton_ID, chrom, start, end
# (the SAME targeton_regions.tsv used for extract_clinvar_vcfs.sh / the
# ClinVar intersect script — NOT the exon map, which only gives exon
# boundaries and can miss flanking intronic/splice sequence that the assay
# actually covers.)
#
# Note: no local .tbi needed for the remote files — htslib fetches the
# remote index (<url>.tbi) automatically via HTTP range requests, the same
# way it fetches the requested region of the VCF itself. Nothing is
# downloaded in full.

set -euo pipefail

# ---------- defaults ----------
REGIONS_TSV=""
EXOMES_URL=""
GENOMES_URL=""
OUTDIR="gnomad_vcfs"

# ---------- argument parsing ----------
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --regions)      REGIONS_TSV="$2"; shift 2 ;;
        --exomes-url)   EXOMES_URL="$2";  shift 2 ;;
        --genomes-url)  GENOMES_URL="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2";      shift 2 ;;
        -h|--help)      usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ---------- validation ----------
if [[ -z "$REGIONS_TSV" ]]; then
    echo "ERROR: --regions is required."
    usage
fi

if [[ -z "$EXOMES_URL" && -z "$GENOMES_URL" ]]; then
    echo "ERROR: at least one of --exomes-url / --genomes-url is required."
    usage
fi

if [[ ! -f "$REGIONS_TSV" ]]; then
    echo "ERROR: regions file not found: $REGIONS_TSV"; exit 1
fi

for tool in bcftools; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found on PATH"; exit 1
    fi
done

if ! bcftools --version 2>/dev/null | head -1 | grep -q .; then
    echo "ERROR: bcftools did not report a version — check installation"; exit 1
fi

mkdir -p "$OUTDIR"

# ---------- helper: query one source (exomes or genomes) for one region ----------
# Args: label(exomes|genomes) url targeton chrom start end
query_source() {
    local label="$1" url="$2" targeton="$3" chrom="$4" start="$5" end="$6"
    local region="${chrom}:${start}-${end}"
    local outfile="${OUTDIR}/${targeton}_gnomad_${label}.vcf.gz"

    echo -n "  ${targeton}  [${label}]  ${region}  ->  ${outfile} ... "

    if bcftools view \
            --regions "$region" \
            --output-type z \
            --output "$outfile" \
            "$url" 2>/dev/null; then

        local n_vars
        n_vars=$(bcftools view -H "$outfile" 2>/dev/null | wc -l)

        if [[ "$n_vars" -eq 0 ]]; then
            echo "OK (0 variants)"
            EMPTY=$((EMPTY+1))
        else
            bcftools index --tbi "$outfile"
            echo "OK (${n_vars} variant(s))"
            SUCCESS=$((SUCCESS+1))
        fi
    else
        echo "FAILED"
        echo "    (remote query failed — check network access, and that the URL/contig" \
             "naming e.g. 'chrX' is correct for this release)"
        FAIL=$((FAIL+1))
    fi
}

# ---------- main loop ----------
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
[[ -n "$EXOMES_URL"  ]] && echo "Exomes source:  $EXOMES_URL"
[[ -n "$GENOMES_URL" ]] && echo "Genomes source: $GENOMES_URL"
echo ""

SUCCESS=0
EMPTY=0
FAIL=0

while IFS=$'\t' read -r -a FIELDS; do
    [[ "${FIELDS[0]}" == "Targeton_ID" ]] && continue

    TARGETON="${FIELDS[$((ID_COL-1))]}"
    CHROM="${FIELDS[$((CHR_COL-1))]}"
    START="${FIELDS[$((START_COL-1))]}"
    END="${FIELDS[$((END_COL-1))]}"

    if [[ -z "$TARGETON" || -z "$CHROM" || -z "$START" || -z "$END" ]]; then
        echo "WARNING: skipping incomplete row: ${FIELDS[*]}"
        continue
    fi

    # gnomAD v4.1 GRCh38 VCFs use 'chr'-prefixed contigs — add prefix if
    # the regions TSV doesn't already have one (keeps this robust either way).
    CHROM_FOR_QUERY="$CHROM"
    [[ "$CHROM_FOR_QUERY" != chr* ]] && CHROM_FOR_QUERY="chr${CHROM_FOR_QUERY}"

    [[ -n "$EXOMES_URL"  ]] && query_source "exomes"  "$EXOMES_URL"  "$TARGETON" "$CHROM_FOR_QUERY" "$START" "$END"
    [[ -n "$GENOMES_URL" ]] && query_source "genomes" "$GENOMES_URL" "$TARGETON" "$CHROM_FOR_QUERY" "$START" "$END"

done < "$REGIONS_TSV"

# ---------- summary ----------
echo ""
echo "Done."
echo "  Queries with variants (indexed): $SUCCESS"
echo "  Empty (no gnomAD hits):          $EMPTY"
echo "  Failed:                          $FAIL"

if [[ "$FAIL" -gt 0 ]]; then
    echo ""
    echo "NOTE: failures usually mean either (a) no network access to"
    echo "      storage.googleapis.com from this environment, or (b) a"
    echo "      contig-naming mismatch. Confirm with e.g.:"
    echo "        bcftools view -h <url> | grep '^##contig' | head"
fi
