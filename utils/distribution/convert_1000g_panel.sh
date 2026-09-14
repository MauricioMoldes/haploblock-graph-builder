#!/bin/sh
# Converts the IGSR/1000 Genomes panel file (sample, pop, super_pop,
# gender -- tab-separated, with header) into our phenotypes.csv schema
# (individual_id,phenotype,value,source), so it can be used directly
# with usecase1_snp_distribution.py and usecase2_background_enrichment.py.
#
# Download the real panel file first:
#   wget https://ftp.1000genomes.org/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#
# Usage: ./convert_1000g_panel.sh <panel_file> <output_phenotypes.csv>

set -e

PANEL="$1"
OUT="${2:-phenotypes_real.csv}"

if [ -z "$PANEL" ]; then
    echo "Usage: $0 <panel_file> <output_phenotypes.csv>" >&2
    exit 1
fi

AWK=awk
command -v mawk >/dev/null 2>&1 && AWK=mawk

echo "individual_id,phenotype,value,source" > "$OUT"

"$AWK" -F'\t' -v out="$OUT" '
NR == 1 { next }
{
    sample = $1; pop = $2; super_pop = $3; gender = $4
    print sample",ancestry,"super_pop",1000G_panel" >> out
    print sample",population,"pop",1000G_panel" >> out
    if (gender != "") print sample",sex,"gender",1000G_panel" >> out
}
' "$PANEL"

echo "Done."
wc -l "$OUT"
echo "Ancestry (super_pop) value counts:"
awk -F, '$2=="ancestry"{print $3}' "$OUT" | sort | uniq -c
