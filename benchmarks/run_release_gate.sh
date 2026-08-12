#!/bin/sh
set -eu

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "usage: $0 FIXTURE_DIRECTORY RESULT_TSV [ROUNDS]" >&2
    exit 2
fi

fixture_dir=$1
result=$2
rounds=${3:-12}
rsomics=${RSOMICS_CNV_BIN:-${CARGO_TARGET_DIR:?CARGO_TARGET_DIR is required}/release/rsomics-cnv}
bcftools=${BCFTOOLS:-/opt/homebrew/bin/bcftools}

case $(uname -s) in
    Darwin) ;;
    *) echo "release gate currently records Darwin /usr/bin/time -lp RSS bytes" >&2; exit 2 ;;
esac

for path in "$fixture_dir/call.vcf" "$fixture_dir/polysomy.vcf" "$rsomics" "$bcftools"; do
    test -f "$path" || { echo "missing $path" >&2; exit 2; }
done
case $rounds in
    *[!0-9]*|'') echo "rounds must be a positive integer" >&2; exit 2 ;;
esac
test "$rounds" -gt 0 || { echo "rounds must be a positive integer" >&2; exit 2; }
test ! -e "$result" || { echo "result already exists: $result" >&2; exit 2; }
test ! -e "$result.meta" || { echo "metadata already exists: $result.meta" >&2; exit 2; }

work=$(mktemp -d "${TMPDIR:?TMPDIR is required}/rsomics-cnv-release-gate.XXXXXX")
cleanup() {
    case $work in
        "$TMPDIR"/rsomics-cnv-release-gate.*) rm -rf -- "$work" ;;
        *) echo "refusing to remove unexpected scratch path: $work" >&2 ;;
    esac
}
trap cleanup EXIT HUP INT TERM

mkdir -p "$(dirname "$result")"
measurements="$work/result.tsv"
metadata="$work/result.meta"
printf 'workflow\tround\torder\ttool\treal_seconds\tuser_seconds\tsystem_seconds\tmax_rss_bytes\n' > "$measurements"

record() {
    workflow=$1
    round=$2
    order=$3
    tool=$4
    timing=$5
    real=$(awk '$1 == "real" { print $2 }' "$timing")
    user=$(awk '$1 == "user" { print $2 }' "$timing")
    system=$(awk '$1 == "sys" { print $2 }' "$timing")
    rss=$(awk '$2 == "maximum" && $3 == "resident" { print $1 }' "$timing")
    test -n "$real" && test -n "$user" && test -n "$system" && test -n "$rss"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$workflow" "$round" "$order" "$tool" "$real" "$user" "$system" "$rss" >> "$measurements"
}

run_call() {
    round=$1
    tool=$2
    order=$3
    output="$work/call-$round-$tool"
    timing="$work/call-$round-$tool.time"
    if [ "$tool" = bcftools ]; then
        /usr/bin/time -lp "$bcftools" cnv -s QUERY -c CONTROL -O 0.2 -o "$output" \
            "$fixture_dir/call.vcf" > "$work/call-$round-$tool.stdout" 2> "$timing"
    else
        /usr/bin/time -lp "$rsomics" call -s QUERY -c CONTROL -O 0.2 -o "$output" \
            "$fixture_dir/call.vcf" > "$work/call-$round-$tool.stdout" 2> "$timing"
    fi
    record call "$round" "$order" "$tool" "$timing"
}

verify_call() {
    upstream=$1
    ours=$2
    for file in dat.QUERY.tab dat.CONTROL.tab cn.QUERY.tab cn.CONTROL.tab; do
        cmp "$upstream/$file" "$ours/$file"
    done
    awk -F '\t' 'BEGIN { OFS="\t" } $1 == "CF" { print } $1 == "RG" { print $1,$2,$3,$4,$5,$6 }' \
        "$upstream/summary.tab" > "$work/call-upstream.summary"
    awk -F '\t' 'BEGIN { OFS="\t" } $1 == "CF" { print } $1 == "RG" { print $1,$2,$3,$4,$5,$6 }' \
        "$ours/summary.tab" > "$work/call-rsomics.summary"
    cmp "$work/call-upstream.summary" "$work/call-rsomics.summary"
}

run_polysomy() {
    round=$1
    tool=$2
    order=$3
    output="$work/polysomy-$round-$tool"
    timing="$work/polysomy-$round-$tool.time"
    if [ "$tool" = bcftools ]; then
        /usr/bin/time -lp "$bcftools" polysomy -s SAMPLE -o "$output" \
            "$fixture_dir/polysomy.vcf" > "$work/polysomy-$round-$tool.stdout" 2> "$timing"
    else
        /usr/bin/time -lp "$rsomics" polysomy -s SAMPLE -o "$output" \
            "$fixture_dir/polysomy.vcf" > "$work/polysomy-$round-$tool.stdout" 2> "$timing"
    fi
    record polysomy "$round" "$order" "$tool" "$timing"
}

verify_polysomy() {
    upstream=$1
    ours=$2
    grep '^DIST' "$upstream/dist.dat" > "$work/polysomy-upstream.dist"
    grep '^DIST' "$ours/dist.dat" > "$work/polysomy-rsomics.dist"
    cmp "$work/polysomy-upstream.dist" "$work/polysomy-rsomics.dist"
    awk -F '\t' '
        function abs(value) { return value < 0 ? -value : value }
        NR == FNR && $1 == "FIT" { fit[$2 SUBSEP $4 SUBSEP $5]=$3; fits++; next }
        NR == FNR && $1 == "CN" { cn[$2]=$3; deviation[$2]=$4; calls++; next }
        $1 == "FIT" {
            key=$2 SUBSEP $4 SUBSEP $5
            if (!(key in fit) || abs(fit[key] - $3) > 0.00001) exit 1
            seen_fits++
        }
        $1 == "CN" {
            if (!($2 in cn) || cn[$2] != $3 || abs(deviation[$2] - $4) > 0.00001) exit 1
            seen_calls++
        }
        END { if (fits != seen_fits || calls != seen_calls) exit 1 }
    ' "$upstream/dist.dat" "$ours/dist.dat"
}

round=1
while [ "$round" -le "$rounds" ]; do
    if [ $((round % 2)) -eq 1 ]; then
        run_call "$round" bcftools 1
        run_call "$round" rsomics 2
    else
        run_call "$round" rsomics 1
        run_call "$round" bcftools 2
    fi
    verify_call "$work/call-$round-bcftools" "$work/call-$round-rsomics"
    round=$((round + 1))
done

round=1
while [ "$round" -le "$rounds" ]; do
    if [ $((round % 2)) -eq 1 ]; then
        run_polysomy "$round" bcftools 1
        run_polysomy "$round" rsomics 2
    else
        run_polysomy "$round" rsomics 1
        run_polysomy "$round" bcftools 2
    fi
    verify_polysomy "$work/polysomy-$round-bcftools" "$work/polysomy-$round-rsomics"
    round=$((round + 1))
done

{
    uname -a
    sw_vers
    sysctl -n machdep.cpu.brand_string
    "$bcftools" --version | sed -n '1,2p'
    "$rsomics" --version
    shasum -a 256 "$rsomics" "$fixture_dir/call.vcf" "$fixture_dir/polysomy.vcf"
    printf 'call records\t'; grep -vc '^#' "$fixture_dir/call.vcf"
    printf 'polysomy records\t'; grep -vc '^#' "$fixture_dir/polysomy.vcf"
    printf 'rounds\t%s\n' "$rounds"
} > "$metadata"

cp "$metadata" "$result.meta"
cp "$measurements" "$result"
