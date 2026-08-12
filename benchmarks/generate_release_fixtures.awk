BEGIN {
    if (call_path == "" || polysomy_path == "") {
        print "call_path and polysomy_path are required" > "/dev/stderr"
        exit 2
    }
    if (records_per_chrom == 0) records_per_chrom = 25000
    chromosomes = 12
    pi = atan2(0, -1)

    print "##fileformat=VCFv4.3" > call_path
    print "##FILTER=<ID=PASS,Description=\"All filters passed\">" > call_path
    for (chrom = 1; chrom <= chromosomes; chrom++) {
        printf "##contig=<ID=chr%d,length=%d>\n", chrom, records_per_chrom * 1000 + 1000 > call_path
    }
    print "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">" > call_path
    print "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">" > call_path
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY\tCONTROL" > call_path
    for (chrom = 1; chrom <= chromosomes; chrom++) {
        for (record = 1; record <= records_per_chrom; record++) {
            segment = int((record - 1) * 4 / records_per_chrom)
            if (segment == 0 || segment == 3) {
                query_baf = (record % 3 == 0) ? 0 : ((record % 3 == 1) ? 0.5 : 1)
                query_lrr = 0
            } else if (segment == 1) {
                query_baf = (record % 2 == 0) ? 0 : 1
                query_lrr = -0.45
            } else {
                query_baf = (record % 4 == 0) ? 0 : ((record % 4 == 1) ? 0.4 : ((record % 4 == 2) ? 0.6 : 1))
                query_lrr = 0.3
            }
            control_baf = (record % 3 == 0) ? 0 : ((record % 3 == 1) ? 0.5 : 1)
            printf "chr%d\t%d\t.\tA\tG\t.\tPASS\t.\tBAF:LRR\t%.7f:%.7f\t%.7f:0.0000000\n", chrom, record * 1000, query_baf, query_lrr, control_baf > call_path
        }
    }
    close(call_path)

    print "##fileformat=VCFv4.3" > polysomy_path
    print "##FILTER=<ID=PASS,Description=\"All filters passed\">" > polysomy_path
    for (chrom = 1; chrom <= chromosomes; chrom++) {
        printf "##contig=<ID=chr%d,length=%d>\n", chrom, records_per_chrom * 1000 + 1000 > polysomy_path
    }
    print "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">" > polysomy_path
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" > polysomy_path
    for (chrom = 1; chrom <= chromosomes; chrom++) {
        copy_number = (chrom <= 4) ? 2 : ((chrom <= 8) ? 3 : 4)
        peaks = copy_number + 1
        for (record = 1; record <= records_per_chrom; record++) {
            peak = (record - 1) % peaks
            mean = peak / copy_number
            u1 = (((record * 73 + chrom * 37) % 9973) + 0.5) / 9973
            u2 = (((record * 151 + chrom * 71) % 9967) + 0.5) / 9967
            z = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
            if (z < -3) z = -3
            if (z > 3) z = 3
            baf = mean + 0.025 * z
            if (baf < 0) baf = 0
            if (baf > 1) baf = 1
            printf "chr%d\t%d\t.\tA\tG\t.\tPASS\t.\tBAF\t%.7f\n", chrom, record * 1000, baf > polysomy_path
        }
    }
    close(polysomy_path)
}
