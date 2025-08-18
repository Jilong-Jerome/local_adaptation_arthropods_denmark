import sys
# input parsing
input_vcf = sys.argv[1]
maf = float(sys.argv[2])
spid = sys.argv[3]
vcf = open(input_vcf,"r")
# output files
out_sf2 = open(f"{spid}_maf{maf}_sweepfinder_input.tsv","w")
out_sf2_pop_header = open(f"{spid}_maf{maf}_sweepfinder_pop_header.tsv","w")

for line in vcf:
    if ("##" in line): # test whether it is a header line
        pass
    elif ("CHROM" in line): # field info header
        pop_list = line.strip("\n").split("\t")[9:]
        out_sf2_pop_header.write('\t'.join(pop_list)+"\n")
    else:
        line_info = line.strip("\n").split("\t")
        chrom = line_info[0]
        pos = line_info[1]
        pops = line_info[9:]
        total_0_count = []
        total_1_count = []
        for pop in pops:
            total_0_count.append(pop.split(":")[0].count("0"))
            total_1_count.append(pop.split(":")[0].count("1"))
        freq = sum(total_1_count)/(sum(total_0_count)+sum(total_1_count))
        if (freq > maf) and (freq < 1-maf):
            info_line = '\t'.join([str(count1) for count1 in total_1_count])
            out_sf2.write(f"{chrom}\t{pos}\t{info_line}\n")





