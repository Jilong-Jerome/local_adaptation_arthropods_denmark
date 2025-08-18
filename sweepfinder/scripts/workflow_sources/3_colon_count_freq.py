import sys
input_geno_info = sys.argv[1]
output_count1 = sys.argv[2]
output_freq = sys.argv[3]

out_count = open(output_count1,"w")
out_freq = open(output_freq,"w")

with open(input_geno_info, 'r') as file:
        for line in file:
            line=line.strip()
            columns=line.split('\t')
            chr_pos=columns[:2]
            nocolon=[col.split(':')[0] for col in columns[2:]]
            count_ones = [col.count('1') for col in nocolon]
            out_count.write('\t'.join(chr_pos + list(map(str, count_ones))) + '\n')
            frequency = [str(float(col) / 100) for col in count_ones]
            out_freq.write('\t'.join(chr_pos + [str(x) for x in frequency]) + '\n')

