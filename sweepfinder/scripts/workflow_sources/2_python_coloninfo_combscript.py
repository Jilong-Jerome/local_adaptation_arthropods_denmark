import sys

for line in sys.stdin:
    line = line.strip()
    columns = line.split('\t')
    nocolon = [col.split(':')[0] for col in columns]
    print('\t'.join(nocolon))
