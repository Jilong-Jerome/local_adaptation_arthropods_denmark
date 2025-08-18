import sys

def generate_sliding_windows(chrom_length, window_size, step_size, output_file):
    with open(output_file, 'w') as f:
        start = 1
        while start + window_size - 1 <= chrom_length:
            end = start + window_size - 1
            f.write(f"{start}\t{end}\n")
            start += step_size

        # Handle the last window if it wasn't included already
        if start <= chrom_length:
            f.write(f"{start}\t{chrom_length}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <chrom_length> <window_size> <step_size> <output_file>")
        sys.exit(1)

    chrom_length = int(sys.argv[1])
    window_size = int(sys.argv[2])
    step_size = int(sys.argv[3])
    output_file = sys.argv[4]

    generate_sliding_windows(chrom_length, window_size, step_size, output_file)
    print(f"Sliding windows written to {output_file}")

