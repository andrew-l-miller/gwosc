# convert_segments.py
import sys

def convert_sciseg_file(infile, outfile):
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        for i, line in enumerate(fin, start=1):
            parts = line.strip().split()
            if len(parts) != 2:
                continue  # skip malformed lines
            start, end = map(int, parts)
            duration = end - start
            fout.write(f"{i} {start} {end} {duration}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 convert_segments.py <input_file> <output_file>")
        sys.exit(1)
    infile = sys.argv[1]
    outfile = sys.argv[2]
    convert_sciseg_file(infile, outfile)
