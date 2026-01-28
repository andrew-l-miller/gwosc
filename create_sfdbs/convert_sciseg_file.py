# convert_segments.py
import sys

def load_sciseg_file(infile):
    scisegs = []
    with open(infile, "r") as fin:
        for line in fin:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            start, end = map(int, parts)
            scisegs.append((start, end))
    return scisegs

def convert_sciseg_file(infile, outfile):
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        for i, line in enumerate(fin, start=1):
            parts = line.strip().split()
            if len(parts) != 2:
                continue  # skip malformed lines
            start, end = map(int, parts)
            duration = end - start
            fout.write(f"{i} {start} {end} {duration}\n")

def time_in_science(t0, t1, scisegs):
    """
    Returns:
      inside_t0 (bool): t0 is within any sciseg
      inside_t1 (bool): t1 is within any sciseg
      t_in (float): total seconds of [t0,t1] that overlap science segments
      t_out (float): total seconds of [t0,t1] that are outside science segments
    Intervals treated as half-open: [start, end)
    """
    if t1 < t0:
        t0, t1 = t1, t0

    def contains(t, seg):
        s, e = seg
        return (s <= t) and (t < e)

    inside_t0 = any(contains(t0, seg) for seg in scisegs)
    inside_t1 = any(contains(t1, seg) for seg in scisegs)

    # overlap of [t0,t1] with union of scisegs
    t_in = 0.0
    for s, e in scisegs:
        # overlap length of [t0,t1] and [s,e]
        a = max(t0, s)
        b = min(t1, e)
        if b > a:
            t_in += (b - a)

    total = float(t1 - t0)
    t_out = total - t_in
    return inside_t0, inside_t1, t_in, t_out


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 convert_segments.py <input_file> <output_file>")
        sys.exit(1)
    infile = sys.argv[1]
    outfile = sys.argv[2]
    convert_sciseg_file(infile, outfile)
