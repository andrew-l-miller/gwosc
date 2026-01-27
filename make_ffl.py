#!/usr/bin/env python3
import os
import argparse

def make_ffl(input_dir, output_file, detector="H1", absolute=False):
    detector = detector.upper()
    if detector not in ("H1", "L1"):
        raise ValueError("Detector must be H1 or L1")

    # GWOSC frame prefix conventions
    prefix = "H-H1" if detector == "H1" else "L-L1"

    files = [
        f for f in os.listdir(input_dir)
        if f.startswith(prefix) and f.endswith(".gwf")
    ]
    files.sort()

    with open(output_file, "w") as ffl:
        for fname in files:
            path = os.path.join(input_dir, fname)
            ffl.write(os.path.abspath(path) + "\n" if absolute else fname + "\n")

    print(f"Wrote {len(files)} entries to {output_file} (prefix='{prefix}')")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create .ffl list of GWOSC .gwf files for a given detector"
    )
    parser.add_argument("input_dir", help="Directory containing .gwf files")
    parser.add_argument("output_file", help="Output .ffl file")
    parser.add_argument(
        "detector",
        nargs="?",
        choices=["H1", "L1"],
        default="H1",
        help="Detector to select files for (default: H1)"
    )
    parser.add_argument(
        "--absolute", action="store_true",
        help="Include absolute paths instead of just filenames"
    )

    args = parser.parse_args()
    make_ffl(args.input_dir, args.output_file, args.detector, args.absolute)

