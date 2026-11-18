import sys
import gzip
import argparse


def _open_auto(path: str, mode: str):
    # mode: "rt" / "wt"
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def filter_ficture_pixels(_args):
    parser = argparse.ArgumentParser(
        prog="filter_ficture_pixels",
        description="Filter Ficture pixel file by P1 threshold, preserving original header/meta format."
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Input ficture pixel file (.tsv or .tsv.gz).")
    parser.add_argument("--output", type=str, required=True,
                        help="Output filtered file (recommend same suffix as input).")
    parser.add_argument("--pcut", type=float, default=0.8,
                        help="Keep rows with P1 >= pcut. Default: 0.8")
    args = parser.parse_args(_args)

    if len(_args) == 0:
        parser.print_help()
        return

    p1_idx = None

    with _open_auto(args.input, "rt") as fin, _open_auto(args.output, "wt") as fout:
        for line in fin:
            # 1) meta lines: keep exactly
            if line.startswith("##"):
                fout.write(line)
                continue
                
            if line.startswith("#"):
                fout.write(line)
                header = line.lstrip("#").rstrip("\n").split("\t")
                p1_idx = header.index("P1")
                continue

            # 3) data lines: keep original line text if pass threshold
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if float(parts[p1_idx]) >= args.pcut:
                fout.write(line)


if __name__ == "__main__":
    filter_ficture_pixels(sys.argv[1:])
