import argparse
from pathlib import Path

def main(input_file: Path, output_file: Path | None ) -> None:
    if not output_file:
        output_file = input_file.with_suffix(".csv")

    if not input_file.exists():
        raise ValueError(f"Cannot find input file {input_file}")


    column_indexer = []
    rows = []
    blosum_matrix = {}
    with open(input_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            elif line.startswith(" "):
                column_indexer = line.strip().split("  ")
            else:
                rows.append(line.strip())


    for row in rows:
        idx = 0
        for cell in row.split(" "):
            if cell != "":
                if idx == 0:
                    from_letter = cell
                    blosum_matrix[from_letter] = []

                else:
                    blosum_matrix[from_letter].append(cell)
                idx += 1

    with open(output_file, "w") as outfh:
        outfh.write(" ,")
        outfh.write(",".join(column_indexer))
        outfh.write(",-")
        outfh.write("\n")

        for from_aa, row in blosum_matrix.items():
            outfh.write(f"{from_aa},")
            outfh.write(",".join(row))
            outfh.write(",0")
            outfh.write("\n")

        outfh.write(f"-,{','.join(['0' for _ in row])},0")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Path to the input file")
    parser.add_argument("--output", type=str, default=None, help="Output path")

    args = parser.parse_args()

    main(Path(args.input), args.output)


