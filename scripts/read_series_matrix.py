import argparse
import csv
from pathlib import Path

import pandas as pd


VALID_HEADER_NAMES = {"ID_REF", "IDREF"}


def find_series_matrix_files(directory: Path) -> list[Path]:
    """Return all series matrix CSV files in a directory tree."""
    return sorted(directory.rglob("*_series_matrix.csv"))


def detect_table_start(file_path: Path) -> tuple[int, int]:
    """
    Find the 1-based line where the data table starts.

    Returns:
        table_start_line: Line number containing the header row.
        metadata_end_line: Last metadata line before the table starts.
    """
    with file_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue

            first_cell = row[0].strip().strip('"')
            if first_cell in VALID_HEADER_NAMES:
                return line_number, line_number - 1

    raise ValueError("Could not find a header row starting with ID_REF or IDREF.")


def read_table_shape(file_path: Path, table_start_line: int) -> tuple[int, int]:
    """Read the actual table and return its shape."""
    dataframe = pd.read_csv(file_path, skiprows=table_start_line - 1)
    return dataframe.shape


def parse_directory(group_name: str, directory: Path) -> None:
    """Parse all series matrix files within one directory."""
    if not directory.exists():
        raise FileNotFoundError(f"{group_name} directory does not exist: {directory}")

    if not directory.is_dir():
        raise NotADirectoryError(f"{group_name} path is not a directory: {directory}")

    files = find_series_matrix_files(directory)
    print(f"\n{group_name} directory: {directory}")

    if not files:
        print("Files found: none")
        return

    print("Files found:")
    for file_path in files:
        print(file_path.as_posix())

    for file_path in files:
        table_start_line, metadata_end_line = detect_table_start(file_path)
        n_rows, n_cols = read_table_shape(file_path, table_start_line)

        print(f"\nFile: {file_path.as_posix()}")
        print(f"Shape: ({n_rows}, {n_cols})")
        print(f"Metadata ends at line: {metadata_end_line}")
        print(f"Table starts at line: {table_start_line}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse GEO series matrix CSV files from AD and PD directories."
    )
    parser.add_argument(
        "--ad-dir",
        required=True,
        help="Directory containing AD series matrix CSV files.",
    )
    parser.add_argument(
        "--pd-dir",
        required=True,
        help="Directory containing PD series matrix CSV files.",
    )
    args = parser.parse_args()

    parse_directory("AD", Path(args.ad_dir))
    parse_directory("PD", Path(args.pd_dir))


if __name__ == "__main__":
    main()
