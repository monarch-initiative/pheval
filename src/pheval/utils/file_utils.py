import difflib
from pathlib import Path


def files_with_suffix(directory: Path, suffix: str):
    """Obtains all files ending in a specified suffix from a given directory."""
    files = [path for path in directory.iterdir() if path.suffix == suffix]
    files.sort()
    return files


def all_files(directory: Path) -> list[Path]:
    """Obtains all files from a given directory."""
    files = [path for path in directory.iterdir()]
    files.sort()
    return files


def is_gzipped(path: Path) -> bool:
    """Confirms whether a file is gzipped."""
    return path.name.endswith(".gz")


def obtain_closest_file_name(file_to_be_queried: Path, file_paths: list[Path]) -> Path:
    """Obtains the closest file name when given a template file name and a list of full path of files to be queried."""
    closest_file_match = Path(
        str(
            difflib.get_close_matches(
                str(file_to_be_queried.name),
                [str(file_path.name) for file_path in file_paths],
            )[0]
        )
    )
    return [
        file_path for file_path in file_paths if Path(closest_file_match) == Path(file_path.name)
    ][0]
