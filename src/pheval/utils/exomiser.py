from pathlib import Path

from pheval.infra.exomiserdb import ExomiserDB


def semsim_to_exomiserdb(input_path: Path, object_prefix: str, subject_prefix: str, db_path: Path):
    """ingests semsim file into exomiser phenotypic database

    Args:
        input_path (Path): semsim input file. e.g phenio-plus-hp-mp.0.semsimian.tsv
        object_prefix (str): object prefix. e.g. MP
        subject_prefix (str): subject prefix e.g HP
        db_path (Path): Exomiser Phenotypic Database Folder Path. (e.g. /exomiser_folder/2209_phenotype/2209_phenotype/)
    """
    exomiserdb = ExomiserDB(db_path)
    exomiserdb.import_from_semsim_file(input_path, object_prefix, subject_prefix)
