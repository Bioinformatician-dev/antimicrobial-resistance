from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def mini_db_fasta() -> Path:
    return DATA_DIR / "mini_resfinder.fasta"


@pytest.fixture(scope="session")
def test_assembly_fasta() -> Path:
    return DATA_DIR / "test_assembly.fasta"


@pytest.fixture()
def blast_db(tmp_path_factory, mini_db_fasta):
    from amr_consensus.db_manager import build_local_database

    db_dir = tmp_path_factory.mktemp("db")
    return build_local_database(mini_db_fasta, db_dir, name="resfinder")
