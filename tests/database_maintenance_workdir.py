import pathlib
import shutil

class DatabaseMaintenanceWorkdir:
    def prepare_workdir(self, work_dir_path: pathlib.Path, test_files_path: pathlib.Path, protein_data_file_path: pathlib.Path):
        # Prepare work directory for test
        if work_dir_path.is_dir():
            shutil.rmtree(str(work_dir_path))
        work_dir_path.mkdir(parents=True, exist_ok=True)
        ## Add protein data
        protein_data_path = work_dir_path.joinpath('protein_data/')
        protein_data_path.mkdir()
        shutil.copy(str(protein_data_file_path), str(protein_data_path))
        ## Add taxonomy data
        taxonomy_data_path = work_dir_path.joinpath('taxonomy_data/')
        taxonomy_data_path.mkdir()
        for dmp_file_path in test_files_path.glob('*.dmp'):
            shutil.copy(str(dmp_file_path), str(taxonomy_data_path))
