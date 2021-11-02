# std imports
import pathlib
import os
from typing import List

# internal imports
from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.models.maintenance_information import MaintenanceInformation
from macpepdb.models.peptide import Peptide
from macpepdb.models.peptide_metadata import PeptideMetadata
from macpepdb.models.protein import Protein
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.proteomics.file_reader.uniprot_text_reader import UniprotTextReader
from macpepdb.tasks.database_maintenance.database_maintenance import DatabaseMaintenance, DatabaseStatus

# test imports
from tests.abstract_database_test_case import AbstractDatabaseTestCase
from tests.database_maintenance_workdir import DatabaseMaintenanceWorkdir

TRYPSIN_MAX_MISSED_CLEAVAGES = 2
TRYPSIN_MIN_PEPTIDE_LENGTH = 5
TRYPSIN_MAX_PEPTIDE_LENGTH = 40

class DigestionToDatabaseTestCase(AbstractDatabaseTestCase, DatabaseMaintenanceWorkdir):
    """[summary]
    Contains test to test the full database digest.
    """
    def test_maintenance(self):
        """
        1. Digests `text_files/proteins.txt` and verify the database
        2. Digest `test_files/B0FIH3_merge.txt` which will merge B0FIH3 
           with a slightly updated version of the protein creating new peptides.
           And verify again.
        """
        initial_digest_file_proteins = self.initial_digestion()
        initial_and_merged_digest_file_proteins = self.merge_digestion(initial_digest_file_proteins)
        self.update_digestion(initial_and_merged_digest_file_proteins)
        

    def initial_digestion(self) -> List[Protein]:
        """
        Digests `text_files/proteins.txt` and verify the database

        Returns
        -------
        The peptides from the file.
        """
        work_dir = pathlib.Path(f"./tmp/{self.id()}_digest")
        test_files_path = pathlib.Path('./test_files')
        protein_data_test_file_path = test_files_path.joinpath('proteins.txt')
        self.prepare_workdir(work_dir, test_files_path, protein_data_test_file_path)

        maintenance = DatabaseMaintenance(
            os.getenv("TEST_MACPEPDB_URL"),
            work_dir,
            4,
            5,
            'Trypsin',
            TRYPSIN_MAX_MISSED_CLEAVAGES,
            TRYPSIN_MIN_PEPTIDE_LENGTH,
            TRYPSIN_MAX_PEPTIDE_LENGTH
        )

        maintenance.start()

        EnzymeClass = DigestEnzyme.get_enzyme_by_name("Trypsin")
        trypsin = EnzymeClass(TRYPSIN_MAX_MISSED_CLEAVAGES, TRYPSIN_MIN_PEPTIDE_LENGTH, TRYPSIN_MAX_PEPTIDE_LENGTH)

        # Read proteins from file
        file_proteins = []
        with protein_data_test_file_path.open("r") as protein_data_test_file:
            protein_file_reader = UniprotTextReader(protein_data_test_file)
            file_proteins = [file_protein for file_protein in protein_file_reader]

        self.verify_database_integrity(file_proteins, trypsin)

        return file_proteins


    def merge_digestion(self, initial_digest_file_proteins: List[Protein]) -> List[Protein]:
        """
        Digests `test_files/B0FIH3_merge.txt`
        which will merge B0FIH3 with a slightly updates version of the protein
        creating new peptides.
        Than the database will be verified.

        Returns
        -------
        List of proteins containing the proteins from the inital digest with the applied merges.
        """
        # Run digest with updated B0FIH3.
        # The old peptide is merged with the new one
        # which has the accesseion 'NEWACC'
        work_dir = pathlib.Path(f"./tmp/{self.id()}_merge")
        test_files_path = pathlib.Path('./test_files')
        protein_data_test_file_path = test_files_path.joinpath('B0FIH3_merge.txt')
        self.prepare_workdir(work_dir, test_files_path, protein_data_test_file_path)

        maintenance = DatabaseMaintenance(
            os.getenv("TEST_MACPEPDB_URL"),
            work_dir,
            1,
            5,
            'Trypsin',
            TRYPSIN_MAX_MISSED_CLEAVAGES,
            TRYPSIN_MIN_PEPTIDE_LENGTH,
            TRYPSIN_MAX_PEPTIDE_LENGTH
        )

        maintenance.start()

        old_file_proteins_len = len(initial_digest_file_proteins)

        merged_file_proteins = []
        with protein_data_test_file_path.open("r") as protein_data_test_file:
            protein_file_reader = UniprotTextReader(protein_data_test_file)
            merged_file_proteins = [merged_file_protein for merged_file_protein in protein_file_reader]


        # Remove all file_proteins which are merged with the proteins in the merge file.
        # In this case only B0FIH3 should be remove and 
        for merged_file_protein in merged_file_proteins:
            for secondary_accession in merged_file_protein.secondary_accessions:
                file_proteins_to_remove = []
                for file_protein in initial_digest_file_proteins:
                    if secondary_accession == file_protein.accession or secondary_accession in file_protein.secondary_accessions:
                        file_proteins_to_remove.append(file_protein)
                for file_protein in file_proteins_to_remove:
                    initial_digest_file_proteins.remove(file_protein)

        # One protein should be removed
        self.assertEqual(len(initial_digest_file_proteins), old_file_proteins_len - 1)

        # Than add the merged proteins to the file_proteins
        # and verfy the database again
        initial_digest_file_proteins = initial_digest_file_proteins + merged_file_proteins
        self.assertEqual(len(initial_digest_file_proteins), old_file_proteins_len)

        EnzymeClass = DigestEnzyme.get_enzyme_by_name("Trypsin")
        trypsin = EnzymeClass(TRYPSIN_MAX_MISSED_CLEAVAGES, TRYPSIN_MIN_PEPTIDE_LENGTH, TRYPSIN_MAX_PEPTIDE_LENGTH)

        self.verify_database_integrity(initial_digest_file_proteins, trypsin)

        return initial_digest_file_proteins

    def update_digestion(self, initial_and_merged_digest_file_proteins: List[Protein]) -> List[Protein]:
        """
        Digests `test_files/NEWACC_updated.txt`
        which will add  B0FIJ1 to NEWACCs secondary accession (which ultimatly deletes and merges B0FIJ1 into NEWACC) 
        and also updating NEWACCs sequence without creating new peptides.
        Than the database will be verified.

        Returns
        -------
        List of proteins containing the proteins from the merge with the applied updates.
        """
        # Run digest with updated B0FIH3.
        # The old peptide is merged with the new one
        # which has the accesseion 'NEWACC'
        work_dir = pathlib.Path(f"./tmp/{self.id()}_merge")
        test_files_path = pathlib.Path('./test_files')
        protein_data_test_file_path = test_files_path.joinpath('NEWACC_updated.txt')
        self.prepare_workdir(work_dir, test_files_path, protein_data_test_file_path)

        maintenance = DatabaseMaintenance(
            os.getenv("TEST_MACPEPDB_URL"),
            work_dir,
            1,
            5,
            'Trypsin',
            TRYPSIN_MAX_MISSED_CLEAVAGES,
            TRYPSIN_MIN_PEPTIDE_LENGTH,
            TRYPSIN_MAX_PEPTIDE_LENGTH
        )

        maintenance.start()

        old_file_proteins_len = len(initial_and_merged_digest_file_proteins)

        merged_file_proteins = []
        with protein_data_test_file_path.open("r") as protein_data_test_file:
            protein_file_reader = UniprotTextReader(protein_data_test_file)
            merged_file_proteins = [merged_file_protein for merged_file_protein in protein_file_reader]

        self.assertEqual(len(merged_file_proteins), 1)
        newacc_upgraded_protein = merged_file_proteins.pop()

        # NEWACC will be updated, so replace it.
        newacc_idx = -1
        for file_protein_index, file_protein in enumerate(initial_and_merged_digest_file_proteins):
            if file_protein.accession == newacc_upgraded_protein.accession:
                newacc_idx = file_protein_index
                break

        initial_and_merged_digest_file_proteins[newacc_idx] = newacc_upgraded_protein

        # Remove B0FIJ1 from file proteins
        file_proteins_to_remove = []
        for file_protein in initial_and_merged_digest_file_proteins:
            if file_protein.accession in newacc_upgraded_protein.secondary_accessions:
                file_proteins_to_remove.append(file_protein)

        for file_protein in file_proteins_to_remove:
            initial_and_merged_digest_file_proteins.remove(file_protein)

        # Length should be decreased by one
        self.assertEqual(len(initial_and_merged_digest_file_proteins), old_file_proteins_len - 1)

        EnzymeClass = DigestEnzyme.get_enzyme_by_name("Trypsin")
        trypsin = EnzymeClass(TRYPSIN_MAX_MISSED_CLEAVAGES, TRYPSIN_MIN_PEPTIDE_LENGTH, TRYPSIN_MAX_PEPTIDE_LENGTH)

        self.verify_database_integrity(initial_and_merged_digest_file_proteins, trypsin)

        return initial_and_merged_digest_file_proteins


    def verify_database_integrity(self, proteins_from_file: List[Protein], enzym: DigestEnzyme):
        """
        Verifies the database by:
        1. Check if all protein from file exists and their attributes are matching
        2. Digest the given proteins and check if:
            2.1 The peptides are found in the database (by primary key)
            2.2 The values which are generated on the fly and not send from the database, e.g. amino acid counts, matches the on in the database.
        3. Check if all proteins and their peptides have association and if the association count matches the actual protein peptides relationships
        4. Check if all peptides have a related metadata record

        Parameters
        ----------
        proteins_from_file : List[Protein]
            Proteins read from the protein file
        enzym : DigestEnzyme
            Enzym for digesting. Shoud match the one which is used for the database creation.
        """
        peptides_from_file_proteins = set()
        for file_protein in proteins_from_file:
            for new_peptide in enzym.digest(file_protein):
                peptides_from_file_proteins.add(new_peptide)

        with self.database_connection.cursor() as database_cursor:
            # Check if protein count in database are equals to set
            database_cursor.execute(f"SELECT count(*) FROM {Protein.TABLE_NAME};")
            self.assertEqual(len(proteins_from_file), database_cursor.fetchone()[0])

            # Check if all proteins are correct proteins 
            for file_protein in proteins_from_file:
                db_protein = Protein.select(
                    database_cursor,
                    WhereCondition(
                        "accession = %s",
                        [file_protein.accession]
                    ),
                )
                self.assertIsNotNone(db_protein)
                self.assertEqual(
                    db_protein.accession,
                    file_protein.accession
                )
                self.assertEqual(
                    db_protein.secondary_accessions,
                    file_protein.secondary_accessions
                )
                self.assertEqual(
                    db_protein.entry_name,
                    file_protein.entry_name
                )
                self.assertEqual(
                    db_protein.name,
                    file_protein.name
                )
                self.assertEqual(
                    db_protein.sequence,
                    file_protein.sequence
                )
                self.assertEqual(
                    db_protein.taxonomy_id,
                    file_protein.taxonomy_id
                )
                self.assertEqual(
                    db_protein.proteome_id,
                    file_protein.proteome_id
                )
                self.assertEqual(
                    db_protein.is_reviewed,
                    file_protein.is_reviewed
                )

            # Check if set count is equals db count
            # Because peptides are not removed from the database it is possible to have more peptides
            # in the database after protein updates than in the file.
            database_cursor.execute(f"SELECT count(*) FROM {Peptide.TABLE_NAME};")
            self.assertLessEqual(len(peptides_from_file_proteins), database_cursor.fetchone()[0])

            for file_peptide in peptides_from_file_proteins:
                db_peptide = Peptide.select(
                    database_cursor,
                    WhereCondition(
                        "partition = %s AND mass = %s AND sequence = %s",
                        [file_peptide.partition, file_peptide.mass, file_peptide.sequence]
                    )
                )
                self.assertIsNotNone(db_peptide)
                self.assertEqual(
                    db_peptide.sequence,
                    file_peptide.sequence
                )
                self.assertEqual(
                    db_peptide.mass,
                    file_peptide.mass
                )
                self.assertEqual(
                    db_peptide.partition,
                    file_peptide.partition
                )
                self.assertEqual(
                    db_peptide.number_of_missed_cleavages,
                    file_peptide.number_of_missed_cleavages
                )

                # Because the amino acid counts are counted on the fly to save I/O and bandwidth, lets check the values in the database
                database_cursor.execute(
                    (
                        "SELECT "
                        "a_count, "
                        "b_count, "
                        "c_count, "
                        "d_count, "
                        "e_count, "
                        "f_count, "
                        "g_count, "
                        "h_count, "
                        "i_count, "
                        "j_count, "
                        "k_count, "
                        "l_count, "
                        "m_count, "
                        "n_count, "
                        "o_count, "
                        "p_count, "
                        "q_count, "
                        "r_count, "
                        "s_count, "
                        "t_count, "
                        "u_count, "
                        "v_count, "
                        "w_count, "
                        "y_count, "
                        "z_count, "
                        "n_terminus, "
                        "c_terminus "
                        f"FROM {Peptide.TABLE_NAME} "
                        "WHERE partition = %s AND mass = %s AND sequence = %s"
                    ),
                    (file_peptide.partition, file_peptide.mass, file_peptide.sequence)
                )
                db_peptide_record = database_cursor.fetchone()
                self.assertIsNotNone(db_peptide_record)
                # file_peptide attributes in the array below have the same order as in the query
                for value_idx, file_peptide_value in enumerate([
                    file_peptide.a_count,
                    file_peptide.b_count,
                    file_peptide.c_count,
                    file_peptide.d_count,
                    file_peptide.e_count,
                    file_peptide.f_count,
                    file_peptide.g_count,
                    file_peptide.h_count,
                    file_peptide.i_count,
                    file_peptide.j_count,
                    file_peptide.k_count,
                    file_peptide.l_count,
                    file_peptide.m_count,
                    file_peptide.n_count,
                    file_peptide.o_count,
                    file_peptide.p_count,
                    file_peptide.q_count,
                    file_peptide.r_count,
                    file_peptide.s_count,
                    file_peptide.t_count,
                    file_peptide.u_count,
                    file_peptide.v_count,
                    file_peptide.w_count,
                    file_peptide.y_count,
                    file_peptide.z_count,
                    file_peptide.get_n_terminus_ascii_dec(),
                    file_peptide.get_c_terminus_ascii_dec()
                ]):
                    self.assertEqual(file_peptide_value, db_peptide_record[value_idx])

            # Check protein/peptide-associations from both directions
            desired_number_of_associations = 0
            for file_protein in proteins_from_file:
                for file_peptide in enzym.digest(file_protein):
                    # Increase association counter
                    desired_number_of_associations += 1
                    database_cursor.execute(
                        (
                            "SELECT true "
                            f"FROM {ProteinPeptideAssociation.TABLE_NAME} "
                            "WHERE protein_accession = %s AND partition = %s AND peptide_mass = %s AND peptide_sequence = %s;"
                        ),
                        (
                            file_protein.accession,
                            file_peptide.partition,
                            file_peptide.mass,
                            file_peptide.sequence,
                        )
                    )
                    is_association_found = database_cursor.fetchone()[0]
                    self.assertIsNotNone(is_association_found)
                    self.assertTrue(is_association_found)

            # Check association counter. Must be equals even after updates.
            database_cursor.execute(f"SELECT count(*) FROM {ProteinPeptideAssociation.TABLE_NAME};") 
            self.assertEqual(desired_number_of_associations, database_cursor.fetchone()[0])

            # Check if peptide metadata equals peptides
            database_cursor.execute(f"SELECT count(*) FROM {PeptideMetadata.TABLE_NAME};")
            metadata_count = database_cursor.fetchone()[0]
            database_cursor.execute(f"SELECT count(*) FROM {Peptide.TABLE_NAME};")
            peptide_count = database_cursor.fetchone()[0]
            self.assertEqual(metadata_count, peptide_count)
            # Check if the current peptides have updated metadata
            for file_peptide in peptides_from_file_proteins:
                file_peptide.fetch_metadata_from_proteins(database_cursor)
                db_metadata = PeptideMetadata.select(database_cursor, file_peptide)
                self.assertIsNotNone(
                    db_metadata,
                    f"metadata for peptide '{file_peptide.sequence}' is missing"
                )
                if db_metadata:
                    self.assertEqual(
                        db_metadata.is_swiss_prot,
                        file_peptide.metadata.is_swiss_prot
                    )
                    self.assertEqual(
                        db_metadata.is_trembl,
                        file_peptide.metadata.is_trembl
                    )
                    self.assertEqual(
                        sorted(db_metadata.taxonomy_ids),
                        sorted(file_peptide.metadata.taxonomy_ids)
                    )
                    self.assertEqual(
                        sorted(db_metadata.unique_taxonomy_ids),
                        sorted(file_peptide.metadata.unique_taxonomy_ids)
                    )
                    self.assertEqual(
                        sorted(db_metadata.proteome_ids),
                        sorted(file_peptide.metadata.proteome_ids)
                    )

            # Check if maintenance mode is false and update timestamp is greater zero
            database_status = MaintenanceInformation.select(
                database_cursor,
                MaintenanceInformation.DATABASE_STATUS_KEY
            )
            self.assertNotEqual(database_status, None)
            self.assertGreater(database_status.values['last_update'], 0)
            self.assertEqual(database_status.values['status'], DatabaseStatus.READY.value)
            self.assertFalse(database_status.values['maintenance_mode'])
