# std imports
import unittest

# internal imports
from macpepdb.models.protein import Protein
from macpepdb.proteomics.enzymes.trypsin import Trypsin

# Peptides for Leptin (UniProt accession Q257X2) digested with 3 missed cleavages, length 0 - 60
# Tested with https://web.expasy.org/peptide_mass/
# Leucine and Isoleucine are replaced with J already! 
DESIRED_RESULTS = {
    "VTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASK",
    "NVIQISNDLENLRDLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSR",
    "DLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLR",
    "QRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLR",
    "INDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSR",
    "SCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC",
    "MDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASKSCPLPQVR",
    "VTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLR",
    "SCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLR",
    "ALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC",
    "DLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSR",
    "MDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASK",
    "QRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSR",
    "TIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSK",
    "VTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSR",
    "ALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLR",
    "CGPLYRFLWLWPYLSYVEAVPIRKVQDDTK",
    "SCPLPQVRALESLESLGVVLEASLYSTEVVALSR",
    "INDISHTQSVSSKQRVTGLDFIPGLHPLLSLSK",
    "MDQTLAIYQQILASLPSRNVIQISNDLENLR",
    "NVIQISNDLENLRDLLHLLAASKSCPLPQVR",
    "FLWLWPYLSYVEAVPIRKVQDDTKTLIK",
    "MRCGPLYRFLWLWPYLSYVEAVPIRK",
    "MRCGPLYRFLWLWPYLSYVEAVPIR",
    "VQDDTKTLIKTIVTRINDISHTQSVSSK",
    "CGPLYRFLWLWPYLSYVEAVPIRK",
    "FLWLWPYLSYVEAVPIRKVQDDTK",
    "CGPLYRFLWLWPYLSYVEAVPIR",
    "ALESLESLGVVLEASLYSTEVVALSR",
    "TLIKTIVTRINDISHTQSVSSKQR",
    "NVIQISNDLENLRDLLHLLAASK",
    "TLIKTIVTRINDISHTQSVSSK",
    "FLWLWPYLSYVEAVPIRK",
    "TIVTRINDISHTQSVSSKQR",
    "QRVTGLDFIPGLHPLLSLSK",
    "FLWLWPYLSYVEAVPIR",
    "MDQTLAIYQQILASLPSR",
    "TIVTRINDISHTQSVSSK",
    "LQGSLQDMLRQLDLSPGC",
    "DLLHLLAASKSCPLPQVR",
    "VTGLDFIPGLHPLLSLSK",
    "KVQDDTKTLIKTIVTR",
    "VQDDTKTLIKTIVTR",
    "INDISHTQSVSSKQR",
    "NVIQISNDLENLR",
    "INDISHTQSVSSK",
    "KVQDDTKTLIK",
    "VQDDTKTLIK",
    "LQGSLQDMLR",
    "DLLHLLAASK",
    "TLIKTIVTR",
    "MRCGPLYR",
    "SCPLPQVR",
    "KVQDDTK",
    "QLDLSPGC",
    "CGPLYR",
    "VQDDTK",
    "TIVTR",
    "TLIK",
	"MR",
	"QR",
	"K"
}

# Peptides for Leptin (UniProt accession Q257X2) digested with 3 missed cleavages, length 0 - 60
# Tested with https://web.expasy.org/peptide_mass/ 
DESIRED_RESULTS_FOR_UNKNOWN_AMINO_ACID_FILTERING = {
    "QRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLR",
    "INDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSR",
    "SCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC",
    "VTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLR",
    "SCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLR",
    "ALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC",
    "QRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSR",
    "VTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSR",
    "ALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLR",
    "CGPLYRFLWLWPYLSYVEAVPIRKVQDDTK",
    "SCPLPQVRALESLESLGVVLEASLYSTEVVALSR",
    "INDISHTQSVSSKQRVTGLDFIPGLHPLLSLSK",
    "MDQTLAIYQQILASLPSRNVIQISNDLENLR",
    "FLWLWPYLSYVEAVPIRKVQDDTKTLIK",
    "MRCGPLYRFLWLWPYLSYVEAVPIRK",
    "MRCGPLYRFLWLWPYLSYVEAVPIR",
    "CGPLYRFLWLWPYLSYVEAVPIRK",
    "FLWLWPYLSYVEAVPIRKVQDDTK",
    "CGPLYRFLWLWPYLSYVEAVPIR",
    "ALESLESLGVVLEASLYSTEVVALSR",
    "FLWLWPYLSYVEAVPIRK",
    "QRVTGLDFIPGLHPLLSLSK",
    "FLWLWPYLSYVEAVPIR",
    "MDQTLAIYQQILASLPSR",
    "LQGSLQDMLRQLDLSPGC",
    "VTGLDFIPGLHPLLSLSK",
    "INDISHTQSVSSKQR",
    "NVIQISNDLENLR",
    "INDISHTQSVSSK",
    "KVQDDTKTLIK",
    "VQDDTKTLIK",
    "LQGSLQDMLR",
    "MRCGPLYR",
    "SCPLPQVR",
    "KVQDDTK",
    "QLDLSPGC",
    "CGPLYR",
    "VQDDTK",
    "TLIK",
	"MR",
	"QR",
	"K"
}


class TrypsinTestCase(unittest.TestCase):
    def test_digest(self):
        # Using Leptin (UniProt  accession: Q257X2)
        leptin = Protein("Q257X2", ["TESTACC"], "LEP_CAPHI", "Leptin", "MRCGPLYRFLWLWPYLSYVEAVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC", 9925, "UP000291000", True, 1145311200)
        trypsin = Trypsin(3, 0, 60)
        peptides = trypsin.digest(leptin)

        self.assertEqual(len(DESIRED_RESULTS), len(peptides))

        for peptide in peptides:
            self.assertTrue(peptide.sequence in DESIRED_RESULTS)
        

    def test_removal_of_sequences_with_unknown_amino_acid(self):
        # Using Leptin (UniProt  accession: Q257X2) with a X in the peptides `DLLHLLAASK` and `TIVTR`
        leptin = Protein("Q257X2", ["TESTACC"], "LEP_CAPHI", "Leptin", "MRCGPLYRFLWLWPYLSYVEAVPIRKVQDDTKTLIKTIXVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHXLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC", 9925, "UP000291000", True, 1145311200)
        trypsin = Trypsin(3, 0, 60)
        peptides = trypsin.digest(leptin)

        self.assertEqual(len(DESIRED_RESULTS_FOR_UNKNOWN_AMINO_ACID_FILTERING), len(peptides))

        for peptide in peptides:
            self.assertEqual(peptide.sequence.count('X'), 0)
            self.assertIn(peptide.sequence, DESIRED_RESULTS_FOR_UNKNOWN_AMINO_ACID_FILTERING)
