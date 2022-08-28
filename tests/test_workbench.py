import unittest
from bioworkbench import biobaseproc

class TestDNAValidation(unittest.TestCase):

    def setUp(self) -> None:
        self.seq1 = "abcdefghijk"
        self.seq2 = "1234567"
        self.seq3 = "actgacgtagctagca"
        self.seq4 = "ACTGGGTG"

    def test_dna_validation(self):
        """
        Test that it can validate sequence as DNA sequence or not.
        """
        self.assertFalse(biobaseproc.validate_dna(self.seq1))
        self.assertFalse(biobaseproc.validate_dna(self.seq2))
        self.assertTrue(biobaseproc.validate_dna(self.seq3))
        self.assertTrue(biobaseproc.validate_dna(self.seq4))

class TestFrequency(unittest.TestCase):
    def setUp(self) -> None:
        self.seq1 = "aattggcc"
        self.seq2 = "bioworkbench"
        self.return1 = biobaseproc.frequency(self.seq1)
        self.return2 = biobaseproc.frequency(self.seq2)

    def test_frequency(self):
        self.assertEqual(self.return1["A"], 2, "Should be 2")
        self.assertEqual(self.return1["T"], 2, "Should be 2")
        self.assertEqual(self.return1["C"], 2, "Should be 2")
        self.assertEqual(self.return1["G"], 2, "Should be 2")
        self.assertEqual(self.return2["H"], 1, "Should be 1")
        self.assertNotEqual(self.return2["W"], 3, "Should be 1")

class TestGCContent(unittest.TestCase):
    def setUp(self) -> None:
        self.seq1 = "aattggcc"
        self.seq2 = "attaatt"

    def test_gc_content (self):
        self.assertEqual(biobaseproc.gc_content(self.seq1), 0.5)
        self.assertEqual(biobaseproc.gc_content(self.seq2), 0.0)

class TestGCContentSubseq(unittest.TestCase):
    def setUp(self) -> None:
        self.seq1 = "atgcgtcaagct"
        self.size1 = 4

    def test_gc_content_subseq(self):
        self.assertEqual(biobaseproc.gc_content_subseq(self.seq1, self.size1), [0.5, 0.5, 0.5])

# TESTS FOR TRANSCRIPTION
# TESTS FOR REVERSE_COMPLEMENT
# TESTS FOR TRANSLATE CODON
# TESTS FOR TRANSLATE SEQUENCE


############# patternfinding.py ################
# ALL TESTS MISSING

if __name__ == "__main__":
    unittest.main()