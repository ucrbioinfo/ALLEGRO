#!/usr/bin/python
import unittest
from chopchop import *


class FunctionsTest(unittest.TestCase):

    def test_get_mismatch_pos(self):
        self.assertEqual(get_mismatch_pos("23"), [])
        self.assertEqual(get_mismatch_pos("0T20A1"), [0,21])
        self.assertEqual(get_mismatch_pos("21T1"), [21])
        self.assertEqual([0], get_mismatch_pos("0A22"))
        self.assertEqual([1], get_mismatch_pos("1A21"))
        self.assertEqual([21], get_mismatch_pos("21A1"))
        self.assertEqual([22], get_mismatch_pos("22A"))
        self.assertEqual([0,22], get_mismatch_pos("A21T"))

    ### FIX
    def test_coordToFasta(self):
        region = ['Cap1_berg-2068c03.p1k:610-1673']
        fastaFile = "temporary/test/sequence.fa"
        outputDir = "temporary/test"
        #evalSeq = lambda name, guideSize, dna, num, fastaFile: eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, allowed=allowed, PAM=args.PAM)
        targetSize = 23
        indexDir = "genomes"
        genome = "Pberghei"

    def test_getMismatchVectors(self):
        pass

class HitTestCase(unittest.TestCase):

    def test_calc_mismatches_wNoMismatch(self):
        line = ["C:chr3:17181253-17181317_4-27", "16", "chr3", "17181258", "255", "23M", "*", "0", "0", "GGGGACAGAAGCTAAACTCATGG", "IIIIIIIIIIIIIIIIIIIIIII", "XA:i:0", "MD:Z:23", "NM:i:0"]
        hit = Hit(line)

        before = hit.matchSeq
        hit.calc_mismatchPos()
        after = hit.matchSeq

        self.assertEqual(before, after)


    def test_calc_mismatches_wMismatch(self):
        line = ["C:chr3:17184885-17185303_87-110", "0", "chr21", "16899387", "255", "23M", "*", "0", "0", "CCTCCTGCTGCGCGCGCGCTCCC", "IIIIIIIIIIIIIIIIIIIIIII", "XA:i:2", "MD:Z:0T20A1", "NM:i:2"]
        result = "tCTCCTGCTGCGCGCGCGCTCaC"

        hit = Hit(line)
        
        before = hit.matchSeq
        hit.calc_mismatchPos()
        after = hit.matchSeq

        self.assertNotEqual(before, after)
        self.assertEqual(result, after)
    

class GuideTestCase(unittest.TestCase):
    """ Tests the Guide class """



    def test_Guide(self):
        self.assertTrue(True)


    def test_addOffTarget_PAMwGG_mismatch(self):
        allowed, count = getMismatchVectors("NGG", False, False)

        correct = ["C:chr1:24415023-24415213_87-110","16","chr12", "4415111255","23M","*","0","0","AGCTTCAACACTGATCTGAGAGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:23","NM:i:0"]
        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)

        pamMut = ["C:chr1:24415023-24415213_87-110","16","chr10","549557255","23M","*","0","0","AGCTTCAACACTGATCTGAGAGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:1","MD:Z:21T1","NM:i:1"]
        guide.addOffTarget(Hit(pamMut), True, 50, allowed, count)

        self.assertEqual(guide.numOffTargets(), 0, msg="Added off target w/mismatch in GG PAM part")


    def test_addOffTarget_PAMwGG_mismatch_noCheck(self):
        allowed, count = getMismatchVectors("NGG", False, False)

        correct = ["C:chr1:24415023-24415213_87-110","16","chr12", "4415111255","23M","*","0","0","AGCTTCAACACTGATCTGAGAGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:23","NM:i:0"]

        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)

        pamMut = ["C:chr1:24415023-24415213_87-110","16","chr10","549557255","23M","*","0","0","AGCTTCAACACTGATCTGAGAGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:1","MD:Z:21T1","NM:i:1"]
        guide.addOffTarget(Hit(pamMut), False, 50, defAllowed, defCount)

        self.assertEqual(guide.numOffTargets(), 1, msg="Did not add off target w/PAM mutation when not checking for PAM mutations")




    def test_Hsu_mismatchInPAMGG(self):        

        # Example from pitx2 - correct w/flag 0
        allowed, count = getMismatchVectors("NGG", True, False)
        correct = ["C:chr14:37446948-37447175_41-64","0","chr14", "37446990","23M","*","0","0","CCGTTGAGTTTGGACCACCATCA","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:23","NM:i:0"]
        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)
        # Different Flag 16
        mmInPAMR = ["C:chr14:37446948-37447175_41-64","16","chr2", "20890061","23M","*","0","0","TGATGGTGGTCCAAACTCAACGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:1", "MD:Z:10G10A1","NM:i:1"]
        guide.addOffTarget(Hit(mmInPAMR), True, 50, allowed, count)
        self.assertEqual(guide.numOffTargets(), 0, msg="Added off target w/PAM GG mismatch in Hsu mode: Different flags 0 vs 16")


        # pitx2 - Correct w/flag 16
        allowed, count = getMismatchVectors("NGG", True, False)
        correct = ["C:chr14:37448936-37449546_21-44","16","chr14", "37448958","23M","*","0","0","AGGTGTGGTTCAAGAATCGACGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:23","NM:i:0"]
        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)
        # Different Flag 0
        mmInPAMF = ["C:chr14:37448936-37449546_21-44","0","chr15", "28360584","23M","*","0","0","CCGTCGATTCTTGAACCACACCT","IIIIIIIIIIIIIIIIIIIIIII","XA:i:1","MD:Z:0A22","NM:i:1"]
        guide.addOffTarget(Hit(mmInPAMF), True, 50, allowed, count)
        self.assertEqual(guide.numOffTargets(), 0, msg="Added off target w/PAM GG mismatch in Hsu mode: Different flags 16 vs 0")


        # Example from spaw - correct w/flag 16
        allowed, count = getMismatchVectors("NGG", True, False)
        correct = ["C:chr5:71786941-71787781_683-706", "16", "chr5", "71787625", "255", "23M", "*", "0", "0", "GGGTGGATTTTGATCAGATTGGG", "IIIIIIIIIIIIIIIIIIIIIII", "XA:i:0", "MD:Z:23", "NM:i:0"]
        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)
        # Same flag 16
        mmInPAMS =["C:chr5:71786941-71787781_683-706", "16", "chr2", "119840356", "255", "23M", "*", "0", "0", "GGGTGGATTTTGATCAGATTGGG", "IIIIIIIIIIIIIIIIIIIIIII", "XA:i:2", "MD:Z:7C14T0", "NM:i:2"]
        guide.addOffTarget(Hit(mmInPAMS), True, 50, allowed, count)
        self.assertEqual(guide.numOffTargets(), 0, msg="Added off target w/PAM GG mismatch in Hsu mode: Same flags 16")


        ### MISSING A 0 vs 0 TEST


    def test_Hsu_mismatchInPAMN(self):        
        allowed, count = getMismatchVectors("NGG", True, False)

        # pitx2 - Correct w/flag 0
        correct = ["C:chr14:37446948-37447175_41-64","0","chr14", "37446990","23M","*","0","0","CCGTTGAGTTTGGACCACCATCA","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:23","NM:i:0"]
        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)
        # Flag 0
        mmInPAM = ["C:chr14:37446948-37447175_41-64","0","chr2", "20890061","23M","*","0","0","CCGTTGAGTTTGGACCACCATCA","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:10G9A2","NM:i:0"]
        guide.addOffTarget(Hit(mmInPAM), True, 50, allowed, count)
        # Flag 16
        mmInPAM = ["C:chr14:37446948-37447175_41-64","16","chr2", "1000","23M","*","0","0","TGATGGTGGTCCAAACTCAACGG","IIIIIIIIIIIIIIIIIIIIIII","XA:i:0","MD:Z:10G9A2","NM:i:0"]
        guide.addOffTarget(Hit(mmInPAM), True, 50, allowed, count)


        self.assertEqual(guide.numOffTargets(), 2, msg="Did not add off target w/PAM N mismatch in Hsu mode")



    def test_Cong_mismatchIn11nearest_PAM(self):        

        # example camk2g1
        allowed, count = getMismatchVectors("NGG", False, True)


        correct = ["C:chr12:35924407-35924556_104-127", "16", "chr12", "35924512", "255", "23M", "*", "0", "0", "GGAGATCAGCAGGCCTGGTTTGG", "IIIIIIIIIIIIIIIIIIIIIII", "XA:i:0", "MD:Z:23", "NM:i:0"]
        guide = Guide(correct[0], correct[1], len(correct[9]), correct[9], True)
        hit = ["C:chr12:35924407-35924556_104-127", "0", "chr1", "13080466", "255", "23M", "*", "0", "0", "CCAAACCAGGCCTGCTGATCTCC", "IIIIIIIIIIIIIIIIIIIIIII", "XA:i:2", "MD:Z:8A11A2", "NM:i:2"]
        guide.addOffTarget(Hit(hit), True, 50, allowed, count)
        self.assertEqual(guide.numOffTargets(), 0, msg="Added off target w/mutation in Cong forbidden zone")




if __name__ == '__main__':
    unittest.main()
