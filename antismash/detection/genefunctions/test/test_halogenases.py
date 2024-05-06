# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

""" KtzQ: Trp-7
    KtzR: Trp-6
    mibH: Trp-5
    bmp2: Pyrrole
    ChlB4: Orsellinic-like
    End30: Hpg
    BhaA: Tyrosine
"""

import json
import unittest
from unittest.mock import patch

from antismash.common import utils, secmet, subprocessing
from antismash.common.subprocessing import hmmscan
from antismash.common.secmet.test.helpers import (
    DummyRecord,
    DummyFeature,
    DummyCDS,
)
from antismash.common.signature import HmmSignature
from antismash.common.test.helpers import FakeHSPHit, FakeHit
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
    phenolic,
    pyrrolic,
)

from antismash.detection.genefunctions.halogenases import specific_analysis
from antismash.detection.genefunctions.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as FDH,
    Match,
)

from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    run_halogenase_phmms,
    extract_residues,
    categorize_on_substrate_level,
    fdh_specific_analysis,
)

from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import FDH_SUBGROUPS


TRANSLATIONS = {
    "ktzQ": (
        "MDDNRIRSILVLGGGTAGWMSACYLSKALGPGVEVTVLEAPSISRIRVGEATIPNLHKVF"
        "FDFLGIAEDEWMRECNASYKAAVRFVNWRTPGDGQATPRRRPDGRPDHFDHLFGQLPEHE"
        "NLPLSQYWAHRRLNGLTDEPFDRSCYVQPELLDRKLSPRLMDGTKLASYAWHFDADLVAD"
        "FLCRFAVQKLNVTHVQDVFTHADLDQRGHITAVNTESGRTLAADLFIDCSGFRSVLMGKV"
        "MQEPFLDMSKHLLNDRAVALMLPHDDEKVGIEPYTSSLAMRSGWSWKIPLLGRFGSGYVY"
        "SSQFTSQDEAAEELCRMWDVDPAEQTFNNVRFRVGRSRRAWVRNCVAIGVSAMFVEPLES"
        "TGLYFSYASLYQLVKHFPDKRFRPILADRFNREVATMYDDTRDFLQAHFSLSPRDDSEFW"
        "RACKELPFADGFAEKVEMYRAGLPVELPVTIDDGHYYGNFEAEFRNFWTNSNYYCIFAGL"
        "GFLPEHPLPVLEFRPEAVDRAEPVFAAVRRRTEELVATAPTMQAYLRRLHQGT"
        ),
    "ktzR": (
        "MTAAYLKTAFGDRLSITVVESSRIGTIGVGEATFSDIQHFFQFL"
        "NLREQDWMPACNATYKLGIRFENWRHVGHHFYQPFEQIRPVYGFPLTDWWLHDAPTDR"
        "FDTDCFVMPNLCEAGRSPRHLDGTLADEDFVEEGDELANRTMSEHQGKSQFPYAYHFE"
        "AALLAKFLTGYAVDRGVEHVVDDVLDVRLDQRGWIEHVVTAEHGEIHGDLFVDCTGFR"
        "GLLLNKALGVPFVSYQDTLPNDSAVALQVPLDMQRRGIVPNTTATAREAGWIWTIPLF"
        "GRVGTGYVYAKDYLSPEEAERTLREFVGPAAADVEANHIRMRIGRSQESWRNNCVAIG"
        "LSSGFVEPLESTGIFFIHHAIEQLVKHFPAADWNPKSRDMYNSAVAHVMDGIREFLVI"
        "HYRGAARADNQYWRDTKTRPLPDGLAERIECWQTQLPDTETIYPYYHGLPPYSYMCIL"
        "MGGGAIRTPASAALALTDQGAAQKEFAAVRDRAAQLRDTLPSHYEYLARMRGLDV"
        ),
    "mibH": (
        "MLKNVVVVGGGTAGWMTASYLTAAFGDRIGVTLVESKRVGSIGVGEATFSTVRHFFEYLG"
        "LEEKEWMPACNATYKLAIRFENWREPGHHFYHPFERQRVVDGFPLTDWWLREPRSDRFDK"
        "DCFLVGTLCDDLKSPRQLNGELFEGGLGGRSAYRTTLAEQTTQFPYAYHFDATLVANYLR"
        "DYAVARGVKHVLDDVQDVALDDRGWISHVVTGESGNLTGDLFIDCTGFRSLLLGKALAEP"
        "FQSYQDSLPNDSAVALRVPQDMENRGLRPCTTATAQEAGWIWTIPLFDRIGTGYVYAGDY"
        "ISPEEAERTLRAFVGPAAEHADANHIKMRIGRSNRHWVNNCVAVGLSSGFVEPLESTGIF"
        "FIQHAIEQLVKHFPDERWDDGLRTAYNKLVNNVMDGVREFLVVHYYAAKRQDNQYWKDAK"
        "TRPLPDGLAERLERWQTRLPDNESVFPHYHGFESYSYVCMLLGLGGLDLKSSPALGLMDA"
        "APARHEFKLVGEQAAELARTLPTQYEYFAQLHRAR"
    ),
    "bmp2": (
        "MDQFKSYDVVIIGSGPAGSLCGIECRKKGLSVLCIEKDEFPRFHIGESLTGNAGQIIRDL"
        "GLADEMNAAGFPDKPGVNVIGSLSKNEFFIPILAPTWQVRRSDFDNMLKRRALEHGVEYQ"
        "QGLVKDVIKHEEKVVGAIYKADGVEHQVRSKVLVDASGQNTFLSRKGIAGKREIEFFSQQ"
        "IASFAHYKNVERDLPPFSTNTTILYSKQYHWSWIIPISPDTDSLGIVIPKDLYYKECKNP"
        "DDAIEWGMEHISPEIRRRFKNAERVGESQSMADFSYRIEPFVGDGWLCIGDAHRFLDPIF"
        "SYGVSFAMKEGIKAADAIKRAIDGNDWKTPFYEYRDWSNGGQQIAADLIRYFWIYPIFFG"
        "YQMQNPDLRDEVIRLLGGCCFDCEGWKAPTIFRNAIEEYDRKQMAG"
    ),
    "ChlB4": (
        "MQPDFDAAIVGGGPAGSAMASYLAEAGLSVAVFESEMFPRPHIGESLVPATMPVLDEIGV"
        "MPDIEAAGFPKKYGAAWTSAESRDVPHNGFTGLDHDFKAAEVMFVERDQPGVHRDYTFHV"
        "DRGKFDLILLKHAESRGAQVFQKTRVLKADFDTDPDLVTLNCRLGPRTLDFTTRMVIDAS"
        "GRQTMLGNQLKVKVPDPVFNQYAIHAWFEGLDRTAMALDPAKRDYIYVHFLPLEDTWMWQ"
        "IPITDTITSVGVVTQKHRFKAASADREKFFWDIVSSRKDIYDALQKAERIRPFKAEGDYS"
        "YAMRQICGDRFLLIGDAARFVDPIFSSGVSVALNSARLAAKDVIAAHRAGDFRKESFATY"
        "EEKLRRAVRNWYEFISVYYRLNILFTAFVQDPRYRIDVLKMLQGDFYDGEEPKALKAMRD"
        "LVTKVENDPEHLWHPYLGTLRAPSAAPTF"
    ),
    "End30": (
        "STLSTLVAMQGHSVLLLEKETFPRYQIGESLLPSTIHGICHLLGVTDELAAAGFPHKRGG"
        "TFRWGASPKPWNFSFSVSSKVSGPTSFAYQVERSKFDKILLDNAARKGVVVRQDRTVTDV"
        "VDDADGRARGLRYTDPDGTEHEVSARYVVDASGNTSRIHKRVGGSRTYSDFFKSLALFGY"
        "FENGKRMPAPYAGNILCVAFGSGWFWYIPLSSTLTSVGAVVRREDAAKVQGDPESALRGL"
        "IDECPMIKEYLADATRVTTGQYGQLRVRKDYSYHHTTFWRPGMVLVGDAACFVDPVFSSG"
        "VHLATYSALLAARSLNSVLAGRIDERRAFDEFEARYRREYGVFYEFLTSFYDMHVDEDSY"
        "FWTAKK"
    ),
    "BhaA": (
        "STVATLVAMQGHRVLLLEKEVFPRYQIGESLLPATVHGVCRMLGISDELANAGFPIKRGG"
        "TFRWGARPEPWTFHFGISAKMAGSTSHAYQVERARFDEMLLNNAKRKGVVVREGCAVTDV"
        "VEDGERVTGARYTDPDGTEREVSARFVIDASGNKSRLYTKVGGSRNYSEFFRSLALFGYF"
        "EGGKRLPEPVSGNILSVAFDSGWFWYIPLSDTLTSVGAVVRREDAEKIQGDREKALNTLI"
        "AECPLISEYLADATRVTTGRYGELRVRKDYSYQQETYWRPGMILVGDAACFVDPVFSSGV"
        "HLATYSALLAARSINSVLAGDLDEKTALNEFELRYRREYGVFYEFLVSFYQMNVNEESYF"
        "WQAKK"
    ),
    "VatD": (
        "MIEVCIIGFGFSAVPLIRELQRTGTEFKIISEESNSVWDALSQSNRLDFDLVSSYLTSFY"
        "SFDLVKDFVEDYYPTSKQFYEMHQRWRKVYENEIIRDRVTRIDNFKEHSVIFTKSGKTLN"
        "AKHVICSTGFSRAIHTHINDIDYSVSNKTFVFDTMGDSANLIISKLIPNNNKIIIRTNGF"
        "NARDKVVPGAGAIYTLDQLEGHNFRYMSHEHYGSVIYGLPIGSKKPILMGDQFPVTVRDD"
        "NYITSKSRPASGTIAIKYWPIDQYADKFGNNLEESISQGYLLNDIAMWLHTGKAIVVPKD"
        "TAINFEKKTITYAGIERAFHQYIKGDPEQPRLPKIMIDGNTPYEYQYRDNFMGVIPRTLN"
        "NVYFIGYTRPYTGGLANIIEMQGLFVHKMITQSEFHQKIHHNLDERIVAYNNHYYGTTKP"
        "RSADHLVYFGFYTDDLARLIGIDYKPSEINSIKDMVFYYAFPNNALKYRLKGEYAVDGVE"
        "DLIKKINEKYYDFIDVFAYLWGTSKMDSVELTEELEQFIRQYFNDMRHKEPYTKFLENYI"
        "QVYRRVKNTRVDETDDYEWSLMVKKASETRDRVLQEFKESGDYQLEENFRNFISNEIELI"
        "QSLMNYKILSVKDGQLKIVPEKIGGHSILENILCKITKILGFQGIIESVLGKNEMLPKIE"
        "AQDLQSLLSLTKPKEYELLYLKP"
    ),
    "CtoA": (
        "MEANPTAGTEVVVIGAGIVGVHSAIQFAKRGLKVVLIDNIVGQKKSFKVGESLLVFS"
        "NMFLRTISELDEFNQKCFPKHGVWFTYGMEGTTSFEEKAEWALESTLPQAMRDAFANKAL"
        "LRAMADDVQIVRPEAEELMQQTARAHPNITFLDTAKVTNVVIAEGGGPHEVTWECKATQR"
        "TGVVRTTWLIDCSGRNRLLAKKLKHAAEDVELNDGFKTTAVWGQFSGIKDEMFGENWVNR"
        "TSDGARSKRDLNTLHLWGDGYWIWVIRLSEGRISVGATYDQRRPPAGAGYREKFWDIIRR"
        "YPLFDGMLSDDNMLEFHVFKNCQHITDTFVSEKRYGMIGDAASVIDAYYSQGVSLALVTS"
        "WHITNIMERDLRERRLDKEYIARVNRHTRQDWHIMRNMVIEKYTSAMADGRFFVMTHLLD"
        "MIIFVGAAFPRYLLVRWLVETQGSTARETPVMREMRRYLEENLYYSKIGSLAPEKVQKVQ"
        "RGLQASLSERARWRVENGVKVTRLKAIVHAPSGLLKFWKLPLSGQREFEDISPKPVKQIP"
        "KWLAMTGEETNPRMLKMARPLMASTFFLMYGYDGLSTAVTKVRQRLERLPGAAATAETTA"
        "AGRRGEAPEPAMNGAAPVRNVLREVSA"
    ),
}


def create_flavin_match(profile, confidence=0., consensus_residues="", substrates=None,
                        target_positions=None, number_of_decorations=None,
                        ):
    return Match(
        profile,
        "flavin",
        "FDH",
        confidence=confidence,
        consensus_residues=consensus_residues,
        substrates=substrates,
        target_positions=target_positions,
        number_of_decorations=number_of_decorations,
    )


class PyrrolicBase(unittest.TestCase):
    def setUp(self):
        self.pyrrole_enzyme = FDH("bmp2")

        self.pyrrole_hmm_result = HalogenaseHmmResult(
            hit_id="pyrrole_FDH",
            bitscore=1000,
            query_id="pyrrole_FDH",
            enzyme_type="FDH",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )

        self.negative_pyrrole_hmm_result = HalogenaseHmmResult(
            hit_id="pyrrole_FDH",
            bitscore=1,
            query_id="pyrrole_FDH",
            enzyme_type="FDH",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )


class IndolicBase(unittest.TestCase):
    def setUp(self):
        self.test_trp_5_match = create_flavin_match(
            "trp_5_FDH", confidence=1, substrates=("tryptophan",), target_positions=5,
            number_of_decorations="mono",
        )
        self.test_trp_6_7_match = create_flavin_match(
            "trp_6_7_FDH", confidence=1, substrates=("tryptophan",),
            target_positions=6, number_of_decorations="mono",
        )

        self.trp_5_hmm_result = HalogenaseHmmResult(
            hit_id="trp_5_FDH",
            bitscore=1000,
            query_id="trp_5_FDH",
            enzyme_type="Flavin-dependent",
            profile=indolic.SPECIFIC_PROFILES[0].path,
        )
        self.trp_6_7_hmm_result = HalogenaseHmmResult(
            hit_id="trp_6_7_FDH",
            bitscore=1000,
            query_id="trp_6_7_FDH",
            enzyme_type="Flavin-dependent",
            profile=indolic.SPECIFIC_PROFILES[1].path,
        )

        tryptophan_single_matches = [self.test_trp_6_7_match]
        tryptophan_matches = [self.test_trp_5_match, self.test_trp_6_7_match]

        # Trp-5 halogenase
        self.trp_5_enzyme_with_matches = FDH("mibH", potential_matches=tryptophan_matches)
        # Trp-6 halogenase
        self.trp_enzyme_with_matches = FDH("ktzR", potential_matches=tryptophan_single_matches)
        # Trp-7 halogenase
        self.trp_empty_enzyme = FDH("ktzQ")

    def specific_analysis_test(self, name, fake_hit,
                               fake_hsps, categorize_return_value,
                               get_best_matches_return_value,
                               analysis_function):
        with patch.object(hmmscan, "run_hmmscan",
                          return_value=[fake_hsps]) as _patched_hmmscan:
            for value in _patched_hmmscan.return_value:
                value.id = name
                value.hits = fake_hit
            with patch.object(subprocessing, "run_hmmsearch",
                              return_value=fake_hit) as _patched_run_hmmsearch:
                for value in _patched_run_hmmsearch.return_value:
                    value.hsps = [fake_hsps]
                with patch.object(substrate_analysis, "categorize_on_substrate_level",
                                  return_value=categorize_return_value):
                    with patch.object(FDH, "get_best_matches",
                                      return_value=get_best_matches_return_value):
                        cds = DummyCDS(locus_tag=name,
                                       translation=TRANSLATIONS[name])
                        record = DummyRecord(features=[cds])
                        with patch.object(secmet.Record, "get_cds_features_within_regions",
                                          return_value=[cds]):
                            return analysis_function(record)


class TestPhenolic(unittest.TestCase):
    def setUp(self):
        self.tyrosine_hmm_result = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=1000,
            query_id="tyrosine-like_hpg_FDH",
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[0].path,
        )
        self.less_confident_tyrosine_hmm_result = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=310,
            query_id="tyrosine-like_hpg_FDH",
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[0].path
        )
        self.hpg_hmm_result = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=600,
            query_id="tyrosine-like_hpg_FDH",
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[0].path,
        )
        self.cycline_orsellinic_hmm_result = HalogenaseHmmResult(
            hit_id="cycline_orsellinic_FDH",
            bitscore=600,
            query_id="cycline_orsellinic_FDH",
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[1].path,
        )

        # tyrosine
        self.tyr_enzyme = FDH("BhaA")

        # Hpg
        self.hpg_enzyme = FDH("End30")

        # other phenolic-substrate halogenase (orsellinic-like)
        self.orsellinic_empty_enzyme = FDH("ChlB4")

        self.tyrosine_match_not_hpg = {
            "Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
            "Hpg": "VALAMI",
        }

    def test_invalid_profile(self):
        invalid_hit = HalogenaseHmmResult("wrong_name", 400, "wrong_name", "foo", "wrong_name")
        with self.assertRaisesRegex(ValueError, "unknown profile"):
            categorize_on_substrate_level(DummyCDS(), FDH("ChlB4"), [invalid_hit])

    def test_categorising_orsellinic(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=phenolic.ORSELLINIC.motif_residues):
            categorize_on_substrate_level(DummyCDS(), self.orsellinic_empty_enzyme,
                                          [self.cycline_orsellinic_hmm_result])
        match = self.orsellinic_empty_enzyme.potential_matches[0]
        assert match.profile == "cycline_orsellinic_FDH"
        assert match.confidence == 1
        assert match.substrates == ("cycline_orsellinic-like",)
        assert match.target_positions == [6, 8]

    def test_categorising_hpg(self):
        with patch.object(substrate_analysis, "extract_residues",
                          return_value=""):
            assert not categorize_on_substrate_level(DummyCDS(), self.orsellinic_empty_enzyme,
                                                     [self.cycline_orsellinic_hmm_result])

        # Tyr and Hpg halogenases
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=phenolic.TYR_HPG.motif_residues):
            categorize_on_substrate_level(DummyCDS(translation=TRANSLATIONS["End30"]), self.hpg_enzyme,
                                          [self.hpg_hmm_result])
        match = self.hpg_enzyme.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.confidence == 1
        assert match.substrates == ("Hpg",)
        assert match.target_positions == [6, 8]

    def test_categorising_tyrosine(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=self.tyrosine_match_not_hpg):
            assert categorize_on_substrate_level(DummyCDS(translation=TRANSLATIONS["BhaA"]), self.tyr_enzyme,
                                                 [self.less_confident_tyrosine_hmm_result])
        match = self.tyr_enzyme.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.confidence == 0.5
        assert match.substrates == ("Tyr",)
        assert match.target_positions == [6, 8]

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={None}):
            assert categorize_on_substrate_level(DummyCDS(), self.tyr_enzyme,
                                                 self.tyrosine_hmm_result) is None

    def test_retrieve_fdh_signature_residues(self):
        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = [FakeHit(1, 2, 1000, "foo")]
            # checking if hit_id != query_id it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("BhaA", hit_id="tyrosine-like_hpg_FDH")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = TRANSLATIONS["BhaA"]
                    hit_profile.seq = TRANSLATIONS["BhaA"]
                    hit.aln = [hit_profile, hit_query]
            cds = DummyCDS(locus_tag="BhaA", translation=TRANSLATIONS["BhaA"])
            assert " " not in cds.translation
            residues = substrate_analysis.retrieve_fdh_signature_residues(
                cds.translation, self.hpg_hmm_result,
                [phenolic.TYROSINE_LIKE_MOTIF.positions, phenolic.HPG_MOTIF.positions],
                enzyme_substrates=["Tyr", "Hpg"],
            )
        assert isinstance(residues, dict)
        assert residues["Tyr"] is not None and residues["Hpg"] is not None


class TestPyrrolic(PyrrolicBase):
    def test_get_consensus_signature(self):
        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = [FakeHit(1, 2, 1000, "foo")]
            # checking if hit_id != query_id it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("bmp2", hit_id="pyrrole_FDH")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = TRANSLATIONS["bmp2"]
                    hit_profile.seq = TRANSLATIONS["bmp2"]
                    hit.aln = [hit_profile, hit_query]
            sigs = pyrrolic.get_consensus_signature(DummyCDS(translation=TRANSLATIONS["bmp2"]),
                                                    self.pyrrole_hmm_result)
            assert sigs == {"pyrrole_FDH": {"mono_di": "RAKDIM", "tetra": "RAKDIM", "unconv_mono_di": "RAKDIM"}}

    def test_score_below_cutoff(self):
        cds = DummyCDS()
        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": pyrrolic.PYRROLE.motif_residues}) as patched:
            assert not categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                     [self.negative_pyrrole_hmm_result])
            patched.assert_called_once_with(cds, self.negative_pyrrole_hmm_result)

    def test_no_motif(self):
        cds = DummyCDS()
        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": None}) as patched:
            assert not categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                     [self.pyrrole_hmm_result])
            patched.assert_called_once_with(cds, self.pyrrole_hmm_result)

    def test_conventional_mono_di(self):
        cds = DummyCDS()
        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": {"mono_di": "DRSVFW"}}) as patched:
            assert categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                 [self.pyrrole_hmm_result])
            match = self.pyrrole_enzyme.potential_matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == "mono_di"
            assert match.profile == "pyrrole_FDH"
            assert match.substrates == ("pyrrole",)
            patched.assert_called_once_with(cds, self.pyrrole_hmm_result)

    def test_unconventional_mono_di(self):
        cds = DummyCDS()
        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": {"unconv_mono_di": "YRRNFN"}}) as patched:
            assert categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                 [self.pyrrole_hmm_result])
            match = self.pyrrole_enzyme.potential_matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == "unconv_mono_di"
            assert match.profile == "pyrrole_FDH"
            assert match.substrates == ("pyrrole",)
            patched.assert_called_once_with(cds, self.pyrrole_hmm_result)

    def test_tetra(self):
        cds = DummyCDS()
        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": {"tetra": "RRYFFA"}}) as patched:
            assert categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                 [self.pyrrole_hmm_result])
            match = self.pyrrole_enzyme.potential_matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == "tetra"
            assert match.profile == "pyrrole_FDH"
            assert match.substrates == ("pyrrole",)
            patched.assert_called_once_with(cds, self.pyrrole_hmm_result)

    def test_residue_mismatch(self):
        residues = {"pyrrole_FDH": {"pyrrole_FDH": "ABCDEF"}}
        assert not pyrrolic.PYRROLE.get_matches_from_hit(residues, self.pyrrole_hmm_result)


class TestIndolic(IndolicBase):
    def test_get_consensus_signature(self):
        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = [FakeHit(1, 2, 1000, "foo")]
            # checking if hit_id != query_id it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("ktzR", hit_id="trp_6_7_FDH")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = TRANSLATIONS["ktzR"]
                    hit_profile.seq = TRANSLATIONS["ktzR"]
                    hit.aln = [hit_profile, hit_query]
            sigs = indolic.get_consensus_signature(DummyCDS(translation=TRANSLATIONS["ktzR"]),
                                                   self.trp_6_7_hmm_result)
            assert sigs == {"trp_6_7_FDH": {}}

    def test_get_best_matches(self):
        assert not self.trp_empty_enzyme.get_best_matches()

        positive_test_best_match = self.trp_enzyme_with_matches.get_best_matches()

        assert len(positive_test_best_match) == 1 and isinstance(positive_test_best_match[0], Match)

        self.trp_enzyme_with_matches.add_potential_match(self.test_trp_6_7_match)
        assert len(self.trp_enzyme_with_matches.potential_matches) == 2

        multiple_matches = self.trp_enzyme_with_matches.get_best_matches()
        assert len(multiple_matches) == 2 and isinstance(positive_test_best_match[0], Match)

        one_potential_match = FDH("test_enzyme", potential_matches=[self.test_trp_5_match])
        multiple_matches = one_potential_match.get_best_matches()
        assert len(multiple_matches) == 1 and isinstance(positive_test_best_match[0], Match)

    # halogenases analysis
    def test_conversion_methods(self):
        original = self.trp_enzyme_with_matches.to_json()
        assert isinstance(original, dict)
        raw = json.loads(json.dumps(original))

        rebuilt = FDH.from_json(raw)
        assert isinstance(rebuilt, FDH)
        assert rebuilt.to_json() == original

    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit(1, 2, 1000, "foo")])
    def test_run_halogenase_phmms(self, run_hmmsearch):
        signature = HmmSignature("foo", "description", 300, "dummy_path")
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=250)]

        negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
        assert not negative_test_halogenase_hmms_by_id

        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=1000)]

        positive_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
        assert positive_test_halogenase_hmms_by_id
        for hit in positive_test_halogenase_hmms_by_id["foo"]:
            assert isinstance(hit, HalogenaseHmmResult)

    def test_search_signature_residues(self):
        target_positions = indolic.TRP_6_MOTIF.positions

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = []
            signature_residues = extract_residues(TRANSLATIONS["ktzR"],
                                                  target_positions, self.trp_6_7_hmm_result)
            assert signature_residues is None

            run_hmmpfam2.return_value = [FakeHit(1, 2, 1000, "foo")]
            # checking if hit_id != query_id it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("foo", hit_id="foo")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = "ghjvghjkbln"
                    hit_profile.seq = "xfdhcgkbjlnkml"
                    hit.aln = [hit_profile, hit_query]

            signature_residues = extract_residues(list(TRANSLATIONS.values())[0],
                                                  target_positions, self.trp_6_7_hmm_result)
            assert signature_residues is None

            # checking if hit_id == query_id it runs til the reference target_positions extraction
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("trp_6_7_FDH", hit_id="trp_6_7_FDH")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = TRANSLATIONS["ktzR"]
                    hit_profile.seq = "xfdhcgkbjlnkml"
                    hit.aln = [hit_profile, hit_query]

            with patch.object(utils, "extract_by_reference_positions", return_value="dummy") as patched:
                signature_residues = extract_residues(TRANSLATIONS["ktzR"], target_positions, self.trp_6_7_hmm_result)
                assert signature_residues == "dummy"
                patched.assert_called_once_with("xfdhcgkbjlnkml", TRANSLATIONS["ktzR"], list(target_positions))

    def test_false_check_for_match(self):
        false_test = FDH("fake_name")
        false_test_residues = "FAKESEQUENCE"
        hit = HalogenaseHmmResult(
            hit_id="mismatch",
            bitscore=1000,
            query_id="query_name",
            enzyme_type="Flavin-dependent",
            profile=indolic.SPECIFIC_PROFILES[1].path,
        )
        with self.assertRaisesRegex(ValueError, "unhandled profile"):
            FDH_SUBGROUPS["trp_5_FDH"].update_match(false_test_residues, false_test, hit)

    def test_strong_trp_5(self):
        cds = DummyCDS()
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"trp_5_FDH": indolic.TRP_5_MOTIF.residues}) as patched:
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [self.trp_5_hmm_result])
            patched.assert_called_once_with(
                cds.translation, self.trp_5_hmm_result,
                (indolic.TRP_5_MOTIF.positions,), [indolic.TRP_5_MOTIF.name],
            )

        match = self.trp_empty_enzyme.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 1
        assert match.substrates == ("tryptophan",)
        assert match.number_of_decorations == "mono"
        assert match.target_positions == [5]

    def test_weak_trp_5(self):
        with patch.object(indolic, "get_consensus_signature",
                          return_value={"trp_5_FDH": indolic.TRP_5.motif_residues}):
            low_quality_hit = HalogenaseHmmResult("trp_5_FDH", 380, "trp_5_FDH",
                                                  "foo", "trp_5_FDH")
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [low_quality_hit])
        match = self.trp_empty_enzyme.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 0.5
        assert match.substrates == ("tryptophan",)
        assert match.number_of_decorations == "mono"
        assert match.target_positions == [5]

    def test_trp_6(self):
        with patch.object(indolic, "get_consensus_signature",
                          return_value={"trp_6_7_FDH": indolic.TRP_6.motif_residues}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme,
                                          [self.trp_6_7_hmm_result])
        match = self.trp_empty_enzyme.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.0
        assert match.substrates == ("tryptophan",)
        assert match.number_of_decorations == "mono"
        assert match.target_positions == [6]

    def test_trp_7(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"trp_6_7_FDH": indolic.TRP_6.motif_residues}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [self.trp_6_7_hmm_result])
        match = self.trp_empty_enzyme.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.
        assert match.target_positions == [7]
        assert match.substrates == ("tryptophan",)
        assert match.number_of_decorations == "mono"

    def test_no_match(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={}):
            assert categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme,
                                                 [self.trp_6_7_hmm_result]) is None

    def test_invalid_profile(self):
        invalid_hit = HalogenaseHmmResult("wrong_name", 400, "wrong_name",
                                          "foo", "wrong_name")
        with self.assertRaisesRegex(ValueError, "unknown profile"):
            categorize_on_substrate_level(DummyCDS(), FDH("mibH"), [invalid_hit])

    @patch.object(indolic, "get_consensus_signature",
                  return_value={"trp_5_FDH": indolic.TRP_5.motif_residues})
    def test_categorise_substrate_no_match(self, _patched_get_consensus_signature):
        result = categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [])
        assert result is None

    @patch.object(indolic, "get_consensus_signature",
                  return_value={"trp_5_FDH": indolic.TRP_5.motif_residues})
    def test_categorise_substrate_good_match(self, _patched_consensus_sig):
        cds = DummyCDS(locus_tag="mibH", translation=TRANSLATIONS["mibH"])
        result = categorize_on_substrate_level(cds, FDH("mibH"), [self.trp_5_hmm_result])
        assert isinstance(result, FDH)
        assert result.cds_name == "mibH"
        assert result.potential_matches == [
            create_flavin_match(profile="trp_5_FDH", confidence=1.0,
                                consensus_residues=indolic.TRP_5_MOTIF.residues,
                                substrates=("tryptophan",), target_positions=[5],
                                number_of_decorations="mono")
        ]


class TestSpecificAnalysis(IndolicBase):
    def test_one_best_match(self):
        positive_test = self.specific_analysis_test("ktzR",
                                                    FakeHit(1, 2, 900, "ktzR"),
                                                    FakeHSPHit("trp_6_7_FDH", "ktzR", bitscore=1500),
                                                    self.trp_enzyme_with_matches,
                                                    [self.test_trp_6_7_match],
                                                    fdh_specific_analysis)
        assert positive_test[0].family == "FDH"
        assert positive_test[0].cds_name == "ktzR"
        assert positive_test[0].cofactor == "flavin"
        assert positive_test[0].substrates == ("tryptophan",)
        assert positive_test[0].number_of_decorations == "mono"
        assert positive_test[0].potential_matches
        assert positive_test and positive_test[0].target_positions == 6

    def test_more_best_match(self):
        positive_test = self.specific_analysis_test("mibH", FakeHit(1, 2, 500, "mibH"),
                                                    FakeHSPHit("trp_5_FDH", "mibH", bitscore=800),
                                                    self.trp_5_enzyme_with_matches,
                                                    [self.test_trp_5_match, self.test_trp_6_7_match],
                                                    fdh_specific_analysis)
        assert len(positive_test) == 1
        assert positive_test[0].potential_matches
        assert positive_test[0].consensus_residues is None
        assert not positive_test[0].number_of_decorations
        assert positive_test[0].target_positions is None

    @patch.object(hmmscan, "run_hmmscan", return_value=[])
    def test_no_hits(self, _patched_hmmscan):
        cds = DummyCDS(locus_tag="mibH",
                       translation=TRANSLATIONS["mibH"])
        record = DummyRecord(features=[cds])
        with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
            results = fdh_specific_analysis(record)

        assert results == []


class TestGeneralEnzymes(IndolicBase):
    def setUp(self):
        super().setUp()

        self.general_empty_enzyme = FDH("CtoA")
        self.unconventional_empty_enzyme = FDH("VatD")

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="VALAMIVALAMI")
    def test_unconventional_fdh_specific_analysis(self, _patched_extract_residues):
        positive_test = self.specific_analysis_test("VatD", FakeHit(1, 2, 200, "unconventional_FDH"),
                                                    FakeHSPHit("unconventional_FDH", "VatD", bitscore=200),
                                                    self.unconventional_empty_enzyme, [],
                                                    fdh_specific_analysis)
        assert isinstance(positive_test[0], FDH)
        assert positive_test[0].substrates is None
        assert not positive_test[0].consensus_residues
        assert not positive_test[0].potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional_fdh_specific_analysis(self, _patched_extract_residues):
        positive_test = self.specific_analysis_test("CtoA",
                                                    FakeHit(1, 2, 200, "all_general_FDH"),
                                                    FakeHSPHit("all_general_FDH", "CtoA", bitscore=200),
                                                    self.general_empty_enzyme, [],
                                                    fdh_specific_analysis)

        assert not positive_test[0].substrates
        assert not positive_test[0].target_positions

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_random(self, _patched_extract_residues):
        result = substrate_analysis.categorize_on_consensus_level(
            DummyCDS(locus_tag="CtoA", translation=TRANSLATIONS["CtoA"]),
            {},
            [HalogenaseHmmResult("CtoA", 200, "all_conventional_FDH", "FDH", "dummy_path")],
        )
        assert result.consensus_residues == {"W.W.I.": "WIWVIR"}
        assert result.substrates is None
        assert not result.potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_specific_analysis(self, _patched_extract_residues):
        categorized_halogenase = self.specific_analysis_test("CtoA", FakeHit(1, 2, 200, "all_general_FDH"),
                                                             FakeHSPHit("all_general_FDH", "CtoA", bitscore=200),
                                                             self.general_empty_enzyme, [],
                                                             specific_analysis)
        assert categorized_halogenase[0].consensus_residues == {"W.W.I.": "WIWVIR"}

    @patch.object(substrate_analysis, "extract_residues", return_value="VALAMI")
    def test_negative_get_conserved_motifs(self, _patched_extract_residues):
        categorized_halogenase = self.specific_analysis_test("VatD", FakeHit(1, 2, 200, "unconventional_FDH"),
                                                             FakeHSPHit("unconventional_FDH", "VatD", bitscore=200),
                                                             self.unconventional_empty_enzyme, [],
                                                             specific_analysis)
        assert categorized_halogenase
        assert categorized_halogenase[0].consensus_residues == {}

    @patch.object(substrate_analysis, "extract_residues", return_value="")
    def test_negative_enzyme_hits_fasta(self, _patched_extract_residues):
        no_potential_matches = self.specific_analysis_test("VatD",
                                                           FakeHit(1, 2, 200, "unconventional_FDH"),
                                                           FakeHSPHit("unconventional_FDH", "VatD", bitscore=200),
                                                           self.unconventional_empty_enzyme, [], fdh_specific_analysis)
        assert not no_potential_matches[0].potential_matches
