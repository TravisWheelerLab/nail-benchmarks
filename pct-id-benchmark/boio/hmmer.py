from dataclasses import dataclass, field
from typing import List
import re

# # hmmsearch :: search profile(s) against a sequence database
# # HMMER 3.4 (Aug 2023); http://hmmer.org/
# # Copyright (C) 2023 Howard Hughes Medical Institute.
# # Freely distributed under the BSD open source license.
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # query HMM file:                  tmp/tmp.hmm
# # target sequence database:        tmp/tmp.fa
# # max ASCII text line length:      unlimited
# # Max sensitivity mode:            on [all heuristic filters off]
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Query:       2HCT  [M=415]
# Accession:   PF03390.19
# Description: 2-hydroxycarboxylate transporter family
# Scores for complete sequences (score includes all domains):
#    --- full sequence ---   --- best 1 domain ---    -#dom-
#     E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
#     ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
#    5.6e-188  610.7  33.9   6.3e-188  610.5  33.9    1.0  1  2HCT-consensus


# Domain annotation for each sequence (and alignments):
# >> 2HCT-consensus
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#    1 !  610.5  33.9  6.3e-188  6.3e-188       1     415 []       1     416 []       1     416 [] 0.99

#   Alignments for each domain:
#   == domain 1  score: 610.5 bits;  conditional E-value: 6.3e-188
#             2HCT   1 iggiplplfllllavlllavltgklpkdligalavllvlGillgeiGkriPilkkylGggailallvpsalvylgllpeevvkavttfmkksnfldlyiaalivgSiLgmdrklLikalvrylvpiligvvaalllgilvGlllGlsvkeallyivlPimaGGvgeGaiPlseiysevlgkdqeellsqlipavvlgnivAillagllnklgkkkpsltGnGkllkkkeeeellkeeeekekkvdlkklgaglllalalyllgtllekligv.ipalalmiilvvivkllglvpeeleegakklykfvskaltlallvgvGvaytdlkeliaaltlqnvviilvvVlgavlgaflvgklvglypieaaitagLcmanlGGtGdvavLsAanRmeLmpFAqistRlGGaivvilasl 415
#                      igg+plplfl+l+av+l+avl++klpkdl+galavllvlGillgeiG+riPilk+ylGggailallv+++lvy++llpeevvkavtt+mkksnfldlyiaali+gSiLgm+rklLikalvryl++il+gvv+alllgilvGll+G+svkea+lyivlPim+GG+g+Ga+Plseiys+vlg+++eel+sqlipa+++gnivAil+a+ll+klg+kkpsltGnG+l+k+keeeell+ee+ekekkvdlkklgaglllalal+llg+llekl+gv i++lalmiilv+i+k+lgl+p+eleegak+lykfvsk+ltlallvg+Gvaytdlkeliaal+l++vvi+l+vVlgavlga+lvgklvglypieaaitagLcmanlGGtGdvavLsAa+RmeLmpFAqis+RlGGaiv+ilasl
#   2HCT-consensus   1 IGGLPLPLFLVLAAVVLAAVLLEKLPKDLVGALAVLLVLGILLGEIGERIPILKEYLGGGAILALLVAALLVYFKLLPEEVVKAVTTLMKKSNFLDLYIAALITGSILGMNRKLLIKALVRYLPVILVGVVVALLLGILVGLLFGISVKEAVLYIVLPIMGGGIGAGAVPLSEIYSSVLGESSEELVSQLIPALTIGNIVAILVAALLKKLGEKKPSLTGNGELVKSKEEEELLEEEKEKEKKVDLKKLGAGLLLALALFLLGKLLEKLLGVkIHELALMIILVAILKALGLIPKELEEGAKQLYKFVSKNLTLALLVGIGVAYTDLKELIAALSLSYVVIVLAVVLGAVLGAALVGKLVGLYPIEAAITAGLCMANLGGTGDVAVLSAADRMELMPFAQISSRLGGAIVLILASL 416
#                      789***************************************************************************************************************************************************************************************************************************************************************************977********************************************************************************************************************************************985 PP


# Internal pipeline statistics summary:
# -------------------------------------
# Query model(s):                            1  (415 nodes)
# Target sequences:                          1  (416 residues searched)
# Passed MSV filter:                         1  (1); expected 1.0 (1)
# Passed bias filter:                        1  (1); expected 1.0 (1)
# Passed Vit filter:                         1  (1); expected 1.0 (1)
# Passed Fwd filter:                         1  (1); expected 1.0 (1)
# Initial search space (Z):                  1  [actual number of targets]
# Domain search space  (domZ):               1  [number of targets reported over threshold]
# # CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# # Mc/sec: 86.57
# //
# [ok]

AMINO_ALPH = "ACDEFGHIKLMNPQRSTVWY"
AMINO_RE_CLASS = rf"[{AMINO_ALPH}{AMINO_ALPH.lower()}\-\.]"

# Query:       2HCT  [M=415]
query_re = re.compile(r'^Query:\s+(\S+)')
# >> 2HCT-consensus
target_re = re.compile(r'^>>\s+(\S+)')


@dataclass
class HmmerAlignment:
    query: str
    target: str
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    query_line: str
    target_line: str

    def __str__(self) -> str:
        return f"{self.query} {self.query_start}..{self.query_end} | {self.target_start}..{self.target_end} {self.target} "


@dataclass
class HmmerOutput:
    alignments: List[HmmerAlignment] = field(default_factory=list)

    @staticmethod
    def parse(buf: str) -> "HmmerOutput":
        output = HmmerOutput()

        lines = buf.splitlines()

        m = next((query_re.search(s) for s in lines if query_re.search(s)), None)
        if m:
            query = m.group(1)
        else:
            print("no query found")
            exit()

        targets = [m.group(1) for s in lines if (m := target_re.search(s))]

        s = next((i for i, s in enumerate(lines) if "Alignments for each domain" in s), None)
        e = next((i for i, s in enumerate(lines) if "Internal pipeline statistics summary" in s), None)

        if s is None or e is None:
            print("failed to find alignments")
            exit()

        ali_lines = "\n".join([s for s in lines[s + 1:e] if s])

        for target in targets:
            res = re.search(
                rf"\s*{query}\s+(?P<q_start>\d+)\s+(?P<q_line>{AMINO_RE_CLASS}+)\s+(?P<q_end>\d+).*\n"
                rf"(?P<mid>.*)\n"
                rf"\s*{target}\s+(?P<t_start>\d+)\s+(?P<t_line>{AMINO_RE_CLASS}+)\s+(?P<t_end>\d+).*\n",
                ali_lines,
                re.MULTILINE,
            )
            output.alignments.append(
                HmmerAlignment(
                    query,
                    target,
                    int(res.group("q_start")),
                    int(res.group("q_end")),
                    int(res.group("t_start")),
                    int(res.group("t_end")),
                    res.group("q_line"),
                    res.group("t_line"),
                )
            )

        return output
