"""
BioScreen — toxin reference configuration.

Shared source of truth for:
  - functional region bounds (in full precursor coordinates)
  - literature-curated catalytic residues
  - short names and mechanism classes

Used by scripts/, the Studio, and the generation wrapper.
"""

TOXIN_CONFIG = {
    "P00588": {
        "name": "Diphtheria",
        "region": (0, 193),
        # Catalytic: H21 (mature) = H53 (precursor), E148 (mature) = E180
        # UniProt active site: 180(E). Verified session 5.
        "catalytic": [53, 180],
        "mechanism": "ADP-ribosylation",
        "functional_class": "enzymatic",
    },
    "P02879": {
        "name": "Ricin",
        "region": (35, 302),
        # Catalytic: Y80/Y123/E177/R180 (mature) = Y115/Y158/E212/R215
        # UniProt active sites: 115, 158, 212, 215. Verified session 5.
        "catalytic": [115, 158, 212, 215],
        "mechanism": "N-glycosidase",
        "functional_class": "enzymatic",
    },
    "P04958": {
        "name": "Tetanus",
        "region": (0, 450),
        # HExxH motif (H233/E234/H237) + E271 (3rd zinc ligand) + Y375.
        # UniProt: active site 234(E), binding 233(H), 237(H). Verified session 5.
        "catalytic": [233, 234, 237, 271, 375],
        "mechanism": "Zinc protease",
        "functional_class": "enzymatic",
    },
    "P09385": {
        "name": "Shiga 2A",
        "region": (0, 319),
        # Catalytic E167/R170 (mature) = E189/R192 (precursor).
        # UniProt active site: 189(E). Verified session 5.
        "catalytic": [189, 192],
        "mechanism": "rRNA depurination",
        "functional_class": "enzymatic",
    },
    "P0DPI0": {
        "name": "Botulinum E",
        "region": (0, 435),
        # HExxH motif at 223/224/227 in the sequence we have.
        # UniProt: active site 224(E), binding 223(H), 227(H). Verified session 5.
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
        "functional_class": "enzymatic",
    },
    "P0DPI1": {
        "name": "Botulinum A",
        "region": (0, 437),
        # Same HExxH positions as Botulinum E. Verified session 5.
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
        "functional_class": "enzymatic",
    },
    "P10149": {
        "name": "Shiga 1A",
        "region": (0, 315),
        # Catalytic E167/R170 (mature) = E189/R192 (precursor).
        # UniProt active site: 189(E). Verified session 5.
        "catalytic": [189, 192],
        "mechanism": "rRNA depurination",
        "functional_class": "enzymatic",
    },
    "P10844": {
        "name": "Botulinum B",
        "region": (0, 435),
        # HExxH motif: H230/E231/H234 (was incorrectly 223/224/227),
        # plus E268 as third zinc ligand.
        # UniProt: active site 231(E), binding 230(H), 234(H), 268(E).
        # Verified session 5.
        "catalytic": [230, 231, 234, 268],
        "mechanism": "Zinc protease",
        "functional_class": "enzymatic",
    },
    "P11140": {
        "name": "Abrin",
        "region": (0, 260),
        # Catalytic Y74/Y113/E164/R167. UniProt active site: 164(E).
        # Config was already correct. Verified session 5.
        "catalytic": [74, 113, 164, 167],
        "mechanism": "N-glycosidase",
        "functional_class": "enzymatic",
    },
    "P11439": {
        "name": "Pseudomonas ExoA",
        "region": (430, 638),
        # Catalytic domain III. Literature (mature numbering) -> UniProt (+25):
        #   H440 -> H465, Y470 -> Y495, Y481 -> Y506, E553 -> E578
        # E578 is UniProt-annotated as active site.
        "catalytic": [465, 495, 506, 578],
        "mechanism": "ADP-ribosylation",
        "functional_class": "enzymatic",
    },
    "P15917": {
        "name": "Anthrax LF",
        "region": (500, 809),
        # HExxH motif: H719/E720/H723, plus third zinc ligand E768.
        # UniProt: active site 720(E), binding 719(H), 723(H), 761(Y), 768(E).
        # Config was previously [686, 690, 728, 735] — wrong, off by ~33.
        # Verified session 5.
        "catalytic": [719, 720, 723, 768],
        "mechanism": "Zinc protease",
        "functional_class": "enzymatic",
    },
    "P33183": {
        "name": "Modeccin",
        "region": (0, 260),
        # UNRESOLVED. UniProt shows only active site 188(E). The current
        # positions [80, 123, 177, 180] point at residues V/K/R/L — likely
        # wrong, but no reliable literature source gives the precise set.
        # Modeccin is an N-glycosidase (same class as ricin/abrin which
        # use YYER), so the expected residues are probably tyrosine-tyrosine-
        # glutamate-arginine. Do not guess without verification.
        "catalytic": [80, 123, 177, 180],
        "mechanism": "N-glycosidase",
        "functional_class": "enzymatic",
        "unresolved": True,
    },
    "P40136": {
        "name": "Anthrax EF",
        "region": (0, 500),
        # UniProt active site: 351(H). Config was already correct.
        # Verified session 5.
        "catalytic": [351],
        "mechanism": "Adenylate cyclase",
        "functional_class": "enzymatic",
    },
    "P01552": {
        "name": "SEB",
        "region": (27, 266),
        # Functional positions discovered by unsupervised conservation
        # across 13 superantigen sequences (Clustal Omega alignment +
        # surface-residue classification). Literature cross-check:
        # TSST-1 Y144 (literature) = position 144 (matches).
        # Residue pattern (mature numbering): Y-D-D-T-K-T-Q-E-D-R-Y-D-Y-N
        "catalytic": [55, 106, 110, 146, 179, 183, 185, 186, 188, 192, 202, 226, 244, 247],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P06886": {
        "name": "TSST-1",
        "region": (40, 234),
        # Functional positions discovered by unsupervised conservation
        # across 13 superantigen sequences (Clustal Omega alignment +
        # surface-residue classification). Literature cross-check:
        # TSST-1 Y144 (literature) = position 144 (matches).
        # Residue pattern (mature numbering): Y-D-D-T-K-T-Q-E-D-R-Y-D-Y-N
        "catalytic": [53, 81, 85, 127, 161, 165, 167, 168, 170, 174, 184, 207, 215, 218],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L5": {
        "name": "SEC3",
        "region": (27, 266),
        # Functional positions discovered by unsupervised conservation
        # across 13 superantigen sequences (Clustal Omega alignment +
        # surface-residue classification). Literature cross-check:
        # TSST-1 Y144 (literature) = position 144 (matches).
        # Residue pattern (mature numbering): Y-D-D-T-K-T-Q-E-D-R-Y-D-Y-N
        "catalytic": [55, 106, 110, 143, 178, 182, 184, 185, 187, 191, 201, 226, 244, 247],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0C0I5": {
        "name": "SPEA",
        "region": (27, 235),
        # Functional positions discovered by unsupervised conservation
        # across 13 superantigen sequences (Clustal Omega alignment +
        # surface-residue classification). Literature cross-check:
        # TSST-1 Y144 (literature) = position 144 (matches).
        # Residue pattern (mature numbering): Y-D-D-T-K-T-Q-E-D-R-Y-D-Y-N
        "catalytic": [44, 93, 97, 118, 151, 155, 157, 158, 160, 164, 174, 198, 215, 218],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L2": {
        "name": "SEA",
        "region": (24, 257),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YGDTKTQEDRYDYN
        "catalytic": [54, 106, 110, 136, 171, 175, 177, 178, 180, 184, 194, 221, 237, 240],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L3": {
        "name": "SED",
        "region": (27, 266),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [55, 106, 110, 143, 178, 182, 184, 185, 187, 191, 201, 226, 244, 247],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L4": {
        "name": "SEE",
        "region": (27, 266),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [55, 106, 110, 143, 178, 182, 184, 185, 187, 191, 201, 226, 244, 247],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L6": {
        "name": "SEG",
        "region": (25, 258),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [54, 102, 106, 139, 171, 175, 177, 178, 180, 184, 194, 219, 238, 241],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L7": {
        "name": "SEH",
        "region": (25, 258),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [54, 102, 106, 139, 171, 175, 177, 178, 180, 184, 194, 219, 238, 241],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0A0L8": {
        "name": "SEI",
        "region": (25, 258),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [54, 102, 106, 139, 171, 175, 177, 178, 180, 184, 194, 219, 238, 241],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0C0I3": {
        "name": "SPEC",
        "region": (33, 571),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: LDIGDTKEDIFENN
        "catalytic": [116, 162, 166, 210, 317, 321, 323, 324, 408, 412, 422, 478, 519, 522],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0C0I6": {
        "name": "SPEG",
        "region": (32, 236),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [47, 95, 99, 123, 154, 158, 160, 161, 163, 167, 177, 202, 218, 221],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0C0I7": {
        "name": "SPEJ",
        "region": (32, 236),
        # Functional positions from unsupervised conservation across
        # 13 superantigen sequences. Residues: YDDTKTQEDRYDYN
        "catalytic": [47, 95, 99, 123, 154, 158, 160, 161, 163, 167, 177, 202, 218, 221],
        "mechanism": "superantigen-binding",
        "functional_class": "superantigen",
        "source": "unsupervised-conservation",
    },
    "P0C0I3": {
        "name": "SLO",
        "region": (33, 571),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFGKQYTGEPYGRKTSSF...
        "catalytic": [138, 149, 152, 161, 173, 174, 186, 187, 235, 239, 251, 258, 261, 281, 287, 297, 298, 301, 320, 328, 332, 334, 343, 344, 345, 350, 353, 355, 357, 364, 385, 394, 395, 426, 432, 433, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "P0C2E9": {
        "name": "PFO",
        "region": (28, 500),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFNKQYTGEPYGRKTSSF...
        "catalytic": [67, 78, 81, 90, 102, 103, 115, 116, 164, 168, 180, 187, 190, 210, 216, 226, 227, 230, 249, 257, 261, 263, 272, 273, 274, 279, 282, 284, 286, 293, 314, 323, 324, 355, 361, 362, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "P13128": {
        "name": "LLO",
        "region": (24, 529),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFGKQYTGEPYGRKTSSF...
        "catalytic": [90, 101, 104, 113, 125, 126, 138, 139, 188, 192, 205, 212, 215, 235, 241, 251, 252, 255, 274, 282, 286, 288, 297, 298, 299, 304, 307, 309, 311, 318, 339, 348, 349, 380, 386, 387, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "Q04IN8": {
        "name": "Pneumolysin",
        "region": (1, 471),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFGKQYTGEPYGRKTSSF...
        "catalytic": [35, 46, 49, 58, 70, 71, 83, 84, 133, 137, 149, 156, 159, 179, 185, 195, 196, 199, 218, 226, 230, 232, 241, 242, 243, 248, 251, 253, 255, 262, 283, 292, 293, 324, 330, 331, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "P0C2J9": {
        "name": "Pneumolysin-v",
        "region": (1, 471),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFGKQYTGEPYGRKTSSF...
        "catalytic": [35, 46, 49, 58, 70, 71, 83, 84, 133, 137, 149, 156, 159, 179, 185, 195, 196, 199, 218, 226, 230, 232, 241, 242, 243, 248, 251, 253, 255, 262, 283, 292, 293, 324, 330, 331, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "Q8P2T7": {
        "name": "SLO-M18",
        "region": (33, 571),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFGKQYTGEPYGRKTSSF...
        "catalytic": [138, 149, 152, 161, 173, 174, 186, 187, 235, 239, 251, 258, 261, 281, 287, 297, 298, 301, 320, 328, 332, 334, 343, 344, 345, 350, 353, 355, 357, 364, 385, 394, 395, 426, 432, 433, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "C1KYD5": {
        "name": "LLO-4b",
        "region": (24, 529),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 47 positions, residues GEKDYPENWYYSQFGKQYTGEPYGRKTSSF...
        "catalytic": [90, 101, 104, 113, 125, 126, 138, 139, 188, 192, 205, 212, 215, 235, 241, 251, 252, 255, 274, 282, 286, 288, 297, 298, 299, 304, 307, 309, 311, 318, 339, 348, 349, 380, 386, 387, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
    "P19265": {
        "name": "ILY",
        "region": (0, 298),
        # CDC pore-forming family. Functional positions from unsupervised
        # conservation across 8 CDC members (7/8 for undecapeptide).
        # Undecapeptide: ECTGLAWEWWR (membrane-inserting motif).
        # 36 positions, residues GEHRGPADWQILDVGEEPLGGPYGRDGSSY...
        "catalytic": [50, 61, 64, 73, 86, 87, 99, 100, 116, 120, 129, 136, 140, 163, 170, 181, 182, 185, 187, 195, 199, 201, 206, 207, 208, 213, 216, 237, 239, 246, 258, 267, 268, 291, 295, 296],
        "mechanism": "pore-formation",
        "functional_class": "pore_forming",
        "source": "unsupervised-conservation",
    },
}

