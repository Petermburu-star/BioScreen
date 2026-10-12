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
    },
    "P02879": {
        "name": "Ricin",
        "region": (35, 302),
        # Catalytic: Y80/Y123/E177/R180 (mature) = Y115/Y158/E212/R215
        # UniProt active sites: 115, 158, 212, 215. Verified session 5.
        "catalytic": [115, 158, 212, 215],
        "mechanism": "N-glycosidase",
    },
    "P04958": {
        "name": "Tetanus",
        "region": (0, 450),
        # HExxH motif (H233/E234/H237) + E271 (3rd zinc ligand) + Y375.
        # UniProt: active site 234(E), binding 233(H), 237(H). Verified session 5.
        "catalytic": [233, 234, 237, 271, 375],
        "mechanism": "Zinc protease",
    },
    "P09385": {
        "name": "Shiga 2A",
        "region": (0, 319),
        # Catalytic E167/R170 (mature) = E189/R192 (precursor).
        # UniProt active site: 189(E). Verified session 5.
        "catalytic": [189, 192],
        "mechanism": "rRNA depurination",
    },
    "P0DPI0": {
        "name": "Botulinum E",
        "region": (0, 435),
        # HExxH motif at 223/224/227 in the sequence we have.
        # UniProt: active site 224(E), binding 223(H), 227(H). Verified session 5.
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
    },
    "P0DPI1": {
        "name": "Botulinum A",
        "region": (0, 437),
        # Same HExxH positions as Botulinum E. Verified session 5.
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
    },
    "P10149": {
        "name": "Shiga 1A",
        "region": (0, 315),
        # Catalytic E167/R170 (mature) = E189/R192 (precursor).
        # UniProt active site: 189(E). Verified session 5.
        "catalytic": [189, 192],
        "mechanism": "rRNA depurination",
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
    },
    "P11140": {
        "name": "Abrin",
        "region": (0, 260),
        # Catalytic Y74/Y113/E164/R167. UniProt active site: 164(E).
        # Config was already correct. Verified session 5.
        "catalytic": [74, 113, 164, 167],
        "mechanism": "N-glycosidase",
    },
    "P11439": {
        "name": "Pseudomonas ExoA",
        "region": (430, 638),
        # Catalytic domain III. Literature (mature numbering) -> UniProt (+25):
        #   H440 -> H465, Y470 -> Y495, Y481 -> Y506, E553 -> E578
        # E578 is UniProt-annotated as active site.
        "catalytic": [465, 495, 506, 578],
        "mechanism": "ADP-ribosylation",
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
        "unresolved": True,
    },
    "P40136": {
        "name": "Anthrax EF",
        "region": (0, 500),
        # UniProt active site: 351(H). Config was already correct.
        # Verified session 5.
        "catalytic": [351],
        "mechanism": "Adenylate cyclase",
    },
}

