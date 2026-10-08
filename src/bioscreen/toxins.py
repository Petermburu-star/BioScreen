"""
BioScreen — toxin reference configuration.

Shared source of truth for:
  - functional region bounds (in full precursor coordinates)
  - literature-curated catalytic residues
  - short names and mechanism classes

Used by scripts/, the Studio, and the generation wrapper.
"""

TOXIN_CONFIG = {
    "P02879": {
        "name": "Ricin",
        "region": (35, 302),
        "catalytic": [80, 123, 177, 180],
        "mechanism": "N-glycosidase",
    },
    "P0DPI1": {
        "name": "Botulinum A",
        "region": (0, 437),
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
    },
    "P10844": {
        "name": "Botulinum B",
        "region": (0, 435),
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
    },
    "P0DPI0": {
        "name": "Botulinum E",
        "region": (0, 435),
        "catalytic": [223, 224, 227],
        "mechanism": "Zinc protease",
    },
    "P10149": {
        "name": "Shiga 1A",
        "region": (0, 315),
        "catalytic": [167, 170],
        "mechanism": "rRNA depurination",
    },
    "P09385": {
        "name": "Shiga 2A",
        "region": (0, 319),
        "catalytic": [167, 170],
        "mechanism": "rRNA depurination",
    },
    "P00588": {
        "name": "Diphtheria",
        "region": (0, 193),
        "catalytic": [21, 148],
        "mechanism": "ADP-ribosylation",
    },
    "P11140": {
        "name": "Abrin",
        "region": (0, 260),
        "catalytic": [74, 113, 164, 167],
        "mechanism": "N-glycosidase",
    },
    "P33183": {
        "name": "Modeccin",
        "region": (0, 260),
        "catalytic": [80, 123, 177, 180],
        "mechanism": "N-glycosidase",
    },
    "P15917": {
        "name": "Anthrax LF",
        "region": (500, 809),
        "catalytic": [686, 690, 728, 735],
        "mechanism": "Zinc protease",
    },
    "P40136": {
        "name": "Anthrax EF",
        "region": (0, 500),
        "catalytic": [351],
        "mechanism": "Adenylate cyclase",
    },
    "P04958": {
        "name": "Tetanus",
        "region": (0, 450),
        # Best-guess set; UniProt active-site verification pending
        "catalytic": [233, 234, 237, 271, 375],
        "mechanism": "Zinc protease",
    },
}
