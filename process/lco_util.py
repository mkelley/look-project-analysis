from collections import defaultdict
import numpy as np
from astropy.stats import sigma_clip
from sbpy.activity import phase_HalleyMarcus

__all__ = [
    'assumed_gmr',
    'filters',
    'rename_target',
    'locations',
    'target_to_filename',
    'color_corrections'
]

# color for color correction
assumed_gmr = defaultdict(lambda: 0.55)
assumed_gmr['141P'] = 0.39
assumed_gmr['7P'] = 0.58
assumed_gmr['C/2020 R4'] = 0.48
assumed_gmr['C/2014 UN271'] = 0.48

# sigma clipped medians
color_corrections = {
    'gp': -0.0863,  # N=1008, stdev 0.035
    'rp': 0.0212,  # N=1147, stdev 0.018
    'ip': 0.0111,  # N=51, stdev 0.025
}

filters = {
    'rp': 'r',
    'gp': 'g',
    'air': 'r',
    'w': 'r',
    'R': 'r',
    'ip': 'i',
    'zs': 'z',
    'B': 'g',
    'V': 'g'
}

rename_target = {
    'C4AGJ62': 'C/2021 A1',
    'N00ht7m': 'C/2021 A2',
    '14F2601': '2021 HO',
    'C/2020 R4 ( ATLAS )': 'C/2020 R4',
    '7P Pons-Winnecke': '7P',
    '7P/Pons-Winnecke': '7P',
    'Comet C/2020 R4': 'C/2020 R4',
    'C/2020 R4 (ATLAS)': 'C/2020 R4',
    '2014 UN271': 'C/2014 UN271',
    'c/2014 UN271': 'C/2014 UN271',
    'P/2016 J3': '414P',
    'A/2019 U5': 'C/2019 U5'
}

target_names = {
    'C/2019 U5': 'C/2019 U5 (PanSTARRS)',
    'C/2014 UN271': 'C/2014 UN271 (Bernardinelli-Bernstein)',
    'C/2017 K2': 'C/2017 K2 (PanSTARRS)',
    'C/2018 F4': 'C/2018 F4 (PanSTARRS)',
    'C/2019 F1': 'C/2019 F1 (ATLAS-Africano)',
    'C/2019 L3': 'C/2019 L3 (ATLAS)',
    'C/2019 U6': 'C/2019 U6 (Lemmon)',
    'C/2020 O2': 'C/2020 O2 (Amaral)',
    'C/2020 Q1': 'C/2020 Q1 (Borisov)',
    'C/2020 R4': 'C/2020 R4 (ATLAS)',
    'C/2020 R7': 'C/2020 R7 (ATLAS)',
    'C/2020 T2': 'C/2020 T2 (Palomar)',
    'C/2020 U4': 'C/2020 U4 (PanSTARRS)',
    'C/2021 A1': 'C/2021 A1 (Leonard)',
    'C/2021 A2': 'C/2021 A2 (NEOWISE)',
    'C/2021 C4': 'C/2021 C4 (ATLAS)',
    'C/2021 C5': 'C/2021 C5 (PanSTARRS)',
    'C/2021 E3': 'C/2021 E3 (ZTF)',
    'P/2020 X1': 'P/2020 X1 (ATLAS)',
    '7P': '7P/Pons-Winnecke',
    '22P': '22P/Kopff',
    '29P': '29P/Schwassmann-Wachmann 1',
    '67P': '67P/Churyumov-Gerasimenko',
    '99P': '99P/Kowal 1',
    '117P': '117P/Helin-Roman-Alu 1',
    '141P': '141P/Machholz 2',
    '156P': '156P/Russell-LINEAR',
    '414P': '414P/STEREO',
    '2021 HO': '2021 HO'
}

locations = {  # approximate
    'ogg': 'F65',
    'elp': 'V37',
    'lsc': 'W86',
    'cpt': 'K92',
    'coj': 'Q63',
    'tfn': 'Z21',
    'tlv': '097'
}


def target_to_filename(target):
    return target.lower().replace('/', '').replace(' ', '')
