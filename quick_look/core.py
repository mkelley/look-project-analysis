import logging
from collections import defaultdict

rho_arcsec = [2, 5, 10, 12, 20]
rho_km = [5e3, 1e4, 2e4]
rho_labels = [str(r) for r in rho_arcsec] + [
    "{}k".format(str(int(r // 1000))) for r in rho_km
]
color_rho = "5"

# color for color correction
assumed_gmr = defaultdict(lambda: 0.55)
assumed_gmr["141P"] = 0.39
assumed_gmr["7P"] = 0.58
assumed_gmr["C/2020 R4"] = 0.48
assumed_gmr["C/2014 UN271"] = 0.48

# sigma clipped means
color_corrections = {
    "1m0a": {
        "gp": -0.0914,  # N=13933, stdev 0.0333
        "rp": 0.0197,  # N=14249, stdev 0.0151
        "ip": 0.0260,  # N=888, stdev 0.0162
    },
    "2m0a": {
        "gp": -0.0279,  # 418, 0.0206
        "rp": 0.0158,  # 438, 0.0175
        "ip": 0.0159,  # 408, 0.0194
        "zs": -0.0037,  # 364, 0.0170
    },
}

filters = {
    "rp": "r",
    "gp": "g",
    "air": "r",
    "w": "r",
    "R": "r",
    "ip": "i",
    "zs": "z",
    "B": "g",
    "V": "g",
    "solar": "r",
}

# Place PCCP objects here until they get assigned a new designation
skip_targets = {
    "P21D3Eu",
    "CA9MP42",
}

rename_target = {
    "A10X1pX": "C/2023 S2",
    "C4AGJ62": "C/2021 A1",
    "N00ht7m": "C/2021 A2",
    "P21JDKb": "A/2023 R3",
    "14F2601": "2021 HO",
    "C/2020 R4 ( ATLAS )": "C/2020 R4",
    "7P Pons-Winnecke": "7P",
    "7P/Pons-Winnecke": "7P",
    "Comet C/2020 R4": "C/2020 R4",
    "C/2020 R4 (ATLAS)": "C/2020 R4",
    "2014 UN271": "C/2014 UN271",
    "c/2014 UN271": "C/2014 UN271",
    "c/2014 Un271": "C/2014 UN271",
    "c/2014 un271": "C/2014 UN271",
    "C/2014 Un271": "C/2014 UN271",
    "P/2016 J3": "414P",
    "A/2019 U5": "C/2019 U5",
    "C/2018 F4": "C/2018 F4-A",
    "A/2021 G2": "C/2021 G2",
    "73P-JD001": "73P-BU",
    "73P-JD001.5": "73P-BU",
    "73p": "73P",
    "RJ001": "73P-BY",
    "364p": "364P",
    "12p": "12P",
    "12P/Pons-Brooks": "12P",
    "chiron": "2060",
    "Chiron": "2060",
    "144p": "144P",
    "144P/Kushida": "144P",
    "62p": "62P",
}

target_names = {
    "C/2019 U5": "C/2019 U5 (PANSTARRS)",
    "C/2014 UN271": "C/2014 UN271 (Bernardinelli-Bernstein)",
    "C/2017 K2": "C/2017 K2 (PANSTARRS)",
    "C/2018 A6": "C/2018 A6 (Gibbs)",
    "C/2018 F4": "C/2018 F4 (PANSTARRS)",
    "C/2018 F4-A": "C/2018 F4-A (PANSTARRS)",
    "C/2018 F4-B": "C/2018 F4-B (PANSTARRS)",
    "C/2019 F1": "C/2019 F1 (ATLAS-Africano)",
    "C/2019 L3": "C/2019 L3 (ATLAS)",
    "C/2019 U6": "C/2019 U6 (Lemmon)",
    "C/2020 O2": "C/2020 O2 (Amaral)",
    "C/2020 Q1": "C/2020 Q1 (Borisov)",
    "C/2020 R4": "C/2020 R4 (ATLAS)",
    "C/2020 R7": "C/2020 R7 (ATLAS)",
    "C/2020 T2": "C/2020 T2 (Palomar)",
    "C/2020 U4": "C/2020 U4 (PANSTARRS)",
    "C/2021 A1": "C/2021 A1 (Leonard)",
    "C/2021 A2": "C/2021 A2 (NEOWISE)",
    "C/2021 C4": "C/2021 C4 (ATLAS)",
    "C/2021 C5": "C/2021 C5 (PANSTARRS)",
    "C/2021 E3": "C/2021 E3 (ZTF)",
    "C/2021 G2": "C/2021 G2 (ATLAS)",
    "C/2021 K2": "C/2021 K2 (MASTER)",
    "C/2021 K3": "C/2021 K3 (Catalina)",
    "C/2021 O3": "C/2021 O3 (PANSTARRS)",
    "C/2021 P4": "C/2021 P4 (ATLAS)",
    "C/2021 Q4": "C/2021 Q4 (Fuls)",
    "C/2021 R2": "C/2021 R2 (PANSTARRS)",
    "C/2021 S3": "C/2021 S3 (PANSTARRS)",
    "C/2021 T4": "C/2021 T4 (Lemmon)",
    "C/2021 Y1": "C/2021 Y1 (ATLAS)",
    "C/2022 E3": "C/2022 E3 (ZTF)",
    "C/2023 S2": "C/2023 S2 (ATLAS)",
    "P/2020 X1": "P/2020 X1 (ATLAS)",
    "7P": "7P/Pons-Winnecke",
    "12P": "12P/Pons-Brooks",
    "19P": "19P/Borrelly",
    "22P": "22P/Kopff",
    "29P": "29P/Schwassmann-Wachmann 1",
    "44P": "44P/Reinmuth 2",
    "57P": "57P/duToit-Neujmin-Delporte",
    "62P": "62P/Tsuchinshan 1",
    "67P": "67P/Churyumov-Gerasimenko",
    "73P": "73P/Schwassmann Wachmann 3",
    "73P-BU": "73P-BU/Schwassmann Wachmann 3 BU",
    "73P-BV": "73P-BV/Schwassmann Wachmann 3 BV",
    "73P-BY": "73P-BY/Schwassmann Wachmann 3 BY",
    "81P": "81P/Wild 2",
    "97P": "Metcalf-Brewington",
    "99P": "99P/Kowal 1",
    "108P": "108P/Ciffreo",
    "117P": "117P/Helin-Roman-Alu 1",
    "120P": "120P/Mueller 1",
    "141P": "141P/Machholz 2",
    "144P": "144P/Kushida",
    "156P": "156P/Russell-LINEAR",
    "157P": "157P/Tritton",
    "169P": "169P/NEAT",
    "191P": "191P/McNaught",
    "285P": "285P/LINEAR",
    "382P": "382P/Larson",
    "364P": "364P/PANSTARRS",
    "395P": "395P/Catalina-NEAT",
    "414P": "414P/STEREO",
    "445P": "445P/Lemmon-PANSTARRS",
    "2005 QN173": "248370 (2005 QN173)",
    "2021 HO": "2021 HO",
    "2060": "95P/Chiron",
}

locations = {  # approximate
    "ogg": "F65",
    "elp": "V37",
    "lsc": "W86",
    "cpt": "K92",
    "coj": "Q63",
    "tfn": "Z21",
    "tlv": "097",
}


def target_to_filename(target):
    return target.lower().replace("/", "").replace(" ", "").replace("-", "")


def setup_logger(name, level):
    logger = logging.getLogger(name)
    [logger.removeHandler(h) for h in list(logger.handlers)]
    logger.setLevel("DEBUG" if level == "DEBUG" else "INFO")
    logger.addHandler(logging.FileHandler("look-project-analysis.log"))
    h = logging.StreamHandler()
    h.setLevel("INFO" if level is None else level)
    logger.addHandler(h)
    formatter = logging.Formatter(
        "%(name)s : %(levelname)s : %(asctime)s : %(message)s"
    )
    for h in logger.handlers:
        h.setFormatter(formatter)
    return logger
