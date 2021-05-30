#  standard one-card record formats of ENSDF file. For further information
# please read the "Evaluated Nuclear Structure Data File - A Manual for Preparation 
# of Data Sets" by Jagdish K. Tuli (2001)

ID_positions = slice(5,9)

Time_scale = {
    'Y'     :   lambda x: x*3.154e7,
    'D'     :   lambda x: x*86400,
    'H'     :   lambda x: x*3600,
    'M'     :   lambda x: x*60,
    'S'     :   lambda x: x*1,
    'MS'    :   lambda x: x*1e-3,
    'US'    :   lambda x: x*1e-6,
    'NS'    :   lambda x: x*1e-9,
    'PS'    :   lambda x: x*1e-12,
    'FS'    :   lambda x: x*1e-15,
    'AS'    :   lambda x: x*1e-18,
    'EV'    :   lambda x: 6.582119569509e-16/x,
    'KEV'   :   lambda x: 6.582119569509e-19/x,
    'MEV'   :   lambda x: 6.582119569509e-22/x
}

ENSDF_Record = {
    "ID":{
        "ID_pa" :   "    ",
        "NUCID" :   slice(  0   , 5 ),
        "DSID"  :   slice(  9   , 39)
    },
    "Parent": {
        "ID_pa" :   "  P.",
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "J"     :   slice(  21  , 39),
        "T"     :   slice(  39  , 49),
        "DT"    :   slice(  49  , 55),
        "QP"    :   slice(  64  , 74),
        "DQP"   :   slice(  74  , 76),
        "ION"   :   slice(  76  , 80)
    },
    "Norm": {
        "ID_pa" :   "  N.",
        "NUCID" :   slice(  0   , 5 ),
        "NR"    :   slice(  9   , 19),
        "DNR"   :   slice(  19  , 21),
        "NT"    :   slice(  21  , 29),
        "DNT"   :   slice(  29  , 31),
        "BR"    :   slice(  31  , 39),
        "DBR"   :   slice(  39  , 41),
        "NB"    :   slice(  41  , 49),
        "DNB"   :   slice(  49  , 55),
        "NP"    :   slice(  55  , 62),
        "DNP"   :   slice(  62  , 63)
    },
    "Level": {
        "ID_pa" :   "  L ",
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "J"     :   slice(  21  , 39),
        "T"     :   slice(  39  , 49),
        "DT"    :   slice(  49  , 55),
        "L"     :   slice(  55  , 64),
        "S"     :   slice(  64  , 74),
        "DS"    :   slice(  74  , 76),
        "C"     :   slice(  76  , 77),
        "MS"    :   slice(  77  , 79),
        "Q"     :   slice(  79  , 80)
    },
    "Beta": {
        "ID_pa" :   "[ \w] B ",
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "I"    :   slice(  21  , 29),
        "DI"   :   slice(  29  , 31),
        "LOGFT" :   slice(  41  , 49),
        "DFT"   :   slice(  49  , 55),
        "C"     :   slice(  76  , 77),
        "UN"    :   slice(  77  , 79),
        "Q"     :   slice(  79  , 80)
    },
    "EC": {
        "ID_pa" :   "[ \w] E ",
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "IB"    :   slice(  19  , 29),
        "DIB"   :   slice(  29  , 31),
        "IE"    :   slice(  31  , 39),
        "DIE"   :   slice(  39  , 41),
        "LOGFT" :   slice(  41  , 49),
        "DFT"   :   slice(  49  , 55),
        "I"     :   slice(  64  , 74),
        "DI"    :   slice(  74  , 76),
        "C"     :   slice(  76  , 77),
        "UN"    :   slice(  77  , 79),
        "Q"     :   slice(  79  , 80)
    },
    "Alpha":{
        "ID_pa" :   "  A ",
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "I"    :   slice(  21  , 29),
        "DI"   :   slice(  29  , 31),
        "HF"    :   slice(  31  , 39),
        "DHF"   :   slice(  39  , 41),
        "C"     :   slice(  76  , 77),
        "Q"     :   slice(  79  , 80)
    },
    "Delayed":{
        "ID_pa" :   "[ \w] [ D][NPA]",
        "Particle": slice(  8   , 9 ),
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "I"    :   slice(  21  , 29),
        "DI"   :   slice(  29  , 31),
        "EI"    :   slice(  31  , 39),
        "T"     :   slice(  39  , 49),
        "DT"    :   slice(  49  , 55),
        "L"     :   slice(  55  , 64),
        "C"     :   slice(  76  , 77),
        "COIN"  :   slice(  77  , 78),
        "Q"     :   slice(  79  , 80)
    },
    "Gamma": {
        "ID_pa" :   "[ \w] G ",
        "NUCID" :   slice(  0   , 5 ),
        "E"     :   slice(  9   , 19),
        "DE"    :   slice(  19  , 21),
        "RI"    :   slice(  21  , 29),
        "DRI"   :   slice(  29  , 31),
        "M"     :   slice(  31  , 41),
        "MR"    :   slice(  41  , 49),
        "DMR"   :   slice(  49  , 55),
        "CC"    :   slice(  55  , 62),
        "DCC"   :   slice(  62  , 64),
        "TI"    :   slice(  64  , 74),
        "DTI"   :   slice(  74  , 76),
        "C"     :   slice(  76  , 77),
        "COIN"  :   slice(  77  , 78),
        "Q"     :   slice(  79  , 80)
    }
}