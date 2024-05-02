names_table = {
    "methane":                                 "methane",
    "ethane":                                  "ethane",
    "propane":                                 "propane",
    "butane":                                  "butane",
    "n-pentane":                               "pentane",
    "n-hexane":                                "hexane",
    "n-heptane":                               "heptane",
    "n-octane":                                "octane",
    "n-nonane":                                "nonane",
    "n-decane":                                "decane",
    "n-undecane":                              "c11",
    "n-dodecane":                              "c12",
    "n-tridecane":                             "c13",
    "n-tetradecane":                           "c14",
    "n-pentadecane":                           "c15",
    "n-hexadecane":                            "c16",
    "n-heptadecane":                           "c17",
    "n-octadecane":                            "c18",
    "n-nonadecane":                            "c19",
    "n-eicosane":                              "c20",
    "n-octacosane":                            "octacosane",
    "n-hexatriacontane":                       None,
    "isobutane":                               "isobutan",
    "isopentane":                              "ipentane",
    "neopentane":                              "neopentane", 
    "2,2-dimethylbutane":                      "neohexane",
    "2,3-dimethylbutane":                      "tetramethylethan",
    "2-methylpentane":                         "ihexane",
    "3-methylpentane":                         "3-methylpentane",
    "2-methylhexane":                          "isoheptane",
    "ethylene":                                "ethylene",
    "propylene":                               "propylen",
    "1-butene":                                "1butene",
    "1-pentene":                               "1pentene",
    "1-hexene":                                "1hexene",
    "1-octene":                                "1octene",
    "cyclopentane":                            "cyclopen",
    "cyclohexane":                             "cyclohex",
    "methylcyclopentane":                      None,
    "methylcyclohexane":                       None,
    "benzene":                                 "benzene",
    "benzene NP":                              None,
    "methylbenzene":                           "toluene",
    "ethylbenzene":                            "ebenzene",
    "m-xylene":                                "mxylene",
    "n-propylbenzene":                         None,
    "n-butylbenzene":                          None,
    "tetralin":                                None,
    "dichlorodifluoromethane":                 "r12",
    "chlorotrifluoromethane":                  "r13",
    "chlorodifluoromethane":                   "r22",
    "chlorodifluoromethane NP":                None,
    "trifluoromethane":                        "r23",
    "trifluoromethane NP":                     None,
    "difluoromethane":                         "r32",
    "difluoromethane NP":                      None,
    "1,1,2-trichloro-1,2,2-trifluoroethane":  "r113",
    "1,2-dichloro-1,1,2,2-tetrafluoroethane": "r114",
    "1,1-dichloro-2,2,2-trifluoroethane":     "r123",
    "1-chloro-1,2,2,2-tetrafluoroethane":     "r124",
    "pentafluoroethane":                      "r125",
    "pentafluoroethane NP":                   None,
    "1,1,1,2-tetrafluoroethane":              "r134a",
    "1,1,1,2-tetrafluoroethane NP":           None,
    "tetrachloromethane":                     None,
    "chloroethane":                           None, # <-- ended conversion of the names here 04.04.2024
    "isopropyl-chloride":                     None,
    "methyl chloride":                        None,
    "chlorobenzene":                          "cbenzene",
    "brombenzene":                            None,
    "1-chlorobutane":                        "1-chlorobutane",
    "perfluoromethane":                      "r14",
    "perfluoroethane":                       "r116",
    "perfluoropropane":                      "r218",
    "perfluorobutane":                       "C4F10",
    "perfluoropentane":                      "C5F12",
    "perfluorohexane":                       "C6F14",
    "perfluoroheptane":                      "C7F16",
    "perfluorooctane":                       "C8F18",
    "eicosafluorononane":                    "C9F20",
    "dimethyl-ether":                        "dme",
    "methyl-ethyl-ether":                    None,
    "diethyl-ether":                         "dee",
    "methyl-n-propyl-ether":                 "methyl-n-propyl_ether",
    "butanone":                              None, # even though it is available in saft fluid file
    "pentanone-2":                           "2-pentanone",
    "pentanone-3":                           "3-pentanone",
    "argon":                                 "argon",
    "carbon monoxide":                       "co",
    "chlorine":                              "chlorine",
    "carbon dioxide":                        "co2",
    "carbon dioxide NP":                     None,
    "sulfur dioxide":                        "so2",
    "nitrogen":                              "nitrogen",
    "carbon disulfide":                      "carbon_disulfide",
    "oxygen":                                None,
    "xenon":                                 None,
    "methanol NP":                           None,
    "methanol":                              "methanol",
    "ethanol NP":                            None,
    "ethanol":                               "ethanol",
    "1-propanol NP":                         None,
    "1-propanol":                            "1-propanol",
    "1-butanol NP":                          None,
    "1-butanol":                             "1-butanol",
    "1-pentanol NP":                         None,
    "1-pentanol":                            "1-pentanol",
    "1-hexanol":                             "1-hexanol",
    "1-heptanol":                            "1-heptanol",
    "1-octanol":                             "1-octanol",
    "1-nonanol":                             "1-nonanol",
    "2-propanol NP":                         None,
    "2-propanol":                            "2-propanol",
    "2-methyl-2-butanol":                    "2-butanol",
    "water":                                 "water",
    "water4c":                               "water4c",
    "methylamine":                           "methylamine",
    "ethylamine":                            "ethylamine",
    "1-propylamine":                         "1-propylamine",
    "2-propylamine":                         "2-propylamine",
    "aniline":                               "aniline",
    "acetic acid":                           "acetic_acid",
    "butyl ethanoate":                       "butyl_ethanoate",
    "acetone":                               "acetone",
    "sulfur hexafluoride":                   "SF6",
    "HFE-7000AA":                            "hfe-7000aa",
    "HFE-7100AA":                            "hfe-7100aa",
    "HFE-7200AA":                            "hfe-7200aa",
    "HFE-7300AA":                            "hfe-7300aa",
    "HFE-7500AA":                            "hfe-7500aa",
    "HFE-7000VV":                            "hfe-7000vv",
    "HFE-7100VV":                            "hfe-7100vv",
    "HFE-7200VV":                            "hfe-7200vv",
    "HFE-7300VV":                            "hfe-7300vv",
    "HFE-7500VV":                            "hfe-7500vv",
    "HFE-7000VVD":                           "hfe-7000vvd",
    "HFE-7100VVD":                           "hfe-7100vvd",
    "HFE-7200VVD":                           "hfe-7200vvd",
    "HFE-7300VVD":                           "hfe-7300vvd",
    "HFE-7500VVD":                           "hfe-7500vvd",
    "davidium32":                            "davidium32",
    "davidium42":                            "davidium42",
}

inv_names_table = {v: k for k, v in names_table.items()}

def convert2trend_name(name):
    # convert the pc_saft name to trend name
    try:
        # print(f"{name} -> {names_table[name]}") # DEBUG trasnformation of the name
        if name in inv_names_table.keys(): # if it is trend name
            return name 
        else: # else it is pcsaft name
            return names_table[name]
    except KeyError as ke:
        print(ke)
        return None

def convert2pcsaft_name(name):
    # convert the trend name to pc_saft name
    try:
        # print(f"{name} -> {names_table[name]}") # DEBUG trasnformation of the name
        if name in inv_names_table.keys(): # if it is trend name
            return inv_names_table[name] # it is actually a trend name already
        else: # else it is pcsaft name
            return name
    except KeyError as ke:
        print(ke)
        return None

def is_valid_subst_name(names, is_global=False):
    # print(name) #DEBUG
    res = []
    for name in names:
        if name in names_table.keys():
            res.append(True)
        else:
            res.append(False)
    if is_global:
        return all(res)
    else:
        return res