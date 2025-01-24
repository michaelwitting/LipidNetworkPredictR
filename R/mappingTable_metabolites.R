#' @name mappingTable
#' 
#' @title Create a mapping table to map between LipidNetworkPredictR and RHEA metabolite names
#' 
#' @description 
#' The function returns a \code{data.frame} containing information of 
#' corresponding features from \code{LipidNetworkPredictR}, the 
#' \code{RHEA} database and \code{ChEBI} database.
#' 
#' @details 
#' The mapping table is derived from the file 
#' `rhea-reactions.txt`.
#' 
#' @return data.frame
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples 
#' LipidNetworkPredictR:::mappingTable()
mappingTable <- function() {
    
    ## name_LipidNetworkPredictR, name_RHEA, ChEBI_ID, RHEA_ID_undetermined, RHEA_ID_left_to_right, RHEA_ID_right_to_left, RHEA_ID_both_directions
    res <- rbind(
    ## acyldhap_to_alkyldhap
    c("M_AcylDHAP",  "a 1-acylglycerone 3-phosphate",  "CHEBI:57534", "RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174"),
    c("M_FAO",       "a long chain fatty alcohol",     "CHEBI:17135", "RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174"),
    c("M_H+",        "H(+)",                           "CHEBI:15378", "RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174"),
    c("M_FA",        "a long-chain fatty acid",        "CHEBI:57560", "RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174"),
    c("M_AlkylDHAP", "1-O-alkylglycerone 3-phosphate", "CHEBI:73315", "RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174"),
        
    ## alkyldhap_to_lpao
    c("M_H+",        "H(+)",                             "CHEBI:15378", "RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178"),
    c("M_NADPH",     "NADPH",                            "CHEBI:57783", "RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178"),
    c("M_AlkylDHAP", "1-O-alkylglycerone 3-phosphate",   "CHEBI:73315", "RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178"),
    c("M_LPA-O",     "1-O-alkyl-sn-glycero-3-phosphate", "CHEBI:58014", "RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178"),
    c("M_NADP",      "NADP(+)",                          "CHEBI:58349", "RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178"),

    ## cerp_to_cer
    c("M_H2O",  "H2O",                                 "CHEBI:15377", "RHEA:33743", "RHEA:33744", "RHEA:33745", "RHEA:33746"),
    c("M_CerP", "an N-acylsphing-4-enine 1-phosphate", "CHEBI:57674", "RHEA:33743", "RHEA:33744", "RHEA:33745", "RHEA:33746"),
    c("M_Pi",   "phosphate",                           "CHEBI:43474", "RHEA:33743", "RHEA:33744", "RHEA:33745", "RHEA:33746"),
    c("M_Cer",  "an N-acylsphing-4-enine",             "CHEBI:52639", "RHEA:33743", "RHEA:33744", "RHEA:33745", "RHEA:33746"),
    
    ## cdpdg_to_pgp
    c("M_Glycerol-3-P", "sn-glycerol 3-phosphate",                                      "CHEBI:57597", "RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596"),
    c("M_CDP-DG",       "a CDP-1,2-diacyl-sn-glycerol",                                 "CHEBI:58332", "RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596"),
    c("M_H+",           "H(+)",                                                         "CHEBI:15378", "RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596"),
    c("M_PGP",          "1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycero-3'-phosphate)", "CHEBI:60110", "RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596"),
    c("M_CMP",          "CMP",                                                          "CHEBI:60377", "RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596"),

    ## cdpdg_to_pi
    c("M_CDP-DG",       "a CDP-1,2-diacyl-sn-glycerol",                        "CHEBI:58332", "RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583"),
    c("M_myo-Inositol", "myo-inositol",                                        "CHEBI:17268", "RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583"),
    c("M_PI",           "a 1,2-diacyl-sn-glycero-3-phospho-(1D-myo-inositol)", "CHEBI:57880", "RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583"),
    c("M_CMP",          "CMP",                                                 "CHEBI:60377", "RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583"),
    c("M_H+",           "H(+)",                                                "CHEBI:15378", "RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583"),

    ## cer_to_cerp
    c("M_Cer",  "an N-acylsphing-4-enine",             "CHEBI:52639", "RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932"),
    c("M_ATP",  "ATP",                                 "CHEBI:30616", "RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932"),
    c("M_ADP",  "ADP",                                 "CHEBI:456216", "RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932"),
    c("M_CerP", "an N-acylsphing-4-enine 1-phosphate", "CHEBI:57674", "RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932"),
    c("M_H+",   "H(+)",                                "CHEBI:15378", "RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932"),
    
    ## cer_to_glccer
    c("M_Cer",         "an N-acylsphing-4-enine",                         "CHEBI:52639", "RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091"),
    c("M_UDP-Glucose", "UDP-alpha-D-glucose",                             "CHEBI:58885", "RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091"),
    c("M_GlcCer",      "a beta-D-glucosyl-(1<->1')-N-acylsphing-4-enine", "CHEBI:22801", "RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091"),
    c("M_H+",          "H(+)",                                            "CHEBI:15378", "RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091"),
    c("M_UDP",         "UDP",                                             "CHEBI:58223", "RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091"),
    
    ## cer_to_sm
    c("M_PC",  "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:18765", "RHEA:18766", "RHEA:18767", "RHEA:18768"),
    c("M_Cer", "an N-acylsphing-4-enine",                  "CHEBI:52639", "RHEA:18765", "RHEA:18766", "RHEA:18767", "RHEA:18768"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol",              "CHEBI:17815", "RHEA:18765", "RHEA:18766", "RHEA:18767", "RHEA:18768"),
    c("M_SM", "a sphingomyelin",                           "CHEBI:17636", "RHEA:18765", "RHEA:18766", "RHEA:18767", "RHEA:18768"),
        
    ## cl_to_lcl
    c("M_CL",        "a cardiolipin",                                                                  "CHEBI:62237", "RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"),
    c("M_H2O",       "H2O",                                                                            "CHEBI:15377", "RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"),
    c("M_1,2,4-LCL", "1'-[1,2-diacyl-sn-glycero-3-phospho],3'-[1-acyl-sn-glycero-3-phospho]-glycerol", "CHEBI:64743", "RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"),
    c("M_FA",        "a fatty acid",                                                                   "CHEBI:28868", "RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"),
    c("M_H+",        "H(+)",                                                                           "CHEBI:15378", "RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"),
   
    ## coa_to_acdhap
    c("M_AcylCoA",            "an acyl-CoA",                   "CHEBI:58342", "RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660"),
    c("M_Dihydroxyacetone-P", "dihydroxyacetone phosphate",    "CHEBI:57642", "RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660"),
    c("M_AcylDHAP",           "a 1-acylglycerone 3-phosphate", "CHEBI:57534", "RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660"),
    c("M_CoA",                "CoA",                           "CHEBI:57287", "RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660"),
    
    ## coa_to_ce
    c("M_Cholesterol", "cholesterol",         "CHEBI:16113", "RHEA:17729", "RHEA:17730", "RHEA:17731", "RHEA:17732"),
    c("M_AcylCoA",     "an acyl-CoA",         "CHEBI:58342", "RHEA:17729", "RHEA:17730", "RHEA:17731", "RHEA:17732"),
    c("M_CE",          "a cholesterol ester", "CHEBI:17002", "RHEA:17729", "RHEA:17730", "RHEA:17731", "RHEA:17732"),
    c("M_CoA",         "CoA",                 "CHEBI:57287", "RHEA:17729", "RHEA:17730", "RHEA:17731", "RHEA:17732"),
    
    ## coa_to_fao
    c("M_AcylCoA", "a long-chain fatty acyl-CoA",        "CHEBI:83139", "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"),
    c("M_H+",      "H(+)",                               "CHEBI:15378", "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"),
    c("M_NADPH",   "NADPH",                              "CHEBI:57783", "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"),
    c("M_FAO",     "a long-chain primary fatty alcohol", "CHEBI:77396", "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"),
    c("M_CoA",     "CoA",                                "CHEBI:57287", "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"),
    c("M_NADP",    "NADP(+)",                            "CHEBI:58349", "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"),
    
    ## coa_to_lpa
    c("M_AcylCoA",      "an acyl-CoA",                     "CHEBI:58342", "RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328"),
    c("M_Glycerol-3-P", "sn-glycerol 3-phosphate",         "CHEBI:57597", "RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328"),
    c("M_LPA",          "a 1-acyl-sn-glycero-3-phosphate", "CHEBI:57970", "RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328"),
    c("M_CoA",          "CoA",                             "CHEBI:57287", "RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328"),
    
    ## dg_to_sn1mg
    c("M_1,2-DG", "a 1,2-diacylglycerol", "CHEBI:49172", "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"),
    c("M_H2O",    "H2O",                  "CHEBI:15377", "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"),
    c("M_1-MG",   "a 1-acylglycerol",     "CHEBI:35759", "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"),
    c("M_FA",     "a fatty acid",         "CHEBI:28868", "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"),
    c("M_H+",     "H(+)",                 "CHEBI:15378", "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol", "CHEBI:17815", "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"),
    c("M_H2O",    "H2O",                      "CHEBI:15377", "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"),
    c("M_1-MG",   "a 1-acyl-sn-glycerol",     "CHEBI:64683", "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"),
    c("M_FA",     "a fatty acid",             "CHEBI:28868", "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"),
    c("M_H+",     "H(+)",                     "CHEBI:15378", "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"),

    ## dg_to_sn2mg
    c("M_H2O",    "H2O",                      "CHEBI:15377", "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol", "CHEBI:17815", "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"),
    c("M_2-MG",   "a 2-acylglycerol",         "CHEBI:17389", "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"),
    c("M_FA",     "a fatty acid",             "CHEBI:28868", "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"),
    c("M_H+",     "H(+)",                     "CHEBI:15378", "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"),

    ## dg_to_pa
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol",            "CHEBI:17815", "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"),
    c("M_ATP",    "ATP",                                 "CHEBI:30616", "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"),
    c("M_PA",     "a 1,2-diacyl-sn-glycero-3-phosphate", "CHEBI:58608", "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"),
    c("M_ADP",    "ADP",                                 "CHEBI:456216", "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"),
    c("M_H+",     "H(+)",                                "CHEBI:15378", "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"),
    
    ## dg_to_pc
    c("M_1,2-DG",      "a 1,2-diacyl-sn-glycerol",                 "CHEBI:17815", "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942"),
    c("M_CDP-Choline", "CDP-choline",                              "CHEBI:58779", "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942"),
    c("M_PC",          "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942"),
    c("M_CMP",         "CMP",                                      "CHEBI:60377", "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942"),
    c("M_H+",          "H(+)",                                     "CHEBI:15378", "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942"),
    
    ## dg_to_pe
    c("M_1,2-DG",           "a 1,2-diacyl-sn-glycerol",                      "CHEBI:17815", "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946"),
    c("M_CDP-Ethanolamine", "CDP-ethanolamine",                              "CHEBI:57876", "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946"),
    c("M_PE",               "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946"),
    c("M_CMP",              "CMP",                                           "CHEBI:60377", "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946"),
    c("M_H+",               "H(+)",                                          "CHEBI:15378", "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946"),
    
    ## dg_to_tg
    c("M_1,2-DG",  "a 1,2-diacyl-sn-glycerol", "CHEBI:17815", "RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871"),
    c("M_AcylCoA", "an acyl-CoA",              "CHEBI:58342", "RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871"),
    c("M_TG",      "a triacyl-sn-glycerol",    "CHEBI:64615", "RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871"),
    c("M_CoA",     "CoA",                      "CHEBI:57287", "RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871"),
    
    ## dgo_to_pco
    c("M_DG-O",        "1-O-alkyl-2-acyl-sn-glycerol",                 "CHEBI:52595", "RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"),
    c("M_CDP-Choline", "CDP-choline",                                  "CHEBI:58779", "RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"),
    c("M_PC-O",        "1-O-alkyl-2-acyl-sn-glycero-3-phosphocholine", "CHEBI:36702", "RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"),
    c("M_CMP",         "CMP",                                          "CHEBI:60377", "RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"),
    c("M_H+",          "H(+)",                                         "CHEBI:15378", "RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"),
 
    ## dgo_to_peo
    c("M_DG-O",             "1-O-alkyl-2-acyl-sn-glycerol",                      "CHEBI:52595", "RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"),
    c("M_CDP-Ethanolamine", "CDP-ethanolamine",                                  "CHEBI:57876", "RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"),
    c("M_PE-O",             "1-O-alkyl-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:60520", "RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"),
    c("M_CMP",              "CMP",                                               "CHEBI:60377", "RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"),
    c("M_H+",               "H(+)",                                              "CHEBI:15378", "RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"),

    ## dhcer_to_cer
    c("M_DhCer",              "an N-acylsphinganine",    "CHEBI:31488", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),
    c("M_Fe2+-cytochrome_b5", "Fe(II)-[cytochrome b5]",  "CHEBI:29033", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),
    c("M_H+",                 "H(+)",                    "CHEBI:15378", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),
    c("M_O2",                 "O2",                      "CHEBI:15379", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),
    c("M_Cer",                "an N-acylsphing-4-enine", "CHEBI:52639", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),
    c("M_Fe3+-cytochrome_b5", "Fe(III)-[cytochrome b5]", "CHEBI:57540", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),
    c("M_H2O",                "H2O",                     "CHEBI:15377", "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547"),

    ## dhcer_to_dhsm
    c("M_PC",     "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:44620", "RHEA:44621", "RHEA:44622", "RHEA:44623"),
    c("M_DhCer",  "an N-acylsphinganine",                     "CHEBI:31488", "RHEA:44620", "RHEA:44621", "RHEA:44622", "RHEA:44623"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol",                 "CHEBI:17815", "RHEA:44620", "RHEA:44621", "RHEA:44622", "RHEA:44623"),
    c("M_DhSM",   "an N-acylsphinganine-1-phosphocholine",    "CHEBI:67090", "RHEA:44620", "RHEA:44621", "RHEA:44622", "RHEA:44623"),

    ## dhsm_to_dhcer
    c("M_DhSM",           "an N-(acyl)-sphingosylphosphocholine", "CHEBI:64583", "RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303"),
    c("M_H2O",            "H2O",                                  "CHEBI:15377", "RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303"),
    c("M_DhCer",          "an N-acyl-sphingoid base",             "CHEBI:83273", "RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303"),
    c("M_H+",             "H(+)",                                 "CHEBI:15378", "RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303"),
    c("M_Phosphocholine", "phosphocholine",                       "CHEBI:295975", "RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303"),
    
    ## fa_to_coa
    c("M_FA",      "a long-chain fatty acid",     "CHEBI:57560", "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"),
    c("M_ATP",     "ATP",                         "CHEBI:30616", "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"),
    c("M_CoA",     "CoA",                         "CHEBI:57287", "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"),
    c("M_AcylCoA", "a long-chain fatty acyl-CoA", "CHEBI:83139", "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"),
    c("M_AMP",     "AMP",                         "CHEBI:456215", "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"),
    c("M_PPi",     "diphosphate",                 "CHEBI:33019", "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"),
    c("M_FA",      "a fatty acid",                "CHEBI:28868", "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"),
    c("M_ATP",     "ATP",                         "CHEBI:30616", "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"),
    c("M_CoA",     "CoA",                         "CHEBI:57287", "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"),
    c("M_AcylCoA", "a fatty acyl-CoA",            "CHEBI:77636", "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"),
    c("M_AMP",     "AMP",                         "CHEBI:456215", "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"),
    c("M_PPi",     "diphosphate",                 "CHEBI:33019", "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"),

    ## lcl_to_cl
    c("M_1,2,4-LCL", "1'-[1,2-diacyl-sn-glycero-3-phospho],3'-[1-acyl-sn-glycero-3-phospho]-glycerol", "CHEBI:64743", "RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842"),
    c("M_AcylCoA",   "an acyl-CoA",                                                                    "CHEBI:58342", "RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842"),
    c("M_CL",        "a cardiolipin",                                                                  "CHEBI:62237", "RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842"),
    c("M_CoA",       "CoA",                                                                            "CHEBI:57287", "RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842"),
    
    ## lnape_to_gpnae
    c("M_H2O",   "H2O",                                         "CHEBI:15377", "RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423"),
    c("M_LNAPE", "N,1-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:85216", "RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423"),
    c("M_FA",    "a fatty acid",                                "CHEBI:28868", "RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423"),
    c("M_H+",    "H(+)",                                        "CHEBI:15378", "RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423"),
    c("M_GPNAE", "N-acyl-sn-glycero-3-phosphoethanolamine",     "CHEBI:85225", "RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423"),

    ## lpa_to_pa
    c("M_LPA",     "a 1-acyl-sn-glycero-3-phosphate",     "CHEBI:57970", "RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712"),
    c("M_AcylCoA", "an acyl-CoA",                         "CHEBI:58342", "RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712"),
    c("M_PA",      "a 1,2-diacyl-sn-glycero-3-phosphate", "CHEBI:58608", "RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712"),
    c("M_CoA",     "CoA",                                 "CHEBI:57287", "RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712"),
    
    ## lpao_to_pao
    c("M_LPA-O",   "1-O-alkyl-sn-glycero-3-phosphate",        "CHEBI:58014", "RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238"),
    c("M_AcylCoA", "an acyl-CoA",                             "CHEBI:58342", "RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238"),
    c("M_PA-O",    "1-O-alkyl-2-acyl-sn-glycero-3-phosphate", "CHEBI:73332", "RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238"),
    c("M_CoA",     "CoA",                                     "CHEBI:57287", "RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238"),
    
    ## sn1lpc_to_fa
    c("M_1-LPC",                 "a 1-acyl-sn-glycero-3-phosphocholine", "CHEBI:58168", "RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180"),
    c("M_H2O",                   "H2O",                                  "CHEBI:15377", "RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180"),
    c("M_FA",                    "a fatty acid",                         "CHEBI:28868", "RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180"),
    c("M_H+",                    "H(+)",                                 "CHEBI:15378", "RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180"),
    c("M_Glycerophosphocholine", "sn-glycerol 3-phosphocholine",         "CHEBI:16870", "RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180"),
    
    ## sn2lpc_to_fa
    c("M_2-LPC",                 "a 2-acyl-sn-glycero-3-phosphocholine", "CHEBI:57875", "RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"),
    c("M_H2O",                   "H2O",                                  "CHEBI:15377", "RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"),
    c("M_FA",                    "a fatty acid",                         "CHEBI:28868", "RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"),
    c("M_H+",                    "H(+)",                                 "CHEBI:15378", "RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"),
    c("M_Glycerophosphocholine", "sn-glycerol 3-phosphocholine",         "CHEBI:16870", "RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"),
   
    ## sn1lpc_to_pc
    c("M_1-LPC",   "a 1-acyl-sn-glycero-3-phosphocholine",     "CHEBI:58168", "RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940"),
    c("M_AcylCoA", "an acyl-CoA",                              "CHEBI:58342", "RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940"),
    c("M_PC",      "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940"),
    c("M_CoA",     "CoA",                                      "CHEBI:57287", "RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940"),

    ## lpco_to_pco
    c("M_1-LPC-O", "1-O-alkyl-sn-glycero-3-phosphocholine",        "CHEBI:30909", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
    c("M_AcylCoA", "an acyl-CoA",                                  "CHEBI:58342", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
    c("M_PC-O",    "1-O-alkyl-2-acyl-sn-glycero-3-phosphocholine", "CHEBI:36702", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
    c("M_CoA",     "CoA",                                          "CHEBI:57287", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
        
    ## sn1lpe_to_fa
    c("M_1-LPE",                      "a 1-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64381", "RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970"),
    c("M_H2O",                        "H2O",                                       "CHEBI:15377", "RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970"),
    c("M_FA",                         "a fatty acid",                              "CHEBI:28868", "RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970"),
    c("M_H+",                         "H(+)",                                      "CHEBI:15378", "RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970"),
    c("M_Glycerophosphoethanolamine", "sn-glycero-3-phosphoethanolamine",          "CHEBI:143890", "RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970"),
    
    ## sn2lpe_to_fa
    ## Reaction in RHEA not exisiting so far
    c("M_H2O",                        "H2O",                                       "CHEBI:15377",  NA, NA, NA, "sn2lpe_to_fa"),
    c("M_2-LPE",                      "a 2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:65213",  NA, NA, NA, "sn2lpe_to_fa"),
    c("M_H+",                         "H(+)",                                      "CHEBI:15378",  NA, NA, NA, "sn2lpe_to_fa"),
    c("M_FA",                         "a fatty acid",                              "CHEBI:28868",  NA, NA, NA, "sn2lpe_to_fa"),
    c("M_Glycerophosphoethanolamine", "sn-glycero-3-phosphoethanolamine",          "CHEBI:143890", NA, NA, NA, "sn2lpe_to_fa"),
        
    ## sn1lpe_to_pe
    c("M_1-LPE",   "a 1-acyl-sn-glycero-3-phosphoethanolamine",     "CHEBI:64381", "RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998"),
    c("M_AcylCoA", "an acyl-CoA",                                   "CHEBI:58342", "RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998"),
    c("M_PE",      "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998"),
    c("M_CoA",     "CoA",                                           "CHEBI:57287", "RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998"),
    
    ## sn1lpi_to_pi
    c("M_1-LPI",   "a 1-acyl-sn-glycero-3-phospho-(1D-myo-inositol)",     "CHEBI:64771", "RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198"),
    c("M_AcylCoA", "an acyl-CoA",                                         "CHEBI:58342", "RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198"),
    c("M_PI",      "a 1,2-diacyl-sn-glycero-3-phospho-(1D-myo-inositol)", "CHEBI:57880", "RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198"),
    c("M_CoA",     "CoA",                                                 "CHEBI:57287", "RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198"),
    
    ## lpeo_to_peo
    ## same as peo_to_lpeo??
    ## reaction not exiisting in Rhea
    c("M_AcylCoA", "an acyl-CoA",                                                     "CHEBI:58342", NA, NA, NA, "lpeo_to_peo"),
    c("M_LPE-O",   "",                                                                "",            NA, NA, NA, "lpeo_to_peo"),
    c("M_CoA",     "CoA",                                                             "CHEBI:57287", NA, NA, NA, "lpeo_to_peo"),
    c("M_PE-O",    "1-(1,2-saturated alkyl)-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:75028", NA, NA, NA, "lpeo_to_peo"),

    ## lpep_to_pep
    c("M_LPE-P",   "1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine",        "CHEBI:77288", "RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248"),
    c("M_AcylCoA", "an acyl-CoA",                                              "CHEBI:58342", "RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248"),
    c("M_PE-P",    "1-O-(1Z-alkenyl)-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:77290", "RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248"),
    c("M_CoA",     "CoA",                                                      "CHEBI:57287", "RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248"),
    
    ## sn1lpg_to_pg
    c("M_1-LPG",   "1-acyl-sn-glycero-3-phospho-(1'-sn-glycerol)",       "CHEBI:64840", "RHEA:33203", "RHEA:33204", "RHEA:33205", "RHEA:33206"),
    c("M_AcylCoA", "an acyl-CoA",                                        "CHEBI:58342", "RHEA:33203", "RHEA:33204", "RHEA:33205", "RHEA:33206"),
    c("M_PG",      "a 1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycerol)", "CHEBI:64716", "RHEA:33203", "RHEA:33204", "RHEA:33205", "RHEA:33206"),
    c("M_CoA",     "CoA",                                                "CHEBI:57287", "RHEA:33203", "RHEA:33204", "RHEA:33205", "RHEA:33206"),
    
    ## sn1mg_to_dg
    c("M_1-MG",    "a 1-acyl-sn-glycerol",     "CHEBI:64683", "RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466"),
    c("M_AcylCoA", "an acyl-CoA",              "CHEBI:58342", "RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466"),
    c("M_1,2-DG",  "a 1,2-diacyl-sn-glycerol", "CHEBI:17815", "RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466"),
    c("M_CoA",     "CoA",                      "CHEBI:57287", "RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466"),
    c("M_1-MG",    "a 1-acylglycerol",     "CHEBI:35759", "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946"),
    c("M_AcylCoA", "an acyl-CoA",          "CHEBI:58342", "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946"),
    c("M_1,2-DG",  "a 1,2-diacylglycerol", "CHEBI:49172", "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946"),
    c("M_CoA",     "CoA",                  "CHEBI:57287", "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946"),
    
    ## sn2mg_to_dg
    c("M_2-MG",    "a 2-acylglycerol",         "CHEBI:17389", "RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950"),
    c("M_AcylCoA", "an acyl-CoA",              "CHEBI:58342", "RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950"),
    c("M_1,2-DG",  "a 1,2-diacyl-sn-glycerol", "CHEBI:17815", "RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950"),
    c("M_CoA",     "CoA",                      "CHEBI:57287", "RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950"),
    c("M_2-MG",    "a 2-acylglycerol",     "CHEBI:17389", "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744"),
    c("M_AcylCoA", "an acyl-CoA",          "CHEBI:58342", "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744"),
    c("M_1,2-DG",  "a 1,2-diacylglycerol", "CHEBI:49172", "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744"),
    c("M_CoA",     "CoA",                  "CHEBI:57287", "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744"),

    ## sn1mg_to_fa
    c("M_1-MG",     "a 1-acylglycerol", "CHEBI:35759", "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"),
    c("M_H2O",      "H2O",              "CHEBI:15377", "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"),
    c("M_FA",       "a fatty acid",     "CHEBI:28868", "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"),
    c("M_Glycerol", "glycerol",         "CHEBI:17754", "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"),
    c("M_H+",       "H(+)",             "CHEBI:15378", "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"),
    
    ## sn2mg_to_fa
    c("M_2-MG",     "a 2-acylglycerol", "CHEBI:17389", "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"),
    c("M_H2O",      "H2O",              "CHEBI:15377", "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"),
    c("M_FA",       "a carboxylate",    "CHEBI:29067", "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"),
    c("M_Glycerol", "glycerol",         "CHEBI:17754", "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"),
    c("M_H+",       "H(+)",             "CHEBI:15378", "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"),
    
    ## sn1mg_to_lpa
    c("M_1-MG", "a 1-acyl-sn-glycerol",            "CHEBI:64683", "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"),
    c("M_ATP",  "ATP",                             "CHEBI:30616", "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"),
    c("M_LPA",  "a 1-acyl-sn-glycero-3-phosphate", "CHEBI:57970", "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"),
    c("M_ADP",  "ADP",                             "CHEBI:456216", "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"),
    c("M_H+",   "H(+)",                            "CHEBI:15378", "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"),
    
    ## sn2mg_to_sn1mg
    ## Reaction in RHEA not exisiting
    c("M_2-MG", "a 1-acylglycerol", "CHEBI:35759", NA, NA, NA, "sn2mg_to_sn1mg"),
    c("M_1-MG", "a 2-acylglycerol", "CHEBI:17389", NA, NA, NA, "sn2mg_to_sn1mg"),

    ## nae_to_fa
    ## Both are correct, but the second reaction is the more generic one
    c("M_H2O",          "H2O",                                   "CHEBI:15377", "RHEA:17505" , "RHEA:17506", "RHEA:17507", "RHEA:17508"),
    c("M_NAE",          "N-(long-chain fatty acyl)ethanolamine", "CHEBI:15897", "RHEA:17505" , "RHEA:17506", "RHEA:17507", "RHEA:17508"),
    c("M_FA",           "a long-chain fatty acid",               "CHEBI:57560", "RHEA:17505" , "RHEA:17506", "RHEA:17507", "RHEA:17508"),
    c("M_Ethanolamine", "ethanolamine",                          "CHEBI:57603", "RHEA:17505" , "RHEA:17506", "RHEA:17507", "RHEA:17508"),
    c("M_H2O",          "H2O",                   "CHEBI:15377", "RHEA:39995" , "RHEA:39996", "RHEA:39997", "RHEA:39998"),
    c("M_NAE",          "an N-acylethanolamine", "CHEBI:52640", "RHEA:39995" , "RHEA:39996", "RHEA:39997", "RHEA:39998"),
    c("M_FA",           "a fatty acid",          "CHEBI:28868", "RHEA:39995" , "RHEA:39996", "RHEA:39997", "RHEA:39998"),
    c("M_Ethanolamine", "ethanolamine",          "CHEBI:57603", "RHEA:39995" , "RHEA:39996", "RHEA:39997", "RHEA:39998"),

    ## nape_to_lnape
    c("M_NAPE",  "N-acyl-1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:62537", "RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463"),
    c("M_H2O",   "H2O",                                                "CHEBI:15377", "RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463"),
    c("M_LNAPE", "N,1-diacyl-sn-glycero-3-phosphoethanolamine",        "CHEBI:85216", "RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463"),
    c("M_FA",    "a fatty acid",                                       "CHEBI:28868", "RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463"),
    c("M_H+",    "H(+)",                                               "CHEBI:15378", "RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463"),
    
    ## nape_to_nae
    c("M_H2O",  "H2O",                                                "CHEBI:15377", "RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162"),
    c("M_NAPE", "N-acyl-1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:62537", "RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162"),
    c("M_PA",   "a 1,2-diacyl-sn-glycero-3-phosphate",                "CHEBI:58608", "RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162"),
    c("M_NAE",  "an N-acylethanolamine",                              "CHEBI:52640", "RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162"),
    c("M_H+",   "H(+)",                                               "CHEBI:15378", "RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162"),

    ## nape_to_pnae
    ## reaction not existing in Rhea so far
    c("M_NAPE",   "N-acyl-1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:62537",  NA, NA, NA, "nape_to_pnae"),
    c("M_H2O",    "H2O",                                                "CHEBI:15377",  NA, NA, NA, "nape_to_pnae"),
    c("M_PNAE",   "N-acylethanolamine phosphate",                       "CHEBI:145538", NA, NA, NA, "nape_to_pnae"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol",                           "CHEBI:17815",  NA, NA, NA, "nape_to_pnae"),
        
    ## napeo_to_nae
    ## reaction not existing in Rhea so far
    c("M_NAPEO", "",                                        "",            NA, NA, NA, "napeo_to_nae"),
    c("M_H2O",   "H2O",                                     "CHEBI:15377", NA, NA, NA, "napeo_to_nae"),
    c("M_NAE",   "an N-acylethanolamine",                   "CHEBI:52640", NA, NA, NA, "napeo_to_nae"),
    c("M_PA-O",  "1-O-alkyl-2-acyl-sn-glycero-3-phosphate", "CHEBI:73332", NA, NA, NA, "napeo_to_nae"),
        
    ## pa_to_cdpdg
    c("M_PA",     "a 1,2-diacyl-sn-glycero-3-phosphate", "CHEBI:58608", "RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"),
    c("M_CTP",    "CTP",                                 "CHEBI:37563", "RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"),
    c("M_H+",     "H(+)",                                "CHEBI:15378", "RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"),
    c("M_CDP-DG", "a CDP-1,2-diacyl-sn-glycerol",        "CHEBI:58332", "RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"),
    c("M_PPi",    "diphosphate",                         "CHEBI:33019", "RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"),
    
    ## pa_to_dg
    c("M_PA",     "a 1,2-diacyl-sn-glycero-3-phosphate", "CHEBI:58608", "RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"),
    c("M_H2O",    "H2O",                                 "CHEBI:15377", "RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol",            "CHEBI:17815", "RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"),
    c("M_Pi",     "phosphate",                           "CHEBI:43474", "RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"),
        
    ## pao_to_dgo
    c("M_PA-O", "1-O-alkyl-2-acyl-sn-glycero-3-phosphate", "CHEBI:73332", "RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242"),
    c("M_H2O",  "H2O",                                     "CHEBI:15377", "RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242"),
    c("M_DG-O", "1-O-alkyl-2-acyl-sn-glycerol",            "CHEBI:52595", "RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242"),
    c("M_Pi",   "phosphate",                               "CHEBI:43474", "RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242"),
        
    ## pc_to_dg
    c("M_PC",             "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607"),
    c("M_H2O",            "H2O",                                      "CHEBI:15377", "RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607"),
    c("M_1,2-DG",         "a 1,2-diacyl-sn-glycerol",                 "CHEBI:17815", "RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607"),
    c("M_H+",             "H(+)",                                     "CHEBI:15378", "RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607"),
    c("M_Phosphocholine", "phosphocholine",                           "CHEBI:295975", "RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607"),

    ## pc_to_sn1lpc
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804"),
    c("M_H2O",   "H2O",                                      "CHEBI:15377", "RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804"),
    c("M_1-LPC", "a 1-acyl-sn-glycero-3-phosphocholine",     "CHEBI:58168", "RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804"),
    c("M_FA",    "a fatty acid",                             "CHEBI:28868", "RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804"),
    c("M_H+",    "H(+)",                                     "CHEBI:15378", "RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804"),

    ## pc_to_sn2lpc
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:18689", "RHEA:18690", "RHEA:18691", "RHEA:18692"),
    c("M_H2O",   "H2O",                                      "CHEBI:15377", "RHEA:18689", "RHEA:18690", "RHEA:18691", "RHEA:18692"),
    c("M_2-LPC", "a 2-acyl-sn-glycero-3-phosphocholine",     "CHEBI:57875", "RHEA:18689", "RHEA:18690", "RHEA:18691", "RHEA:18692"),
    c("M_FA",    "a fatty acid",                             "CHEBI:28868", "RHEA:18689", "RHEA:18690", "RHEA:18691", "RHEA:18692"),
    c("M_H+",    "H(+)",                                     "CHEBI:15378", "RHEA:18689", "RHEA:18690", "RHEA:18691", "RHEA:18692"),
        
    ## pc_to_pa
    c("M_PC",      "a 1,2-diacyl-sn-glycero-3-phosphocholine", "CHEBI:57643", "RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"),
    c("M_H2O",     "H2O",                                      "CHEBI:15377", "RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"),
    c("M_PA",      "a 1,2-diacyl-sn-glycero-3-phosphate",      "CHEBI:58608", "RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"),
    c("M_Choline", "choline",                                  "CHEBI:15354", "RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"),
    c("M_H+",      "H(+)",                                     "CHEBI:15378", "RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"),
        
    ## pc_to_ps
    c("M_PC",       "a 1,2-diacyl-sn-glycero-3-phosphocholine",   "CHEBI:57643", "RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091"),
    c("M_L-Serine", "L-serine",                                   "CHEBI:33384", "RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091"),
    c("M_PS",       "a 1,2-diacyl-sn-glycero-3-phospho-L-serine", "CHEBI:57262", "RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091"),
    c("M_Choline",  "choline",                                    "CHEBI:15354", "RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091"),
    
    ## pco_to_lpco
    c("M_PC-O",  "1-O-alkyl-2-acyl-sn-glycero-3-phosphocholine", "CHEBI:36702", "RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234"),
    c("M_H2O",   "H2O",                                          "CHEBI:15377", "RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234"),
    c("M_LPC-O", "1-O-alkyl-sn-glycero-3-phosphocholine",        "CHEBI:30909", "RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234"),
    c("M_FA",    "a fatty acid",                                 "CHEBI:28868", "RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234"),
    c("M_H+",    "H(+)",                                         "CHEBI:15378", "RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234"),
    
    ## lpco_to_lpao
    c("M_1-LPC-O", "1-O-alkyl-sn-glycero-3-phosphocholine", "CHEBI:30909", "RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930"),
    c("M_H2O",     "H2O",                                   "CHEBI:15377", "RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930"),
    c("M_1-LPA-O", "1-O-alkyl-sn-glycero-3-phosphate",      "CHEBI:58014", "RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930"),
    c("M_Choline", "choline",                               "CHEBI:15354", "RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930"),
    c("M_H+",      "H(+)",                                  "CHEBI:15378", "RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930"),

    ## lpco_to_mgo
    c("M_1-LPC-O",        "1-O-alkyl-sn-glycero-3-phosphocholine", "CHEBI:30909", "RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086"),
    c("M_H2O",            "H2O",                                   "CHEBI:15377", "RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086"),
    c("M_1-MG-O",         "1-O-alkyl-sn-glycerol",                 "CHEBI:15850", "RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086"),
    c("M_H+",             "H(+)",                                  "CHEBI:15378", "RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086"),
    c("M_Phosphocholine", "phosphocholine",                        "CHEBI:295975", "RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086"),
        
    ## lpco_to_pco
    c("M_1-LPC-O", "1-O-alkyl-sn-glycero-3-phosphocholine",        "CHEBI:30909", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
    c("M_AcylCoA", "an acyl-CoA",                                  "CHEBI:58342", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
    c("M_PC-O",    "1-O-alkyl-2-acyl-sn-glycero-3-phosphocholine", "CHEBI:36702", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),
    c("M_CoA",     "CoA",                                          "CHEBI:57287", "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"),

    ## pe_to_dg
    c("M_H2O",            "H2O",                                           "CHEBI:15377", "RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954"),
    c("M_PE",             "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954"),
    c("M_P-Ethanolamine", "phosphoethanolamine",                           "CHEBI:58190", "RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954"),
    c("M_1,2-DG",         "a 1,2-diacyl-sn-glycerol",                      "CHEBI:17815", "RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954"),
    c("M_H+",             "H(+)",                                          "CHEBI:15378", "RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954"),

    ## pe_to_sn1lpe
    c("M_PE",    "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"),
    c("M_H2O",   "H2O",                                           "CHEBI:15377", "RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"),
    c("M_1-LPE", "a 1-acyl-sn-glycero-3-phosphoethanolamine",     "CHEBI:64381", "RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"),
    c("M_FA",    "a fatty acid",                                  "CHEBI:28868", "RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"),
    c("M_H+",    "H(+)",                                          "CHEBI:15378", "RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"),
        
    ## pe_to_sn2lpe
    c("M_PE",    "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"),
    c("M_H2O",   "H2O",                                           "CHEBI:15377", "RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"),
    c("M_2-LPE", "a 2-acyl-sn-glycero-3-phosphoethanolamine",     "CHEBI:65213", "RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"),
    c("M_FA",    "a fatty acid",                                  "CHEBI:28868", "RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"),
    c("M_H+",    "H(+)",                                          "CHEBI:15378", "RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"),
    
    ## pe_to_nape_sn1
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine",           "CHEBI:57643", "RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191"),
    c("M_PE",    "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine",      "CHEBI:64612", "RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191"),
    c("M_2-LPC", "a 2-acyl-sn-glycero-3-phosphocholine",               "CHEBI:57875", "RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191"),
    c("M_H+",    "H(+)",                                               "CHEBI:15378", "RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191"),
    c("M_NAPE",  "N-acyl-1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:62537", "RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191"),
        
    ## pe_to_nape_sn2
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine",           "CHEBI:57643", "RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195"),
    c("M_PE",    "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine",      "CHEBI:64612", "RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195"),
    c("M_1-LPC", "a 1-acyl-sn-glycero-3-phosphocholine",               "CHEBI:58168", "RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195"),
    c("M_H+",    "H(+)",                                               "CHEBI:15378", "RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195"),
    c("M_NAPE",  "N-acyl-1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:62537", "RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195"),
  
    ## pe_to_pa
    ## reaction not present in Rhea so far
    c("M_H2O",          "H2O",                                           "CHEBI:15377", NA, NA, NA, "pe_to_pa"),
    c("M_PE",           "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", NA, NA, NA, "pe_to_pa"),
    c("M_Ethanolamine", "ethanolamine",                                  "CHEBI:57603", NA, NA, NA, "pe_to_pa"),
    c("M_H+",           "H(+)",                                          "CHEBI:15378", NA, NA, NA, "pe_to_pa"),
    c("M_PA",           "a 1,2-diacyl-sn-glycero-3-phosphate",           "CHEBI:58608", NA, NA, NA, "pe_to_pa"),

    ## pe_to_ps
    c("M_PE",           "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609"),
    c("M_L-Serine",     "L-serine",                                      "CHEBI:33384", "RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609"),
    c("M_PS",           "a 1,2-diacyl-sn-glycero-3-phospho-L-serine",    "CHEBI:57262", "RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609"),
    c("M_Ethanolamine", "ethanolamine",                                  "CHEBI:57603", "RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609"),

    ## peo_to_lpeo
    ## reaction not present in Rhea so far
    c("M_H2O",   "H2O",                                                             "CHEBI:15377", NA, NA, NA, "peo_to_lpeo"),
    c("M_PE-O",  "1-(1,2-saturated alkyl)-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:75028", NA, NA, NA, "peo_to_lpeo"),
    c("M_H+",    "H(+)",                                                            "CHEBI:15378", NA, NA, NA, "peo_to_lpeo"),
    c("M_LPE-O", "",                                                                "",            NA, NA, NA, "peo_to_lpeo"),
    c("M_FA",    "a fatty acid",                                                    "CHEBI:28868", NA, NA, NA, "peo_to_lpeo"),
        
    ## peo_to_napeo_sn1
    ## reaction not present in Rhea so far
    c("M_PE-O",  "1-(1,2-saturated alkyl)-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:75028", NA, NA, NA, "peo_to_napeo_sn1"),
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine",                        "CHEBI:57643", NA, NA, NA, "peo_to_napeo_sn1"),
    c("M_NAPEO", "",                                                                "",            NA, NA, NA, "peo_to_napeo_sn1"),
    c("M_2-LPC", "a 2-acyl-sn-glycero-3-phosphocholine",                            "CHEBI:58168", NA, NA, NA, "peo_to_napeo_sn1"),
        
    ## peo_to_napeo_sn2
    ## reaction not prsent in Rhea so far
    c("M_PE-O",  "1-(1,2-saturated alkyl)-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:75028", NA, NA, NA, "peo_to_napeo_sn2"),
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine",                        "CHEBI:57643", NA, NA, NA, "peo_to_napeo_sn2"),
    c("M_NAPEO", "",                                                                "",            NA, NA, NA, "peo_to_napeo_sn2"),
    c("M_1-LPC", "a 1-acyl-sn-glycero-3-phosphocholine",                            "CHEBI:58168", NA, NA, NA, "peo_to_napeo_sn2"),
    
    ## peo_to_pep
    c("M_PE-O",               "1-(1,2-saturated alkyl)-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:75028", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),
    c("M_Fe2+-cytochrome_b5", "Fe(II)-[cytochrome b5]",                                          "CHEBI:29033", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),
    c("M_H+",                 "H(+)",                                                            "CHEBI:15378", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),
    c("M_O2",                 "O2",                                                              "CHEBI:15379", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),
    c("M_PE-P",               "1-O-(1Z-alkenyl)-2-acyl-sn-glycero-3-phosphoethanolamine",        "CHEBI:77290", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),
    c("M_Fe3+-cytochrome_b5", "Fe(III)-[cytochrome b5]",                                         "CHEBI:29034", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),
    c("M_H2O",                "H2O",                                                             "CHEBI:15377", "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"),

    ## pep_to_lpep
    c("M_PE-P",   "1-O-(1Z-alkenyl)-2-acyl-sn-glycero-3-phosphoethanolamine",  "CHEBI:77290", "RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"),
    c("M_H2O",    "H2O",                                                       "CHEBI:15377", "RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"),
    c("M_LPE-P", "1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine",          "CHEBI:77288", "RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"),
    c("M_FA",    "a carboxylate",                                              "CHEBI:29067", "RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"),
    c("M_H+",    "H(+)",                                                       "CHEBI:15378", "RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"),

    ## lpep_to_fal
    c("M_1-LPE-P",                    "1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine", "CHEBI:77288", "RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908"),
    c("M_H2O",                        "H2O",                                               "CHEBI:15377", "RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908"),
    c("M_FAL",                        "a 2,3-saturated aldehyde",                          "CHEBI:73359", "RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908"),
    c("M_Glycerophosphoethanolamine", "sn-glycero-3-phosphoethanolamine",                  "CHEBI:143890", "RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908"),
    
    ## lpep_to_lpap
    c("M_1-LPE-P",      "1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine", "CHEBI:77288", "RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206"),
    c("M_H2O",          "H2O",                                               "CHEBI:15377", "RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206"),
    c("M_LPA-P",        "1-O-(1Z-alkenyl)-sn-glycero-3-phosphate",           "CHEBI:77283", "RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206"),
    c("M_Ethanolamine", "ethanolamine",                                      "CHEBI:57603", "RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206"),
    c("M_H+",           "H(+)",                                              "CHEBI:15378", "RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206"),

    ## lpep_to_mgp
    c("M_1-LPE-P",             "1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine", "CHEBI:77288", "RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202"),
    c("M_H2O",                 "H2O",                                               "CHEBI:15377", "RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202"),
    c("M_1-MG-P",              "1-O-(1Z-alkenyl)-sn-glycerol",                      "CHEBI:77297", "RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202"),
    c("M_H+",                  "H(+)",                                              "CHEBI:15378", "RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202"),
    c("M_Phosphoethanolamine", "phosphoethanolamine",                               "CHEBI:58190", "RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202"),

    ## pep_to_napep_sn1
    c("M_PE-P",  "1-O-(1Z-alkenyl)-2-acyl-sn-glycero-3-phosphoethanolamine",       "CHEBI:77290",  "RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599"),
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine",                       "CHEBI:57643",  "RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599"),
    c("M_2-LPC", "a 2-acyl-sn-glycero-3-phosphocholine",                           "CHEBI:57875",  "RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599"),
    c("M_H+",    "H(+)",                                                           "CHEBI:15378",  "RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599"),
    c("M_NAPEP", "N-acyl-1-[(1Z)-alkenyl]-2-acyl-sn-glycero-3-phosphoethanolamine","CHEBI:140451", "RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599"),
    

    ## pep_to_napep_sn2
    ## reaction currently not existing in Rhea
    c("M_PE-P",  "1-O-(1Z-alkenyl)-2-acyl-sn-glycero-3-phosphoethanolamine",        "CHEBI:77290",  NA, NA, NA, "pep_to_napep_sn2"),
    c("M_PC",    "a 1,2-diacyl-sn-glycero-3-phosphocholine",                        "CHEBI:57643",  NA, NA, NA, "pep_to_napep_sn2"),
    c("M_NAPEP", "N-acyl-1-[(1Z)-alkenyl]-2-acyl-sn-glycero-3-phosphoethanolamine", "CHEBI:140451", NA, NA, NA, "pep_to_napep_sn2"),
    c("M_1-LPC", "a 1-acyl-sn-glycero-3-phosphocholine",                            "CHEBI:58168", NA, NA, NA, "pep_to_napep_sn2"),
        
    ## pg_to_cl
    c("M_PG",     "1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycerol)", "CHEBI:64716", "RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"),
    c("M_CDP-DG", "a CDP-1,2-diacyl-sn-glycerol",                     "CHEBI:58332", "RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"),
    c("M_CL",     "a cardiolipin",                                    "CHEBI:62237", "RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"),
    c("M_CMP",    "CMP",                                              "CHEBI:60377", "RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"),
    c("M_H+",     "H(+)",                                             "CHEBI:15378", "RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"),
    
    ## pgp_to_pg
    c("M_PGP", "1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycero-3'-phosphate)", "CHEBI:60110", "RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754"),
    c("M_H2O", "H2O",                                                          "CHEBI:15377", "RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754"),
    c("M_PG",  "1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycerol)",             "CHEBI:64716", "RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754"),
    c("M_Pi",  "phosphate",                                                    "CHEBI:43474", "RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754"),
        
    ## pi_to_dg
    c("M_PI",               "a 1,2-diacyl-sn-glycero-3-phospho-(1D-myo-inositol)", "CHEBI:57880", "RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"),
    c("M_H2O",              "H2O",                                                 "CHEBI:15377", "RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"),
    c("M_myo-Inositol-1-P", "1D-myo-inositol 1-phosphate",                         "CHEBI:58433", "RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"),
    c("M_1,2-DG",           "a 1,2-diacyl-sn-glycerol",                            "CHEBI:17815", "RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"),
    c("M_H+",               "H(+)",                                                "CHEBI:15378", "RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"),

    ## pi_to_sn1lpi
    c("M_PI",    "a 1,2-diacyl-sn-glycero-3-phospho-(1D-myo-inositol)", "CHEBI:57880", "RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"),
    c("M_H2O",   "H2O",                                                 "CHEBI:15377", "RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"),
    c("M_1-LPI", "a 1-acyl-sn-glycero-3-phospho-(1D-myo-inositol)",     "CHEBI:64771", "RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"),
    c("M_FA",    "a fatty acid",                                        "CHEBI:28868", "RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"),
    c("M_H+",    "H(+)",                                                "CHEBI:15378", "RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"),
    
    ## ps_to_pe
    c("M_PS",  "a 1,2-diacyl-sn-glycero-3-phospho-L-serine",    "CHEBI:57262", "RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"),
    c("M_H+",  "H(+)",                                          "CHEBI:15378", "RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"),
    c("M_PE",  "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", "CHEBI:64612", "RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"),
    c("M_CO2", "CO2",                                           "CHEBI:16526", "RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"),

    ## sm_to_cer    
    c("M_H2O",            "H2O",                                              "CHEBI:15377", "RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647"),
    c("M_SM",             "N-(hexadecanoyl)-sphing-4-enine-1-phosphocholine", "CHEBI:78646", "RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647"),
    c("M_H+",             "H(+)",                                             "CHEBI:15378", "RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647"),
    c("M_Cer",            "N-hexadecanoylsphing-4-enine",                     "CHEBI:72959", "RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647"),
    c("M_Phosphocholine", "phosphocholine",                                   "CHEBI:295975", "RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647"),
        
    ## sphinga_to_dhcer
    c("M_AcylCoA",     "a fatty acyl-CoA",         "CHEBI:77636", "RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427"),
    c("M_Sphinganine", "a sphingoid base",         "CHEBI:84410", "RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427"),
    c("M_DhCer",       "an N-acyl-sphingoid base", "CHEBI:83273", "RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427"),
    c("M_CoA",         "CoA",                      "CHEBI:57287", "RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427"),
    c("M_H+",          "H(+)",                     "CHEBI:15378", "RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427"),
    
    ## tg_to_dg
    c("M_TG",     "a triacyl-sn-glycerol",    "CHEBI:64615", "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"),
    c("M_H2O",    "H2O",                      "CHEBI:15377", "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"),
    c("M_1,2-DG", "a 1,2-diacyl-sn-glycerol", "CHEBI:17815", "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"),
    c("M_FA",     "a fatty acid",             "CHEBI:28868", "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"),
    c("M_H+",     "H(+)",                     "CHEBI:15378", "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"),
    c("M_TG",     "a triacyl-sn-glycerol", "CHEBI:64615", "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"),
    c("M_H2O",    "H2O",                   "CHEBI:15377", "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"),
    c("M_1,2-DG", "a 1,2-diacylglycerol",  "CHEBI:49172", "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"),
    c("M_FA",     "a fatty acid",          "CHEBI:28868", "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"),
    c("M_H+",     "H(+)",                  "CHEBI:15378", "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"))
    
    ## set colnames 
    colnames(res) <- c("name_LipidNetworkPredictR", "name_RHEA", "ChEBI_ID", 
        "RHEA_ID_undetermined", "RHEA_ID_left_to_right", 
        "RHEA_ID_right_to_left", "RHEA_ID_both_directions")
    
    ## return the object
    res
}
