# prepare Data for visualization in shiny web app
library(data.table)
source("tracesMethods.R")
load("~/sonas/html/HEK293_SEC_SWATH_data_analysis/DownstreamAnalysis/SWATH_MS2_Corum_BioPlex_StringDB/StringDB/workspaceItems/string.up.09.netDist.degree1groups")
string_subnetworks <- string.up.09.netDist.degree1groups
rm(string.up.09.netDist.degree1groups)
calibration_functions <- readRDS("SECMWcalibration.rda")

up <- fread("~/sonas/databases/Uniprot/uniprot-all-human9606-20161130_extended.tab")
up[, Mass:=as.numeric(gsub(",", ".", Mass))]
names(up) <- gsub(" ", "_", names(up))


prot_int_r1 <- readRDS("../../SWATH_Extractions/SECfractions_Rep1/Interphase/protTraces_raw_a_fdr1pc_top2.rda")
prot_mit_r1 <- readRDS("../../SWATH_Extractions/SECfractions_Rep1/Mitosis/protTraces_raw_a_fdr1pc_top2.rda")
prot_int_r2 <- readRDS("../../SWATH_Extractions/SECfractions_Rep2/Interphase/protTraces_raw_a_fdr1pc_top2.rda")
prot_mit_r2 <- readRDS("../../SWATH_Extractions/SECfractions_Rep2/Mitosis/protTraces_raw_a_fdr1pc_top2.rda")
prot_int_r3 <- readRDS("../../SWATH_Extractions/SECfractions_Rep3/Interphase/protTraces_raw_a_fdr1pc_top2.rda")
prot_mit_r3 <- readRDS("../../SWATH_Extractions/SECfractions_Rep3/Mitosis/protTraces_raw_a_fdr1pc_top2.rda")

trace_annotation_cum <- unique(rbind(prot_int_r1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw)],
                                     prot_int_r2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw)],
                                     prot_int_r3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw)],
                                     prot_mit_r1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw)],
                                     prot_mit_r2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw)],
                                     prot_mit_r3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw)]))


pep_int_r1 <- readRDS("../../SWATH_Extractions/SECfractions_Rep1/Interphase/pepTraces_raw_a_cs3_fdr1pc.rda")
pep_mit_r1 <- readRDS("../../SWATH_Extractions/SECfractions_Rep1/Mitosis/pepTraces_raw_a_cs3_fdr1pc.rda")
pep_int_r2 <- readRDS("../../SWATH_Extractions/SECfractions_Rep2/Interphase/pepTraces_raw_a_cs3_fdr1pc.rda")
pep_mit_r2 <- readRDS("../../SWATH_Extractions/SECfractions_Rep2/Mitosis/pepTraces_raw_a_cs3_fdr1pc.rda")
pep_int_r3 <- readRDS("../../SWATH_Extractions/SECfractions_Rep3/Interphase/pepTraces_raw_a_cs3_fdr1pc.rda")
pep_mit_r3 <- readRDS("../../SWATH_Extractions/SECfractions_Rep3/Mitosis/pepTraces_raw_a_cs3_fdr1pc.rda")

trace_annotation_cum_peptides <- unique(rbind(pep_int_r1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw, id)],
                                              pep_int_r2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw, id)],
                                              pep_int_r3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw, id)],
                                              pep_mit_r1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw, id)],
                                              pep_mit_r2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw, id)],
                                              pep_mit_r3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Protein_names, protein_mw, id)]))

DiffCorrTable <- readRDS("../../CompareTraces/Hela_CCSEC_IvsM_DiffCorrTable_combined.rda")
trace_annotation_cum <- merge(trace_annotation_cum, unique(DiffCorrTable[, .(protein_id, pearson_cor_avg,
                                                                             pearson_cor_avg_rank,
                                                                             difference_normalized_avg,
                                                                             difference_normalized_avg_rank)]), by = "protein_id", all.x = TRUE)

save.image("data.rda")
