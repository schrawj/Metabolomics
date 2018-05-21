#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.03.06.
#' 
#' Volcano plots depicting the univariate associations of each compound 
#' with MRD and relapse.  Figure 1 panels A and B.
#' 
#' 2018.04.26.
#' 
#' Explanation of approach to identifying top metabolites has changed.
#' 
#' Exclude currency metabolites, then retain top 20 by K-W p-value.
#' 
#' Points to plot and label must change to reflect this.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

require(stringr); require(ggplot2); require(ggrepel); require(dplyr)

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load("./Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata")

#' List of currency metabolites to exclude from plots.
met.signif <- filter(met.signif, !(compound %in% c('phosphate',"cytidine 5'-diphosphocholine","inosine 5'-monophosphate (IMP)", "adenosine 3',5'-cyclic monophosphate (cAMP)", "adenosine 3'-monophosphate (3'-AMP)")))

#' Cosmetic changes for brevity and formatting.
met.signif$super.pathway <- ifelse(met.signif$super.pathway == 'Cofactors and Vitamins', 'Cofactor/vitamin', met.signif$super.pathway)
met.signif$super.pathway <- ifelse(met.signif$super.pathway == 'Amino Acid', 'Amino acid', met.signif$super.pathway)

met.signif$compound <- str_replace(met.signif$compound, '[*]$', '') # strips trailing asterisks
met.signif$compound <- str_replace(met.signif$compound, '[(].+[)]$', '') # trailing parentheses and contents



# MRD plot ----------------------------------------------------------------

print(ggplot(data = met.signif, aes(x = log2(met.signif$fold.change.mrd), y = -log10(met.signif$mrd.pvalue.kruskal))) + 
        geom_point(aes(color = met.signif$super.pathway)) +
        scale_color_manual(values = c('blue2','red','gold','violetred','springgreen3','purple','cyan','tan4')) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() +
        theme(text = element_text(size = 17.5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank()) +
        scale_x_continuous(limits = c(-3,5)) +
        guides(colour = guide_legend(title = 'Super Pathway', override.aes = list(size=5))) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-keto-3-deoxy-gluconate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.2, nudge_y = 0.6, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'gamma-carboxyglutamate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.9, nudge_y = 0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound %in% c('malate','fumarate'), compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'alanine', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'histidine', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.4, nudge_y = -0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dihydroorotate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.75, nudge_y = -1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '4-hydroxy-2-oxoglutaric acid', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 1.85, nudge_y = -0.75, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-hydroxyoctanoate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 1.75, nudge_y = -0.25, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'succinate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.35, nudge_y = 0.35, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pyruvate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pyrraline', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "3-phosphoglycerate", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.125, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-oleoyl-GPI ", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-linoleoyl-GPI ", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-palmitoyl-GPI ", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -1, nudge_y = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-arachidonoyl-GPI ", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.1, nudge_y = 0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-stearoyl-GPI ", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'acetoacetate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-arachidonylglycerol ", compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.5, show.legend = FALSE))



# Relapse plot ------------------------------------------------------------

print(ggplot(data = met.signif, aes(x = log2(met.signif$fold.change.relapse), y = -log10(met.signif$relapse.pvalue.kruskal))) + 
        geom_point(aes(color = met.signif$super.pathway)) +
        scale_color_manual(values = c('blue2','red','gold','violetred','springgreen3','purple','cyan','tan4')) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() +
        theme(text = element_text(size = 17.5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank()) +
        scale_x_continuous(limits = c(-3,5)) +
        guides(colour = guide_legend(title = 'Super Pathway', override.aes = list(size=5))) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dimethylglycine', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '10-heptadecenoate ', compound, ''),
                            color = super.pathway), cex = 7.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'myristate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -1, nudge_y = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '15-methylpalmitate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'palmitoleate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.6, nudge_y = 0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '13-methylmyristate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.25, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'margarate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'valine', compound, ''),
                            color = super.pathway), cex = 7.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '10-nonadecenoate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -2, nudge_y = 0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pentadecanoate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 1.6, nudge_y = 0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'picolinate', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -1.25, nudge_y = 0.05, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'palmitate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.75, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'oleate/vaccenate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -2, nudge_y = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '1-palmitoyl-2-linoleoyl-GPI ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 2, nudge_y = -0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'nonadecanoate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.25, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'S-adenosylhomocysteine ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.9, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-stearoyl-GPE ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 1.75, nudge_y = -0.45, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '1-linoleoyl-GPI ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -1.5, nudge_y = -0.05, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'linoleate ', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = -0.4, nudge_y = -0.15, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'urea', compound, ''),
                            color = super.pathway), cex = 7.5, nudge_x = 0.6, nudge_y = -0.3, show.legend = FALSE))
        