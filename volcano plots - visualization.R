#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2018.03.06.
#' 
#' Volcano plots depicting the univariate associations of each compound with MRD and 
#' relapse.
#' 
#' Likely these will constitute figure 1 and figure 2. 
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load("./Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata")

require(stringr)

#' Function that condenses the code required to generate a volcano plot somewhat.
volcano.plot <- function(a, x, y, z){
  ggplot(data = a, aes(x = log2(x), y = -log10(y))) + 
    geom_point(aes(color = z))
}

#' Cosmetic changes for brevity and formatting.
met.signif$super.pathway <- ifelse(met.signif$super.pathway == 'Cofactors and Vitamins', 'Cofactor/vitamin', met.signif$super.pathway)
met.signif$super.pathway <- ifelse(met.signif$super.pathway == 'Amino Acid', 'Amino acid', met.signif$super.pathway)

met.signif$compound <- str_replace(met.signif$compound, '[*]$', '') # strips trailing asterisks
met.signif$compound <- str_replace(met.signif$compound, '[(].+[)]$', '') # trailing parentheses and contents

#' Flags that allow easy simultaneous plotting of compounds that are significant by at least one test.
met.signif$mrd.sig.flag <- ifelse(met.signif$mrd.pvalue.t.test <= 0.001 | met.signif$mrd.pvalue.kruskal <= 0.001, 1, 0)
met.signif$rel.sig.flag <- ifelse(-log10(met.signif$relapse.pvalue.kruskal) >= 2.5 | -log10(met.signif$relapse.pvalue.t.test) >= 2.5, 1, 0)



# Plots from t-tests ------------------------------------------------------

require(ggplot2); require(ggrepel)

print(volcano.plot(met.signif, met.signif$fold.change.mrd, met.signif$mrd.pvalue.t.test, met.signif$super.pathway) +
        scale_color_discrete('Super Pathway') +
        geom_text_repel(aes(label=ifelse(-log10(mrd.pvalue.t.test) > 3, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        ggtitle('Compounds Associated With End-Induction Minimal Residual Disease') +
        xlab('Log2 Fold Change in MRD Positive Plasma') +
        ylab('-Log10 p-value'))

print(volcano.plot(met.signif, met.signif$fold.change.relapse, met.signif$relapse.pvalue.t.test, met.signif$super.pathway) +
        scale_color_discrete('Super Pathway') +
        geom_text_repel(aes(label=ifelse(-log10(relapse.pvalue.t.test) > 2.5, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        ggtitle('Metabolites Associated with Relapse') +
        xlab('Log2 Fold Change in Relapse Plasma') +
        ylab('-Log10 p-value'))



# Plots from K-W tests ----------------------------------------------------

require(ggplot2); require(ggrepel)

print(volcano.plot(met.signif, met.signif$fold.change.mrd, met.signif$mrd.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super Pathway') +
        geom_text_repel(aes(label=ifelse(-log10(mrd.pvalue.kruskal) > 3, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        ggtitle('Compounds Associated With End-Induction Minimal Residual Disease') +
        xlab('Log2 Fold Change in MRD Positive Plasma') +
        ylab('-Log10 p-value'))

print(volcano.plot(met.signif, met.signif$fold.change.relapse, met.signif$relapse.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super Pathway') +
        geom_text_repel(aes(label=ifelse(-log10(relapse.pvalue.kruskal) > 2.5, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        ggtitle('Metabolites Associated with Relapse') +
        xlab('Log2 Fold Change in Relapse Plasma') +
        ylab('-Log10 p-value'))



# Plots from both tests ---------------------------------------------------

require(ggplot2); require(ggrepel)

print(volcano.plot(met.signif, met.signif$fold.change.mrd, met.signif$mrd.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super pathway') +
        geom_text_repel(aes(label=ifelse(met.signif$mrd.sig.flag == 1, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        ggtitle('Compounds Associated With End-Induction Minimal Residual Disease') +
        xlab('Log2 fold change in plasma from children with subsequent MRD positivity') +
        ylab('-Log10 p-value (Kruskal-Wallis test)'))

print(volcano.plot(met.signif, met.signif$fold.change.relapse, met.signif$relapse.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super pathway') +
        geom_text_repel(aes(label=ifelse(met.signif$rel.sig.flag == 1, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        ggtitle('Metabolites Associated with Relapse') +
        xlab('Log2 fold change in plasma from children with subsequent relapse') +
        ylab('-Log10 p-value (Kruskall-Wallist test)'))



# Plots optimized for tiling ----------------------------------------------

require(ggplot2); require(ggrepel)

print(volcano.plot(met.signif, met.signif$fold.change.mrd, met.signif$mrd.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super pathway') +
        geom_text_repel(aes(label=ifelse(met.signif$mrd.sig.flag == 1, compound, ''),
                            color = super.pathway), cex = 5, force = 30) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
#        ggtitle('Compounds Associated With End-Induction Minimal Residual Disease') +
        xlab('Log2 fold change in plasma from children with subsequent MRD positivity') +
        theme(axis.title.x = element_text(hjust = 0)) +
        ylab('-Log10 p-value (Kruskal-Wallis test)'))

print(volcano.plot(met.signif, met.signif$fold.change.relapse, met.signif$relapse.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super pathway') +
        geom_text_repel(aes(label=ifelse(met.signif$rel.sig.flag == 1, compound, ''),
                            color = super.pathway), cex = 3) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 2.5), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
#        ggtitle('Metabolites Associated with Relapse') +
        xlab('Log2 fold change in plasma from children with subsequent relapse') +
        theme(axis.title.x = element_text(hjust = 0)) +
        theme(axis.title.y = element_blank())) 



# Tiled plots with optimized text jittering -------------------------------

require(ggplot2); require(ggrepel)

print(volcano.plot(met.signif, met.signif$fold.change.mrd, met.signif$mrd.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super pathway') +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-keto-3-deoxy-gluconate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, nudge_y = 0.6, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'gamma-carboxyglutamate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.9, nudge_y = 0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound %in% c('malate','fumarate','phosphate'), compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'alanine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'histidine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.4, nudge_y = -0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dihydroorotate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.75, nudge_y = -1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '4-hydroxy-2-oxoglutaric acid', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1.7, nudge_y = -0.75, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-hydroxyoctanoate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1.75, nudge_y = -0.25, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'succinate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.35, nudge_y = 0.35, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pyruvate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pyrraline', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "cytidine 5'-diphosphocholine", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.075, nudge_y = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-oleoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-linoleoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-palmitoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1, nudge_y = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-arachidonoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.1, nudge_y = 0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-stearoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'acetoacetate', compound, ''),
                           color = super.pathway), cex =8.5, nudge_x = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-arachidonylglycerol ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.5, show.legend = FALSE) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank()))

print(volcano.plot(met.signif, met.signif$fold.change.relapse, met.signif$relapse.pvalue.kruskal, met.signif$super.pathway) +
        scale_color_discrete('Super pathway') +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'valine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, nudge_y = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "13-methylmyristate", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1.5, nudge_y = -0.1, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'gamma-CEHC', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1.25, nudge_y = 0.1, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'betaine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.5, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'picolinate', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == '1-linoleoyl-GPI ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1, nudge_y = -0.1, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dopamine sulfate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_y = -0.2, show.legend = FALSE) +    
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'N-carboxyethylphenylalanine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, show.legend = FALSE) +   
        geom_text_repel(aes(label=ifelse(met.signif$compound == '10-nonadecenoate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, nudge_y = -0.2, show.legend = FALSE) +  
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'margarate ', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +  
        geom_text_repel(aes(label=ifelse(met.signif$compound == "15-methylpalmitate", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'palmitoleate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1.25, nudge_y = 0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'myristate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, nudge_y = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '10-heptadecenoate ', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dimethylglycine', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "inosine 5'-monophosphate ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.5, show.legend = FALSE) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 2.5), linetype = 'dotdash') +
        theme_bw() + 
        theme(text = element_text(size = 17.5)) +
        theme(axis.title.x = element_text(hjust = 0)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank())) 



# Tiled, optimized jittering, high contrast colors ------------------------

require(ggplot2); require(ggrepel)

print(ggplot(data = met.signif, aes(x = log2(met.signif$fold.change.mrd), y = -log10(met.signif$mrd.pvalue.kruskal))) + 
        geom_point(aes(color = met.signif$super.pathway)) +
        scale_color_manual(values = c('blue2','red','gold','violetred','springgreen3','purple','cyan','deepskyblue')) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() +
        theme(text = element_text(size = 17.5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank()) +
        scale_x_continuous(limits = c(-3,5)) +
        guides(colour = guide_legend(title = 'Super Pathway', override.aes = list(size=5))) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-keto-3-deoxy-gluconate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, nudge_y = 0.6, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'gamma-carboxyglutamate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.9, nudge_y = 0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound %in% c('malate','fumarate','phosphate'), compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'alanine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'histidine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.4, nudge_y = -0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dihydroorotate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.75, nudge_y = -1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '4-hydroxy-2-oxoglutaric acid', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1.7, nudge_y = -0.75, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '2-hydroxyoctanoate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1.75, nudge_y = -0.25, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'succinate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.35, nudge_y = 0.35, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pyruvate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'pyrraline', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "cytidine 5'-diphosphocholine", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.075, nudge_y = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-oleoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-linoleoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-palmitoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1, nudge_y = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-arachidonoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.1, nudge_y = 0.3, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-stearoyl-GPI ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'acetoacetate', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "1-arachidonylglycerol ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -0.5, show.legend = FALSE))

print(ggplot(data = met.signif, aes(x = log2(met.signif$fold.change.relapse), y = -log10(met.signif$relapse.pvalue.kruskal))) + 
        geom_point(aes(color = met.signif$super.pathway)) +
        scale_color_manual(values = c('blue2','red','gold','violetred','springgreen3','purple','cyan','deepskyblue')) +
        geom_vline(aes(xintercept = 0), linetype = 'dotdash') +
        geom_hline(aes(yintercept = 3), linetype = 'dotdash') +
        theme_bw() +
        theme(text = element_text(size = 17.5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank()) +
        scale_x_continuous(limits = c(-3,5)) +
        guides(colour = guide_legend(title = 'Super Pathway', override.aes = list(size=5))) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'valine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.2, nudge_y = -0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "13-methylmyristate", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1.5, nudge_y = -0.1, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'gamma-CEHC', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1.25, nudge_y = 0.1, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'betaine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.5, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'picolinate', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == '1-linoleoyl-GPI ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -1, nudge_y = -0.1, show.legend = FALSE) +        
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dopamine sulfate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_y = -0.2, show.legend = FALSE) +    
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'N-carboxyethylphenylalanine', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, show.legend = FALSE) +   
        geom_text_repel(aes(label=ifelse(met.signif$compound == '10-nonadecenoate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, nudge_y = -0.2, show.legend = FALSE) +  
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'margarate ', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +  
        geom_text_repel(aes(label=ifelse(met.signif$compound == "15-methylpalmitate", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'palmitoleate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 1.25, nudge_y = 0.1, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'myristate ', compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = -2, nudge_y = 0.2, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == '10-heptadecenoate ', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == 'dimethylglycine', compound, ''),
                            color = super.pathway), cex =8.5, show.legend = FALSE) +
        geom_text_repel(aes(label=ifelse(met.signif$compound == "inosine 5'-monophosphate ", compound, ''),
                            color = super.pathway), cex =8.5, nudge_x = 0.5, show.legend = FALSE))