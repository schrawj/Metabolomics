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

#' Function that condenses the code required to generate a volcano plot somewhat.
volcano.plot <- function(a, x, y, z){
  ggplot(data = a, aes(x = log2(x), y = -log10(y))) + 
    geom_point(aes(color = z))
}

require(ggplot2); require(ggrepel)

met.signif$super.pathway <- ifelse(met.signif$super.pathway == 'Cofactors and Vitamins', 'Cofactor/Vitamin', met.signif$super.pathway)



# Generate plots: from t-tests --------------------------------------------

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



# Generate plots: from Kruskal-Wallist tests ------------------------------

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
