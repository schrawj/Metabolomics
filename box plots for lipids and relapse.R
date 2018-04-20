#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2018.02.07.
#' 
#' Generate some boxplots showcasing decreased lipid concentrations among relapsed 
#' patients.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20180306.1.rdata")

par(mfrow=c(3,3),         #' Sets number of rows and columns
    cex.main = 2.25,      #' Specify print title at 2.25x magnification for subsequent plots
    cex.axis = 2)         #' Specify print axis test at 2x magnification for subsequent plots

boxplot(met$`arachidate (20:0)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Arachidate (20:0); p=0.01',
        ylim = c(0, 3))
points(factor(met$relapse), met$`arachidate (20:0)`)


boxplot(met$`linoleate (18:2n6)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Linoleate (18:2n6); p=0.01',
        ylim = c(0, 3))
points(factor(met$relapse), met$`linoleate (18:2n6)`)


boxplot(met$`linolenate [alpha or gamma; (18:3n3 or 6)]` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Linolenate (18:3n3 or n6); p=0.02',
        ylim = c(0, 5))
points(factor(met$relapse), met$`linolenate [alpha or gamma; (18:3n3 or 6)]`)


boxplot(met$`stearate (18:0)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Stearate (18:0); p=0.02')
points(factor(met$relapse), met$`stearate (18:0)`)


boxplot(met$`dihomo-linoleate (20:2n6)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Dihomo-linoleate (20:2n6); p=0.04',
        ylim = c(0, 5))
points(factor(met$relapse), met$`dihomo-linoleate (20:2n6)`)


boxplot(met$`palmitate (16:0)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Palmitate (16:0), p=0.03')
points(factor(met$relapse), met$`palmitate (16:0)`)


boxplot(met$`oleate/vaccenate (18:1)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Oleate/Vaccenate (18:1); p=0.04')
points(factor(met$relapse), met$`oleate/vaccenate (18:1)`)


boxplot(met$`eicosenoate (20:1)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Eicosanoate (20:1); p=0.04',
        ylim = c(0, 5))
points(factor(met$relapse), met$`eicosenoate (20:1)`)


boxplot(met$`dihomo-linolenate (20:3n3 or n6)` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Dihomo-linolenate (20:3n3 or n6); p=0.04',
        ylim = c(0, 3))
points(factor(met$relapse), met$`dihomo-linolenate (20:3n3 or n6)`)







boxplot(met$`behenate (22:0)*` ~ met$relapse, 
        col = c('gold', 'darkgreen'), 
        main = 'Behenate (22:0); p=0.05',
        ylim = c(0, 4))
points(factor(met$relapse), met$`behenate (22:0)*`)

