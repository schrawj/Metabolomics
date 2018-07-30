#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.30.
#' 
#' Some code that will generate synthetic datasets that resemble our 
#' metabolomics discovery set, and perform power analysis via simulation
#' based on user-supplied parameters.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------


# Power analysis for K99 aim 1 --------------------------------------------

repetitions <- 5000
significant <- matrix(nrow=repetitions, ncol=1)
fold.change <- log(1.25)
n <- 50
prop.mrd.neg.nci.hr <- 0.33
prop.mrd.pos.nci.hr <- 0.67
prop.mrd.neg.unfav.cyto <- 0.1
prop.mrd.pos.unfav.cyto <- 0.25

for (i in 1:repetitions){
  mrd.neg <- data.frame(sim.metab = rlnorm(n, meanlog = 0, sdlog= 0.5),
                        sim.nci = rbinom(n, 1, 0.33), 
                        sim.cyto = rbinom(n, 1, 0.1),
                        sim.immuno = rbinom(n, 1, 0.08),
                        mrd = 0)
  mrd.pos <- data.frame(sim.metab = rlnorm(n, meanlog = fold.change, sdlog = 0.75),
                        sim.nci = rbinom(n, 1, 0.66), 
                        sim.cyto = rbinom(n, 1, 0.25),
                        sim.immuno = rbinom(n, 1, 0.12),
                        mrd = 1)
  sim.dat <- rbind(mrd.neg, mrd.pos)
  
  fit <- glm(mrd ~ sim.metab + sim.nci + sim.cyto + sim.immuno, data = sim.dat, family = binomial (link = 'logit'))

  significant[i, ] <- summary(fit)$coefficients[1,4]

}

(table(significant[,1] <= 0.05)[2]/5000)*100

