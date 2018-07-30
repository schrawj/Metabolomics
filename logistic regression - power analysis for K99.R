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


# Scratch paper -----------------------------------------------------------

repetitions <- 5000
significant <- matrix(nrow=repetitions, ncol=1)
n <- 45
fold.change <- log(1.5)
prop.mrd.neg.nci.hr <- 0.33
prop.mrd.pos.nci.hr <- 0.67
prop.mrd.neg.unfav.cyto <- 0.1
prop.mrd.pos.unfav.cyto <- 0.25

for (i in 1:repetitions){
  mrd.neg <- data.frame(sim.metab = rlnorm(n, meanlog = 0, sdlog= 0.5),
                        sim.nci = rbinom(n, 1, prop.mrd.neg.nci.hr), 
                        sim.cyto = rbinom(n, 1, 0.1),
                        sim.immuno = rbinom(n, 1, 0.08),
                        mrd = 0)
  mrd.pos <- data.frame(sim.metab = rlnorm(n, meanlog = fold.change, sdlog = 0.75),
                        sim.nci = rbinom(n, 1, prop.mrd.pos.nci.hr), 
                        sim.cyto = rbinom(n, 1, 0.25),
                        sim.immuno = rbinom(n, 1, 0.12),
                        mrd = 1)
  sim.dat <- rbind(mrd.neg, mrd.pos)
  
  fit <- glm(mrd ~ sim.metab + sim.nci + sim.cyto + sim.immuno, data = sim.dat, family = binomial (link = 'logit'))
  
  significant[i, ] <- summary(fit)$coefficients[1,4]
  
}

(table(significant[,1] <= 0.05)[2]/5000)*100

mydat <- data.frame( v1 = rep( c(3,6,9), each=2 ),
                     v2 = rep( 0:1, 3 ), 
                     resp=c(0.0025, 0.00395, 0.003, 0.0042, 0.0035, 0.002) )

fit0 <- glm( resp ~ poly(v1, 2, raw=TRUE)*v2, data=mydat,
             weight=rep(100000,6), family=binomial)
b0 <- coef(fit0)


simfunc <- function( beta=b0, n=10000 ) {
  w <- sample(1:6, n, replace=TRUE, prob=c(3, rep(1,5))) #' Tells R to sample of size n from rows 1:6, with replacement, with 3x the weight for the first element.  The author of this problem specified that row 1 in mydat represented a condition that should be overweighted 3:1.
  mydat2 <- mydat[w, 1:2] #' creates a synthetic dataset with 10,000 rows with the combinations of v1 and v2 weighted as in the statement above.
  eta <- with(mydat2,  cbind( 1, v1, 
                              v1^2, v2,
                              v1*v2,
                              v1^2*v2 ) %*% beta ) #' a matrix with one column.  Appears to be log-scale response rate calculated by multiplying the combinations of variables by their betas?
  p <- exp(eta)/(1+exp(eta)) #' converts from log-scale to odds?
  mydat2$resp <- rbinom(n, 1, p) #' randomly generates 10,000 binomial random variables with probability of taking 1 equal to p.
  
  fit1 <- glm( resp ~ poly(v1, 2)*v2, data=mydat2,
               family=binomial)
  fit2 <- update(fit1, .~ poly(v1,2) )
  anova(fit1,fit2, test='Chisq')[2,5]
}

out <- replicate(100, simfunc(b0, 10000))
mean( out <= 0.05 )
hist(out)
abline(v=0.05, col='lightgrey')



# Example of post-hoc power analysis for log reg --------------------------

#' Taken from: https://stats.stackexchange.com/questions/35940/simulation-of-logistic-regression-power-analysis-designed-experiments/36040#36040

set.seed(1)

repetitions = 1000
N = 10000
n = N/8
var1  = c(   .03,    .03,    .03,    .03,    .06,    .06,    .09,   .09)
var2  = c(     0,      0,      0,      1,      0,      1,      0,     1)
rates = c(0.0025, 0.0025, 0.0025, 0.00395, 0.003, 0.0042, 0.0035, 0.002)

var1    = rep(var1, times=n)
var2    = rep(var2, times=n)
var12   = var1**2
var1x2  = var1 *var2
var12x2 = var12*var2

significant = matrix(nrow=repetitions, ncol=7)

startT = proc.time()[3]
for(i in 1:repetitions){
  responses          = rbinom(n=N, size=1, prob=rates)
  model              = glm(responses~var1+var2+var12+var1x2+var12x2, 
                           family=binomial(link="logit"))
  significant[i,1:5] = (summary(model)$coefficients[2:6,4]<.05)
  significant[i,6]   = sum(significant[i,1:5])
  modelDev           = model$null.deviance-model$deviance
  significant[i,7]   = (1-pchisq(modelDev, 5))<.05
}
endT = proc.time()[3]
endT-startT

sum(significant[,1])/repetitions      # pre-specified effect power for var1
sum(significant[,2])/repetitions      # pre-specified effect power for var2
sum(significant[,3])/repetitions      # pre-specified effect power for var12
sum(significant[,4])/repetitions      # pre-specified effect power for var1X2
sum(significant[,5])/repetitions      # pre-specified effect power for var12X2
sum(significant[,7])/repetitions      # power for likelihood ratio test of model
sum(significant[,6]==5)/repetitions   # all effects power
sum(significant[,6]>0)/repetitions    # any effect power
sum(significant[,4]&significant[,5])/repetitions   # power for interaction terms



# Scratch paper -----------------------------------------------------------

mydat2 <- mydat[sample(1:6, 10000, replace=TRUE, prob=c(3, rep(1,5))), 1:2]

eta <- with(mydat2,  cbind( 1, v1, 
                            v1^2, v2,
                            v1*v2,
                            v1^2*v2 ) %*% b0 )
p <- exp(eta)/(1+exp(eta)) 
mydat2$resp <- rbinom(10000, 1, p)



# Generate synthetic datasets based on metabolite distributions -----------

load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20180306.1.rdata")

require(ggplot2)

p <- ggplot(data = met) + geom_density(aes(x = pyruvate)) + facet_wrap(~mrd)
print(p)

rlnorm()

aggregate(pyruvate ~ mrd, data = met, mean)
aggregate(pyruvate ~ mrd, data = met, sd)

tmp <- data.frame(sim.data = c(rlnorm(50, meanlog = 0, sdlog= 0.5),
                               rlnorm(50, meanlog = 0.75, sdlog = 0.6)),
                  mrd = rep(0:1, each = 50))
print(ggplot(data = tmp) + geom_density(aes(x=sim.data)) + facet_wrap(~mrd) + xlim(0, 5))

fit <- glm(mrd ~ sim.data, data = tmp, family = binomial (link = 'logit'))
summary(fit)

mrd.neg <- data.frame(sim.metab = rlnorm(45, meanlog = 0, sdlog= 0.5),
                      sim.nci = rbinom(45, 1, 0.33), 
                      sim.cyto = rbinom(45, 1, 0.1),
                      sim.immuno = rbinom(45, 1, 0.08),
                      mrd = 0)
mrd.pos <- data.frame(sim.metab = rlnorm(45, meanlog = 0.5, sdlog = 0.75),
                      sim.nci = rbinom(45, 1, 0.66), 
                      sim.cyto = rbinom(45, 1, 0.25),
                      sim.immuno = rbinom(45, 1, 0.12),
                      mrd = 1)
sim.dat <- rbind(mrd.neg, mrd.pos)

print(ggplot(data = sim.dat) + geom_density(aes(x=sim.metab)) + facet_wrap(~mrd) + xlim(0, 5))

fit <- glm(mrd ~ sim.metab + sim.nci + sim.cyto + sim.immuno, data = sim.dat, family = binomial (link = 'logit'))
summary(fit)
exp(fit$coefficients)

