library(ape)          # -- provides rcoal, rtree
library(apTreeshape)  # -- provides rtreeshape
library(geiger)       # -- provides sim.char
library(lme4)         # -- provides lmer
library(afex)         # -- provides mixed, for p-values

## -- ensure reproducibility
set.seed(2009)

## -- process CLI args
ca <- commandArgs()[-(1:5)]
ca <- as.numeric(ca)

pair.var     <- ca[1]
species.var  <- ca[2]
residual.var <- ca[3]
effect.size  <- ca[4]
nsim         <- as.integer(ca[5])

pval.method <- "KR" # -- use Kenward-Roger for mixed models

## -- simulation parameters
npairs      <- 30
nspecies    <- 2 * npairs
nreplicates <- 5

pair.sd     <- sqrt(pair.var)
pair.sd     <- 1.0
species.sd  <- sqrt(species.var)
residual.sd <- sqrt(residual.var)

## -- origin offset vector
vorigin <- c(0.0, effect.size)

## -- Extract same parameters from lmer model as summary()$coefs does from lm
lmer.summary <- function(fitted) {
    fitted.est <- fixef(fitted)[2]
    fitted.sd  <- sqrt(diag(vcov(fitted)))[2]
    fitted.t   <- fitted.est / fitted.sd
    fitted.p   <- 2 * pnorm(-abs(fitted.t))
    c(fitted.est, fitted.sd, fitted.t, fitted.p)
}

## -- Output matrix containing estimates, variances, tests and p-values for the three methods
outmat <- matrix(0, nsim, 12)

for (k in 1:nsim)
{
    ## -- generate Yule phylogenetic tree
    tree <- as.phylo(rtreeshape(1, npairs, model="yule")[[1]])

    ## -- pair offset vector, via evolutionary tree simulation
    vpair <- sim.char(tree, matrix(c(pair.sd ** 2)))[1:npairs,1,1]
    vpair <- vpair * sqrt(pair.var) / sd(vpair)

    ## -- species offset vector
    vspecies <- rnorm(nspecies, sd=species.sd)

    ## -- allocate vectors to hold the simulated data
    value     <- c()
    origin    <- c()
    pair      <- c()
    replicate <- c()
    species   <- c()

    ## -- simulate dataset!
    for (.origin in 1:2) {
        for (.pair in 1:npairs) {
            for (.replicate in 1:nreplicates) {
                species.index <- 1 + 2 * (.pair - 1) + (.origin - 1)
                components <- c(vorigin[.origin],        # -- fixed effect of origin
                                vpair[.pair],            # -- random effect of pair
                                vspecies[species.index], # -- random effect of species
                                rnorm(1, sd=residual.sd) # -- residual (random) effect
                                )

                ## -- augment vectors
                value     <- c(value, sum(components))
                origin    <- c(origin, .origin)
                pair      <- c(pair, .pair)
                replicate <- c(replicate, .replicate)
                species   <- c(species, species.index)
            }
        }
    }

    origin <- factor(origin) # -- or else afex::mixed will complain

    model1 <- lm(value ~ origin)
    data <- data.frame(value=value, origin=origin, pair=pair, species=species)
    model2 <- mixed(value ~ origin + (1|pair),
                    data,
                    check.contrasts=FALSE,
                    method=pval.method, progress=FALSE)
    model3 <- mixed(value ~ origin + (1|pair) + (1|species),
                    data,
                    check.contrasts=FALSE,
                    method=pval.method, progress=FALSE)

    outmat[k,1:4]  <- summary(model1)$coef[2,]
    outmat[k,5:8]  <- c(0,0,0,model2$anova.table$p.value)
    outmat[k,9:12] <- c(0,0,0,model3$anova.table$p.value)

}

# -- Compute Type I error when effect size is zero, power when it's not
m1 <- sum(outmat[,4] < 0.05) / nsim
m2 <- sum(outmat[,8] < 0.05) / nsim
m3 <- sum(outmat[,12] < 0.05) / nsim

# -- Summarize
c(pair.var, species.var, residual.var, effect.size, m1, m2, m3)
