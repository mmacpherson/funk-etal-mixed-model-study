## -- Setup.
library(afex) # -- provides mixed, wrapper around lme4
pval.method <- "KR" # - ie, Kenward-Roger


herb <- read.table("Funk_2010_herbivory_study.csv", h=T, sep=",")

traits <- c("phenolics", "toughness", "thickness", "leafN")
results <- c()
for (trait in traits) {

    herb["trait"] <- herb[trait]

    model1.trait <- lm(trait ~ origin, herb)
    model2.trait <- mixed(trait ~ origin + (1|pair),
                              herb,
                              check.contrasts=FALSE,
                              method=pval.method,
                              progress=FALSE)
    model3.trait <- mixed(trait ~ origin + (1|pair) + (1|species),
                              herb,
                              check.contrasts=FALSE,
                              method=pval.method,
                              progress=FALSE)

    vals <- c(

        ## -- pull out fixed effect due to origin
        as.numeric(model1.trait$coefficients[2]),
        model2.trait$full.model@beta[2],
        model3.trait$full.model@beta[2],

        ## -- pull out p.value
        summary(model1.trait)$coef[2,4],
        model2.trait$anova.table$p.value,
        model3.trait$anova.table$p.value

    )
    ## -- append results
    results <- rbind(results, t(matrix(vals, nrow=2, byrow=T)))
}

## -- rearrange into coherent data frame & write to disk
results <- as.data.frame(results)
names(results) <- c("origin", "p.value")

results$trait <- rep(traits, each=3)
results$model <- rep(c("model.one", "model.two", "model.three"), 4)

results <- results[,c(3,4,1,2)]
write.table(results, "funk_traits_model_fits.csv",quote=F,sep=",",row.names=F)
