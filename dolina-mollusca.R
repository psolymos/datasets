## rbind egyeles es avarminta, method valtozo + random hatas -- glmm
## az adattabla adat a csomagban

library(mefa4)
source("~/repos/opticut/R/opticut.R")
setwd("~/Dropbox/mollusca/dolina-zita/minden/")
mc <- read.csv("dol-moll-count.csv")
ms <- read.csv("dol-moll-samp.csv")
mt <- read.csv("dol-moll-taxa.csv")
md <- read.csv("dolina-tobor.csv")

#mc <- data.frame(mc, mt[match(mc$species, mt$codes),])
ms2 <- nonDuplicated(mc, subsample, TRUE)
mss <- data.frame(ms2, ms[match(ms2$sample, ms$sample.id),])
mss <- data.frame(mss, md[match(mss$dolina, md$Tobor.szama),])
tmp <- as.character(mss$subsample)
mss$method <- as.factor(substr(tmp, nchar(tmp), nchar(tmp)))
tmp <- as.character(mss$subsample)
mss$method <- as.factor(substr(tmp, nchar(tmp), nchar(tmp)))


mss2 <- nonDuplicated(mss, sample, TRUE)

excl <- c("NONE", levels(mc$species)[grep("indet", levels(mc$species))])

ym <- as.matrix(Xtab(catch ~ sample + species, mc,
    subset=mc$segment %in% c("a","j","f"), cdrop=excl))
ya <- as.matrix(Xtab(catch ~ sample + species, mc,
    subset=mc$segment %in% c("a","j","f") & mc$method == "Q", cdrop=excl))
ye <- as.matrix(Xtab(catch ~ sample + species, mc,
    subset=mc$segment %in% c("a","j","f") & mc$method == "T", cdrop=excl))

yyy <- as.matrix(Xtab(catch ~ subsample + species, mc,
    subset=mc$segment %in% c("a","j","f"), cdrop=excl))
mmm <- Mefa(yyy, mss, mt)

ckeep <- c("sample","dolina","microhab", "mhab", "method","aspect","stratum",
    "lmoist","lthick")
dolina <- list(
    xtab=as.matrix(xtab(mmm)),
    samp=samp(mmm)[,ckeep],
    taxa=taxa(mmm)[,c("scientific.name","family")])
table(dolina$samp$stratum, as.integer(dolina$samp$stratum))
levels(dolina$samp$stratum) <- c("1bottom", "3edge", "2middle", "4outside")
dolina$samp$stratum <- as.factor(as.character(dolina$samp$stratum))
table(dolina$samp$stratum, as.integer(dolina$samp$stratum))

table(dolina$samp$mhab, dolina$samp$microhab)
# "A" "H" "L" "R"
levels(dolina$samp$mhab) <- c("LI", "DW", "TL", "RO")
table(dolina$samp$mhab, dolina$samp$microhab)

## data set for opticut package
save(dolina, file="~/repos/opticut/data/dolina.rda")

## standard example

library(opticut)
data(dolina)
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
Y <- dolina$xtab[,colSums(dolina$xtab > 0) >= 5]
## opticut results
dol <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb", comb="all")
summary(dol)
plot(dol)

dol2 <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb", comb="rank")
summary(dol2)
plot(dol2,which=4)

