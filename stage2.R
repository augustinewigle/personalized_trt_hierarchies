# Stage 2 analysis ---------------------------------------------------

library(reshape2)
library(tidyverse)

# Function ----------------------------------------
rank_indiv <- function(imat, trtmat, blippars, post_lm) {
  
  # Each row represents the same patient under diff trt
  X_indiv <- matrix(nrow = ncol(trtmat), ncol = length(blippars), dimnames = list(NULL, blippars))
  
  for(p in 1:ncol(imat)) {
    
    for(t in 1:ncol(trtmat)) {
      
      ix <- ((p-1)*(5)+1):((p)*5)
      
      X_indiv[t,ix] <- imat[,p]*trtmat[,t]
      
    }
    
  }
  
  # Add main effects to matrix
  X_indiv[1:ncol(trtmat), (max(ix)+1):ncol(X_indiv)] <- t(trtmat)
  
  # Create posterior prediction for relative outcome
  
  postpreds_lm <- matrix(nrow = nrow(post_lm), ncol = ncol(trtmat), dimnames = list(NULL, colnames(trtmat)))
  
  for(tr in 1:ncol(trtmat)) {
    
    postpreds_lm[,tr] <- apply(post_lm, 1, function(b, x) b%*%x, x = X_indiv[tr,])
    
  }
  
  # reshape for ggplot
  postdf_lm <- melt(postpreds_lm)
  
  # Make plot
  # Overlay poseriors
  plt <- ggplot(subset(postdf_lm, Var2 != "Sertraline"), aes(x = value, fill = Var2)) +
    geom_density(alpha = 0.25) +
    labs(fill = "Treatment vs. Sertraline", x = "Predicted outcome vs. Sertraline", title = "Q-learning, Patient X") +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0, lty = "dashed") 
  
  # Rank with SUCRA
  rankdf_lm <- postdf_lm %>% group_by(Var1) %>% mutate(ranking = rank(-value)) %>% group_by(Var2) %>% 
    summarise(Erank = mean(ranking)) %>% mutate(sucra = (6-Erank)/5) %>% ungroup() %>% mutate(sucra_rank = rank(-sucra))
  
  return(list(plot = plt, sucra = rankdf_lm))
  
}

# Data set-up ------------------------------------------------
#' Data needed to set up:
#' ns = number of studies
#' P = number of tailoring variables (be careful when it comes to multi-category variables)
#' ntrt = vector of the number of treatment arms in each study
#' G = number of treatments in network
#' nvar = (ntrt-1) x (P+1)
#' y = ns x max(nvar) matrix containing vectors of study-level estimates of the blip model params
#' Sigmahat = ns x max(nvar) x max(nvar) array containing study-level covariance matrices for blip model pars
#' X = ns x max(nvar) x (P+1)(G-1) array containing design matrix describing treatment structure of study i 
#' NOTE THIS IS NOT YET ADAPTED TO THE CASE WHERE ANY STUDIES DON'T INCLUDE THE OVERALL REFERENCE TREATMENT

ns <- 3
P <- 11
ntrt <- c(5, 2, 3) # stard, embarc, revamp
G <- 6
nvar <- (ntrt-1)*(P+1)

# Define stage 1 estimates (found by fitting linear models to confidential data)

stard_ests <- read.csv("coefficients_stard.csv")
revamp_ests <- read.csv("coefficients_revamp.csv")
embarc_ests <- read.csv("coefficients_embarc.csv")

stard_varobj <- read.csv("vcov_stard.csv")
stard_vcov <- as.matrix(stard_varobj[,-c(1)])
rownames(stard_vcov) <- colnames(stard_vcov)

embarc_varobj <- read.csv("vcov_embarc.csv")
embarc_vcov <- as.matrix(embarc_varobj[,-c(1)])
rownames(embarc_vcov) <- colnames(embarc_vcov)

revamp_varobj <- read.csv("vcov_revamp.csv")
revamp_vcov <- as.matrix(revamp_varobj[,-c(1)])
rownames(revamp_vcov) <- colnames(revamp_vcov)


# Define order of Phi first
allvarnames <- unique(c(stard_ests$param, revamp_ests$param, embarc_ests$param))
blippars <- sort(allvarnames[str_detect(allvarnames, "treatment")]) # Order of vars in Phi

# Create X based on Phi ordering
Xflm <- array(data = 0, c(ns, max(nvar), length(blippars)), dimnames = list(studies = c("STARD", "EMBARC", "REVAMP"), 
                                                                            arm = paste0("par", 1:max(nvar)), 
                                                                            var = blippars))
yflm <- matrix(NA, nrow = ns, ncol = max(nvar))
Sigmahatflm <- array(NA, dim = c(ns, max(nvar), max(nvar)))
lm_studies <- list(stard_ests,
                   embarc_ests, 
                   revamp_ests)
lm_vcov <- list(stard_vcov,
                embarc_vcov,
                revamp_vcov)

for(i in 1:ns) {
  
  blipix <- which(str_detect(lm_studies[[i]]$param, "treatment"))
  
  studyparslm <- lm_studies[[i]]$param[blipix]
  
  # Create X
  for(j in 1:length(blippars)) {
    
    if(blippars[j] %in% studyparslm) {
      
      ix <- which(studyparslm == blippars[j])
      
      Xflm[i,ix,j] <- 1
      
    }
    
  }
  
  # Create y
  yflm[i,1:nvar[i]] <- lm_studies[[i]]$value[blipix]
  
  # Create Sigmahat
  Sigmahatflm[i,1:nvar[i],1:nvar[i]] <- lm_vcov[[i]][blipix, blipix]
  
}

# Stage 2 - Fit models using Stan --------------------------------
library(rstan)
set.seed(2025)
initsf <- list(list(Phi = rnorm((P+1)*(G-1), 0, 100)),
               list(Phi = runif((P+1)*(G-1), -10, 10)),
               list(Phi = rnorm((P+1)*(G-1), 0, 10)),
               list(Phi = runif((P+1)*(G-1), -20, 20)))

dat_lmf <- list(ns = ns,
                nvar = nvar,
                P = P, 
                G = G,
                y1 = yflm[1,1:nvar[1]],
                Sigmahat1 = Sigmahatflm[1, 1:nvar[1], 1:nvar[1]],
                X1 = Xflm[1, 1:nvar[1], 1:((G-1)*(P+1))],
                y2 = yflm[2,1:nvar[2]],
                Sigmahat2 = Sigmahatflm[2, 1:nvar[2], 1:nvar[2]],
                X2 = Xflm[2, 1:nvar[2], 1:((G-1)*(P+1))],
                y3 = yflm[3,1:nvar[3]],
                Sigmahat3 = Sigmahatflm[3, 1:nvar[3], 1:nvar[3]],
                X3 = Xflm[3, 1:nvar[3], 1:((G-1)*(P+1))],
                prior_sd = sqrt(1/0.0001))

# Synthesize q-learning results
fit_lm <- stan(file = "fullcovarmatrix.stan", pars = "phi", 
               iter = 5000,
               data = dat_lmf, seed = 26)

traceplot(fit_lm, pars = "phi")

post_lm <- extract(fit_lm, pars = "phi")$phi
colnames(post_lm) <- blippars

# Do individualized ranking stuff ------------------------------------------------------

# Look at coefficients to help pick covariate profiles
lm_s <- summary(fit_lm, pars = "phi")$summary

res <- data.frame(lm_s)
res$param <- blippars

res <- res %>% mutate(pargroup = str_extract(param, "[^:]*"))
res$pargroup[startsWith(res$pargroup, "treatment")] <- "maineff"
res$treatment <- str_extract(res$param, "treatment.*")
res$trt2 <- str_replace(res$treatment, "treatment", "")
res$trt3 <- str_replace_all(res$trt2, c("BUP" = "Bupropion",
                                        "ESCIT" = "Escitalopram",
                                             "CIT" = "Citalopram",
                                             "BUS" = "Buspirone",
                                             "VEN" = "Venlafaxine"))
res$pargroup2 <- str_replace_all(res$pargroup,
                                 c("age" = "Age",
                                   "baseline_hrsd" = "Baseline HRSD-17",
                                   "chronicTRUE" = "Chronic episode",
                                   "edutot" = "Years education",
                                   "emp_statusUnemployed/Retired" = "Unemployed/Retired",
                                   "marital_statusDivorced/Widowed" = "Divorced/Widowed",
                                   "marital_statusSingle" = "Single",
                                   "morethan3epi" = "Had > 3 episodes",
                                   "onsetAge" = "Age at onset",
                                   "sexM" = "Male",
                                   "thous" = "Household size",
                                   "maineff" = "Main effect"))

ggplot(res, aes(x = trt3, y = mean, colour = trt3)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.3) +
  geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width = 0.75, linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~pargroup2) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "Treatment", color = "Treatment", y = "Estimate")

# create a vector representing a "reference" individual
individ <- data.frame(age = 0, 
                      baseline_hrsd = 0,
                      chronic = 0,
                      edutot = 0,
                      employment_statusUnemployedRetired = 0,
                      martial_statusDivorcedWidowed = 0,
                      martial_statusSingle = 0,
                      morethan3epi = 0,
                      onsetage = 0,
                      sexM = 1,
                      thous = 1)

# create a vector representing someone who is older and has somewhat more severe case
ind2 <- data.frame(age = 0, 
                   baseline_hrsd = 0,
                   chronic = 0,
                   edutot = 0,
                   employment_statusUnemployedRetired = 1,
                   martial_statusDivorcedWidowed = 0,
                   martial_statusSingle = 0,
                   morethan3epi = 1,
                   onsetage = 0,
                   sexM = 0,
                   thous = -1)

# treatments
trts <- data.frame(Bupropion    = c(1,0,0,0,0),
                   `Citalopram+Bupropion` = c(0,1,0,0,0),
                   `Citalopram+Buspirone` = c(0,0,1,0,0),
                   Escitalopram  = c(0,0,0,1,0),
                   Venlafaxine    = c(0,0,0,0,1),
                   Sertraline    = c(0,0,0,0,0))


pat0 <- rank_indiv(imat = individ, 
                   trtmat = trts, blippars = blippars, 
                   post_lm = post_lm)

pat1 <- rank_indiv(imat = ind2, 
                   trtmat = trts, blippars = blippars, post_lm = post_lm)

pA <-pat0$plot+coord_cartesian(ylim = c(0,0.2), xlim = c(-18,18))+
  labs(title = "Patient A", x = "Expected Outcome vs. Sertraline")+
  theme_bw() +
  theme(legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
pB<- pat1$plot+coord_cartesian(ylim = c(0,0.2), xlim = c(-18,18)) + 
  labs(title = "Patient B", x = "Expected Outcome vs. Sertraline") +
  theme_bw()

library(ggpubr)

ggarrange(pA, pB, nrow = 2, common.legend = TRUE, legend = "right")
ggsave("pApB.png", width = 7, height = 8, units = "in", dpi = 330)

# Do it naively - just rank psi zeros from postpred
naive_post <- post_lm[,c("treatmentBUP", 
                          "treatmentCIT+BUP", 
                          "treatmentCIT+BUS",
                          "treatmentESCIT",
                          "treatmentVEN")]
naive_post <- cbind(naive_post, treatmentSER = rep(0, nrow(naive_post)))
naive_ranks <- t(apply(-naive_post, 1, rank))
naive_sucras <- (6-apply(naive_ranks, 2, mean))/5

cbind(pat0$sucra, pat1$sucra)


print(xtable::xtable(data.frame(pat0$sucra$Var2, naive_sucras, rank(-naive_sucras),pat0$sucra$sucra, pat0$sucra$sucra_rank,
                                pat1$sucra$sucra, pat1$sucra$sucra_rank)), include.rownames = F)
