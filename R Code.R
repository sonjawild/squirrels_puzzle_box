#### Behavioral type predicts three stages of anthropogenic resource discovery and exploitation in California ground squirrels #####

# two notes for better understanding:
 # More disturbed population = CROW
 # Less disturbed population = PARADISE

# The variable 'mobility' in the model is called 'exploration' in the data set

# 1) Data prep ------------------------------------------------------------

# 1a) Install packages and load data ---------------------------------------


# packages 
#install.packages("survival")
library(survival)
library(car)
library(tidyverse)

library(broom)
library(officer)
library(flextable)


# set working directory
setwd("Data/")


# load latency data - solving
latency.data.solve <- read.csv("latency.data.solve.csv", row.names = 1)

# latency until discovery
latency.data.discovery <- read.csv("latency.data.discovery.csv", row.names = 1)

names(latency.data.solve)

# Variables of interest:
# id: squirrel identity
# timestamp: date and time of third solve for those with at least three solve (or time of discovery for discovery df)
# time: latency to third solve experimental hours 
# training.solves: how many solves did the squirrel perform during the training phase (with additional peanut butter on the lever). controls for puzzle exposure prior to trial phase
# solve.rate.day: solves/day detected on puzzle box AFTER learning to solve 
# total solves: total number of solves for knowledgeable individuals
# total.arrivals: number of visits to the puzzle box AFTER learning to solve
# n_observation_opportunities: number of times a focal individual had opportunity to watch another squirrel perform a solve. Opportunity = either if present at the same box at the same time, or focal arrived within 10s of a solve occuring.
# eigenvector: measure of sociability (eigenvector centrality based on co-feeding)
# exploration: median number of RFID feeders visited during an experimental day (out of a total of 18 possible) - called 'mobility' in MS
# n_days_net: on how many days with RFID feeders was this squirrel registered
# log_num_trapped_per_day: log number of trapping events per trapping day (average)
# any_beh_prop: proportion of trapping events the squirrel showed a fear response (chatter, struggle, call)
# age: p for juveniles, a for adults
# sex: M for male, F for female
# colony: Crow or Paradise
# prop_days_trapped: proportion of trapping days the squirrel was trapped
# bold: boldness score (from PCA on trappability and fear responses in trap)


# 1b) Calculate boldness --------------------------------------------------

# We provide the boldness score in the above data frame, but here present the raw data and results from the PCA.

# trap.behav.sum <- read.csv("trap.behav.sum.csv")
# 
# 
# # run PCA
# 
# trap_data_scaled <- scale(trap.behav.sum[,c("any_beh_prop", "log_num_trapped_per_day")])
# 
# # load library
# library(psych)
# # Run PCA
# pca_bh <- principal(trap_data_scaled, nfactors = 1, rotate = "none")
# summary(pca_bh)
# 
# # Factor analysis with Call: principal(r = trap_data_scaled, nfactors = 1, rotate = "none")
# # 
# # Test of the hypothesis that 1 factor is sufficient.
# # The degrees of freedom for the model is -1  and the objective function was  0.33 
# # The number of observations was  78  with Chi Square =  24.47  with prob <  NA 
# # 
# # The root mean square of the residuals (RMSA) is  0.32 
# 
# pca_bh$loadings
# 
# # Loadings:
# #   PC1   
# # any_beh_prop            -0.827
# # log_num_trapped_per_day  0.827
# # 
# # PC1
# # SS loadings    1.368
# # Proportion Var 0.684
# 
# # we can see that trap behavior and trappability load in opposite directions and account for 68.4% of the variance.
# # Higher PC scores means more trappable and fewer fear responses in the trap
# 
# # add the scores to the data frame
# trap.behav.sum$bold <- pca_bh$scores[,1]
# 
# # this score corresponds to the 'bold' column in the latency data frames (for both discovery and solving)
# 
# # note that 12 individuals who were never trapped during summer 2024 received the minimum boldness score. We here provide a list of the 12:
# 
# # id log_num_trapped_per_day any_beh_prop age sex   colony prop_days_trapped      bold
# # 1  BF648B1                      NA           NA   A   F Paradise                 0 -2.141668
# # 2  BF64922                      NA           NA   A   F Paradise                 0 -2.141668
# # 3  D47CB1A                      NA           NA   A   M Paradise                 0 -2.141668
# # 4  D47CAFC                      NA           NA   A   F     Crow                 0 -2.141668
# # 5  D47CAE5                      NA           NA   A   F Paradise                 0 -2.141668
# # 6  D47CAFF                      NA           NA   A   M     Crow                 0 -2.141668
# # 7  D47CAD8                      NA           NA   A   M Paradise                 0 -2.141668
# # 8  D47CBAD                      NA           NA   A   M Paradise                 0 -2.141668
# # 9  D47CB5C                      NA           NA   A   F Paradise                 0 -2.141668
# # 10 D47CB7B                      NA           NA   A   F     Crow                 0 -2.141668
# # 11 D47CB54                      NA           NA   A   F Paradise                 0 -2.141668
# # 12 D47CB3B                      NA           NA   A   F     Crow                 0 -2.141668

# 1c) Sociability (eigenvector centrality) --------------------------------

# Since computation of social networks is computationally intense, we have already added the eigenvector centrality as a column to the above data frames. We here provide the code and raw objects needed to run gaussian mixture models, create networks and extract network position. 

# # Read network data:
# 
# Crow.net.data <- read.csv("Crow.net.data.csv")
# 
# Paradise.net.data <- read.csv("Paradise.net.data.csv")
# 
# # load libraries
# library(asnipe)
# library(igraph)
# 
# # this fuction generates the co-feeding groups and extracts the eigenvector centralities
# get_network_data <- function(net.data, IDs){
#   
#   # subset to indidividuals
#   net.data <- net.data[net.data$PIT  %in% names(IDs),]
#   net.data$date_time <- as.numeric(net.data$date_time)
#   
#   gmm <- gmmevents(
#     time = net.data$date_time,
#     identity = net.data$PIT,
#     location = net.data$date_logger_antenna,
#     splitGroups = TRUE,
#     verbose = TRUE
#   )
#   
#   
#   net <- get_network(
#     association_data = gmm$gbi,
#     data_format = "GBI",
#     association_index = "SRI",
#     times = gmm$metadata$Start,
#     identities = colnames(gmm$gbi)
#   )
#   
#   # how to extract degree from net object
#   g.net <- graph_from_adjacency_matrix(net, mode = "undirected",
#                                        weighted = TRUE, diag = TRUE)
#   
#   eigen <- eigen_centrality(g.net, directed=F, weights=E(g.net))
#   
#   
#   object <- NULL
#   object$gmm <- gmm
#   object$net <- net
#   object$eigen <- eigen
#   
#   return(object)
# }
# 
# Crow.discovery.net <- get_network_data(net.data = Crow.net.data, IDs = Crow.ids)
# #save(Crow.discovery.net, file="Crow.discovery.net.RData")
# load("Crow.discovery.net.RData")
# 
# Paradise.discovery.net <- get_network_data(net.data = Paradise.net.data, IDs = Paradise.ids)
# save(Paradise.discovery.net.RData")
#load("Paradise.discovery.net.RData")
# 
# # these objects contain three slots:
# # gmm: gmm object resulting from the gmm_events function (gaussian mixture model)
# # net: the social network based on the simple ratio index
# # eigen: eigenvector centralities for all individuals in the network (this value corresponds to the column 'eigen' in the data frames loaded above)


# 2) VIFs discovery -----------------------------------------------------------------

# calculate variance inflation factors to test for multi-collinearity among predictors

# for sex and age, we convert to a numeric (0/1) instead of factor to be able to calculate correlations
# juveniles = 1; adults = 0 (baseline)
# males = 1; females = 0 (baseline)

latency.data.discovery$age_num <- ifelse(
  latency.data.discovery$age =="P", 1, 0
)

latency.data.discovery$sex_num <- ifelse(
  latency.data.discovery$sex =="M", 1, 0
)

latency.data.discovery$colony_num <- ifelse(
  latency.data.discovery$colony =="Crow", 1, 0
)


vif.discovery <- lm(time ~ mobility + bold + eigenvector + age_num + sex_num + other_present + colony_num, data = latency.data.discovery)

vif(vif.discovery)

# mobility          bold   eigenvector       age_num       sex_num other_present    colony_num 
# 2.038491      1.103433      2.045622      1.028957      1.021955      1.047715      1.111295

cor(latency.data.discovery %>%
      dplyr::select(mobility, eigenvector, bold, age_num, sex_num, other_present, colony_num))

#                 mobility eigenvector        bold    age_num     sex_num other_present  colony_num
# mobility      1.00000000  0.72823836  0.20566795 0.09179018  0.02643648    0.10492018  0.04829580
# eigenvector   0.72823836  1.00000000  0.14128935 0.16846391  0.04315227    0.06241689  0.06280325
# bold          0.20566795  0.14128935  1.00000000 0.07944433 -0.06077849    0.13524148  0.31710518
# age_num       0.09179018  0.16846391  0.07944433 1.00000000  0.02866440    0.02806675  0.01972676
# sex_num       0.02643648  0.04315227 -0.06077849 0.02866440  1.00000000   -0.04809861 -0.02471211
# other_present 0.10492018  0.06241689  0.13524148 0.02806675 -0.04809861    1.00000000  0.20731879
# colony_num    0.04829580  0.06280325  0.31710518 0.01972676 -0.02471211    0.20731879  1.00000000


# 3) Latency - discovery --------------------------------------------------

# create a survival object
# those that did not discover the resource are right censored (+) sign
latency.data.discovery$time_discover <- latency.data.discovery$time
latency.data.discovery$surv <- with(latency.data.discovery, Surv(time_discover, !is.na(timestamp)))

# scale some variables
latency.data.discovery[c("mobility", "eigenvector", "bold")] <- lapply(latency.data.discovery[c("mobility", "eigenvector", "bold")],scale)

# Fit the Cox proportional hazards model
fit_discover <- coxph(surv ~ mobility + bold + eigenvector + age_num + sex_num + other_present + colony_num, data = latency.data.discovery)


# check for proportional hazards assumption
# we need to make sure that the effects of a predictor does not change over time
# significant p values indicate a violation of the proportional hazards assumption

cox.zph(fit_discover)
# chisq df       p
# mobility       5.2919  1 0.02142
# bold           0.3937  1 0.53034
# eigenvector    5.0045  1 0.02528
# age_num        8.7349  1 0.00312
# sex_num        0.1075  1 0.74297
# other_present  8.7929  1 0.00302
# colony_num     0.0653  1 0.79825
# GLOBAL        25.5240  7 0.00061



# several covariates do not meet the assumption - we re-run the model and allow them change with time

fit_discover <- coxph(surv ~ tt(mobility) + bold + tt(eigenvector) + tt(age_num) + sex_num + tt(other_present) + colony_num, data = latency.data.discovery,  tt = function(x, t, ...) x * log(t))


summary(fit_discover)


# Call:
#   coxph(formula = surv ~ tt(mobility) + bold + tt(eigenvector) + 
#           tt(age_num) + sex_num + tt(other_present) + colony_num, data = latency.data.discovery, 
#         tt = function(x, t, ...) x * log(t))
# 
# n= 92, number of events= 74 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)   
# tt(mobility)       0.07120   1.07379  0.08732  0.815  0.41486   
# bold               0.40819   1.50409  0.15150  2.694  0.00705 **
#   tt(eigenvector)   -0.01788   0.98228  0.07902 -0.226  0.82102   
# tt(age_num)       -0.09811   0.90655  0.09691 -1.012  0.31136   
# sex_num           -0.53604   0.58506  0.24911 -2.152  0.03142 * 
#   tt(other_present)  0.19965   1.22098  0.12538  1.592  0.11130   
# colony_num         0.97847   2.66037  0.30078  3.253  0.00114 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# tt(mobility)         1.0738     0.9313    0.9049    1.2742
# bold                 1.5041     0.6649    1.1177    2.0241
# tt(eigenvector)      0.9823     1.0180    0.8413    1.1468
# tt(age_num)          0.9066     1.1031    0.7497    1.0962
# sex_num              0.5851     1.7092    0.3591    0.9533
# tt(other_present)    1.2210     0.8190    0.9550    1.5611
# colony_num           2.6604     0.3759    1.4754    4.7970
# 
# Concordance= 0.706  (se = 0.036 )
# Likelihood ratio test= 37.47  on 7 df,   p=4e-06
# Wald test            = 35.5  on 7 df,   p=9e-06
# Score (logrank) test = 38.94  on 7 df,   p=2e-06



# how to interpret the tt. E.g. tt(age_num): age difference matters early but then diminishes over time (-0.09). Juveniles (coded as 1) are slower at first but then converge later. The effect is non-significant though since it spans 0.

# Bolder individuals, females and Crow individuals faster at discovering.
# Chance of discovery increases by 50% (1.50) per unit increase in boldness.
# Males had a lower hazard of discovery than females (HR = 0.58).
# Crow had a higher hazard of discovery compared to Paradise (HR = 2.66)


# Concordance: model's abilty to correctly rank the latencies (0.5=chance)
# p values - overall statistical significance

# extract results as a table:


table_discovery <- tidy(fit_discover, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
                  "tt(mobility)" = "Mobility (time-varying)",
                  "bold" = "Boldness",
                  "tt(eigenvector)" = "Eigenvector centrality (time-varying)",
                  "tt(age_num)" = "Age [J:A] (time-varying)",
                  "sex_num" = "Sex [M:F]",
                  "tt(other_present)" = "Conspecific presence (time-varying)",
                  "colony_num" = "Population [more:less disturbed]")
  ) %>%
 dplyr:: select(
    Predictor = term,
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    z = statistic,
    p = p.value
  )

# round to two digits
table_discovery <- table_discovery %>%
  mutate(
    across(c(HR, CI_low, CI_high), ~ round(., 2)),
    z = round(z, 2),
    p = signif(p, 2)
  )

# save as word

table_discovery_fmt <- table_discovery %>%
  dplyr::mutate(
    HR = round(HR, 2),
    `95% CI (lower)` = round(`CI_low`, 2),
    `95% CI (upper)` = round(`CI_high`, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )%>%
  dplyr::select(-CI_low, -CI_high)

# change order
table_discovery_fmt <- table_discovery %>%
  transmute(
    Predictor,
    HR = round(HR, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(table_discovery_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_flextable(ft)

setwd("..")
print(doc, target = "Tables/table_discovery.docx")



# 4) VIFs- solving ----------------------------------------------------


latency.data.solve$age_num <- ifelse(
  latency.data.solve$age == "P", 1, 0
)

latency.data.solve$sex_num <- ifelse(
  latency.data.solve$sex == "M", 1, 0
)

latency.data.solve$colony_num <- ifelse(
  latency.data.solve$colony == "Crow", 1, 0
)

vif <- lm(latency.data.solve$time ~ mobility + bold + eigenvector + age_num + sex_num + n_observation_opportunities + log(training.solves+1) + colony_num, data = latency.data.solve)

vif(vif)
# mobility                        bold                 eigenvector                     age_num 
# 1.896118                    1.409199                    1.970772                    1.344322 
# sex_num n_observation_opportunities    log(training.solves + 1)                  colony_num 
# 1.053205                    1.333734                    1.722095                    1.197973 

cor(latency.data.solve %>%
      dplyr::select(mobility, eigenvector, bold, n_observation_opportunities, sex_num, age_num, colony_num))

#                              mobility eigenvector        bold n_observation_opportunities     sex_num
# mobility                     1.00000000  0.65786264  0.06005912                  0.39904068  0.06272376
# eigenvector                  0.65786264  1.00000000 -0.01355535                  0.41359721  0.05567223
# bold                         0.06005912 -0.01355535  1.00000000                  0.02976622 -0.16324285
# n_observation_opportunities  0.39904068  0.41359721  0.02976622                  1.00000000  0.06442346
# sex_num                      0.06272376  0.05567223 -0.16324285                  0.06442346  1.00000000
# age_num                      0.06514797  0.17133787  0.07076381                 -0.08352743 -0.09137478
# colony_num                  -0.06915099 -0.02475148  0.21880003                  0.03378395 -0.09004924
#                               age_num  colony_num
# mobility                     0.06514797 -0.06915099
# eigenvector                  0.17133787 -0.02475148
# bold                         0.07076381  0.21880003
# n_observation_opportunities -0.08352743  0.03378395
# sex_num                     -0.09137478 -0.09004924
# age_num                      1.00000000  0.15059735
# colony_num                   0.15059735  1.00000000




# 5) Latency solving -----------------------------------------------------

# create a survival object
latency.data.solve$time_solve <- latency.data.solve$time
latency.data.solve$surv <- with(latency.data.solve, Surv(time_solve, !is.na(timestamp)))


# scale some variables
latency.data.solve[c("mobility", "eigenvector", "bold", "n_observation_opportunities")] <- lapply(latency.data.solve[c("mobility", "eigenvector", "bold", "n_observation_opportunities")],scale)
latency.data.solve$log_training_solves <- scale(log(latency.data.solve$training.solves+1))


# Fit the Cox proportional hazards model
fit_solve <- coxph(surv ~   mobility + bold + eigenvector + age + sex + log_training_solves + n_observation_opportunities + colony_num, data = latency.data.solve)

# check for proportional hazards assumption

cox.zph(fit_solve)

# chisq df     p
# mobility                     3.238  1 0.072
# bold                         2.022  1 0.155
# eigenvector                  1.215  1 0.270
# age                          0.245  1 0.621
# sex                          0.522  1 0.470
# log_training_solves          0.295  1 0.587
# n_observation_opportunities  2.378  1 0.123
# colony_num                   0.493  1 0.483
# GLOBAL                      12.852  8 0.117

# assumption met for all covaraites

summary(fit_solve)

# Call:
#   coxph(formula = surv ~ mobility + bold + eigenvector + age + 
#           sex + log_training_solves + n_observation_opportunities + 
#           colony_num, data = latency.data.solve)
# 
# n= 58, number of events= 27 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)    
# mobility                     0.51265   1.66972  0.27516  1.863 0.062443 .  
# bold                         0.92565   2.52350  0.32355  2.861 0.004224 ** 
#   eigenvector                  0.69403   2.00177  0.30737  2.258 0.023949 *  
#   ageP                         0.05902   1.06079  0.59718  0.099 0.921276    
# sexM                        -0.34814   0.70600  0.46842 -0.743 0.457348    
# log_training_solves          0.28098   1.32443  0.25508  1.102 0.270653    
# n_observation_opportunities  0.59585   1.81457  0.17898  3.329 0.000871 ***
#   colony_num                   0.03877   1.03953  0.50320  0.077 0.938592    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# mobility                        1.670     0.5989    0.9737     2.863
# bold                            2.524     0.3963    1.3384     4.758
# eigenvector                     2.002     0.4996    1.0959     3.656
# ageP                            1.061     0.9427    0.3291     3.419
# sexM                            0.706     1.4164    0.2819     1.768
# log_training_solves             1.324     0.7550    0.8034     2.183
# n_observation_opportunities     1.815     0.5511    1.2777     2.577
# colony_num                      1.040     0.9620    0.3877     2.787
# 
# Concordance= 0.841  (se = 0.036 )
# Likelihood ratio test= 45.88  on 8 df,   p=3e-07
# Wald test            = 31.99  on 8 df,   p=9e-05
# Score (logrank) test = 48.76  on 8 df,   p=7e-08


# extract results as a table:


table_solve <- tidy(fit_solve, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
                  "mobility" = "Mobility",
                  "bold" = "Boldness",
                  "eigenvector" = "Eigenvector centrality",
                  "ageP" = "Age [J:A]",
                  "sexM" = "Sex [M:F]",
                  "log_training_solves" = "Log # training solves",
                  "n_observation_opportunities" = "Observation opportunities (solving)",
                  "colony_num" = "Population [more:less disturbed]")
  ) %>%
  dplyr:: select(
    Predictor = term,
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    z = statistic,
    p = p.value
  )

# round to two digits
table_solve <- table_solve %>%
  mutate(
    across(c(HR, CI_low, CI_high), ~ round(., 2)),
    z = round(z, 2),
    p = signif(p, 2)
  )

# save as word

table_solve_fmt <- table_solve %>%
  dplyr::mutate(
    HR = round(HR, 2),
    `95% CI (lower)` = round(`CI_low`, 2),
    `95% CI (upper)` = round(`CI_high`, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )%>%
  dplyr::select(-CI_low, -CI_high)

# change order
table_solve_fmt <- table_solve %>%
  transmute(
    Predictor,
    HR = round(HR, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(table_solve_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Tables/table_solve.docx")


# 6) Exploitation success ------------------------------------------------------

# subset the data to those who learned to solve
solving <- latency.data.solve[!is.na(latency.data.solve$timestamp),]

# VIFs
vif <- lm(total.solves ~ mobility + bold + eigenvector + age_num + sex_num + colony_num, data = solving)

vif(vif)
# mobility        bold eigenvector     age_num     sex_num  colony_num 
# 2.030468    1.238846    1.920832    1.280135    1.450145    1.328935


# create a log_days column

solving$log_days <- log(solving$total.solve.days)

library(MASS)
# Fit negative binomial GLM to account for overdispersion
fit_success_glm <- glm.nb( total.solves ~ mobility + bold + eigenvector + age + sex + colony_num + offset(log_days), data = solving )

# View summary
summary(fit_success_glm)

# Call:
#   glm.nb(formula = total.solves ~ mobility + bold + eigenvector + 
#            age + sex + colony_num + offset(log_days), data = solving, 
#          init.theta = 3.288238154, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.8659     0.3336   5.593 2.24e-08 ***
#   mobility      0.5812     0.1499   3.878 0.000105 ***
#   bold          0.4267     0.1493   2.858 0.004264 ** 
#   eigenvector   0.1079     0.1554   0.694 0.487496    
# ageP         -0.5973     0.2843  -2.101 0.035654 *  
#   sexM          0.7485     0.2866   2.612 0.009009 ** 
#   colony_num    0.2370     0.2903   0.816 0.414292    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(3.2882) family taken to be 1)
# 
# Null deviance: 81.526  on 26  degrees of freedom
# Residual deviance: 27.755  on 20  degrees of freedom
# AIC: 252.2
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  3.288 
# Std. Err.:  0.997 
# 
# 2 x log-likelihood:  -236.196 




# extract results as table:

glm_table <- tidy(
  fit_success_glm,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  dplyr::mutate(
    Predictor = recode(term,
                       "(Intercept)" = "Intercept",
                       "mobility" = "Mobility",
                       "bold" = "Boldness",
                       "eigenvector" = "Eigenvector centrality",
                       "ageP" = "Age [J:A]",
                       "sexM" = "Sex [M:F]",
                       "colony_num" = "Population [more:less disturbed]"
    )
  ) %>%
  dplyr::select(
    Predictor,
    Estimate = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    z = statistic,
    p = p.value
  )


glm_table_fmt <- glm_table %>%
  transmute(
    Predictor,
    Estimate = round(Estimate, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(glm_table_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
#  body_add_par("Negative binomial GLM predicting total solves", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "Tables/GLM_solving_success.docx")




library(DHARMa)
res <- simulateResiduals(fit_success_glm)
plot(res)

# no obvious overdispersion or model misfit - means we can trust the output


# 9) Figure ---------------------------------------------------------------

# plot boldness, mobility and eigenvector centrality in the two populations

bold.plot <- ggplot(aes(x=colony, y=bold), data=latency.data.discovery)+
  geom_boxplot()+
  geom_jitter(width=0.15, height=0.07, alpha=0.4)+
  theme_bw()+
  ylab("Boldness score [PCA]")+
  xlab("Population")+
  scale_x_discrete(
    labels = c(
      "Crow" = "More disturbed",
      "Paradise" = "Less disturbed"
    )
  )

# expl.plot <- ggplot(aes(x=population, y=mobility), data=plot)+
#   geom_boxplot()+
#   geom_jitter(width=0.15, height=0.07, alpha=0.4)+
#   theme_bw()+
#   ylab("Mobility")+
#   xlab("Population")+
#   scale_x_discrete(
#     labels = c(
#       "Crow" = "More disturbed",
#       "Paradise" = "Less disturbed"
#     )
#   )

soc.plot <- ggplot(aes(x=colony, y=eigenvector), data=latency.data.discovery)+
  geom_boxplot()+
  geom_jitter(width=0.15, height=0.07, alpha=0.4)+
  theme_bw()+
  ylab("Eigenvector centrality")+
  xlab("Population")+
  scale_x_discrete(
    labels = c(
      "Crow" = "More disturbed",
      "Paradise" = "Less disturbed"
    )
  )

library(ggpubr)

ggarrange(bold.plot,  soc.plot, labels = c("a", "b"), ncol=2)


ggsave("Supplementary Figures/Figure_beh_type_colony.tiff", units="in", width=7, height=3, dpi=300, compression = 'lzw')



