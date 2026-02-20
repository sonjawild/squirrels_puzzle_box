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
Crow.latency.data.solve <- read.csv("Crow.latency.data.solve.csv", row.names = 1)
Paradise.latency.data.solve <- read.csv("Paradise.latency.data.solve.csv", row.names = 1)

# latency until discovery
Crow.latency.data.discovery <- read.csv("Crow.latency.data.discovery.csv", row.names = 1)
Paradise.latency.data.discovery <- read.csv("Paradise.latency.data.discovery.csv", row.names = 1)

names(Crow.latency.data.solve)

# Variables of interest:
# id: squirrel identity
# timestamp: date and time of third solve for those with at least three solve (or time of discovery for discovery df)
# time: latency to third solve experimental hours 
# training.solves: how many solves did the squirrel perform during the training phase (with additional peanut butter on the lever). controls for puzzle exposure prior to trial phase
# solve.rate: solves/visit AFTER learning to solve 
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
Crow.latency.data.discovery$age_num <- ifelse(
  Crow.latency.data.discovery$age == "P", 1, 0
)

Crow.latency.data.discovery$sex_num <- ifelse(
  Crow.latency.data.discovery$sex == "M", 1, 0
)

vif.Crow <- lm(rep(1, nrow(Crow.latency.data.discovery)) ~ exploration + bold + eigenvector + age_num + sex_num + other_present, data = Crow.latency.data.discovery)

vif(vif.Crow)

# exploration          bold   eigenvector       age_num       sex_num other_present 
# 2.000947      1.110828      2.084664      1.326707      1.075397      1.180914 

cor(Crow.latency.data.discovery %>%
  dplyr::select(exploration, eigenvector, bold, other_present, age_num, sex_num))

#              exploration eigenvector        bold other_present     age_num     sex_num
# exploration    1.00000000  0.68694410  0.18559710   -0.16605971  0.33598919 -0.07208602
# eigenvector    0.68694410  1.00000000  0.11115082   -0.15261967  0.41722881 -0.06836602
# bold           0.18559710  0.11115082  1.00000000    0.18400789 -0.02656312 -0.08284678
# other_present -0.16605971 -0.15261967  0.18400789    1.00000000  0.09128709 -0.16158054
# age_num        0.33598919  0.41722881 -0.02656312    0.09128709  1.00000000  0.10114435
# sex_num       -0.07208602 -0.06836602 -0.08284678   -0.16158054  0.10114435  1.00000000


Paradise.latency.data.discovery$age_num <- ifelse(
  Paradise.latency.data.discovery$age == "P", 1, 0
)

Paradise.latency.data.discovery$sex_num <- ifelse(
  Paradise.latency.data.discovery$sex == "M", 1, 0
)



vif.Paradise <- lm(rep(1, nrow(Paradise.latency.data.discovery)) ~ exploration + bold + eigenvector + age_num + sex_num + other_present, data = Paradise.latency.data.discovery)

vif(vif.Paradise)

# exploration          bold   eigenvector       age_num       sex_num other_present 
# 3.844647      1.650093      3.244103      1.235529      1.126134      3.053111 

cor(Paradise.latency.data.discovery %>%
      dplyr::select(exploration, eigenvector, bold, age_num, sex_num, other_present))

#                  exploration eigenvector        bold     age_num     sex_num other_present
# exploration     1.0000000   0.8023467  0.22096012 -0.21193649  0.14953535     0.6883407
# eigenvector     0.8023467   1.0000000  0.15602390 -0.26577424  0.23668082     0.6659147
# bold            0.2209601   0.1560239  1.00000000  0.23738565 -0.01523406    -0.2050971
# age_num        -0.2119365  -0.2657742  0.23738565  1.00000000 -0.08883363    -0.1666667
# sex_num         0.1495354   0.2366808 -0.01523406 -0.08883363  1.00000000     0.2842676
# other_present   0.6883407   0.6659147 -0.20509706 -0.16666667  0.28426762     1.0000000


# exploration is relatively high but below the used criterior of 5


# 3) Latency - discovery --------------------------------------------------


# 3.1) Crow -----------------------------------------------------

# create a survival object
# those that did not discover the resource are right censored (+) sign
Crow.latency.data.discovery$time_discover <- Crow.latency.data.discovery$time
Crow.latency.data.discovery$surv <- with(Crow.latency.data.discovery, Surv(time_discover, !is.na(timestamp)))

# scale some variables
Crow.latency.data.discovery[c("exploration", "eigenvector", "bold")] <- lapply(Crow.latency.data.discovery[c("exploration", "eigenvector", "bold")],scale)

# Fit the Cox proportional hazards model
fit_Crow_discover <- coxph(surv ~ exploration + bold + eigenvector + age_num + sex_num + other_present, data = Crow.latency.data.discovery)

# check for proportional hazards assumption
# we need to make sure that the effects of a predictor does not change over time
# significant p values indicate a violation of the proportional hazards assumption

cox.zph(fit_Crow_discover)

# chisq df       p
# exploration    7.371  1 0.00663
# bold           0.940  1 0.33215
# eigenvector    3.329  1 0.06808
# age_num        6.952  1 0.00837
# sex_num        0.285  1 0.59331
# other_present  7.680  1 0.00559
# GLOBAL        22.732  6 0.00089

# several covariates do not meet the assumption - we re-run the model and allow them change with time


fit_Crow_discover <- coxph(surv ~ tt(exploration) + bold + eigenvector + tt(age_num) + sex_num + tt(other_present), data = Crow.latency.data.discovery,  tt = function(x, t, ...) x * log(t))


summary(fit_Crow_discover)

# Call:
#   coxph(formula = surv ~ tt(exploration) + bold + eigenvector + 
#           tt(age_num) + sex_num + tt(other_present), data = Crow.latency.data.discovery, 
#         tt = function(x, t, ...) x * log(t))
# 
# n= 57, number of events= 50 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# tt(exploration)    0.13257   1.14176  0.09321  1.422   0.1549  
# bold               0.43912   1.55134  0.18449  2.380   0.0173 *
#   eigenvector       -0.25154   0.77760  0.18433 -1.365   0.1724  
# tt(age_num)       -0.13261   0.87580  0.14046 -0.944   0.3451  
# sex_num           -0.72563   0.48402  0.32477 -2.234   0.0255 *
#   tt(other_present)  0.14053   1.15088  0.14589  0.963   0.3354  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
#                     exp(coef) exp(-coef) lower .95 upper .95
# tt(exploration)      1.1418     0.8758    0.9511    1.3706
# bold                 1.5513     0.6446    1.0806    2.2271
# eigenvector          0.7776     1.2860    0.5418    1.1160
# tt(age_num)          0.8758     1.1418    0.6650    1.1534
# sex_num              0.4840     2.0660    0.2561    0.9148
# tt(other_present)    1.1509     0.8689    0.8647    1.5318
# 
# Concordance= 0.671  (se = 0.043 )
# Likelihood ratio test= 18.47  on 6 df,   p=0.005
# Wald test            = 17.08  on 6 df,   p=0.009
# Score (logrank) test = 18.07  on 6 df,   p=0.006

# how to interpret the tt. E.g. tt(age_num): age difference matters early but then diminishes over time (-0.13). Juveniles (coded as 1) are slower at first but then converge later. The effect is non-significant though since it spans 0.

# Bolder individuals and females faster at discovering.
# Chance of discovery increases by 55% (1.55) per unit increase in boldness.
# Males had a lower hazard of discovery than females (HR = 0.48).


# Concordance: model's abilty to correctly rank the latencies (0.5=chance)
# p values - overall statistical significance

# extract results as a table:


Crow_table_discovery <- tidy(fit_Crow_discover, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
                  "tt(exploration)" = "Mobility (time-varying)",
                  "bold" = "Boldness",
                  "eigenvector" = "Eigenvector centrality",
                  "tt(age_num)" = "Age [J:A] (time-varying)",
                  "sex_num" = "Sex [M:F]",
                  "tt(other_present)" = "Conspecific presence (time-varying)")
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
Crow_table_discovery <- Crow_table_discovery %>%
  mutate(
    across(c(HR, CI_low, CI_high), ~ round(., 2)),
    z = round(z, 2),
    p = signif(p, 2)
  )

# save as word

Crow_table_discovery_fmt <- Crow_table_discovery %>%
  dplyr::mutate(
    HR = round(HR, 2),
    `95% CI (lower)` = round(`CI_low`, 2),
    `95% CI (upper)` = round(`CI_high`, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )%>%
  dplyr::select(-CI_low, -CI_high)

# change order
Crow_table_discovery_fmt <- Crow_table_discovery %>%
  transmute(
    Predictor,
    HR = round(HR, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(Crow_table_discovery_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_flextable(ft)

setwd("..")
print(doc, target = "Tables/More disturbed_table_discovery.docx")


# 3.2) Paradise -----------------------------------------------------

# create a survival object
# those that did not discover the resource are right censored (+) sign
Paradise.latency.data.discovery$time_discover <- Paradise.latency.data.discovery$time
Paradise.latency.data.discovery$surv <- with(Paradise.latency.data.discovery, Surv(time_discover, !is.na(timestamp)))


# scale some variables
Paradise.latency.data.discovery[c("exploration", "eigenvector", "bold")] <- lapply(Paradise.latency.data.discovery[c("exploration", "eigenvector", "bold")],scale)

# Fit the Cox proportional hazards model
fit_Paradise_discover <- coxph(surv ~ exploration + bold + eigenvector + age_num + sex_num + other_present, data = Paradise.latency.data.discovery)

# check for proportional hazards assumption

cox.zph(fit_Paradise_discover)

# chisq df     p
# exploration    1.2622  1 0.261
# bold           0.1945  1 0.659
# eigenvector    3.6540  1 0.056
# age_num        2.8146  1 0.093
# sex_num        0.0497  1 0.824
# other_present  2.2512  1 0.134
# GLOBAL        10.2778  6 0.113

# no significant p-values, assumption is met

summary(fit_Paradise_discover)

# Call:
#   coxph(formula = surv ~ exploration + bold + eigenvector + age_num + 
#           sex_num + other_present, data = Paradise.latency.data.discovery)
# 
# n= 35, number of events= 24 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)
# exploration    0.1288    1.1375   0.5558  0.232    0.817
# bold           0.4544    1.5753   0.2937  1.548    0.122
# eigenvector   -0.2819    0.7543   0.4256 -0.662    0.508
# age_num       -0.5433    0.5808   0.4885 -1.112    0.266
# sex_num       -0.3700    0.6907   0.4534 -0.816    0.414
# other_present  1.6392    5.1509   1.6820  0.975    0.330
# 
# exp(coef) exp(-coef) lower .95 upper .95
# exploration      1.1375     0.8791    0.3827     3.381
# bold             1.5753     0.6348    0.8859     2.801
# eigenvector      0.7543     1.3257    0.3275     1.737
# age_num          0.5808     1.7217    0.2229     1.513
# sex_num          0.6907     1.4477    0.2841     1.680
# other_present    5.1509     0.1941    0.1906   139.201
# 
# Concordance= 0.599  (se = 0.069 )
# Likelihood ratio test= 4.36  on 6 df,   p=0.6
# Wald test            = 4.41  on 6 df,   p=0.6
# Score (logrank) test = 4.6  on 6 df,   p=0.6

# extract results as a table:

Paradise_table_discovery <- tidy(fit_Paradise_discover, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
                  "exploration" = "Mobility",
                  "bold" = "Boldness",
                  "eigenvector" = "Eigenvector centrality",
                  "age_num" = "Age [J:A]",
                  "sex_num" = "Sex [M:F]",
                  "other_present" = "Conspecific presence")
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
Paradise_table_discovery <- Paradise_table_discovery %>%
  mutate(
    across(c(HR, CI_low, CI_high), ~ round(., 2)),
    z = round(z, 2),
    p = signif(p, 2)
  )

# save as word

Paradise_table_discovery_fmt <- Paradise_table_discovery %>%
  dplyr::mutate(
    HR = round(HR, 2),
    `95% CI (lower)` = round(`CI_low`, 2),
    `95% CI (upper)` = round(`CI_high`, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )%>%
  dplyr::select(-CI_low, -CI_high)

# change order
Paradise_table_discovery_fmt <- Paradise_table_discovery %>%
  transmute(
    Predictor,
    HR = round(HR, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(Paradise_table_discovery_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Tables/Less disturbed_table_discovery.docx")


# we find no evidence for an effect of any of the covariates

# Plot survival curves
Paradise_surv_fit_discover <- survfit(fit_Paradise_discover)
plot(Paradise_surv_fit_discover, xlab = "Time (hours)", ylab = "Probability discovery [Paradise]")


# 4) VIFs- solving ----------------------------------------------------

# we calculate the VIFs for the data set on solving latency
# VIFS

Crow.latency.data.solve$age_num <- ifelse(
  Crow.latency.data.solve$age == "P", 1, 0
)

Crow.latency.data.solve$sex_num <- ifelse(
  Crow.latency.data.solve$sex == "M", 1, 0
)

vif.Crow <- lm(rep(1, nrow(Crow.latency.data.solve)) ~ exploration + bold + eigenvector + age_num + sex_num + n_observation_opportunities + log(training.solves+1), data = Crow.latency.data.solve)

vif(vif.Crow)
# exploration                        bold                 eigenvector                     age_num 
# 1.676171                    1.515596                    1.798417                    1.516741 
# sex_num n_observation_opportunities    log(training.solves + 1) 
# 1.047316                    1.114400                    1.967497 

cor(Crow.latency.data.solve %>%
      dplyr::select(exploration, eigenvector, bold, n_observation_opportunities, sex_num, age_num))

# exploration eigenvector        bold n_observation_opportunities     sex_num     age_num
# exploration                  1.00000000  0.60814182  0.15059993                  0.17726939 -0.04838431  0.27294573
# eigenvector                  0.60814182  1.00000000  0.04451408                  0.19242244 -0.07742678  0.36130861
# bold                         0.15059993  0.04451408  1.00000000                  0.06739589 -0.05786852 -0.02730807
# n_observation_opportunities  0.17726939  0.19242244  0.06739589                  1.00000000 -0.12690456 -0.10621433
# sex_num                     -0.04838431 -0.07742678 -0.05786852                 -0.12690456  1.00000000 -0.02706660
# age_num                      0.27294573  0.36130861 -0.02730807                 -0.10621433 -0.02706660  1.00000000


Paradise.latency.data.solve$age_num <- ifelse(
  Paradise.latency.data.solve$age == "P", 1, 0
)

Paradise.latency.data.solve$sex_num <- ifelse(
  Paradise.latency.data.solve$sex == "M", 1, 0
)

vif.Paradise <- lm(rep(1, nrow(Paradise.latency.data.solve)) ~ exploration + bold + eigenvector + age_num + sex_num + n_observation_opportunities + log(training.solves+1), data = Paradise.latency.data.solve)


vif(vif.Paradise)
# exploration                        bold                 eigenvector                     age_num 
# 2.460066                    1.368221                    4.131665                    1.561001 
# sex_num n_observation_opportunities    log(training.solves + 1) 
# 1.403829                    3.568174                    1.608025


cor(Paradise.latency.data.solve %>%
      dplyr::select(exploration, eigenvector, bold, n_observation_opportunities, sex_num, age_num))


# exploration eigenvector          bold n_observation_opportunities    sex_num     age_num
# exploration                  1.0000000000  0.75702575  0.0002108029                  0.64667051  0.1947589 -0.21690996
# eigenvector                  0.7570257516  1.00000000 -0.0974497154                  0.78411721  0.3159323 -0.26682683
# bold                         0.0002108029 -0.09744972  1.0000000000                 -0.04003202 -0.3039990  0.18010418
# n_observation_opportunities  0.6466705105  0.78411721 -0.0400320178                  1.00000000  0.3960989 -0.06055484
# sex_num                      0.1947589137  0.31593233 -0.3039990070                  0.39609892  1.0000000 -0.20916501
# age_num                     -0.2169099648 -0.26682683  0.1801041781                 -0.06055484 -0.2091650  1.00000000


# 5) Latency solving -----------------------------------------------------


# 5.1) Crow ---------------------------------------------------------------


# create a survival object
Crow.latency.data.solve$time_solve <- Crow.latency.data.solve$time
Crow.latency.data.solve$surv <- with(Crow.latency.data.solve, Surv(time_solve, !is.na(timestamp)))


# scale some variables
Crow.latency.data.solve[c("exploration", "eigenvector", "bold", "n_observation_opportunities")] <- lapply(Crow.latency.data.solve[c("exploration", "eigenvector", "bold", "n_observation_opportunities")],scale)
Crow.latency.data.solve$log_training_solves <- scale(log(Crow.latency.data.solve$training.solves+1))


# Fit the Cox proportional hazards model
fit_Crow_solve <- coxph(surv ~   exploration + bold + eigenvector + age + sex + log_training_solves + n_observation_opportunities, data = Crow.latency.data.solve)

# check for proportional hazards assumption

cox.zph(fit_Crow_solve)

# chisq df    p
# exploration                 2.453628  1 0.12
# bold                        0.566501  1 0.45
# eigenvector                 0.070559  1 0.79
# age                         0.000981  1 0.98
# sex                         0.154638  1 0.69
# log_training_solves         0.016496  1 0.90
# n_observation_opportunities 1.203429  1 0.27
# GLOBAL                      8.364942  7 0.30

# assumption met for all covaraites

summary(fit_Crow_solve)

# Call:
#   coxph(formula = surv ~ exploration + bold + eigenvector + age + 
#           sex + log_training_solves + n_observation_opportunities, 
#         data = Crow.latency.data.solve)
# 
# n= 40, number of events= 19 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# exploration                  0.01063   1.01069  0.34409  0.031   0.9754  
# bold                         0.89766   2.45385  0.38953  2.304   0.0212 *
#   eigenvector                  0.87599   2.40126  0.40851  2.144   0.0320 *
#   ageP                        -0.29689   0.74312  0.76724 -0.387   0.6988  
# sexM                        -0.74598   0.47427  0.63083 -1.183   0.2370  
# log_training_solves          0.17338   1.18932  0.34180  0.507   0.6120  
# n_observation_opportunities  0.42722   1.53299  0.19895  2.147   0.0318 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# exploration                    1.0107     0.9894    0.5149     1.984
# bold                           2.4539     0.4075    1.1436     5.265
# eigenvector                    2.4013     0.4164    1.0782     5.348
# ageP                           0.7431     1.3457    0.1652     3.343
# sexM                           0.4743     2.1085    0.1377     1.633
# log_training_solves            1.1893     0.8408    0.6086     2.324
# n_observation_opportunities    1.5330     0.6523    1.0380     2.264
# 
# Concordance= 0.8  (se = 0.05 )
# Likelihood ratio test= 26.27  on 7 df,   p=5e-04
# Wald test            = 20.75  on 7 df,   p=0.004
# Score (logrank) test = 27.62  on 7 df,   p=3e-04


# extract results as a table:


Crow_table_solve <- tidy(fit_Crow_solve, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
                  "exploration" = "Mobility",
                  "bold" = "Boldness",
                  "eigenvector" = "Eigenvector centrality",
                  "ageP" = "Age [J:A]",
                  "sexM" = "Sex [M:F]",
                  "log_training_solves" = "Log # training solves",
                  "n_observation_opportunities" = "Observation opportunities (solving)")
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
Crow_table_solve <- Crow_table_solve %>%
  mutate(
    across(c(HR, CI_low, CI_high), ~ round(., 2)),
    z = round(z, 2),
    p = signif(p, 2)
  )

# save as word

Crow_table_solve_fmt <- Crow_table_solve %>%
  dplyr::mutate(
    HR = round(HR, 2),
    `95% CI (lower)` = round(`CI_low`, 2),
    `95% CI (upper)` = round(`CI_high`, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )%>%
  dplyr::select(-CI_low, -CI_high)

# change order
Crow_table_solve_fmt <- Crow_table_solve %>%
  transmute(
    Predictor,
    HR = round(HR, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(Crow_table_solve_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Tables/More disturbed_table_solve.docx")


# Plot survival curves
Crow_surv_fit_solve <- survfit(fit_Crow_solve)
plot(Crow_surv_fit_solve, xlab = "Time (hours)", ylab = "Probability solve [Crow]")

# bolder individuals faster, more sociable individuals are faster
# faster with more observation opportunities, suggestive of a potential effect of social learning

# 5.2) Paradise -----------------------------------------------------

# create a survival object
Paradise.latency.data.solve$time_solve <- Paradise.latency.data.solve$time
Paradise.latency.data.solve$surv <- with(Paradise.latency.data.solve, Surv(time_solve, !is.na(timestamp)))


# scale some variables
Paradise.latency.data.solve[c("exploration", "eigenvector", "bold", "n_observation_opportunities")] <- lapply(Paradise.latency.data.solve[c("exploration", "eigenvector", "bold", "n_observation_opportunities")],scale)
Paradise.latency.data.solve$log_training_solves <- scale(log(Paradise.latency.data.solve$training.solves+1))


# Fit the Cox proportional hazards model
fit_Paradise_solve <- coxph(surv ~   exploration + bold + eigenvector + age + sex + log_training_solves + n_observation_opportunities, data = Paradise.latency.data.solve)

# Warning message:
#   In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#                  Ran out of iterations and did not converge

# too many covariates for the data - we reduce by excluding eigenvector since it is highly correlated with observation opportunities
# also remove age, sex and training solves
fit_Paradise_solve <- coxph(surv ~   exploration + bold   +  n_observation_opportunities, data = Paradise.latency.data.solve)

# check for proportional hazards assumption

cox.zph(fit_Paradise_solve)

# chisq df    p
# exploration                 0.669  1 0.41
# bold                        0.566  1 0.45
# n_observation_opportunities 2.393  1 0.12
# GLOBAL                      2.552  3 0.47

# assumption met for all covaraites

summary(fit_Paradise_solve)

# Call:
#   coxph(formula = surv ~ exploration + bold + n_observation_opportunities, 
#         data = Paradise.latency.data.solve)
# 
# n= 18, number of events= 8 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# exploration                  3.202    24.578    1.322 2.421   0.0155 *
#   bold                         3.100    22.187    1.372 2.259   0.0239 *
#   n_observation_opportunities  2.839    17.092    1.191 2.384   0.0171 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# exploration                     24.58    0.04069     1.840     328.3
# bold                            22.19    0.04507     1.507     326.6
# n_observation_opportunities     17.09    0.05851     1.656     176.4
# 
# Concordance= 0.972  (se = 0.027 )
# Likelihood ratio test= 28.34  on 3 df,   p=3e-06
# Wald test            = 6.69  on 3 df,   p=0.08
# Score (logrank) test = 25.42  on 3 df,   p=1e-05


# extract results as a table:


Paradise_table_solve <- tidy(fit_Paradise_solve, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
                  "exploration" = "Mobility",
                  "bold" = "Boldness",
               #   "eigenvector" = "Eigenvector centrality",
              #    "ageP" = "Age [J:A]",
              #    "sexM" = "Sex [M:F]",
              #    "log_training_solves" = "Log # training solves",
                  "n_observation_opportunities" = "Observation opportunities (solving)")
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
Paradise_table_solve <- Paradise_table_solve %>%
  mutate(
    across(c(HR, CI_low, CI_high), ~ round(., 2)),
    z = round(z, 2),
    p = signif(p, 2)
  )

# save as word

Paradise_table_solve_fmt <- Paradise_table_solve %>%
  dplyr::mutate(
    HR = round(HR, 2),
    `95% CI (lower)` = round(`CI_low`, 2),
    `95% CI (upper)` = round(`CI_high`, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )%>%
  dplyr::select(-CI_low, -CI_high)

# change order
Paradise_table_solve_fmt <- Paradise_table_solve %>%
  transmute(
    Predictor,
    HR = round(HR, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(Paradise_table_solve_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Tables/Less disturbed_table_solve.docx")



# Plot survival curves
Paradise_surv_fit_solve <- survfit(fit_Paradise_solve)
plot(Paradise_surv_fit_solve, xlab = "Time (hours)", ylab = "Probability solve [Paradise]")


# bolder individuals faster, more exploratory individuals are faster
# faster with more observation opportunities, suggestive of a potential effect of social learning


# 6) Solving success ------------------------------------------------------


# 6.1) Crow -----------------------------------------------------
# subset the data to those who learned to solve

Crow.solving <- Crow.latency.data.solve[!is.na(Crow.latency.data.solve$timestamp),]

# VIFs
vif.Crow <- lm(rep(1, nrow(Crow.solving)) ~ exploration + bold + eigenvector + age_num + sex_num, data = Crow.solving)

vif(vif.Crow)

# exploration        bold eigenvector     age_num     sex_num 
# 2.532429    1.079777    2.444341    1.619359    1.299356 


# create an offset variables (total arrivals) - these are the arrivals AFTER the squirrel learned to solve (disregards the ones before)
# so we get the solving rate after learning (performance)

Crow.solving$log_arrivals <- log(Crow.solving$total.arrivals)


library(MASS)
# Fit negative binomial GLM to account for overdispersion
fit_Crow_success_glm <- glm.nb(
  total.solves ~ exploration + bold + eigenvector + age + sex + offset(log_arrivals),
  data = Crow.solving
)

# View summary
summary(fit_Crow_success_glm)

# Call:
#   glm.nb(formula = total.solves ~ exploration + bold + eigenvector + 
#            age + sex + offset(log_arrivals), data = Crow.solving, init.theta = 2.865965988, 
#          link = log)
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)   
# (Intercept) -0.89650    0.27605  -3.248  0.00116 **
#   exploration  0.26795    0.25978   1.031  0.30232   
# bold         0.67438    0.21704   3.107  0.00189 **
#   eigenvector  0.05973    0.22766   0.262  0.79305   
# ageP        -0.77905    0.39732  -1.961  0.04991 * 
#   sexM         0.22709    0.37953   0.598  0.54960   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(2.866) family taken to be 1)
# 
# Null deviance: 43.077  on 18  degrees of freedom
# Residual deviance: 18.889  on 13  degrees of freedom
# AIC: 171.37
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  2.87 
# Std. Err.:  1.00 
# 
# 2 x log-likelihood:  -157.366 

# extract results as table:

glm_table_Crow <- tidy(
  fit_Crow_success_glm,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  dplyr::mutate(
    Predictor = recode(term,
                       "(Intercept)" = "Intercept",
                       "exploration" = "Mobility",
                       "bold" = "Boldness",
                       "eigenvector" = "Eigenvector centrality",
                       "ageP" = "Age [J:A]",
                       "sexM" = "Sex [M:F]"
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


glm_table_Crow_fmt <- glm_table_Crow %>%
  transmute(
    Predictor,
    Estimate = round(Estimate, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(glm_table_Crow_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
#  body_add_par("Negative binomial GLM predicting total solves", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "Tables/More disturbed_GLM_solving_success.docx")




library(DHARMa)
res <- simulateResiduals(fit_Crow_success_glm)
plot(res)

# no obvious overdispersion or model misfit - means we can trust the output


# 6.2) Paradise -----------------------------------------------------
# subset the data to those who learned to solve

Paradise.solving <- Paradise.latency.data.solve[!is.na(Paradise.latency.data.solve$timestamp),]

# VIFs
vif.Paradise <- lm(rep(1, nrow(Paradise.solving)) ~ exploration + bold + eigenvector, data = Paradise.solving)

# exploration        bold eigenvector 
# 1.472094    3.510242    3.663553 

vif(vif.Paradise)


# create an offset variables (total arrivals) - these are the arrivals AFTER the squirrel learned to solve (disregards the ones before)
# so we get the solving rate after learning (performance)
Paradise.solving$log_arrivals <- log(Paradise.solving$total.arrivals)

# Fit negative binomial GLM
# we again have to fit a skinny model since we have very limited data
fit_Paradise_success_glm <- glm.nb(
  total.solves ~ exploration + bold + eigenvector + offset(log_arrivals),
  data = Paradise.solving
)

# View summary
summary(fit_Paradise_success_glm)

# 
# Call:
#   glm.nb(formula = total.solves ~ exploration + bold + eigenvector + 
#            offset(log_arrivals), data = Paradise.solving, init.theta = 4.335512139, 
#          link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -1.3506     0.4842  -2.789  0.00528 **
#   exploration   0.5215     0.2015   2.589  0.00964 **
#   bold          0.1980     0.3611   0.548  0.58346   
# eigenvector   0.1425     0.4446   0.321  0.74859   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(4.3355) family taken to be 1)
# 
# Null deviance: 16.0068  on 7  degrees of freedom
# Residual deviance:  8.1271  on 4  degrees of freedom
# AIC: 88.284
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  4.34 
# Std. Err.:  2.34 
# 
# 2 x log-likelihood:  -78.284 


# extract results as table:

glm_table_Paradise <- tidy(
  fit_Paradise_success_glm,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  dplyr::mutate(
    Predictor = recode(term,
                       "(Intercept)" = "Intercept",
                       "exploration" = "Mobility",
                       "bold" = "Boldness",
                       "eigenvector" = "Eigenvector centrality"
                #       "ageP" = "Age [J:A]",
                #       "sexM" = "Sex [M:F]"
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


glm_table_Paradise_fmt <- glm_table_Paradise %>%
  transmute(
    Predictor,
    Estimate = round(Estimate, 2),
    `95% CI (lower)` = round(CI_low, 2),
    `95% CI (upper)` = round(CI_high, 2),
    z = round(z, 2),
    p = ifelse(p < 0.001, "<0.001", signif(p, 2))
  )

ft <- flextable(glm_table_Paradise_fmt) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  #  body_add_par("Negative binomial GLM predicting total solves", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "Tables/Less disturbed_GLM_solving_success.docx")



res <- simulateResiduals(fit_Paradise_success_glm)
plot(res)


# Does not seem to converge


# 7) Chi-square comparisons -----------------------------------------------


# 7.1) Discovery ----------------------------------------------------------

data.d <- matrix(c(50, 7,   # Crow: 50 discovered, 7 did not (57 - 50)
                 24, 11), # Paradise: 24 discovered, 11 did not (35 - 24)
               nrow = 2,
               byrow = TRUE)

rownames(data.d) <- c("Crow", "Paradise")
colnames(data.d) <- c("Discovered", "Not_discovered")


# Chi-square test
chisq.test(data.d)


# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  data.d
# X-squared = 3.9086, df = 1, p-value = 0.04804


# 7.2) Solving ------------------------------------------------------------



data.s <- matrix(c(19, 21,   # Crow: 19 learned, 21 did not 
                   8, 10), # Paradise: 8 learned, 10 did not 
                 nrow = 2,
                 byrow = TRUE)

rownames(data.s) <- c("Crow", "Paradise")
colnames(data.s) <- c("Learned", "Learned")


# Chi-square test
chisq.test(data.s)

# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  data.s
# X-squared = 3.1705e-31, df = 1, p-value = 1

# 9) Figure ---------------------------------------------------------------

# plot boldness, exploration and eigenvector centrality in the two populations

# we need to reload the raw data since we have scaled it within populations

Crow.latency.data.discovery.plot <- read.csv("Data/Crow.latency.data.discovery.csv", row.names = 1)
Paradise.latency.data.discovery.plot <- read.csv("Data/Paradise.latency.data.discovery.csv", row.names = 1)

Crow.latency.data.discovery.plot$population <- "Crow"
Paradise.latency.data.discovery.plot$population <- "Paradise"

plot.comb <- rbind(Crow.latency.data.discovery.plot, Paradise.latency.data.discovery.plot)

bold.plot <- ggplot(aes(x=population, y=bold), data=plot.comb)+
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

# expl.plot <- ggplot(aes(x=population, y=exploration), data=plot.comb)+
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

soc.plot <- ggplot(aes(x=population, y=eigenvector), data=plot.comb)+
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



