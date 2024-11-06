library(tidyverse)
library(ppmSuite)
library(dplyr)
library(salso)
library(coda)
library(plotly)
rm(list = ls()) # Clean the work space
set.seed(4365)

# We have decided hitting percentage, blocking percentage, dig percentage as our 
# covariates in the model

# Read in the conference data (3 conferences, ONLY conference games)
rowstats <- read.csv("womenvb_conf3_data_pass_per.csv")

end_set_stats <- rowstats %>% group_by(matchid, game, teamid) %>% 
  filter(row_number()==n()) %>% 
  mutate(Home.HP = (h_attack_kills - h_attack_errors)/tot_h_attacks) %>% 
  mutate(home_win = case_when(
    wonlost == 1 & hteam == teamid ~ 1, 
    wonlost == 0 & hteam == teamid ~ 0, 
    wonlost == 1 & hteam != teamid ~ 0, 
    wonlost == 0 & hteam != teamid ~ 1
  )) %>% dplyr::ungroup()

# Filter for home team rows only, so we don't get repeated data
end_set_stats <- end_set_stats %>% 
  filter(team == 1)

# Standardize:
std_Home.HP <- (end_set_stats$Home.HP - mean(end_set_stats$Home.HP))/sd(end_set_stats$Home.HP)
std_hblock_per <- (end_set_stats$hblock_per - mean(end_set_stats$hblock_per))/sd(end_set_stats$hblock_per)
std_hd_per <- (end_set_stats$hd_per - mean(end_set_stats$hd_per))/sd(end_set_stats$hd_per)
end_set_stats <- cbind(end_set_stats, std_Home.HP, std_hblock_per, std_hd_per)

end_set_stats$hd_per[end_set_stats$hd_per == "Inf"] <- 0
end_set_stats$hd_per[end_set_stats$hd_per == "NaN"] <- 0
end_set_stats$Home.HP[end_set_stats$Home.HP == "Inf"] <- 0
end_set_stats$Home.HP[end_set_stats$Home.HP == "NaN"] <- 0
end_set_stats$hblock_per[end_set_stats$hblock_per == "Inf"] <- 0
end_set_stats$hblock_per[end_set_stats$hblock_per == "NaN"] <- 0
end_set_stats$Home.HP =  ifelse(is.na(end_set_stats$Home.HP), 0, end_set_stats$Home.HP)
end_set_stats$hd_per = ifelse(is.na(end_set_stats$hd_per), 0, end_set_stats$hd_per)
end_set_stats$hblock_per =  ifelse(is.na(end_set_stats$hblock_per), 0, end_set_stats$hblock_per)

# MODEL SETUP
co <- c(-100000, 0, 100000)
v <- .5
simParms <- c(0.0, 1.0, v, 1.0, 2.0, 0.1)
draws <- 80000 
burn <- 30000 # 30000
thin <- 50 # 50
nout <- (draws - burn)/thin # should be 1000 min

# get the end of game summaries
# end_set_stats <- end_set_stats %>% mutate(unique_set = paste0(matchid, game))
# game_ids <- end_set_stats$unique_set
# 
# # seed 300, 400, 500, 600
# set.seed(12) # change this every time you want to do a different randomized data split
# samp_games <- sample(game_ids, size = 400, replace = FALSE)
# 
# train <- end_set_stats %>% filter(game_ids %in% samp_games)
# test <- end_set_stats %>% filter(!game_ids %in% samp_games)

end_set_stats <- end_set_stats %>% mutate(train = case_when(
  date < "2018-10-14" ~ 1 ,
  TRUE ~ 0))
train = end_set_stats %>% filter(train == 1)
test = end_set_stats %>% filter(train == 0)
xtrain = train %>% dplyr::select(std_Home.HP, std_hd_per, std_hblock_per)
xtest = test %>% dplyr::select(std_Home.HP, std_hd_per, std_hblock_per)
ytrain = train$home_win
ytest = test$home_win

# YOU DO NOT NEED TO RUN MODEL, JUST READ IN RDS FILE BELOW IT
fit_1 <- ordinal_ppmx(y = ytrain,
                      co = co,
                      X = xtrain,
                      M = 1,
                      Xpred = xtest,
                      similarity_function = 1,
                      consim = 1,
                      calibrate = 0,
                      simParms = simParms,
                      modelPriors=c(0, 20, 1, 1),
                      draws = draws,
                      burn = burn,
                      thin = thin,
                      verbose = TRUE)
# Save fit
# saveRDS(fit_1, "final_cofit3.rds")
# saveRDS(fit_1, "2final_cofit3.rds")


# fit_1 <- readRDS("2final_cofit3.rds")
# Check out of sample prediction (.7907801, so 80% accurate out of sample)
preds_ordinal_ppmx = apply(fit_1$ord.ppred, 2, median)
table(preds_ordinal_ppmx, ytest)
ellie_pred <- ifelse(preds_ordinal_ppmx > 0.50, 1, 0)
table(ellie_pred, ytest)
sum(diag(table(ellie_pred, ytest)))/length(ellie_pred)

# clusters range from 8-22, avg: 12.84312
max(fit_1$nclus)
min(fit_1$nclus)
mean(fit_1$nclus)
plot(as.mcmc(fit_1$nclus))

source("NEWnimblePredictionFunction.R")

# THIS IS WHERE WE ARE WORKING
# Trying to find ideal number of clusters
# clusters = salso(fit_1$Si, maxNClusters = 5)
# table <- table(clusters)


# look at means within clusters

SS <- NA
MAD <- NA
# use loss = "binder"
for (j in 1:10){
  clusters = salso(fit_1$Si, loss = "binder", maxNClusters = j)
  nclusters = length(unique(clusters))
  i = 1
  v_var = 1
  ytrain = apply(fit_1$ord.fitted.values, 2, median)
  game_preds <- rep(NA, nrow(test))

  for(i in 1:nrow(test)){
    game_preds[i] = CppmxPrediction(data = as.matrix(xtrain),
                                    new = c(test$std_Home.HP[i],
                                            test$std_hd_per[i],
                                            test$std_hblock_per[i]),
                                    clusters = as.numeric(clusters),
                                    numClusters = nclusters,
                                    nobs_in_clus = table(clusters), 
                                    v = v_var, v0 = 1, mu = 0, m0 = 0,
                                    y = ytrain)
  }
  SS[j] <- sum((game_preds - test$home_win)^2)
  MAD[j] <- median(abs(game_preds - test$home_win))
}


write.csv(SS,'SS_final.csv')
write.csv(MAD, 'MAD_final.csv')
MAD.5 <- read.csv('MAD_final.csv')
SS.5 <- read.csv('SS_final.csv')
SS
min(SS) 
range(SS)
plot(SS)

MAD
min(MAD)
plot(MAD)
range(MAD)


SS.5
plot(SS.5)
min(SS.5$x) # 6 clusters
range(SS.5)

MAD.5
plot(MAD.5)
range(MAD.5)
min(MAD.5$x) # 2 clusters

# plot(xtrain$std_Home.HP, ytrain, col = salso(fit_1$Si, loss = "binder"))
# 
# 
# clusters
# clusters = salso(fit_1$Si, loss = "binder", maxNClusters = 5)
# table(clusters)





# Run 10 times with different splits

# Run 1 (seed(1)) size = 375
# SS min <- 2
# MAD min <- 2

# Run 2 (seed(2)) size = 250
# SS min <- 2
# MAD min <- 2

# Run 3 (seed(3)), size = 500
# SS min <- 8
# MAD min <- 15

# Run 4 (seed(4)), size = 450
# SS min <- 2
# MAD min <- 2

# Run 5 (seed(5)), size = 400
# SS min <- 14
# MAD min <- 2

# Run 6 (seed(6)), size = 350 
# SS min <- 2
# MAD min <- 2

# Run 7 (seed(7)), size = 475
# SS min <- 2
# MAD min <- 2

# Run 8 (seed(8)), size = 450
# SS min <- 8
# MAD min <- 2

# Run 9 (seed(9)), size = 450
# SS min <- 2
# MAD min <- 2

# Run 10 (seed(10)), size = 450
# SS min <- 2
# MAD min <- 2

# Run 11 (seed(11)), size = how we normally split
# SS min <- 4
# MAD min <- 2

# Run 12 (seed(12)), size = 
# SS min <- 2
# MAD min <- 4


# Mean SS min
# SS <- c(2,2,8,2,14,2,2,8,2,2,4,2)
# mean(SS)
# # Mean MAD min
# MAD <- c(2,2,15,2,2,2,2,2,2,2,2,4)
# mean(MAD)



# 
pe.pin2 = salso(fit_1$Si,loss = "binder", maxNClusters = 2)
pe.pin3 = salso(fit_1$Si, loss = "binder", maxNClusters = 3)
pe.pin4 = salso(fit_1$Si, loss = "binder", maxNClusters = 4)
pe.pin5 = salso(fit_1$Si, loss = "binder", maxNClusters = 5)
pe.pin6 = salso(fit_1$Si, loss = "binder", maxNClusters = 6)
pe.pin7 = salso(fit_1$Si, loss = "binder", maxNClusters = 7)
pe.pin8 = salso(fit_1$Si, loss = "binder", maxNClusters = 8)
pe.pin10 = salso(fit_1$Si, loss = "binder", maxNClusters = 10)


loc2 <- as.vector(pe.pin2)
passdat2 <- as.data.frame(cbind(loc2,train))

fig21 <- plot_ly(x=passdat2$std_Home.HP, y=passdat2$std_hd_per, z=passdat2$std_hblock_per, type="scatter3d", mode="markers", color = ~passdat2$loc2,
                 colors = c("lightblue", "forestgreen"),
                 symbol = ~as.factor(ytrain),  symbols = c('circle', 'x'), marker = list(size = 5, line = list(
                   width = 0.5
                 )))
fig21 <- fig21 %>% layout(scene = list(xaxis = list(title = 'Hits'),
                                       yaxis = list(title = 'Digs'),
                                       zaxis = list(title = 'Blocks')))
fig21

########################

loc3 <- as.vector(pe.pin3)
passdat3 <- as.data.frame(cbind(loc3,train))

fig31 <- plot_ly(x=passdat3$std_Home.HP, y=passdat3$std_hd_per, z=passdat3$std_hblock_per, type="scatter3d", mode="markers", color = ~passdat3$loc3,
                 colors = c("lightblue", "forestgreen"),
                 symbol = ~as.factor(ytrain),  symbols = c('circle', 'x'), marker = list(size = 5, line = list(
                   width = 0.5
                 )))
fig31 <- fig31 %>% layout(scene = list(xaxis = list(title = 'Hits'),
                                       yaxis = list(title = 'Digs'),
                                       zaxis = list(title = 'Blocks')))
fig31

#######################################################

loc4 <- as.vector(pe.pin4)
passdat4 <- as.data.frame(cbind(loc4,train))

fig41 <- plot_ly(x=passdat4$std_Home.HP, y=passdat4$std_hd_per, z=passdat4$std_hblock_per, type="scatter3d", mode="markers", color = ~passdat4$loc4,
                 symbol = ~as.factor(ytrain),  symbols = c('circle', 'x'), marker = list(size = 5, line = list(
                   width = 0.5
                 )))
fig41 <- fig41 %>% layout(scene = list(xaxis = list(title = 'Hits'),
                                       yaxis = list(title = 'Digs'),
                                       zaxis = list(title = 'Blocks')))
fig41

################################################

loc6 <- as.vector(pe.pin6)
passdat6 <- as.data.frame(cbind(loc6,train))

fig61 <- plot_ly(x=passdat6$std_Home.HP, y=passdat6$std_hd_per, z=passdat6$std_hblock_per, type="scatter3d", mode="markers", color = ~passdat6$loc6,
                 symbol = ~as.factor(ytrain),  symbols = c('circle', 'x'), marker = list(size = 5, line = list(
                   width = 0.5
                 )))
fig61 <- fig61 %>% layout(scene = list(xaxis = list(title = 'Hits'),
                                       yaxis = list(title = 'Digs'),
                                       zaxis = list(title = 'Blocks')))
fig61

############################
loc8 <- as.vector(pe.pin8)
passdat8 <- as.data.frame(cbind(loc8,train))

fig81 <- plot_ly(x=passdat8$std_Home.HP, y=passdat8$std_hd_per, z=passdat8$std_hblock_per, type="scatter3d", mode="markers", color = ~passdat8$loc8,
                 symbol = ~as.factor(ytrain),  symbols = c('circle', 'x'), marker = list(size = 5, line = list(
                   width = 0.5
                 )))
fig81 <- fig81 %>% layout(scene = list(xaxis = list(title = 'Hits'),
                                       yaxis = list(title = 'Digs'),
                                       zaxis = list(title = 'Blocks')))
fig81

###################################
loc10 <- as.vector(pe.pin10)
passdat10 <- as.data.frame(cbind(loc10,train))

fig101 <- plot_ly(x=passdat10$std_Home.HP, y=passdat10$std_hd_per, z=passdat10$std_hblock_per, type="scatter3d", mode="markers", color = ~passdat10$loc10,
                 symbol = ~as.factor(ytrain),  symbols = c('circle', 'x'), marker = list(size = 5, line = list(
                   width = 0.5
                 )))
fig101 <- fig101 %>% layout(scene = list(xaxis = list(title = 'Hits'),
                                       yaxis = list(title = 'Digs'),
                                       zaxis = list(title = 'Blocks')))
fig101



####################################
colnames(passdat2)[1] = 'clusters'
colnames(passdat4)[1] = 'clusters'
colnames(passdat6)[1] = 'clusters'
colnames(passdat8)[1] = 'clusters'
colnames(passdat10)[1] = 'clusters'


# Cluster 2
new21 <- passdat2 %>% 
  filter(clusters == 1)
new21$home_win
new21 %>% group_by(home_win) %>% summarise(count = n())
new21 <- new21 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new21)

new22 <- passdat2 %>% 
  filter(clusters == 2)
new22$home_win
new22 %>% group_by(home_win) %>% summarise(count = n())
new22 <- new22 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new22)

# Cluster 4
new41 <- passdat4 %>% 
  filter(clusters == 1)
new41$home_win
new41 %>% group_by(home_win) %>% summarise(count = n())
new41 <- new41 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new41)

new42 <- passdat4 %>% 
  filter(clusters == 2)
new42$home_win
new42 %>% group_by(home_win) %>% summarise(count = n())
new42 <- new42 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new42)

new43 <- passdat4 %>% 
  filter(clusters == 3)
new43$home_win
new43 %>% group_by(home_win) %>% summarise(count = n())
new43 <- new43 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new43)

new44 <- passdat4 %>% 
  filter(clusters == 4)
new44$home_win
new44 %>% group_by(home_win) %>% summarise(count = n())
new44 <- new44 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new44)

# Cluster 6
new61 <- passdat6 %>% 
  filter(clusters == 1)
new61$home_win
new61 %>% group_by(home_win) %>% summarise(count = n())
new61 <- new61 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new61)

new62 <- passdat6 %>% 
  filter(clusters == 2)
new62$home_win
new62 %>% group_by(home_win) %>% summarise(count = n())
new62 <- new62 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new62)

new63 <- passdat6 %>% 
  filter(clusters == 3)
new63$home_win
new63 %>% group_by(home_win) %>% summarise(count = n())
new63 <- new63 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new63)

new64 <- passdat6 %>% 
  filter(clusters == 4)
new64$home_win
new64 %>% group_by(home_win) %>% summarise(count = n())
new64 <- new64 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new64)

new65 <- passdat65 %>% 
  filter(clusters == 5)
new65$home_win
new65 %>% group_by(home_win) %>% summarise(count = n())
new65 <- new65 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new65)

new66 <- passdat6 %>% 
  filter(clusters == 6)
new66$home_win
new66 %>% group_by(home_win) %>% summarise(count = n())
new66 <- new66 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new66)

# Cluster 8
new81 <- passdat8 %>% 
  filter(clusters == 1)
new81$home_win
new81 %>% group_by(home_win) %>% summarise(count = n())
new81 <- new81 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new81)

new82 <- passdat8 %>% 
  filter(clusters == 2)
new82$home_win
new82 %>% group_by(home_win) %>% summarise(count = n())
new82 <- new82 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new82)

new83 <- passdat8 %>% 
  filter(clusters == 3)
new83$home_win
new83 %>% group_by(home_win) %>% summarise(count = n())
new83 <- new83 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new83)

new84 <- passdat8 %>% 
  filter(clusters == 4)
new84$home_win
new84 %>% group_by(home_win) %>% summarise(count = n())
new84 <- new84 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new84)

new85 <- passdat8 %>% 
  filter(clusters == 5)
new85$home_win
new85 %>% group_by(home_win) %>% summarise(count = n())
new85 <- new85 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new85)

new86 <- passdat8 %>% 
  filter(clusters == 6)
new86$home_win
new86 %>% group_by(home_win) %>% summarise(count = n())
new86 <- new86 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new86)

new87 <- passdat8 %>% 
  filter(clusters == 7)
new87$home_win
new87 %>% group_by(home_win) %>% summarise(count = n())
new87 <- new87 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new87)

new88 <- passdat8 %>% 
  filter(clusters == 8)
new88$home_win
new88 %>% group_by(home_win) %>% summarise(count = n())
new88 <- new88 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new88)

# Cluster 10
new10.1 <- passdat10 %>% 
  filter(clusters == 1)
new10.1$home_win
new10.1 %>% group_by(home_win) %>% summarise(count = n())
new10.1 <- new10.1 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.1)

new10.2 <- passdat10 %>% 
  filter(clusters == 2)
new10.2$home_win
new10.2 %>% group_by(home_win) %>% summarise(count = n())
new10.2 <- new10.2 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.2)

new10.3 <- passdat10 %>% 
  filter(clusters == 3)
new10.3$home_win
new10.3 %>% group_by(home_win) %>% summarise(count = n())
new10.3 <- new10.3 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.3)

new10.4 <- passdat10 %>% 
  filter(clusters == 4)
new10.4$home_win
new10.4 %>% group_by(home_win) %>% summarise(count = n())
new10.4 <- new10 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.4)

new10.5 <- passdat10 %>% 
  filter(clusters == 5)
new10.5$home_win
new10.5 %>% group_by(home_win) %>% summarise(count = n())
new10.5 <- new10.5 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.5)

new10.6 <- passdat10 %>% 
  filter(clusters == 6)
new10.6$home_win
new10.6 %>% group_by(home_win) %>% summarise(count = n())
new10.6 <- new10.6 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.6)

new10.7 <- passdat10 %>% 
  filter(clusters == 7)
new10.7$home_win
new10.7 %>% group_by(home_win) %>% summarise(count = n())
new10.7 <- new10.7 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.7)

new10.8 <- passdat10 %>% 
  filter(clusters == 8)
new10.8$home_win
new10.8 %>% group_by(home_win) %>% summarise(count = n())
new10.8 <- new10.8 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.8)

new10.9 <- passdat10 %>% 
  filter(clusters == 9)
new10.9$home_win
new10.9 %>% group_by(home_win) %>% summarise(count = n())
new10.9 <- new10.9 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.9)

new10.10 <- passdat10 %>% 
  filter(clusters == 10)
new10.10$home_win
new10.10 %>% group_by(home_win) %>% summarise(count = n())
new10.10 <- new10.10 %>% 
  select(clusters, std_Home.HP, std_hblock_per, std_hd_per)
summary(new10.10)

