# Just running the default values so I can go through the function and debug all
#the things

library(mlbench)
n_ <- 101
sim_train <- mlbench.friedman1(n = n_)  |> as.data.frame()
x_train <- sim_train |> dplyr::select(dplyr::starts_with("x"))
y_train <- sim_train$y
x_test <- mlbench.friedman1(n = n_) |> as.data.frame() |> dplyr::select(dplyr::starts_with("x"))
n_tree <- 1
node_min_size = 5
n_mcmc = 2000
n_burn = 500
alpha = 0.95
beta = 2
df = 3
sigquant = 0.9
kappa = 2
tau = 100
scale_bool = TRUE
stump = FALSE
no_rotation_bool = FALSE
numcut = 100L # Defining the grid of split rules
usequants = FALSE
source("R/other_functions.R")

# Splines parameters
nIknots = 3
dif_order = 2
