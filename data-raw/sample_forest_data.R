library(data.table)
library(mvtnorm)
library(linpk)
source("../R/prepare_forest_data.R")

covariate.values <- list(
    sex = c("Male", "Female"),
    wt = c(40, 70, 120),
    crcl = c(50, 100, 150))

reference.values <- list(
    sex = "Male",
    wt = 70,      # kg
    crcl = 100)   # mL/min

model <- function(param, covar) {
    with(c(as.list(param), as.list(covar)), {
        cli <- cl *
            exp(cl_sex * (sex == "Female")) *
            ((wt/70)^cl_wt) *
            ((crcl/90)^cl_crcl)

        vci <- v * ((wt/70)^v_wt)

        amt <- 10 * wt  # 10 units/kg

        tmax <- linpk:::Tmax.oral1cpt(cli, vci, ka)
        t.obs <- sort(unique(c(tmax, seq(0, 8, 0.1))))  # 0-8h

        y <- pkprofile(t.obs, cl=cli, vc=vci, ka=ka, dose=list(amt=amt))
        s <- secondary(y)

        list(cl=cli, vc=vci, cmax=s$Cmax, auc=s$AUC)
    })
}

theta <- fread(text="
       cl         v        ka     cl_wt    cl_sex   cl_crcl      v_wt 
10.303500 12.201500  4.243550  0.583708 -0.255737  0.891369  0.857272 
")


# Sample from the variance-covariance matrix of theta
covmat <- read.table(header=T, text="
            theta1       theta2       theta3       theta4       theta5       theta6       theta7
theta1  0.09795280  0.010994600  1.33304e-03 -0.015920800 -1.07415e-02  0.008765400  0.001717700
theta2  0.01099460  0.089091200 -1.93659e-04  0.002583560 -8.04885e-04  0.006237480  0.007393160
theta3  0.00133304 -0.000193659  5.01075e-03 -0.001949220 -4.14722e-05  0.001420430 -0.001423240
theta4 -0.01592080  0.002583560 -1.94922e-03  0.011589700  2.87094e-03 -0.000247892  0.000493268
theta5 -0.01074150 -0.000804885 -4.14722e-05  0.002870940  2.45663e-03 -0.000131245  0.000453010
theta6  0.00876540  0.006237480  1.42043e-03 -0.000247892 -1.31245e-04  0.014281200  0.000416080
theta7  0.00171770  0.007393160 -1.42324e-03  0.000493268  4.53010e-04  0.000416080  0.015352000
")

params <- as.data.table(rmvnorm(1000, as.numeric(theta), as.matrix(covmat)))
names(params) <- names(theta)
params <- rbind(theta, params) # First row is the final point estimates

sample_forest_data <- prepare_forest_data(model, params, covariate.values, reference.values)

dir.create("../data", F, T)
save(sample_forest_data, file="../data/sample_forest_data.rda")

