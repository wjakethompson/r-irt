### Read in data
dat <- read.csv("Project_data.csv", header = TRUE, stringsAsFactors = FALSE)


### Which items need to be recoded (i.e., score 1, 2, or 1, 2, 3, instead of 0, 1 or 0, 1, 2)?
poly <- 31:35


### Libraries
needed_packages <- c("mirt", "plyr", "ggplot2", "tidyr", "dplyr")
for(i in 1:length(needed_packages)){
    haspackage <- suppressMessages(require(needed_packages[i], character.only = TRUE))
    if(!haspackage){
        install.packages(needed_packages[i])
        suppressMessages(library(needed_packages[i], character.only = TRUE))
    }
}


# Send item data to own object and recode
x <- as.matrix(dat[, seq(2, 36)])
for(i in poly){
    x[,i] <- x[,i] - 1
}


# Calculate total scores
X <- rowSums(x, na.rm = TRUE)

# Get sample characteristics
N <- nrow(x) ##people
n <- ncol(x) ##items
ncat <- rep(0, n)

for (j in 1:n) {
    ncat[j] <- sum(table(x[,j]) >= 0)
}
NC <- max(ncat)


### Classical test theory statistics
#Item Stats
item.means <- rep(0, n)
item.disc <- rep(0, n)
item.freq <- matrix(0, n, NC+1)

item.means <- colMeans(x, na.rm=TRUE)
for (j in 1:n) {
    item.disc[j] <- cor(x[, j], X, use = "pairwise.complete.obs")
    temp <- as.data.frame(table(x[,j]))
    for (l in 1:(NC+1)) {
        for (k in 1:ncat[j]) {
            if (temp$Var1[k]==(l-1)) {
                item.freq[j,l] <- temp$Freq[k]	
            }
        }
    }
}

CTT.items <- cbind(seq(1, n), item.means, item.disc, item.freq)
write.csv(CTT.items, file = "CTT.csv")

### Remove Items 10, 17
x <- x[,-c(10, 17)]
n <- ncol(x)
ncat <- ncat[-c(10,17)]

### IRT Analysis
item.lab <- colnames(x)

# Fit the IRT model to the data
SPECS <- mirt.model('F = 1-33
                    PRIOR = (1-33, a1, lnorm, 0, 1),
                    (1-28, g, norm, -1.39, 1),
                    (1-28, d, norm, 0, 2),
                    PRIOR = (29-33, d1, norm, 0, 2),
                    (29-33, d2, norm, 0, 2),
                    (29-33, d3, norm, 0, 2),
                    (29-33, d4, norm, 0, 2)')

item.list <- c(rep('3PL', 28), rep('graded', 5))

mod <- mirt(data = x, model = SPECS, itemtype = item.list)
RAW <- coef(mod, IRTpars = TRUE)

PARM <- do.call("rbind", lapply(RAW[-length(RAW)], function(y){vec <- as.vector(y);length(vec) <- NC; vec}))
PARM[which(is.na(PARM))] <- 0

a.est <- rep(0, n)
c.est <- rep(0, n)
b.est <- matrix(0, n, NC-1)

for (j in 1:n) {
    a.est[j] <- PARM[j,1]
    b.est[j,1] <- PARM[j,2]
    if (ncat[j] == 2){ 
        c.est[j] <- PARM[j,3]
    } else if (ncat[j] >2){
        b.est[j,2:(ncat[j]-1)] <- PARM[j,3:ncat[j]]
    }
}

theta.est <- fscores(mod, method = "EAP", full.scores = TRUE, scores.only = TRUE)
write.csv(theta.est, file = "theta_est.csv")

# Calculate model probabilities, expected scores
expect <- matrix(0, N, n)
for (i in 1:N){
    for (j in 1:n){
        pstar <- rep(0, ncat[j]+1)
        pstar[1] <- 1
        for (l in 1:(ncat[j]-1)){
            pstar[l+1] <- c.est[j] + (1-c.est[j]) * (1 / (1 + exp(-a.est[j]*(theta.est[i]-b.est[j,l]))))
        }
        for (k in 1:ncat[j]){
            expect[i,j] <- expect[i,j] + (k-1)*(pstar[k] - pstar[k+1])
        }
    }
}

# Evaluation of Local Item Dependence (LID)
d <- x - expect
Q3 <- as.matrix(cor(d, use = "pairwise.complete.obs"))

# Save Q3 matrix to a CSV file, then write unique values to a vector
write.csv(Q3, file = "Q3.csv")
Q3_vec <- unique(Q3[lower.tri(Q3)])

### Residual Aanlysis, Evaluation of Item, Test Fit
NQPTS <- 10                                             #number of quadrature points
S <- (max(theta.est) - min(theta.est)) / NQPTS          #width of quadrature points
B <- vector(mode = "numeric", NQPTS + 1)                #quadrature node cut points
id <- vector(mode = "numeric", N)                       #quadrature node of each person
Num_q <- matrix(data = 0, nrow = NQPTS, ncol = n)       #number of people per node per item
Theta_q <- matrix(data = 0, nrow = NQPTS, ncol = n)     #average theta per node per item
Correct_q <- matrix(data = 0, nrow = NQPTS, ncol = n)   #sum of item scores per node
Prop_q <- matrix(data = 0, nrow = NQPTS, ncol = n)      #average item scores per node
Prob_q <- matrix(data = 0, nrow = NQPTS, ncol = n)      #expected score for each node per item
Var <- matrix(data = 0, nrow = NQPTS, ncol = n)         #variance of residuals
SE <- matrix(data = 0, nrow = NQPTS, ncol = n)          #standard error of residuals
DF <- rep(0, n)                                         #degrees of freedom for chi2 per item
SR <- matrix(data = 0, nrow = NQPTS, ncol = n)          #standardized residuals
Bar_top <- matrix(data = 0, nrow = NQPTS, ncol = n)     #upper 95% CI of expected proportion correct
Bar_bottom <- matrix(data = 0, nrow = NQPTS, ncol = n)  #lower 95% CI of expected proportion correct
chi2 <- vector(mode = "numeric", n)                     #chi2 value for each item
pval <- vector(mode = "numeric", n)                     #pvalues for the chi2 of each item

# Set Quadrature cut points
B[1] <- min(theta.est)
for(k in 2:(NQPTS + 1)){
    B[k] <- B[k-1] + S
}

# Calculate average score and theta per node per item
for(i in 1:N){
    for(k in 1:NQPTS){
        if(theta.est[i] >= B[k]){
            id[i] <- k
        }
    }
    for(k in 1:NQPTS){
        if(id[i] == k){
            for(j in 1:n){
                if(!is.na(x[i,j])){
                    Num_q[k,j] <- Num_q[k,j] + 1
                    Theta_q[k,j] <- Theta_q[k,j] + theta.est[i]
                    Correct_q[k,j] <- Correct_q[k,j] + x[i,j]
                }
            }
        }
    }
}
for(k in 1:NQPTS){
    for(j in 1:n){
        if(Num_q[k,j] >= 5){
            Prop_q[k,j] <- Correct_q[k,j] / Num_q[k,j]
            Theta_q[k,j] <- Theta_q[k,j] / Num_q[k,j]
        } else if(Num_q[k,j] < 5){
            Prop_q[k,j] <- NA
            Theta_q[k,j] <- NA
        }
    }
}

# Calculate standardized residuals
for(j in 1:n){
    for(k in 1:NQPTS){
        if(Num_q[k,j] >= 5){
            pstar <- rep(0, ncat[j]+1)
            pstar[1] <- 1
            for(l in 1:(ncat[j]-1)){
                pstar[l+1] <- c.est[j] + (1 - c.est[j]) * (1 / (1 + exp(-a.est[j] * (Theta_q[k,j] - b.est[j,l]))))
            }
            for(l in 1:ncat[j]){
                Prob_q[k,j] <- Prob_q[k,j] + (l-1) * (pstar[l] - pstar[l+1])
            }
            for(l in 1:ncat[j]){
                Var[k,j] <- Var[k,j] + (pstar[l] - pstar[l+1]) * (1 - (pstar[l] - pstar[l+1]))
            }
            if(Var[k,j] > 0){
                SE[k,j] <- sqrt(Var[k,j] / Num_q[k,j])
                DF[j] <- DF[j] + 1
                SR[k,j] <- (Prop_q[k,j] - Prob_q[k,j]) / SE[k,j]
                Bar_top[k,j] <- Prob_q[k,j] + (1.96 * SE[k,j])
                Bar_bottom[k,j] <- Prob_q[k,j] - (1.96 * SE[k,j])
                chi2[j] <- chi2[j] + (SR[k,j])^2
            }
        } else if(Num_q[k,j] < 5){
            SE[k,j] <- NA
            SR[k,j] <- NA
            Bar_top[k,j] <- NA
            Bar_bottom[k,j] <- NA
        }
    }
    pval[j] <- pchisq(chi2[j], DF[j], lower.tail = FALSE)
}

# Calculate expected scores over theta for each item
tt <- seq(-3, 3, 0.01)
ex <- matrix(data = 0, nrow = length(tt), ncol = n)
for(i in 1:length(tt)){
    for(j in 1:n){
        pstar <- rep(0, ncat[j] + 1)
        pstar[1] <- 1
        for(l in 1:(ncat[j]-1)){
            pstar[l+1] <- c.est[j] + (1-c.est[j]) * (1 / (1 + exp(-a.est[j] * (tt[i] - b.est[j,l]))))
        }
        for(l in 1:ncat[j]){
            ex[i,j] <- ex[i,j] + (l-1) * (pstar[l] - pstar[l+1])
        }
    }
}

### Save item parameter estimates
item.est <- cbind(item.lab, a.est, c.est, b.est, chi2, pval, DF)
write.csv(item.est, file = "IRT_par_est.csv", row.names = FALSE)


### Item Information
INFO <- matrix(data = 0, nrow = length(tt), ncol = n)
for(i in 1:length(tt)){
    for(j in 1:n){
        pstar <- rep(0, ncat[j] + 1)
        pstar[1] <- 1
        for(l in 1:(ncat[j]-1)){
            pstar[l+1] <- c.est[j] + (1 - c.est[j]) * (1 / (1 + exp(-a.est[j] * (tt[i] - b.est[j,l]))))
        }
        for(l in 1:ncat[j]){
            pp1 <- pstar[l] - pstar[l+1]
            a1 <- ((a.est[j]^2) * (pstar[l] * (1 - pstar[l]) - pstar[l+1] * (1 - pstar[l+1]))^2) / pp1
            INFO[i,j] <- INFO[i,j] + a1
        }
    }
}

TEST_INFO <- rowSums(INFO)

### Correct Errorbars to remove values > 1 and < 0
for(i in 1:nrow(Bar_bottom)){
    for(j in 1:ncol(Bar_bottom)){
        Bar_bottom[i,j] <- max(0, Bar_bottom[i,j])
        Bar_top[i,j] <- min((ncat[j]-1), Bar_top[i,j])
    }
}


###Open up a PDF file to save all the graphs
pdf(file="IRT_graphs.pdf")

# Histogram of total scores
X <- rowSums(x,na.rm=TRUE)
Scores <- as.data.frame(rowSums(x,na.rm=TRUE))
ggplot(Scores, aes(x = Scores[,"rowSums(x, na.rm = TRUE)"])) +
    geom_histogram(binwidth = 4, color = "blue", fill = "blue", alpha = 0.2) +
    labs(x = "Total Score", y = "Count", title = "Total Score Distribution") +
    scale_x_continuous(breaks = seq(0, 60, 10)) +
    scale_y_continuous(breaks = seq(0, 150, 20))

# Histogram of theta estimates
theta.est <- as.vector(theta.est)
theta <- data.frame(theta = theta.est)
ggplot(theta, aes(x = theta)) +
    geom_histogram(binwidth = 0.5, color = "orange", fill = "orange", alpha = 0.2) +
    labs(x = expression(theta), y = "Count", title = "Distribution of Latent Trait Scores") +
    scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
    scale_y_continuous(breaks = seq(0, 200, 20))

# Scatterplot of total score vs. theta
Comparison <- data.frame(TotalScore = X, Theta = theta.est)
ggplot(Comparison, aes(x = Theta, y = TotalScore)) +
    geom_point(alpha = 0.4) +
    labs(x = expression(theta), y = "Total Score", title = "Total Score Comparison to Theta") +
    scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
    scale_y_continuous(lim = c(0, 55), breaks = seq(0, 60, 5))

# Test information
info <- data.frame(tt = tt, TEST_INFO = TEST_INFO, SE = 1/sqrt(TEST_INFO))
ggplot(info, aes(x = tt, y = TEST_INFO)) +
    geom_line(size = 1, color = "red") +
    labs(x = expression(theta), y = "Information", title = "Test Information Across Theta") +
    scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
    scale_y_continuous(lim = c(0, 50), breaks = seq(0, 50, 10))

ggplot(info, aes(x = tt, y = SE)) +
    geom_line(size = 1, color = "green") +
    labs(x = expression(theta), y = "Standard Error", title = "Standard Error Across Theta") +
    scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25))

# Item Plots
ItemCharac <- as.data.frame(cbind(tt, ex))
Obs <- as.data.frame(cbind(Theta_q[,1], Prop_q))
for(j in 1:n){
    ymax <- ncat[j] - 1
    ybreak <- ifelse(ncat[j] == 2, 0.1, 0.5)
    coordmin <- -ybreak / 2
    coordmax <- ybreak / 2
    graphTitle <- bquote(atop('Item'~.(j)~'Fit',
                              chi^2~"="~.(round(chi2[j],2))*","~italic(p)~"="~.(round(pval[j],4))))
    
    print(
    ggplot() + 
        geom_errorbar(aes(x = Obs[,1], ymin = Bar_bottom[,j], ymax = Bar_top[,j]), width = 0.1, color = "red") +
        geom_line(aes(x = ItemCharac[,"tt"], y = ItemCharac[,j+1]), size = 1) +
        scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
        scale_y_continuous(breaks = seq(0, ymax, ybreak)) +
        coord_cartesian(ylim = c(coordmin, ymax+coordmax)) +
        geom_point(aes(x = Obs[,1], y = Obs[,j+1]), color = "blue", size = 3) + 
        labs(x = expression(theta), y = expression(paste("E(x|",theta,")")),
             title = graphTitle)
    )
    
    print(
    ggplot() +
        geom_line(aes(x = tt, y = INFO[,j]), color = "red", size = 1) +
        scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
        expand_limits(y = 0) +
        labs(x = expression(theta), y = "Item Information", title = paste0("Item ", j, " Information"))
    )
}

# Density function of Q3 stats
xlab <- bquote(atop(italic(Q)[3]~'(Residual Correlations)',
                    "mean ="~.(round(mean(Q3_vec, na.rm = TRUE), 3))*','~'sd ='~
                        .(round(sd(Q3_vec, na.rm = TRUE), 3))))
binwidth_q <- (max(Q3_vec) - min(Q3_vec)) / 30

ggplot() +
    geom_histogram(aes(x = Q3_vec, y = ..density..),
                   binwidth = 0.02, color = "purple", fill = "purple", alpha = 0.5) +
    geom_density(aes(x = Q3_vec), color = "black", fill = "purple", alpha = 0.2) +
    scale_x_continuous(lim = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, 0.05)) +
    labs(x = xlab,
         y = "Density", title = "Local Independence Evaluation")

# Density function of standardized residuals
SR_plot <- as.vector(SR)
binwidth_sr <- (max(SR_plot) - min(SR_plot)) / 30
xmin <- round_any(min(SR_plot), 1, f = floor)
xmax <- round_any(max(SR_plot), 1, f = ceiling)
outer <- max(abs(xmin), abs(xmax))
xmin <- -outer
xmax <- outer

ggplot() +
    geom_histogram(aes(x = SR_plot, y = ..density..),
                   binwidth = 0.3, color = "gold", fill = "gold", alpha = 0.5) +
    geom_density(aes(x = SR_plot), color = "black", fill = "gold", alpha = 0.2) +
    scale_x_continuous(lim = c(xmin, xmax), breaks = seq(xmin, xmax, 1)) +
    labs(x = "Standardized Residuals", y = "Density", title = "Test-Level Fit")

dev.off()


### ICCs
pdf(file = "CharacteristicCurves.pdf")

# TCC
ggplot() +
    geom_line(aes(x = tt, y = rowSums(ex, na.rm = TRUE)), size = 1, color = "red") +
    scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
    expand_limits(y = 0) +
    scale_y_continuous(breaks = seq(0, 100, 10)) + 
    labs(x = expression(theta), y = "Expected Score", title = "Test Characteristic Curve")

# ICCs
for(i in 1:n){
    PlotData <- data.frame(Theta = seq(-3, 3, 0.01))
    if(ncat[i] == 2){
        PlotData$Expected <- c.est[i] + ((1-c.est[i]) / (1 + exp(-a.est[i] * (PlotData$Theta - b.est[i,1]))))
        print(
            ggplot(PlotData, aes(x = Theta, y = Expected)) +
                geom_line(size = 1) +
                scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
                scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.1)) +
                labs(x = expression(theta), y = "Probability of Success",
                     title = paste0("ICC: Item ", i))
        )
    } else{
        for(c in 1:(ncat[i]+1)){
            if(c == 1){
                PlotData[,paste0("Cum", c)] <- 1
            } else if(c == (ncat[i]+1)){
                PlotData[,paste0("Cum", c)] <- 0
            } else{
                PlotData[,paste0("Cum", c)] <- 1 / (1 + (exp(-a.est[i] * (PlotData$Theta - b.est[i,c-1]))))
            }
        }
        PlotProb <- PlotData %>% select(Theta)
        for(c in 2:(ncol(PlotData)-1)){
            PlotProb[,paste0("Prob", c-1)] <- PlotData[,c] - PlotData[,c+1]
        }
        PlotData <- gather(PlotData, "Category", "CumProb", 2:ncol(PlotData))
        PlotProb <- gather(PlotProb, "Category", "Probability", 2:ncol(PlotProb))
        
        CumLabels <- paste0("P*", 1:(ncat[i]+1))
        
        print(
            ggplot(PlotData, aes(x = Theta, y = CumProb, color = Category)) +
                geom_line(size = 1) +
                scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
                scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.1)) +
                labs(x = expression(theta), y = expression(paste("P(X >= x|", theta, ")")),
                     title = paste0("Cumulative ICC: Item ", i)) +
                theme(legend.position = "bottom", legend.title = element_blank()) +
                scale_color_discrete(labels = CumLabels)
        )
        
        print(
            ggplot(PlotProb, aes(x = Theta, y = Probability, color = Category)) +
                geom_line(size = 1) +
                scale_x_continuous(lim = c(-3, 3), breaks = seq(-3, 3, 1)) +
                scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.1)) +
                labs(x = expression(theta), y = expression(paste("P(X = x|", theta, ")")),
                     title = paste0("Category ICC: Item ", i)) +
                theme(legend.position = "bottom", legend.title = element_blank()) +
                scale_color_discrete(labels = c(1:ncat[i]))
        )
    }
}
dev.off()
