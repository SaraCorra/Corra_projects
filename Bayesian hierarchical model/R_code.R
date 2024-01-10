library (bookdown)
library (officedown)
library(metadat) # datasets
library(tidyverse)
library(dplyr)
library(mcmcr) # as.mcmc
library(ggplot2)
library(RColorBrewer) # palettes
library(scales) # show_col
library(wesanderson) # yellow
library(patchwork)
library(coda)
library(flextable)
library(latex2exp) # latex in ggplots
library(LaplacesDemon)
library(gridExtra)
library(kableExtra)
library (HDInterval) # hdi's

mypal <- brewer.pal (11, 'RdBu')
mypal [12] <- wes_palette(name = ("Cavalcanti1"), n = 5 ) [1]

# data featuring
infarct <- dat.yusuf1985 
infarct <- na.omit (infarct)
infarct$ntot <- infarct$n1i + infarct$n2i
infarct <- infarct %>%
  filter (ai != 0 & ci != 0 & table == '9')
first.ratio <- infarct$ai / (infarct$n1i - infarct$ai)
second.ratio <- infarct$ci / (infarct$n2i - infarct$ci)
infarct$units <- log (first.ratio) - log (second.ratio)
infarct$sigma2 <- 1/infarct$ai + 1/(infarct$n1i-infarct$ai) +
  1/infarct$ci + 1/(infarct$n2i-infarct$ci)
units <- infarct$units

# dataset 
descrplot_df <- data.frame (study = infarct$trial, logodds = infarct$units, sd = infarct$sigma2, tot = infarct$ntot, study = infarct$trial)

# confidence intervals
descrplot_df$lower_ci <-descrplot_df$logodds - 1.96 * sqrt(1/descrplot_df$tot)
descrplot_df$upper_ci <- descrplot_df$logodds + 1.96 * sqrt(1 /descrplot_df$tot)
descrplot_df_order <- descrplot_df[order(descrplot_df$logodd), ]

plot_forest <- ggplot(descrplot_df_order, aes(x = logodds, y = study)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_segment(aes(x = lower_ci, xend = upper_ci, y = study, yend = study), color = mypal [1], size = 1) +
  geom_point(aes(x = logodds, y = study, size = tot), color = mypal [1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = study), height = 0.2, color = mypal [1]) +
  xlim(c(-1.5, 1.3)) + 
  labs(x = "Empirical logits", y = "", title = "Forest plot of empirical logits", size = "Patients") +
  scale_size(breaks = c(50, 500, 1000, 5000), labels = c("50", "500", "1000", "5000")) +
  theme_bw () +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = c(.98, .98),
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.key.size = unit(0.1, "cm"))
plot_forest

# checking normality
df <- data.frame(units = infarct$units)

norm.hist=ggplot(df, aes(x = units)) +
  geom_histogram(aes(y = ..density..),binwidth = 0.55,
                 color = "lightsteelblue4", fill = "lightsteelblue2", alpha = 0.5)+ 
  labs(title = 'Histogram', y='', x='')+
  theme_bw()+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


norm.box=ggplot(df, aes(y = units)) +
  geom_boxplot(color = "lightsteelblue4", fill = "lightsteelblue2", alpha = 0.5)+ 
  labs(title = 'Boxplot', y='', x='')+
  theme_bw()+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


qqplot=ggplot(df, aes(sample = units)) +
  stat_qq(color = 'grey38', size = 1) +
  stat_qq_line(col='lightsteelblue4', size=1)+
  labs(title = "Normal QQ plot", y='', x='') +
  theme_bw()+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

norm.hist+norm.box+qqplot

# variables needed for the gibbs sampling
sigmaj=infarct$sigma2
m <- nrow (infarct)
n.vec <- infarct$ntot # length 15
ybar.vec <- infarct$units # length 15
ybar.tot <- mean(ybar.vec)
var.between <- var(ybar.vec) # of means, it's a scalar
var.within <- mean(infarct$sigma2) # length 15

### FULL CONDITIONALS
mu_fcl <- function(mu0, gamma20, theta.vec, m, tau2) {
  gamma2n <- (m/tau2 + 1/gamma20)^-1
  mun <- gamma2n*(mu0/gamma20 + mean(theta.vec)*m/tau2)
  draw <- rnorm(1, mean = mun, sd = sqrt(gamma2n))
  return(draw)
}

tau2_fcl <- function(tau20, eta0, theta.vec, m, mu) {
  etan <- eta0 + m
  tau2n <- (eta0*tau20 + sum((theta.vec - mu)^2))/etan
  draw <- 1/rgamma(1, shape = etan/2, rate = etan*tau2n/2)
  return(draw)
}

theta_fcl <- function(ybar.vec,  mu, sigma2, tau2) {
  sigma <- diag((1/sigma2 + 1/tau2)^-1)
  mean <- (ybar.vec/sigma2 + mu/tau2)/(1/sigma2 + 1/tau2)
  draw <- mvtnorm::rmvnorm(1, mean = mean, sigma = sigma)
  return(draw)
}

# Gibbs function 
# It allows to set the burnin, thin the chain and set initial values.
gibbs <- function(G, burnin, thin, inits,
                  data, ybar.vec, m, 
                  mu0, gamma20, eta0, tau20, sigma20) {
  n.iter <- burnin + G*thin
  current_state <- inits
  
  sample <- matrix(0, nrow = G, ncol = 17)
  
  g <- 1
  
  for (i in 1:n.iter) {
    current_state[1] <- mu_fcl(mu0, gamma20, theta.vec = current_state[3:17], m, tau2 = current_state[2])
    current_state[2] <- tau2_fcl(tau20, eta0, theta.vec = current_state[3:17], m, mu = current_state[1])
    current_state[3:17] <- theta_fcl(ybar.vec, mu = current_state[1], 
                                     sigma2 = sigma20, tau2 = current_state[2])
    
    if ((i>burnin) & (i%%thin == 0)) {
      sample[g,] <- current_state
      g <- g + 1
    }
  }
  
  colnames(sample) <- c("mu", "tau2", paste("theta", 1:15, sep = ""))
  return(sample)
}

G <- 5000; burnin <- 1000; thin <- 4
inits <- c(ybar.tot, var.between, ybar.vec)
mu0 <- 0; gamma20 <- 4
tau20 <- 4; eta0 <- 2
sigma20 = infarct$sigma2

set.seed(1)
out <- gibbs(G, burnin, thin, inits,
             infarct, ybar.vec, m, 
             mu0, gamma20, eta0, tau20, sigma20)
#saveRDS (out, "gibbs.RDS")

#out <- readRDS ("gibbs.RDS")
mu <- out [,1] # extract the corresponding column in the mcmc output
tau2 <- out [,2]
theta1 <- out [,3]
theta_hats <- apply (out [, 3:17], 2, mean) 
# I computed the mean for each column
probs <- apply (out [, 3:17], 2, function (x) mean (x < 0))
# average number of times the mean is equal to zero 
# it means the risk is increased

# Convergence 

df.out=as.data.frame(out)
colnames(df.out)=c('mu','tau',paste0('theta',1:15))

g1=ggplot(df.out, aes(x = mu)) +
  geom_histogram(aes(y = ..density..),binwidth = 0.1, color = "white", fill = "lightsteelblue2", alpha = 0.5)+ 
  geom_density(aes(x=mu),color = "lightsteelblue4", lwd = 1.3, adjust=2) +
  geom_errorbar(aes(y=0,xmin=quantile(mu,0.05), xmax=quantile(mu,0.95), width=0.02))+
  labs(x = "Values", y = '', 
       title = TeX("$Marginal\\  Posterior \\ \\mu$"))+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

autocorr <- acf(df.out$mu, plot = FALSE)

# Dataframe with lag values and autocorrelation
df <- data.frame(Lag = autocorr$lag, Autocorrelation = autocorr$acf)

a1=ggplot(df, aes(x = Lag, y = Autocorrelation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Lag") +
  ylab("") +
  labs(title='Autocorrelation')+theme_bw()+
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

g2=ggplot(df.out, aes(x = tau)) +
  geom_histogram(aes(y = ..density..),binwidth = 0.2, color = "white", fill = "lightsteelblue2", alpha = 0.5)+ 
  geom_density(aes(x=tau),color = "lightsteelblue4", lwd = 1.3, adjust=2) +
  geom_errorbar(aes(y=0,xmin=quantile(tau,0.05), xmax=quantile(tau,0.95), width=0.02))+
  labs(x = "Values", y = '', 
       title = TeX("$Marginal \\  Posterior \\ \\tau$"))+
  coord_cartesian(xlim=c(0,3.3))+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

autocorr <- acf(df.out$tau, plot = FALSE)

df <- data.frame(Lag = autocorr$lag, Autocorrelation = autocorr$acf)

a2=ggplot(df, aes(x = Lag, y = Autocorrelation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Lag") +
  ylab("") +
  labs(title='Autocorrelation')+
  theme_bw()+
  theme(plot.title = element_text(size = 12, hjust = 0.5))

g3=ggplot(df.out, aes(x = theta1)) +
  geom_histogram(aes(y = ..density..),binwidth = 0.2, color = "white", fill = "lightsteelblue2", alpha = 0.5)+ 
  geom_density(aes(x=theta1),color = "lightsteelblue4", lwd = 1.3, adjust=2) +
  geom_errorbar(aes(y=0,xmin=quantile(theta1,0.05), xmax=quantile(theta1,0.95), width=0.02))+
  labs(x = "Values", y = '', 
       title = TeX("$Marginal \\  Posterior \\ \\theta_1$"))+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

autocorr <- acf(df.out$theta1, plot = FALSE)

df <- data.frame(Lag = autocorr$lag, Autocorrelation = autocorr$acf)

a3=ggplot(df, aes(x = Lag, y = Autocorrelation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Lag") +
  ylab("") +
  labs(title='Autocorrelation')+
  theme_bw()+
  theme(plot.title = element_text(size = 12, hjust = 0.5))

##################################

stationarity.plot <- function(x, ...) {
  S <- length(x)
  scan <- 1:S
  ng <- min(round(S/100), 10)
  group <- S * ceiling(ng * scan / S) / ng
  
  df <- data.frame(x = x, group = group)
  
  getPalette = colorRampPalette(brewer.pal(9, "PuBu"))
  
  ggplot(df, aes(x = as.factor(group), y = x, fill = as.factor(group))) +
    geom_boxplot(outlier.size=0.8) +
    xlab("iteration") +
    scale_fill_manual(values = getPalette(10)) +  
    guides(fill = FALSE)+
    labs(y='')+
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
    scale_x_discrete(breaks=c(0,2500,5000))
}

box.1=stationarity.plot(out[, 1])+labs(title=TeX('$Stationarity \\ \\mu$'))
box.2=stationarity.plot(out[, 2])+labs(title=TeX('$Stationarity \\ \\tau$'))
box.3=stationarity.plot(out[, 3])+labs(title=TeX('$Stationarity \\ \\Theta_1$'))

pl<- list(g1,a1,box.1,g2,a2,box.2,g3,a3,box.3)

grid.arrange(grobs= lapply(pl, "+", theme(plot.margin=margin(4,0.2,4,0.2))), nrow=3, ncol=3)

# I check the distribution for all parameters
color=hcl.colors(n = 17,palette = "Spectral")
par(mfrow=c(3,5))
for (i in 3:17){
  hist(out[,i], br=20, main=colnames(out)[i], xlab='',
       col=color[i])
}

# I check the overall autocorrelation 
par(mfrow=c(3,5))
for (i in 3:17){
  acf(out[,i],main=colnames(out)[i])
}

stationarity.plot_simple<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }

par(mfrow=c(3,5))
for (i in 3:17){
  stationarity.plot_simple(out[,i],xlab="iteration",
                           ylab='',main=colnames(out)[i], col=color[i])
}

ess=sapply(1:17, function(x) effectiveSize(mcmc(out[,x])) ) 
df.ess=as.data.frame(t(ess))
colnames(df.ess)=c('mu','tau', paste('theta',1:15))
df.ess.first=df.ess[1,1:9]; df.ess.first
df.ess.second=df.ess[1,10:17]; df.ess.second

# Monte Carlo Standard error
MCSE <-  apply(out,2,sd)/sqrt(effectiveSize(out))

df_tails <- data.frame(Evemy = out[, 4], Mueller = out[, 6])
melted_df <- df_tails %>% pivot_longer(everything(), names_to = "Variable", values_to = "Value")

plot_tails <- ggplot(melted_df, aes(x = Value, fill = Variable)) +
  geom_density(alpha = 0.3, adjust = 2) +
  labs(x = "Estimated effect", y = "Density") +
  ggtitle ("Posterior distribution for the effects of Evemy and Mueller")+
  scale_fill_manual(values = c("Evemy" = mypal [2], "Mueller" = mypal [9])) +
  guides(fill = guide_legend(title = "Variables")) +
  theme_bw () +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
plot_tails

hdi_mu <- hdi (mu)
hdi_tau <- hdi (tau2)
sd_between <- sqrt (tau2)
med <- median (sd_between)

plot_shrink1 <- ggplot(infarct, aes(units, theta_hats)) +
  geom_point(size = 2, col = mypal [1], pch = 16) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  labs(x = TeX ("$\\y_{j}$"), y = TeX ("$\\theta_{j}$")) +
  theme_bw ()

varsorted_df <- data.frame (thetas = theta_hats, ys = ybar.vec, sigmas = infarct$sigma2)
varsorted_df <- varsorted_df %>%
  arrange(sigmas)

plot_shrink2 <- ggplot(data = varsorted_df, aes (x = sigmas, y = (ys - thetas))) +
  geom_point(aes(x = sigmas, y = (ys - thetas)), size = 2, shape = 16) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  labs(x = "variance", y = TeX ("$\\y_{j} - \\theta_{j}$")) +
  theme_bw ()

plot_shrink3 <- ggplot(varsorted_df) +
  geom_hline(yintercept = c(2.5, 0.5), linetype = "dashed") +
  geom_segment(aes(x = thetas, xend = ys, y = 0.5, yend = 2.5), color = mypal[2]) +
  geom_point(aes(x = thetas, y = 0.5, color = "Bayesian effects"), size = 2) +
  geom_point(aes(x = ys, y = 2.5, color = "Empirical logits"), size = 2) +
  geom_line(aes(x = thetas, y = 0.5, group = 1), color = mypal[3]) +
  geom_line(aes(x = ys, y = 2.5, group = 1), color = mypal[1]) +
  ylim(0.4, 3.55) +
  xlab("Values") +
  ylab("") +
  labs(color = "") +  
  scale_color_manual(values = c("Bayesian effects" = mypal[3], "Empirical logits" = mypal[1])) + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = c(.95, .97),
    legend.justification = c(1, 1),
    legend.direction = "horizontal")

plot_shrink <-(plot_shrink1 + plot_shrink2 + plot_shrink3) + plot_annotation(title = "Shrinkage plots", theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
plot_shrink

color_values <- c(mypal[2], mypal[10])
labels <- c("Norris", "von Essen")

plot_6vs10 <- ggplot() +
  scale_color_manual(values = c (mypal [4], mypal [8], "black"), labels = c ("Norris", "von Essen", "Population mean")) +
  geom_segment(aes(x = mean(mu), xend = mean(mu), y = 0, yend = 0.96, color = "black"),  lwd = 1.1) +
  stat_density(aes(out[, 8]), color = mypal [10], fill = mypal [10], adjust = 2, alpha = 0.07, lwd = 1.1) +
  stat_density(aes(out[, 12]), color = mypal [2], fill = mypal [2], adjust = 2, alpha = 0.07, lwd = 1.1) +
  #geom_hline (yintercept = 0, lwd = .8, color = 1) +
  geom_errorbarh(aes(xmin = descrplot_df[6, 6], xmax = descrplot_df[6, 7], y = -0.06, color = mypal [8]), height = 0.03, size = 1) +
  geom_errorbarh(aes(xmin = descrplot_df[10, 6], xmax = descrplot_df[10, 7], y = -0.11, color = mypal [4]), height = 0.03, size = 1) +
  geom_point(aes(infarct[6, 9], -0.06), color = mypal [4],  size = 3, shape = 18) +
  geom_point(aes(infarct[10, 9], -0.11), color = mypal [8], size = 3, shape = 18) +
  labs(x = "estimates", y = "Relative frequency") +
  coord_cartesian(xlim = c(-2.2, 2.2), ylim = c(-0.11, 1.01)) +
  guides (color = guide_legend(title = NULL), 
          fill = guide_legend(title = NULL)) +
  ggtitle("Posterior distribution for Norris and von Essen studies' parameters")  +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = c(.95, .95),
    legend.justification = c(1, 1),
    legend.direction = "vertical")
plot_6vs10

mn= function(mu,tau2,units,sigma){
  
  mn.num=mu/tau2+units/sigma
  mn.den=1/tau2+1/sigma
  
  return(mn.num/mn.den)
}

mu_sens <-c (-1, -.5, 0, .5, 1); gamma_sens <- c (1, 2, 4, 7, 10)
par.grid <- expand.grid (mu_sens, gamma_sens)
set.seed (1)
sens_out.m <- sapply(1:nrow(par.grid), function(x) {
  gibbs (G, burnin, thin, inits, infarct, ybar.vec, m, par.grid [x, 1], 
         par.grid [x, 2], eta0, tau20, sigma20) [, 1]
})
saveRDS (sens_out.m, "sens_mugamma.RDS")

tau_sens <- c(1, 1.5, 2, 2.5, 4, 5, 7); eta_sens <- 1:7
par.grid <- expand.grid (tau_sens, eta_sens)
set.seed (1)
sens_out.t <- sapply(1:nrow(par.grid), function(x) {
  gibbs (G, burnin, thin, inits, infarct, ybar.vec, m, mu0, gamma20, par.grid [x, 2], 
         par.grid [x, 1], sigma20) [, 1]
})
saveRDS (sens_out.t, "sens_taueta.RDS")


sens_out.m <- readRDS ("sens_mugamma.RDS")
mu_sens <-c (-1, -.5, 0, .5, 1); gamma_sens <- c (1, 2, 4, 7, 10)
par.grid.m <- expand.grid (mu_sens, gamma_sens)
par.grid.m$prob <- apply (sens_out.m, 2, function (x) mean (x < 0))
plot_sens.m <- ggplot(data=par.grid.m, aes(factor(par.grid.m$Var1), factor(par.grid.m$Var2), fill= par.grid.m$prob)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label = round(prob,3)), color = "white", size = 2.5)+
  scale_fill_gradient(low='sandybrown', high='firebrick4')+
  scale_x_discrete(name = TeX("$\\mu_0$"), expand = c(0, 0))+
  scale_y_discrete(name = TeX("$\\gamma^2_0$"), expand = c(0, 0))+
  guides(fill = guide_colourbar (title = "prob", format = "%0.2f")) + 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.width = unit (.7, "cm"),
        legend.key.height = unit(.3, "cm"))

sens_out.t <- readRDS ("sens_taueta.RDS")
tau_sens <- c(1,1.5, 2,2.5, 3,5,7); eta_sens <- 1:7
par.grid.t <- expand.grid (tau_sens, eta_sens)
par.grid.t$prob <- apply (sens_out.t, 2, function (x) mean (x < 0))
plot_sens.t <- ggplot(data=par.grid.t, aes(factor(par.grid.t$Var2), factor(par.grid.t$Var1), fill= par.grid.t$prob)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label = round(prob,3)), color = "white", size = 2.5)+
  scale_fill_gradient(low='sandybrown', high='firebrick4')+
  scale_x_discrete(name = TeX("$\\eta_0$"), expand = c(0, 0))+
  scale_y_discrete(name = TeX("$\\tau^2_0$"), expand = c(0, 0))+
  guides(fill = guide_colourbar (title = "prob", format = "%0.2f")) + 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.width = unit (.7, "cm"),
        legend.key.height = unit(.3, "cm"),
        legend.justification = "center") # +
#ggtitle('Sensitivity analysis for probabilities')+
#theme(plot.title = element_text(size = 17, face = "bold",hjust = 0.5),
#     axis.title.y=element_text(angle=0, vjust=0.6))
plot_sens <-(plot_sens.m + plot_sens.t)  + plot_annotation(title = "Sensitivity analysis for probabilities", theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
plot_sens


var.within <- mean(infarct$sigma2)
tau=out[,2]
R.post=tau/(var.within+tau)

prior=data.frame(mu.prior=rnorm(5000, mu0, sqrt(gamma20)),
                 tau2.prior=1/rgamma(5000, shape = eta0/2, rate = eta0*tau20/2))

R.prior <- prior$tau2.prior/(var.within + prior$tau2.prior)

R.df <- data.frame(distr = rep(c("prior", "posterior"), each = G),
                   value = c(R.prior, R.post))

prior.df=as.data.frame(R.prior)
post.df=as.data.frame(R.post)
df=cbind(prior.df,post.df)
prior_density=density(df$R.prior)
post_density=density(df$R.post)

y.axis <- TeX('$\\frac{\\tau}{\\tau + \\sigma} $')
ggplot(df) + 
  geom_density(aes(x=R.prior, col='Prior'),
               alpha = 0.1, size=1, adjust=2) + 
  geom_ribbon(data = data.frame(x = prior_density$x, y = prior_density$y),
              aes(x = x, ymin = 0, ymax = y), fill = "lightskyblue3", alpha = 0.2) +
  geom_density(aes(x=R.post, col='Posterior'), adjust=2,
               alpha = 0.1, size=1) +
  geom_ribbon(data = data.frame(x = post_density$x, y = post_density$y),
              aes(x = x, ymin = 0, ymax = y), fill = "salmon", alpha = 0.2) +
  geom_vline(xintercept=0.5, linetype='dashed')+
  coord_cartesian(xlim=c(0.2,1))+
  theme_bw()+
  labs(x = 'values', y = y.axis, title = 'Variance comparison', color = '') +
  theme(legend.position = 'bottom')+
  scale_color_manual(name='Density',
                     breaks=c('Prior', 'Posterior'),
                     values=c('Prior'='blue', 'Posterior'='red'))+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.background = element_blank(), legend.position = 'bottom',
        axis.title.y=element_text(angle=0, vjust=0.6))

quantInv <- function(distr, value) ecdf(distr)(value)
# reference https://stackoverflow.com/questions/9123800/how-do-i-calculate-the-probability-for-a-given-quantile-in-r

prob.prior=1-quantInv(R.prior, 0.5)  # 0.9942
prob.post=1-quantInv(R.post, 0.5)   # 0.669


mu0 <- 0; gamma20 <- 4
tau20 <- c(1,1.5, 2,2.5, 3,5,7); eta0 <- 1:7
sigma20 = infarct$sigma2; 

priors=expand.grid(eta0,tau20)
colnames(priors)=c('eta0', 'tau20')

priors=priors %>% 
  mutate(mu0=mu0) %>% 
  mutate(gamma20=gamma20) 


sens.out=mapply(function(mu0, gamma20, tau20, eta0) {
  
  out=gibbs(G, burnin, thin, inits,
            infarct, ybar.vec, m, 
            mu0, gamma20, eta0, tau20, sigma20)
  
}, priors$mu0, priors$gamma20, priors$tau20, priors$eta0)

tau=sapply(1:nrow(priors), function(x) sens.out[5001:10000,x])

tau.prior=1/rgamma(5000, shape = eta0/2, rate = eta0*tau20/2)

R.post=sapply(1:nrow(priors), function(x) tau[,x]/(var.within+tau[,x]))
R.prior <- tau.prior/(var.within + tau.prior)

probs=sapply(1:nrow(priors), function(x) 1-quantInv(R.post[,x], 0.5))

post.ratio=expand.grid(eta0,tau20) %>% 
  mutate(probs=probs) 
colnames (post.ratio) <- c ("eta0", "tau20", "probs")

saveRDS(post.ratio,'post_ratio_new.RDS')


post.ratio=readRDS('post_ratio_new.RDS')

ggplot(data=post.ratio, aes(factor(eta0), factor(tau20), fill= probs)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label = round(probs,3)), color = "white", size = 4.5)+
  scale_fill_gradient(low='sandybrown', high='firebrick4')+
  scale_x_discrete(name = TeX("$\\eta_0$"), expand = c(0, 0))+
  scale_y_discrete(name = TeX("$\\tau^2_0$"), expand = c(0, 0))+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 10,
                                title = "Prob"))+
  ggtitle('Sensitivity analysis for variance comparison')+
  theme(plot.title = element_text(size = 17, face = "bold",hjust = 0.5),
        axis.title.y=element_text(angle=0, vjust=0.6))

# separate and pooled model 

x <- seq(-4, 8, length.out = 10000)
tau20 = 5

data=data.frame(estimate=infarct$units, 
                sigma=sqrt(infarct$sigma2),
                v=rep(sqrt(tau20),15))

normal.sep.g <- mapply(function(estimate, sigma, v, x) {
  vn <- (1 / sigma^2 + 1 / v^2)^(-1 / 2)
  mn <- (m / v^2 + estimate / sigma^2) * vn^2
  dnorm(x,mn,vn)
}, data$estimate, data$sigma, data$v, MoreArgs = list(x = x)) %>% 
  as.data.frame() %>%  cbind(x) %>% 
  pivot_longer(cols = !x, names_to = "ind", values_to = "p")

# POOLED

sigma.p=weighted.mean(sqrt(infarct$sigma2), infarct$ntot)
estimate.p=weighted.mean(infarct$units, infarct$ntot)

v=sqrt(tau20)

vn.p <- (1 / sigma.p^2 + 1 / v^2)^(-1 / 2)
mn.p <- (m / v^2 + estimate.p / sigma.p^2) * vn.p^2

df_pool <- data.frame(x = x, p = dnorm(x, mn.p, vn.p))

# PLOT

g1=ggplot(data = normal.sep.g) +
  geom_line(aes(x = x, y = p, group = ind, col='Separate'), size=1) +
  geom_line(data=df_pool, aes(x = x, y = p, col='Pooled'), size=1) +
  labs(x = expression(theta), y = '', title = 'Separate and Pooled model', color = '') +
  theme_bw()+
  scale_color_manual(name='Method',
                     breaks=c('Separate', 'Pooled'),
                     values=c('Separate'='blue', 'Pooled'='red'))+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.background = element_blank(), legend.position = 'bottom')

theta.out=out[,3:17]
theta.melt=reshape::melt(theta.out)

g2=ggplot(data = theta.melt) +
  geom_density(aes(x = value, group = X2), col='navyblue', size=1, adjust=2) +
  geom_segment(x=2, xend=8,y=0,yend=0, col='navyblue', size=1)+
  labs(x = expression(theta), y = '', title = 'Hierarchical model', color = '') +
  coord_cartesian(xlim=c(-4,8))+
  theme_bw()+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.background = element_blank(), legend.position = 'none')

g1/g2


normal.sep <- mapply(function(estimate, sigma, v) {
  vn <- (1 / sigma^2 + 1 / v^2)^(-1 / 2)
  mn <- (m / v^2 + estimate / sigma^2) * vn^2
  c(mn, vn)
}, data$estimate, data$sigma, data$v)

normal.sep=t(normal.sep)

conf.sep=mapply(function(mn,vn) {
  q10=qnorm(0.1,mn,vn)
  q50=qnorm(0.5,mn,vn)
  q90=qnorm(0.9,mn,vn)
  c(q10,q50,q90)
}, normal.sep[,1], normal.sep[,2] )

conf.sep=as.data.frame(t(conf.sep))
colnames(conf.sep)=c('q10','q50','q90')

conf.hier=apply(theta.out,2, function(x) {
  quantile(x,probs=c(0.1,0.5,0.9))
})

conf.hier=as.data.frame(t(conf.hier))
colnames(conf.hier)=c('q10','q50','q90')
pool.mean=mn.p

plot_sep_int <- conf.sep %>%
  ggplot(aes(x=1:15,y=q50,ymin=q10,ymax=q90)) +
  geom_pointrange(color="blue", alpha=0.3, linewidth=1.2, size=0.6) +
  geom_hline(yintercept=pool.mean, linetype="dashed")+
  labs(x="study number", y="Probability", title="Separate model")+
  theme_bw()+
  annotate("text", x=0, y=pool.mean, label = "Pooled mean",
           vjust=-0.2, hjust=0.3)+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


plot_hier_int <- conf.hier %>%
  ggplot(aes(x=1:15,y=q50,ymin=q10,ymax=q90)) +
  geom_pointrange(color="blue", alpha=0.3, linewidth=1.2,size=0.6) +
  geom_hline(yintercept=pool.mean, linetype="dashed")+
  labs(x="study number", y="Probability", title="Hierarchical model")+
  theme_bw()+
  annotate("text", x=0, y=pool.mean, label = "Pooled mean",
           vjust=-0.2, hjust=0.3)+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


plot_sep_int/plot_hier_int


# RANKING 

theta.rank=t(apply(theta.out, 1, rank)) # 5000 x 15
res=matrix(NA, 15,15)

for (i in 1:15){
  res[,i]=sapply(1:15, function(x) mean(theta.rank[,i]==x))
}

colnames(res)=infarct$trial
res.melt=reshape::melt(res)

colourCount = length(unique(res.melt$X1))
getPalette = colorRampPalette(brewer.pal(9, "Reds"))

ggplot(data=res.melt, aes(x=value, y=factor(X2), fill=factor(X1))) +
  scale_fill_manual(values = rev(getPalette(colourCount)))+
  guides(fill=guide_legend(title="Ranking", nrow=2))+
  geom_bar(stat="identity")+
  theme_minimal()+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.2, hjust=1, size=13),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        legend.background = element_blank(), legend.position = 'bottom',
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  labs(x = 'probability', y ='trials', 
       title = 'Study ranking', color = '')

yusuf.prob=sum(res[1:3,7])
# 73.06



