library(ape)
library(parallel)
library(geoscale)
library(Rcpp)
source("src/geoscale_axis.R")
source("src/posterior_boxplot.R")
source("src/simulate_epochs.R")
sourceCpp("src/simulate_epochs.cpp")

# read the epochs
timescale = timescales$ICS2015
epochs = read.table("data/epochs.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
epoch_times = epochs$start

new_epochs = timescale[timescale$Start <= max(epoch_times) & timescale$Type == "Epoch",]
new_epoch_times = new_epochs$Start

epoch_map = c(findInterval(new_epoch_times, epoch_times, left.open=TRUE) + 1,5)
epoch_map = rev(abs(epoch_map - 6))

# read the data
obs_data = read.table("data/ingroup/taxa.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")

# compute the statistic
obs_nums = calculateObservedEpochNums(obs_data, epochs$start)

#######################
# constant-rate model #
#######################

constant_rate_samples = read.table("output/tree_model_constant_rate_MCMC_run_1/params.log", header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=FALSE)

constant_speciation_rates = constant_rate_samples[,paste0("speciation_rate[",epoch_map,"]")]
colnames(constant_speciation_rates) = paste0("speciation_rate[",1:length(epoch_map),"]")

constant_extinction_rates = constant_rate_samples[,paste0("extinction_rate[",epoch_map,"]")]
colnames(constant_extinction_rates) = paste0("extinction_rate[",1:length(epoch_map),"]")

constant_fossilization_rates = constant_rate_samples[,paste0("fossilization_rate[",epoch_map,"]")]
colnames(constant_fossilization_rates) = paste0("fossilization_rate[",1:length(epoch_map),"]")

constant_rate_samples = cbind(constant_speciation_rates, constant_extinction_rates, constant_fossilization_rates, stem_age=constant_rate_samples$stem_age)
constant_epoch_rates   = makeEpochRates(constant_rate_samples, new_epoch_times)

# constant_epoch_rates   = makeEpochRates(constant_rate_samples, epoch_times)
constant_lambda_epochs = constant_epoch_rates$speciation_rate
constant_mu_epochs     = constant_epoch_rates$extinction_rate
constant_phi_epochs    = constant_epoch_rates$fossilization_rate
constant_age           = constant_rate_samples$stem_age
constant_rho           = 27 / 111

# simulate the events
nsim = 50000
nsamples = nrow(constant_rate_samples)
pts = round(seq(1, nsim, length.out=100))

constant_sims = mclapply(1:nsim, function(i) {

  # get the sample
  sample_idx = sample.int(nsamples, size=1)
  sim = simulateEpochsConditional(constant_age[sample_idx], new_epoch_times, constant_lambda_epochs[[sample_idx]], constant_mu_epochs[[sample_idx]], constant_phi_epochs[[sample_idx]], constant_rho)

  while ( is.null(sim) ) {
    cat(".", sep="")
    sample_idx = sample.int(nsamples, size=1)
    sim = simulateEpochsConditional(constant_age[sample_idx], new_epoch_times, constant_lambda_epochs[[sample_idx]], constant_mu_epochs[[sample_idx]], constant_phi_epochs[[sample_idx]], constant_rho)
  }

  res = data.frame( samples_per_epoch = I(list(sim$num_taxa_per_epoch)),
                    sampled_ltt       = I(list(epochLTT(sim, 420, num_bins = 1001))))

  if (i %% 100 == 0)
    cat(i, "\n")
  
  return(res)
  
}, mc.cores=6, mc.preschedule=FALSE)

# constant_sims = vector("list", nsim)
# bar = txtProgressBar(style=3, width=40)
# for(i in 1:nsim) {
#   
#   # get the sample
#   sample_idx = sample.int(nsamples, size=1)
#   sim = simulateEpochsConditional(constant_age[sample_idx], new_epoch_times, constant_lambda_epochs[[sample_idx]], constant_mu_epochs[[sample_idx]], constant_phi_epochs[[sample_idx]], constant_rho)
#   
#   while ( is.null(sim) ) {
#     cat(".", sep="")
#     sample_idx = sample.int(nsamples, size=1)
#     sim = simulateEpochsConditional(constant_age[sample_idx], new_epoch_times, constant_lambda_epochs[[sample_idx]], constant_mu_epochs[[sample_idx]], constant_phi_epochs[[sample_idx]], constant_rho)
#   }
#   
#   res = data.frame( samples_per_epoch = I(list(sim$num_taxa_per_epoch)),
#                     sampled_ltt       = I(list(epochLTT(sim, 420, num_bins = 1001))))
#   
#   constant_sims[[i]] = res
#   
#   setTxtProgressBar(bar, i / nsim)
#   
# }

constant_samples_per_epoch  = lapply(constant_sims, function(x) x$samples_per_epoch[[1]] )
constant_sampled_ltt        = lapply(constant_sims, function(x) x$sampled_ltt[[1]] )
constant_sampled_ltt        = do.call(rbind, constant_sampled_ltt)
constant_sampled_ltt_median = apply(constant_sampled_ltt, 2, function(x) median(x,na.rm=TRUE) )
constant_sampled_ltt_mean   = colMeans(constant_sampled_ltt, na.rm=TRUE)
constant_sampled_ltt_CI_95  = apply(constant_sampled_ltt, 2, quantile, prob=c(0.025, 0.975), na.rm=TRUE)
constant_sampled_ltt_CI_50  = apply(constant_sampled_ltt, 2, quantile, prob=c(0.25, 0.75), na.rm=TRUE)
constant_sampled_ltt_times  = seq(0, 420, length.out = length(constant_sampled_ltt_median))

#######################
# variable-rate model #
#######################

variable_rate_samples  = read.table("output/tree_model_phi_variable_MCMC_run_1/params.log", header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=FALSE)

variable_speciation_rates = variable_rate_samples[,paste0("speciation_rate[",epoch_map,"]")]
colnames(variable_speciation_rates) = paste0("speciation_rate[",1:length(epoch_map),"]")

variable_extinction_rates = variable_rate_samples[,paste0("extinction_rate[",epoch_map,"]")]
colnames(variable_extinction_rates) = paste0("extinction_rate[",1:length(epoch_map),"]")

variable_fossilization_rates = variable_rate_samples[,paste0("fossilization_rate[",epoch_map,"]")]
colnames(variable_fossilization_rates) = paste0("fossilization_rate[",1:length(epoch_map),"]")

variable_rate_samples = cbind(variable_speciation_rates, variable_extinction_rates, variable_fossilization_rates, stem_age=variable_rate_samples$stem_age)

variable_epoch_rates   = makeEpochRates(variable_rate_samples, new_epoch_times)

# plot(colMeans(variable_rate_samples[,paste0("fossilization_rate[",1:5,"]")]))
# plot(colMeans(variable_rate_samples[,paste0("fossilization_rate[",1:5,"]")])[epoch_map])
# plot(colMeans(variable_fossilization_rates))
# plot(colMeans(do.call(rbind, variable_phi_epochs)))

# variable_epoch_rates   = makeEpochRates(variable_rate_samples, epoch_times)
variable_lambda_epochs = variable_epoch_rates$speciation_rate
variable_mu_epochs     = variable_epoch_rates$extinction_rate
variable_phi_epochs    = variable_epoch_rates$fossilization_rate
variable_age           = variable_rate_samples$stem_age
variable_rho           = 27 / 111

# simulate the events
nsamples = nrow(variable_rate_samples)
pts = round(seq(1, nsim, length.out=100))

variable_sims = mclapply(1:nsim, function(i) {
  
  # get the sample
  sample_idx = sample.int(nsamples, size=1)
  sim = simulateEpochsConditional(variable_age[sample_idx], new_epoch_times, variable_lambda_epochs[[sample_idx]], variable_mu_epochs[[sample_idx]], variable_phi_epochs[[sample_idx]], variable_rho)
  
  while ( is.null(sim) ) {
    cat(".", sep="")
    sample_idx = sample.int(nsamples, size=1)
    sim = simulateEpochsConditional(variable_age[sample_idx], new_epoch_times, variable_lambda_epochs[[sample_idx]], variable_mu_epochs[[sample_idx]], variable_phi_epochs[[sample_idx]], variable_rho)
  }
  
  res = data.frame( samples_per_epoch = I(list(sim$num_taxa_per_epoch)),
                    sampled_ltt       = I(list(epochLTT(sim, 420, num_bins = 1001))))
  
  if (i %% 100 == 0)
    cat(i, "\n")
  
  return(res)
  
}, mc.cores=6, mc.preschedule=FALSE)

# variable_sims = vector("list", nsim)
# bar = txtProgressBar(style=3, width=40)
# for(i in 1:nsim) {
#   
#   # get the sample
#   sample_idx = sample.int(nsamples, size=1)
#   sim = simulateEpochsConditional(variable_age[sample_idx], new_epoch_times, variable_lambda_epochs[[sample_idx]], variable_mu_epochs[[sample_idx]], variable_phi_epochs[[sample_idx]], variable_rho)
#   
#   while ( is.null(sim) ) {
#     cat(".", sep="")
#     sample_idx = sample.int(nsamples, size=1)
#     sim = simulateEpochsConditional(variable_age[sample_idx], new_epoch_times, variable_lambda_epochs[[sample_idx]], variable_mu_epochs[[sample_idx]], variable_phi_epochs[[sample_idx]], variable_rho)
#   }
#   
#   res = data.frame(samples_per_epoch = I(list(sim$num_taxa_per_epoch)),
#                    sampled_ltt       = I(list(epochLTT(sim, 420, num_bins = 1001))))
#   
#   variable_sims[[i]] = res
#   
#   setTxtProgressBar(bar, i / nsim)
#   
# }

variable_samples_per_epoch  = lapply(variable_sims, function(x) x$samples_per_epoch[[1]] )
variable_sampled_ltt        = lapply(variable_sims, function(x) x$sampled_ltt[[1]] )
variable_sampled_ltt        = do.call(rbind, variable_sampled_ltt)
variable_sampled_ltt_median = apply(variable_sampled_ltt, 2, function(x) median(x,na.rm=TRUE) )
variable_sampled_ltt_mean   = colMeans(variable_sampled_ltt, na.rm=TRUE)
variable_sampled_ltt_CI_95  = apply(variable_sampled_ltt, 2, quantile, prob=c(0.025, 0.975), na.rm=TRUE)
variable_sampled_ltt_CI_50  = apply(variable_sampled_ltt, 2, quantile, prob=c(0.25, 0.75), na.rm=TRUE)
variable_sampled_ltt_times  = seq(0, 420, length.out = length(variable_sampled_ltt_median))

# plot the curves
cols = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#D55E00", "#0072B2", "#999999")[c(2,1)]

# create the epochs
timescale = timescales$ICS2015

timescale$Name = as.character(timescale$Name)
timescale$Name[is.na(timescale$Name)] = ""

timescale$Abbrev = as.character(timescale$Abbrev)
timescale$Abbrev[is.na(timescale$Abbrev)] = ""

# correct the end Permian age
timescale$Start[timescale$Start == 252.170] = 251.9
timescale$End[timescale$End == 252.170] = 251.9
timescale[timescale$Name == "Paleogene",]$Name = "Pg"
epochs = timescale[timescale$Type == "Epoch",]


pdf("figures/num_lineages_mean.pdf", height=3, width=8)
par(mar=c(4,5,0,5)+0.1, lend=2)
plot(constant_sampled_ltt_times, y=constant_sampled_ltt_mean, xlim=c(420, 0), ylim=c(0,1200), type="n", log="", col="orange", xaxt="n", yaxt="n", xlab=NA, ylab=NA, bty="n")
for(i in nrow(epochs):1) {
  if ( i %% 2 == 1 ) {
    this_epoch = epochs[i,]
    polygon( x=c(this_epoch[1:2], this_epoch[2:1]), y=c(-10000, -10000, 10000, 10000), border=NA, col="grey90")
  }
}
abline(v=seq(0, 1000, by=50), lty=2)
abline(h=111, lty=3)
lines(constant_sampled_ltt_times, y=constant_sampled_ltt_mean, type="l", col=cols[1], lwd=2)
lines(variable_sampled_ltt_times, y=variable_sampled_ltt_mean, type="l", col=cols[2], lwd=2)
geoscaleAxis(timescale, fraction=0.95, las.axis=2, cex.lab=0.7, outer=FALSE, lwd.axis=0)
legend("topright", legend=c("constant", "variable"), lty=1, lwd=2, col=cols, bg="white", inset=c(0.05,0.05))
axis(2, lwd=1, lwd.tick=1, las=2)
mtext("expected number of lineages", side=2, line=3.5)
dev.off()

pdf("figures/num_lineages_median.pdf", height=3, width=8)
par(mar=c(4,5,0,5)+0.1, lend=2)
plot(constant_sampled_ltt_times, y=constant_sampled_ltt_mean, xlim=c(420, 0), ylim=c(0,500), type="n", log="", col="orange", xaxt="n", yaxt="n", xlab=NA, ylab=NA, bty="n")
for(i in nrow(epochs):1) {
  if ( i %% 2 == 1 ) {
    this_epoch = epochs[i,]
    polygon( x=c(this_epoch[1:2], this_epoch[2:1]), y=c(-10000, -10000, 10000, 10000), border=NA, col="grey90")
  }
}
abline(v=seq(0, 1000, by=50), lty=2)
abline(h=111, lty=3)
lines(constant_sampled_ltt_times, y=constant_sampled_ltt_median, type="l", col=cols[1], lwd=2)
lines(variable_sampled_ltt_times, y=variable_sampled_ltt_median, type="l", col=cols[2], lwd=2)
geoscaleAxis(timescale, fraction=0.95, las.axis=2, cex.lab=0.7, outer=FALSE, lwd.axis=0)
legend("topright", legend=c("constant", "variable"), lty=1, lwd=2, col=cols, bg="white", inset=c(0.05,0.05))
axis(2, lwd=1, lwd.tick=1, las=2)
mtext("expected number of lineages", side=2, line=3.5)
dev.off()


cat("\nMean: ")
cat(max(constant_sampled_ltt_mean), "\t", max(variable_sampled_ltt_mean), sep="")

cat("\nMedian: ")
cat(max(constant_sampled_ltt_median), "\t", max(variable_sampled_ltt_median), sep="")

cat("\nCI constant: \n")
cat(constant_sampled_ltt_CI_95[,which.max(constant_sampled_ltt_mean)])

cat("\nCI variable: \n")
cat(variable_sampled_ltt_CI_95[,which.max(variable_sampled_ltt_mean)])


# PPS
epoch_colors = rgb(
  epochs$Col_R,
  epochs$Col_G,
  epochs$Col_B,
  maxColorValue = 255
)


obs_nums = calculateObservedEpochNums(obs_data, new_epoch_times)
obs_nums = obs_nums[-length(obs_nums)]

constant_pps = do.call(rbind, constant_samples_per_epoch)
variable_pps = do.call(rbind, variable_samples_per_epoch)

constant_pps = lapply(1:ncol(constant_pps), function(x) constant_pps[,x])
variable_pps = lapply(1:ncol(variable_pps), function(x) variable_pps[,x])

sub_epochs = epochs[1:length(obs_nums),]
epoch_labels = character(nrow(sub_epochs))
for(i in 1:length(epoch_labels)) {
  this_epoch = sub_epochs$Name[i]
  if ( this_epoch %in% c("Upper","Lower","Middle") ) {
    this_epoch = paste0(this_epoch, " ", sub_epochs$Part_of[i] )
  }
  epoch_labels[i] = this_epoch
}


pdf("figures/sample_pps.pdf", height=5, width=8)
par(mar=c(0,0,0,0), oma=c(10,4,2,4)+0.1, mfrow=1:2)
boxplot.HPD(constant_pps, outline=FALSE, col=NA, border=NA, frame=FALSE, xaxt="n", yaxt="n", xlim=c(20,1), ylim=c(0,150))
polygon(x=c(-100,100,100,-100), y=c(-1000,-1000,1000,1000), col="grey90")
boxplot.HPD(constant_pps, outline=FALSE, col=epoch_colors, border=epoch_colors, add=TRUE, las=2, names=epoch_labels, lty=1, ylim=c(0,150))
points(obs_nums)
mtext("number of samples", side=2, line=2.5)
mtext("epoch", side=1, line=8.5)
mtext("constant model", side=3, line=0.5)

boxplot.HPD(variable_pps, outline=FALSE, col=NA, border=NA, frame=FALSE, xaxt="n", yaxt="n", xlim=c(20,1), ylim=c(0,150))
polygon(x=c(-100,100,100,-100), y=c(-1000,-1000,1000,1000), col="grey90")
boxplot.HPD(variable_pps, outline=FALSE, col=epoch_colors, border=epoch_colors, add=TRUE, yaxt="n", las=2, names=epoch_labels, lty=1, ylim=c(0,150))
points(obs_nums)
mtext("epoch", side=1, line=8.5)
mtext("variable model", side=3, line=0.5)
dev.off()


# histograms

quant = 0.9

pdf("figures/sample_pps_hist.pdf", height=2.5, width=6)
par(mar=c(0,0,0,0), oma=c(4,4,2,4) + 0.1, mfrow=c(1,3), lend=2)

this_const = constant_pps[epoch_labels == "Mississippian"][[1]]
this_var   = variable_pps[epoch_labels == "Mississippian"][[1]]
this_obs   = obs_nums[epoch_labels == "Mississippian"][[1]]

xlim = c(0, pmax(quantile(this_const, prob=quant),quantile(this_var, prob=quant)))

const_dens = density(this_const, from=0, to=xlim[2])
var_dens   = density(this_var, from=0, to=xlim[2])
ylim       = c(0, max(c(const_dens$y, var_dens$y)))

plot(const_dens, type="s", xaxt="n", col=cols[1], zero.line=FALSE, xlim=xlim, xlab=NA, ylab=NA, xaxt="n", yaxt="n", main=NA, ylim=ylim)
lines(var_dens, type="s", xaxt="n", col=cols[2])
abline(v=this_obs, lty=2)
legend("topright", legend=c("constant", "variable"), lty=1, lwd=1, col=cols, bg="white", inset=c(0.05,0.05))
abline(v=median(this_const), lty=2, col=cols[1])
abline(v=median(this_var), lty=2, col=cols[2])
axis(1, lwd=0, lwd.tick=1, las=2)
mtext("predictive density", side=2, line=1.5)
mtext("number of samples", side=1, line=2.7)
mtext("Mississippian", side=3, line=0.5)

this_const = constant_pps[epoch_labels == "Pennsylvanian"][[1]]
this_var   = variable_pps[epoch_labels == "Pennsylvanian"][[1]]
this_obs   = obs_nums[epoch_labels == "Pennsylvanian"][[1]]

xlim = c(0, pmax(quantile(this_const, prob=quant),quantile(this_var, prob=quant)))

const_dens = density(this_const, from=0, to=xlim[2])
var_dens   = density(this_var, from=0, to=xlim[2])
ylim       = c(0, max(c(const_dens$y,var_dens$y)))

plot(const_dens, type="s", xaxt="n", col=cols[1], zero.line=FALSE, xlim=xlim, xlab=NA, ylab=NA, xaxt="n", yaxt="n", main=NA, ylim=ylim)
lines(var_dens, type="s", xaxt="n", col=cols[2])
abline(v=this_obs, lty=2)
abline(v=median(this_const), lty=2, col=cols[1])
abline(v=median(this_var), lty=2, col=cols[2])
axis(1, lwd=0, lwd.tick=1, las=2)
mtext("Pennsylvanian", side=3, line=0.5)
mtext("number of samples", side=1, line=2.7)


this_const = constant_pps[epoch_labels == "Cisuralian"][[1]]
this_var   = variable_pps[epoch_labels == "Cisuralian"][[1]]
this_obs   = obs_nums[epoch_labels == "Cisuralian"][[1]]

xlim = c(0, pmax(quantile(this_const, prob=quant),quantile(this_var, prob=quant)))

const_dens = density(this_const, from=0, to=xlim[2])
var_dens   = density(this_var, from=0, to=xlim[2])
ylim       = c(0, max(c(const_dens$y,var_dens$y)))

plot(const_dens, type="s", xaxt="n", col=cols[1], zero.line=FALSE, xlim=xlim, xlab=NA, ylab=NA, xaxt="n", yaxt="n", main=NA, ylim=ylim)
lines(var_dens, type="s", xaxt="n", col=cols[2])
abline(v=this_obs, lty=2)
abline(v=median(this_const), lty=2, col=cols[1])
abline(v=median(this_var), lty=2, col=cols[2])
axis(1, lwd=0, lwd.tick=1, las=2)
mtext("Cisuralian", side=3, line=0.5)
mtext("number of samples", side=1, line=2.7)

dev.off()
















