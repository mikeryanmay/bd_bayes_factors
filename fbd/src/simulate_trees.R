library(TESS)
source("src/simulateGBDP.R")

# simulate 10 trees with 10 taxa under a fossilized birth-death model

############
# settings #
############

time     = 1
shift    = 0.5 * time
num_cats = 2
colors   = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#D55E00", "#0072B2", "#999999")

#######################################
# asynchronous diversification events #
#######################################

# speciation rates
lambda_functions = lapply(1:num_cats, function(x) {
  function(t) {
    ifelse(t > shift, 4, 4)
  }
})

# extinction rates
mu_functions = lapply(1:num_cats, function(x) {
  function(t) {
    2
  }
})

# sampling rates
phi_functions = lapply(1:num_cats, function(x) {
  function(t) {
    ifelse(t > shift, 3, 0.5)
  }
})

# destructive-sampling rates
delta_functions = lapply(1:num_cats, function(x) {
  function(t) {
    1e-16
  }
})


######################################
# synchronous diversification events #
######################################

# mass-speciation events
upsilon_times = c(-1)
upsilon       = matrix(rep(0.5, each=num_cats), nrow=num_cats, ncol=length(upsilon_times))

# mass-extinction events
gamma_times = c(-1)
gamma       = matrix(rep(0.7, each=num_cats), nrow=num_cats, ncol=length(gamma_times))

# mass-sampling events
rho_times = c(0)
rho       = matrix(rep(1.0, each=num_cats), nrow=num_cats, ncol=length(rho_times))

# mass-destructive-sampling events
xi_times = c(-1)
xi       = matrix(rep(0.6, each=num_cats), nrow=num_cats, ncol=length(xi_times))

###################################
# spontaneous state-change events #
###################################

# spontaneous state changes
H = matrix(0.5, num_cats, num_cats) / (num_cats - 1)
diag(H) = 0
diag(H) = -rowSums(H)

##################################################
# diversification-associated state-change events #
##################################################

# speciation-associated state-change
Omega = array(0, dim = c(num_cats, num_cats, num_cats))
for(i in 1:num_cats) {
  Omega[i,,i]  = 0
  Omega[,i,i]  = 0
  Omega[i,i,i] = 1.0
}

# mass-extinction-associated state-change
Z = matrix(0, num_cats, num_cats) / (num_cats - 1)
Z[row(Z) - col(Z) == 1] = 0.4
diag(Z) = 1 - rowSums(Z)


lineages = simulateLHBDFP(lambda_functions,
                          mu_functions,
                          phi_functions,
                          delta_functions,
                          upsilon, upsilon_times,
                          gamma, gamma_times, Z,
                          rho, rho_times,
                          xi, xi_times,
                          H, Omega,
                          2, time)

# par(mar=c(0,0,0,0)+0.1)
# lineages$plot(edge_colors=rep("black", num_cats), sample_colors=rep("black", num_cats),
#               pch=19, cex=0.7, lwd=2, lend=2)

# simulate the trees
stratigraphic_intervals = seq(0, time, length.out=21)

reps = 10
for(i in 1:reps) {

  cat("simulation ", i, "\n", sep="")
    
  # simulate a tree with 10 extants and 5 extincts
  repeat {
    
    # simulate the tree
    lineages = simulateLHBDFP(lambda_functions,
                              mu_functions,
                              phi_functions,
                              delta_functions,
                              upsilon, upsilon_times,
                              gamma, gamma_times, Z,
                              rho, rho_times,
                              xi, xi_times,
                              H, Omega,
                              2, time)
    
    # drop unsampled lineages
    lineages$dropUnsampled()
    
    # write to file (also construct the newick string)
    dir = paste0("data/dataset_", i)
    dir.create(dir, showWarnings = FALSE, recursive=TRUE)
    lineages$writeNewickString( paste0(dir,"/tree.nex"))

    # plot the tree
    pdf(paste0(dir,"/tree.pdf"), height=4, width=4)
    par(mar=c(0,0,0,0)+0.1)
    lineages$plot(edge_colors=rep("black", num_cats), sample_colors=rep("black", num_cats),
                  pch=19, cex=0.7, lwd=2, lend=2)
    dev.off()
    
    # check the number of lineages    
    tree        = read.tree(text=lineages$newick_string)
    num_extant  = sum(lineages$lineages$end_time == 0.0)
    num_extinct = length(tree$tip.label) - num_extant
    
    cat(i, "\t", num_extant,"\t", num_extinct, "\n")
    if ( num_extant == 9 & num_extinct == 3 ) {
      break
    }
    
  }
  
  # write the taxon data
  taxon_file = paste0(dir, "/taxa.tsv")
  cat("taxon\tmin\tmax\n", sep="", file = taxon_file)
  
  sampled_lineages = lineages$lineages[lineages$lineages$status == "sampled",]
  for(j in 1:nrow(sampled_lineages)) {
    
    this_sample = sampled_lineages[j,]
    taxon_label = paste0("t_", this_sample$desc)
    
    if ( this_sample$end_time == 0.0 ) {
      start_time = 0.0
      end_time   = 0.0
    } else {
      start_time = stratigraphic_intervals[findInterval(this_sample$end_time, stratigraphic_intervals)]
      end_time   = stratigraphic_intervals[findInterval(this_sample$end_time, stratigraphic_intervals) + 1]
    }
    
    cat(taxon_label, "\t", start_time, "\t", end_time, "\n", sep="", file=taxon_file, append=TRUE)
    
  }
  
}

# for(n in ntaxa) {
# 
#     # simulate the tree
#     tree = tess.sim.taxa.age(1, n, age, lambda - mu, mu=0, samplingProbability=rho)[[1]]
# 
# 	  for(i in 1:reps) {
# 
#       for(seq in seq_model) {
# 	    
# 	        # make the directories
# 	        dir  = paste0("data/yule/",seq,"/n_",n,"/dataset_", i)
# 	        dir.create(dir, showWarnings = FALSE, recursive=TRUE)
# 	        
# 	        # write the tree
# 	        write.nexus(tree, file=paste0(dir, "/tree.nex"))
# 	        
#     }
#     
#   }
# 
# }

