library(TESS)

# simulate 10 trees with n taxa under a birth-death model
lambda    = 4
mu        = 2
age       = 1
rho       = 1.0
ntaxa     = c(8, 16, 32, 64)
# size      = c(0, 1, 10, 100, Inf)
reps      = 10
seq_model = c("JC")

# bd trees
for(n in ntaxa) {
    
  	for(i in 1:reps) {
  	
  		for(seq in seq_model) {
  	
    		  # simulate the tree
    		  tree = tess.sim.taxa.age(1, n, age, lambda, mu, samplingProbability=rho)[[1]]
  		  
  		    # make the directories
  		    dir  = paste0("data/bd/", seq, "/n_", n,"/dataset_", i)
  		    dir.create(dir, showWarnings = FALSE, recursive=TRUE)
  		    
  		    # write the tree
  		    write.nexus(tree, file=paste0(dir, "/tree.nex"))
  		    
      }
      
    }
  
}
