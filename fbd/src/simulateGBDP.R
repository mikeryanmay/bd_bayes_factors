library(ape)
library(stringr)

heterogeneousBDFPTree = setRefClass(

    Class = "heterogeneousBDFPTree",

    fields = list(

        "lineages" = "data.frame",
        "num_states" = "numeric",

        "has_coords" = "logical",

        "horizontal_segments" = "matrix",
        "vertical_segments"   = "matrix",
        "sampled_points"      = "matrix",

        "lineage_counter" = "numeric",
        "splits_counter"  = "numeric",
        "sample_counter"  = "numeric",

        "mark_drop" = "logical",
        "num_descs" = "numeric",

        "phy_num_ints" = "numeric",
        "phy_num_tips" = "numeric",
        "phy_edge_mat" = "matrix",
        "phy_edge_ln"  = "numeric",

        "newick_string" = "character",
        "tip_labels"    = "character",
        "tip_states"    = "character"

    ),

    methods = list(

        initialize = function(x = NULL, ns = NULL) {

            if ( is.null(x) == FALSE ) {
                lineages <<- x
                num_states <<- ns
                sort()
            }

            has_coords <<- FALSE

        },

        save = function(file) {
            lin = .self
            base::save(lin, file = file)
        },
        
        sort = function() {

            # for each lineage, compute the number of descendants
            # at a split, always put the species with the most descendants
            # on the right

            num_descs <<- numeric( nrow(lineages) )

            # start the recursive algorithm at the root
            init_nodes = lineages$desc[lineages$anc == lineages[1,1]]
            for(node in init_nodes) {
                computeNodeNumDescs(node)
            }

            # now sort
            for(i in 1:nrow(lineages)) {

                # if this is a split
                if ( lineages$status[i] == "split" ) {

                    # get the descendants of this split
                    this_anc = lineages$desc[i]
                    desc_idx = which(lineages$anc == this_anc)

                    # get the number on the left and right
                    left_num  = num_descs[ desc_idx[1] ]
                    right_num = num_descs[ desc_idx[2] ]

                    # put the one with the most descendants on the right
                    if ( left_num > right_num ) {
                        lineages[ desc_idx, ] <<- lineages[ rev(desc_idx), ]
                    }

                }

            }

        },

        writeNewickString = function(file) {

            # sort the tree
            sort()

            # initialize the strings
            newick_string <<- ""
            tip_labels    <<- character()
            tip_states    <<- character()

            # start the recursive algorithm at the root
            init_nodes = lineages[1,2]
            for(node in init_nodes) {
                nodeNewickString(node, 0)
            }

            # finalize the string
            newick_string <<- str_c(newick_string,";")

            # create the file
            cat("#NEXUS\n\n", sep="", file=file)

            # create the taxon block
            cat("Begin taxa;\n", sep="", file=file, append=TRUE)
            cat("Dimensions ntax=", length(tip_labels), ";\n", sep="", file=file, append=TRUE)
            cat("Taxlabels\n", sep="", file=file, append=TRUE)
            cat(paste0("\t", tip_labels), sep="\n", file=file, append=TRUE)
            cat(";\n", sep="", file=file, append=TRUE)
            cat("End;\n\n", sep="", file=file, append=TRUE)

            # create the data block
            cat("Begin characters;\n", sep="", file=file, append=TRUE)
            cat("Dimensions nchar=1;\n", sep="", file=file, append=TRUE)
            cat("Format datatype=Standard gap=- missing=? symbols=\"", 1:num_cats ,"\";\n", sep="", file=file, append=TRUE)
            cat("Matrix\n", sep="", file=file, append=TRUE)
            max_length = max(str_length(tip_labels))
            tip_label_pad = str_pad(tip_labels, width=max_length, side="right")
            cat(paste0("\t",tip_label_pad, "\t", tip_states), sep="\n", file=file, append=TRUE)
            cat(";\n", sep="", file=file, append=TRUE)
            cat("End;\n\n", sep="", file=file, append=TRUE)

            # write the string to path
            cat("Begin trees;\n", sep="", file=file, append=TRUE)
            cat("\ttree TREE1 = [&R]", newick_string, "\n", sep="", file=file, append=TRUE)
            cat("End;\n\n", sep="", file=file, append=TRUE)

            # p = read.tree(text=newick_string)
            # plot.phylo(p)
            #
            # strsplit(newick_string, "(", fixed=TRUE)
            # strsplit(newick_string, ")", fixed=TRUE)

        },

        nodeNewickString = function(node, accum_time) {

            # get row index
            row_idx = which(lineages$desc == node)

            # get descendant index
            desc_idx = which(lineages$anc == lineages$desc[row_idx])

            # get the status
            status = lineages$status[row_idx]

            # if this is a split
            if ( status == "split" ) {

                # add open parens
                newick_string <<- str_c(newick_string, "(")

                # call on left descendant
                nodeNewickString( lineages$desc[desc_idx[1]], 0 )

                # add comma
                newick_string <<- str_c(newick_string, ",")

                # call on right descendant
                nodeNewickString( lineages$desc[desc_idx[2]], 0 )

                # close parens
                newick_string <<- str_c(newick_string, ")")

                # add branch length
                newick_string <<- str_c(newick_string, ":", round(accum_time + lineages$start_time[row_idx] - lineages$end_time[row_idx], 16) )

            } else if ( status == "unsampled" | status == "extinct" | (status == "sampled" & length(desc_idx) == 0) ) {

                # a terminal

                # make the string
                label = paste0("t_",lineages$desc[row_idx])
                state = lineages$state[row_idx]
                bl    = round(accum_time + lineages$start_time[row_idx] - lineages$end_time[row_idx], 16)
                newick_string <<- str_c(newick_string, label, "[&state=",state,"]:", bl)

                tip_labels <<- c(tip_labels, label)
                tip_states <<- c(tip_states, state)

            } else if ( status == "sampled" & length(desc_idx) == 1 ) {

                # a sampled ancestor

                # add open parens
                newick_string <<- str_c(newick_string, "(")

                # accumulate time, pass to descendant
                accum_time = accum_time + lineages$start_time[row_idx] - lineages$end_time[row_idx]
                nodeNewickString( lineages$desc[desc_idx[1]], 0)

                # close parens
                newick_string <<- str_c(newick_string, ")")

                # make the string
                label = paste0("t_",lineages$desc[row_idx])
                state = lineages$state[row_idx]
                bl    = round(accum_time, 16)
                newick_string <<- str_c(newick_string, label, "[&state=",state,"]:", bl)

                tip_labels <<- c(tip_labels, label)
                tip_states <<- c(tip_states, state)

            } else if ( status == "switched" | status == "knuckle" ) {

                # accumulate time, pass to descendant
                accum_time = accum_time + lineages$start_time[row_idx] - lineages$end_time[row_idx]
                nodeNewickString( lineages$desc[desc_idx[1]], accum_time )

            } else {
                stop("How did you get here???")
            }

        },

        computeNodeNumDescs = function(node) {

            # get the row for this node
            row_idx = which(lineages$desc == node)

            # get all the descendants of this node
            desc_idx = which(lineages$anc == node)
            num_desc = length(desc_idx)

            # determine if this node is a terminal
            if ( num_desc == 0 ) {

                # determine if to count this node
                if ( lineages$status[row_idx] == "sampled" | lineages$status[row_idx] == "extinct" ) {
                    num_descs[row_idx] <<- 1
                }

            } else {

                # call this function on each descendant
                desc_nodes = lineages$desc[desc_idx]
                for(desc_node in desc_nodes) {
                    computeNodeNumDescs( desc_node )
                }

                # compute the number of this branch
                num_descs[row_idx] <<- sum(num_descs[desc_idx])

            }

        },

        plot = function(edge_colors, sample_colors, new=TRUE, ...) {

            # get the coordinates
            computeCoordinates()

            # compute the ranges
            age = max(lineages$start_time)

            # make an empty plot
            if ( new == TRUE ) {
                plot.new()
                plot.window(xlim = c(age, 0),
                            ylim = c(1, lineage_counter), ...)
            }

            # plot the horizontal segments
            segments(
                y0  = horizontal_segments[,1],
                x0  = horizontal_segments[,2],
                x1  = horizontal_segments[,3],
                col = edge_colors[horizontal_segments[,4]],
                ...
            )

            # plot the vertical segments
            segments(
                x0  = vertical_segments[,1],
                y0  = vertical_segments[,2],
                y1  = vertical_segments[,3],
                col = edge_colors[vertical_segments[,4]],
                ...
            )

            # plot the samples
            points(
                x   = sampled_points[,1],
                y   = sampled_points[,2],
                col = sample_colors[sampled_points[,3]],
                ...
            )

        },

        computeCoordinates = function() {

            if ( has_coords == FALSE ) {

                # compute the coordinates

                # create the matrix
                num_lineages = nrow(lineages)
                horizontal_segments  <<- matrix(NA, nrow=num_lineages, ncol=4)
                vertical_segments    <<- matrix(NA, nrow=0, ncol=4)
                sampled_points       <<- matrix(NA, nrow=0, ncol=3)

                lineage_counter <<- 0
                splits_counter  <<- 0
                sample_counter  <<- 0

                # start the recursive algorithm at the root
                init_nodes = lineages$desc[lineages$anc == lineages[1,1]]
                for(node in init_nodes) {
                    computeNodeCoordinates(node)
                }

                # compute the vertical segments
                has_coords <<- TRUE

            }

            return()

        },

        computeNodeCoordinates = function(node) {

            # get the row for this node
            row_idx = which(lineages$desc == node)

            # get all the descendants of this node
            desc_idx = which(lineages$anc == node)
            num_desc = length(desc_idx)

            # get the state
            node_state = lineages$state[row_idx]

            # determine if this node is a terminal
            if ( num_desc == 0 ) {

                # it's a terminal!
                lineage_counter <<- lineage_counter + 1
                horizontal_segments[row_idx,] <<- c(lineage_counter, lineages$end_time[row_idx], lineages$start_time[row_idx], node_state)

            } else {

                # call this function on each descendant
                desc_nodes = lineages$desc[desc_idx]
                for(desc_node in desc_nodes) {
                    computeNodeCoordinates( desc_node )
                }

                # now compute the coordinates for this node

                # the y-coordinate is the average x-coordinate of the descendants
                horizontal_segments[row_idx,1] <<- mean(horizontal_segments[desc_idx,1])

                # now the x-coordinates
                horizontal_segments[row_idx,2] <<- lineages$end_time[row_idx]
                horizontal_segments[row_idx,3] <<- lineages$start_time[row_idx]
                horizontal_segments[row_idx,4] <<- node_state

            }

            # if this node is a split, add a new vertical segment
            if ( lineages$status[row_idx] == "split" ) {

                # increment the counter
                splits_counter <<- splits_counter + 1

                # compute the coordinates
                x_coord  = horizontal_segments[row_idx,2]
                y_coords = horizontal_segments[desc_idx,1]

                # bind the new coordinates
                middle_coord = mean(y_coords)
                desc_states  = lineages$state[desc_idx]
                vertical_segments <<- rbind(vertical_segments,
                                            c(x_coord, middle_coord, y_coords[1], desc_states[1]),
                                            c(x_coord, middle_coord, y_coords[2], desc_states[2]))

            } else if ( lineages$status[row_idx] == "sampled" ) {

                # increment the counter
                sample_counter <<- sample_counter + 1

                # compute the coordinates
                x_coord = horizontal_segments[row_idx,2]
                y_coord = horizontal_segments[row_idx,1]

                # bind the new coordinates
                sampled_points <<- rbind(sampled_points, c(x_coord, y_coord, node_state))

            }

        },

        dropUnsampled = function() {

            has_coords <<- FALSE

            # visit each tip and mark it for exclusion if it does not have
            # sampled descendants
            mark_drop <<- logical(nrow(lineages))

            # visit each node
            init_nodes = lineages$desc[lineages$anc == 1]
            for(node in init_nodes) {
                computeInclude(node)
            }

            # now drop the nodes that shouldn't be included
            lineages <<- lineages[mark_drop == FALSE,]

            # re-sort
            sort()

        },

        computeInclude = function(node) {

            # get the row for this node
            row_idx = which(lineages$desc == node)

            # get all the descendants of this node
            desc_idx = which(lineages$anc == node)
            num_desc = length(desc_idx)

            # determine if this node is a terminal
            if ( num_desc == 0 ) {

                # determine if to include this
                mark_drop[row_idx] <<- lineages$status[row_idx] == "unsampled" | lineages$status[row_idx] == "extinct"

            } else {

                # call this function on each descendant
                desc_nodes = lineages$desc[desc_idx]
                for(desc_node in desc_nodes) {
                    computeInclude( desc_node )
                }

                if ( (lineages$status[row_idx] == "sampled") == TRUE ) {
                    mark_drop[row_idx] <<- FALSE
                } else {
                    # check if any descendants are included
                    mark_drop[row_idx] <<- all(mark_drop[desc_idx] == TRUE)
                }

                # if this is a split, and only has one sampled descendant,
                # mark it as a knuckle
                if ( lineages$status[row_idx] == "split" & sum(mark_drop[desc_idx]) == 1 ) {
                    lineages$status[row_idx] <<- "knuckle"
                }

            }

        },

        dropFossils = function() {

            has_coords <<- FALSE

            # visit each tip and mark it for exclusion if it does not have
            # sampled descendants
            mark_drop <<- logical(nrow(lineages))

            # visit each node
            init_nodes = lineages$desc[lineages$anc == 1]
            for(node in init_nodes) {
                computeIncludeFossils(node)
            }

            # now drop the nodes that shouldn't be included
            lineages <<- lineages[mark_drop == FALSE,]

            # re-sort
            sort()

        },

        computeIncludeFossils = function(node) {

            # get the row for this node
            row_idx = which(lineages$desc == node)

            # get all the descendants of this node
            desc_idx = which(lineages$anc == node)
            num_desc = length(desc_idx)

            # determine if this node is a terminal
            if ( num_desc == 0 ) {

                # determine if to include this
                mark_drop[row_idx] <<- lineages$status[row_idx] == "unsampled" | lineages$status[row_idx] == "extinct" | lineages$end_time[row_idx] != 0.0

            } else {

                # call this function on each descendant
                desc_nodes = lineages$desc[desc_idx]
                for(desc_node in desc_nodes) {
                    computeIncludeFossils( desc_node )
                }

                if ( (lineages$status[row_idx] == "sampled") == TRUE & (lineages$end_time[row_idx] == 0.0) == TRUE ) {
                    mark_drop[row_idx] <<- FALSE
                } else {
                    # check if any descendants are included
                    mark_drop[row_idx] <<- all(mark_drop[desc_idx] == TRUE)
                }

                # if this is a split, and only has one sampled descendant,
                # mark it as a knuckle
                if ( lineages$status[row_idx] == "split" & sum(mark_drop[desc_idx]) == 1 ) {
                    lineages$status[row_idx] <<- "knuckle"
                }

            }

        }

    )

)


simulateLHBDFP = function(
    lambda,        # a 1 x k vector of speciation rates
    mu,            # a 1 x k vector of extinction rates
    phi,           # a 1 x k vector of sampling rates
    delta,         # a 1 x k vector of destructive-sampling rates
    upsilon,       # a k x n matrix of mass-speciation split probabilities (per-category probs are in columns)
    upsilon_times, # a 1 x n vector of times for mass-speciation events
    gamma,         # a k x n matrix of mass-extinction survival probabilities (per-category probs are in columns)
    gamma_times,   # a 1 x n vector of times for mass-extinction events
    Z,             # a k x k matrix of change probabilities between each pair of states during a mass-extinction event
    rho,           # a k x n matrix of mass-sampling probabilities (per-category probs are in columns)
    rho_times,     # a 1 x n vector of times for mass-sampling events
    xi,            # a k x n matrix of mass-destructive-sampling probabilities (per-category probs are in columns)
    xi_times,      # a 1 x n vector of times for mass-destructive-sampling events
    H,             # a k x k matrix of rates of change between each pair of categories
    Omega,         # a k x k x k array of cladogenetic probabilities
    init,          # the initial category
    t,             # the time to simulate
    condition = c("survival", "MRCA")
) {

    repeat {

        # simulate a tree
        sim = simLHBDFP(lambda, mu, phi, delta,
                        upsilon, upsilon_times,
                        gamma, gamma_times, Z,
                        rho, rho_times,
                        xi, xi_times,
                        H, Omega,
                        init, t)

        # check conditions
        cond = TRUE

        # check survival
        if ( "survival" %in% condition ) {
            cond = cond & any(sim$lineages$end_time == 0 & sim$lineages$status == "sampled")
        }

        # check MRCA
        if ( "MRCA" %in% condition ) {
            cond = cond & sum(sim$lineages$end_time == 0 & sim$lineages$status == "sampled") >= 2
        }

        # if all of the conditions are met, end
        if ( cond == TRUE ) {
            break
        }

    }

    # return
    return(sim)

}


simLHBDFP = function(
    lambda,        # a 1 x k vector of speciation rates
    mu,            # a 1 x k vector of extinction rates
    phi,           # a 1 x k vector of sampling rates
    delta,         # a 1 x k vector of destructive-sampling rates
    upsilon,       # a k x n matrix of mass-speciation split probabilities (per-category probs are in columns)
    upsilon_times, # a 1 x n vector of times for mass-speciation events
    gamma,         # a k x n matrix of mass-extinction survival probabilities (per-category probs are in columns)
    gamma_times,   # a 1 x n vector of times for mass-extinction events
    Z,             # a k x k matrix of change probabilities between each pair of states during a mass-extinction event
    rho,           # a k x n matrix of mass-sampling probabilities (per-category probs are in columns)
    rho_times,     # a 1 x n vector of times for mass-sampling events
    xi,            # a k x n matrix of mass-destructive-sampling probabilities (per-category probs are in columns)
    xi_times,      # a 1 x n vector of times for mass-destructive-sampling events
    H,             # a k x k matrix of rates of change between each pair of categories
    Omega,         # a k x k x k array of cladogenetic probabilities
    init,          # the initial category
    t              # the time to simulate
) {

    # how many categories?
    num_cats = length(lambda)

    # get the dominant rates for each type of asynchronous event
    dominant_lambdas = sapply(lambda, findDominantRate, t, 0)
    dominant_mus     = sapply(mu, findDominantRate, t, 0)
    dominant_phis    = sapply(phi, findDominantRate, t, 0)
    dominant_deltas  = sapply(delta, findDominantRate, t, 0)

    # make a container for the lineages
    lin = data.frame(anc=1, desc=2, start_time=t, end_time=NA, state=init, status="active", stringsAsFactors=FALSE)
    lin_index  = 2
    num_active = 2

    # simulate forward in time
    current_time = t
    repeat {

        # determine which lineages are active
        active_lineages     = which(lin$status == "active")
        num_active_lineages = length(active_lineages)

        # if all lineages are inactive, end
        if ( num_active_lineages == 0 ) {
            break
        }

        # get the active state for each lineage
        active_states = lin[active_lineages,]$state

        # get the waiting times for each lineage
        this_waiting_time = -Inf
        this_lineage      = 0
        this_event_type   = 0
        for(i in 1:num_active_lineages) {

            # get the state for this lineage
            this_active_state = active_states[i]

            # simulate event times
            speciation_time = simulateWaitingTime(lambda[[this_active_state]],
                                                  current_time,
                                                  dominant_lambdas[this_active_state])
            
            extinction_time = simulateWaitingTime(mu[[this_active_state]],
                                                  current_time,
                                                  dominant_mus[this_active_state])

            fossilization_time = simulateWaitingTime(phi[[this_active_state]],
                                                     current_time,
                                                     dominant_phis[this_active_state])

            destruction_time = simulateWaitingTime(delta[[this_active_state]],
                                                     current_time,
                                                     dominant_deltas[this_active_state])

            change_time = simulateWaitingTime(function(t){-H[this_active_state,this_active_state]},
                                              current_time,
                                              -H[this_active_state,this_active_state])

            # check if these events are younger than the current first event time
            if ( speciation_time > this_waiting_time ) {
                this_waiting_time = speciation_time
                this_lineage      = i
                this_event_type   = 1
            }

            if ( extinction_time > this_waiting_time ) {
                this_waiting_time = extinction_time
                this_lineage      = i
                this_event_type   = 2
            }

            if ( fossilization_time > this_waiting_time ) {
                this_waiting_time = fossilization_time
                this_lineage      = i
                this_event_type   = 3
            }

            if ( destruction_time > this_waiting_time ) {
                this_waiting_time = destruction_time
                this_lineage      = i
                this_event_type   = 4
            }

            if ( change_time > this_waiting_time ) {
                this_waiting_time = change_time
                this_lineage      = i
                this_event_type   = 5
            }

        } # end loop over lineages
        
        this_lineage = active_lineages[this_lineage]
        event_type   = this_event_type

        # compute the next event time
        next_time = this_waiting_time

        # compute the "wall time"
        upsilon_wall_time = upsilon_times[upsilon_times < current_time]
        gamma_wall_time   = gamma_times[gamma_times < current_time]
        rho_wall_time     = rho_times[rho_times < current_time]
        xi_wall_time      = xi_times[xi_times < current_time]
        wall_time         = max(c(upsilon_wall_time,
                                  gamma_wall_time,
                                  rho_wall_time,
                                  xi_wall_time,
                                  0))

        # if next time exceeds wall time, set to wall time and end
        if ( next_time < wall_time ) {

            # this is a synchronous event

            # truncate the current time
            current_time = wall_time

            if ( current_time %in% upsilon_times ) {

                # mass-speciation event
                new_lineages = data.frame(anc=NULL, desc=NULL, start_time=NULL, end_time=NULL, state=NULL, status=NULL, stringsAsFactors=FALSE )

                # get the split probs
                this_mass_event  = which(upsilon_times == current_time)
                these_mass_probs = upsilon[,this_mass_event]

                # check each row
                for(i in active_lineages) {

                    # check for speciation
                    this_prob = these_mass_probs[lin$state[i]]
                    if ( runif(1) < this_prob ) {

                        # update the previous lineage
                        lin$end_time[i] = current_time
                        lin$status[i]   = "split"

                        # get the states for the descendants
                        this_clado_matrix = Omega[,,lin$state[i]]
                        left_state  = sample.int( nrow(this_clado_matrix), size=1, prob=rowSums(this_clado_matrix) )
                        right_state = sample.int( nrow(this_clado_matrix), size=1, prob=this_clado_matrix[left_state,])

                        # create the new lineages
                        new_lin = data.frame(anc=lin$desc[i],
                                             desc=1:2 + lin_index,
                                             start_time=current_time, end_time=NA,
                                             state=c(left_state, right_state), status="active", stringsAsFactors=FALSE)

                        # append the new lineages
                        new_lineages = rbind(new_lineages, new_lin)

                        # increment the lineage index
                        lin_index = lin_index + 2

                    }

                }

                # attach the new lineages
                lin = rbind(lin, new_lineages)

            } else if ( current_time %in% gamma_times ) {

                # mass-extinction event

                # get the extinction probs
                this_mass_event  = which(gamma_times == current_time)
                these_mass_probs = gamma[,this_mass_event]

                for(i in active_lineages) {

                    # check for speciation
                    this_prob = these_mass_probs[lin$state[i]]
                    if ( runif(1) < this_prob ) {

                        # update the previous lineage
                        lin$end_time[i] = current_time
                        lin$status[i]   = "extinct"

                    } else {

                        # check for a state-change event
                        state_probs = Z[lin$state[i],]
                        new_state   = sample.int(num_cats, size=1, prob=state_probs)

                        if (new_state != lin$state[i]) {

                            # execute state-change event
                            lin$end_time[i] = current_time
                            lin$status[i]   = "switched"

                            # create the new lineages
                            new_lin = data.frame(anc=lin$desc[i],
                                                 desc=1 + lin_index,
                                                 start_time=current_time, end_time=NA,
                                                 state=new_state, status="active", stringsAsFactors=FALSE)

                            # append the new lineages
                            lin = rbind(lin, new_lin)

                            # increment the lineage index
                            lin_index = lin_index + 1

                        }

                    }

                }

            } else if ( current_time %in% rho_times ) {

                # mass-sampling event
                rho_probs = rho[,rho_times == current_time]

                # for each active lineage, compute the sampling prob
                rho_prob_per_lineage = rho_probs[lin$state[active_lineages]]

                # determine which lineages get sampled
                sampled_lineages = as.logical(rbinom(num_active_lineages, size=1, prob=rho_prob_per_lineage ))

                # mark these lineages as sampled
                lin$end_time[active_lineages[sampled_lineages]] = current_time
                lin$status[active_lineages[sampled_lineages]]   = "sampled"

                # create the new lineages
                num_sampled = sum(sampled_lineages)

                if ( num_sampled > 0 & current_time > 0 ) {

                    new_lin = data.frame(anc=lin$desc[active_lineages[sampled_lineages]],
                                         desc=1:num_sampled + lin_index,
                                         start_time=current_time, end_time=NA,
                                         state=lin$state[active_lineages[sampled_lineages]], status="active", stringsAsFactors=FALSE)

                    # append the new lineages
                    lin = rbind(lin, new_lin)

                    # increment the lineage index
                    lin_index = lin_index + num_sampled

                }

            } else if ( current_time %in% xi_times ) {

                # mass-destructive-sampling event
                xi_probs = xi[,xi_times == current_time]

                # for each active lineage, compute the extinction prob
                xi_prob_per_lineage = xi_probs[lin$state[active_lineages]]

                # determine which lineages go extinct
                sampled_lineages = as.logical(rbinom( num_active_lineages, size=1, prob=xi_prob_per_lineage ))

                # mark these lineages as extinct
                lin$end_time[active_lineages[sampled_lineages]] = current_time
                lin$status[active_lineages[sampled_lineages]]   = "sampled"

            } else {
                stop("Oops!")
            } # end if/else over synchronous events

            # now, move on to the next step
            if ( current_time <= 0 ) {
                break
            } else {
                next
            }

        } else {

            # this is an asynchronous event
            current_time = next_time

            if ( event_type == 1 ) {

                # speciation event

                # update the previous lineage
                lin$end_time[this_lineage] = current_time
                lin$status[this_lineage]   = "split"

                # get the states for the descendants
                this_clado_matrix = Omega[,,lin$state[this_lineage]]
                left_state  = sample.int( nrow(this_clado_matrix), size=1, prob=rowSums(this_clado_matrix) )
                right_state = sample.int( nrow(this_clado_matrix), size=1, prob=this_clado_matrix[left_state,])

                # create the new lineages
                new_lin = data.frame(anc=lin$desc[this_lineage],
                                     desc=1:2 + lin_index,
                                     start_time=current_time, end_time=NA,
                                     state=c(left_state, right_state), status="active", stringsAsFactors=FALSE)

                # append the new lineages
                lin = rbind(lin, new_lin)

                # increment the lineage index
                lin_index = lin_index + 2

            } else if ( event_type == 2 ) {

                # extinction event

                # update the previous lineage
                lin$end_time[this_lineage] = current_time
                lin$status[this_lineage]   = "extinct"

            } else if ( event_type == 3 ) {

                # sampling event

                # update the previous lineage
                lin$end_time[this_lineage] = current_time
                lin$status[this_lineage]   = "sampled"

                # create the new lineages
                new_lin = data.frame(anc=lin$desc[this_lineage],
                                     desc=1 + lin_index,
                                     start_time=current_time, end_time=NA,
                                     state=lin$state[this_lineage], status="active", stringsAsFactors=FALSE)

                # append the new lineages
                lin = rbind(lin, new_lin)

                # increment the lineage index
                lin_index = lin_index + 1

            } else if ( event_type == 4 ) {

                # destructive-sampling event

                # update the previous lineage
                lin$end_time[this_lineage] = current_time
                lin$status[this_lineage]   = "sampled"

            } else if ( event_type == 5 ) {

                # state-change event

                # update the previous lineage
                lin$end_time[this_lineage] = current_time
                lin$status[this_lineage]   = "switched"

                # choose the new category
                cat_probs = H[lin$state[this_lineage],]
                cat_probs[lin$state[this_lineage]] = 0
                new_cat = sample.int(num_cats, size=1, prob=cat_probs)

                # create the new lineages
                new_lin = data.frame(anc=lin$desc[this_lineage],
                                     desc=1 + lin_index,
                                     start_time=current_time, end_time=NA,
                                     state=new_cat, status="active", stringsAsFactors=FALSE)

                # append the new lineages
                lin = rbind(lin, new_lin)

                # increment the lineage index
                lin_index = lin_index + 1

            } else {
                stop("Oops!")
            }

        }

    } # end repeat until simulation complete

    # mark unsampled lineages
    lin$status[is.na(lin$end_time)]   = "unsampled"
    lin$end_time[is.na(lin$end_time)] = 0

    # create the reference-class object
    sim = heterogeneousBDFPTree(lin, num_cats)

    # return
    return(sim)

}



findDominantRate = function(rate_function, t0, t1) {

    # optimize the rate over the interval
    optim_rate = optimize(rate_function, lower=t1, upper=t0, maximum=TRUE)$objective

    # include the end points
    max_rate   = max(c( optim_rate, rate_function(t0), rate_function(t1) ))

    return(max_rate)

}

simulateWaitingTime = function(rate_function, initial_time, dominant_rate) {

    current_time = initial_time
    repeat {

        # simulate a waiting time
        current_time = current_time - rexp(1, dominant_rate)

        # check if event is too late
        if ( current_time < 0 ) {
            return(-Inf)
        }
        
        # check if this event is retained
        retain_prob = rate_function(current_time) / dominant_rate
        if ( runif(1) < retain_prob ) {
            # end simulation
            break
        }

    }

    # return the waiting time
    return(current_time)

}
































