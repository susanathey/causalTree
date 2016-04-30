# only suitable for full binary tree:
sons <- function(parent, leaf) {
	max <- floor(log(max(leaf), 2))
	result <- sonshelper(parent, leaf, 0, max)
	return (result)
}

sonshelper <- function(parent, leaf, count, maxdepth) {
	if (count > maxdepth) {
		# rarely happens
		stop ("This node is not in the tree.")
	} else {
		if (!is.na(match(parent, leaf))) {
			return (parent)
		} else {
			left_son <- 2 * parent
			right_son <-  2 * parent + 1
			left_sons <- sonshelper(left_son, leaf, count + 1, maxdepth)
			right_sons <- sonshelper(right_son, leaf, count + 1, maxdepth)
			result <- c(left_sons, right_sons)
		}
		return(result)    
	} 
}

estimate.causalTree <- function(object, data, treatment, na.action = na.causalTree)
{
	if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
	data$w <- treatment
	# get the leaf of the object
	leaf <- as.numeric(row.names(object$frame)[which(object$frame$var == "<leaf>")])


	# get the node id for each observation:
	#Terms <- delete.response(object$terms)
	Terms <- object$terms
	#data <- model.frame(Terms, data, na.action = na.action,
	data <- model.frame(Terms, data, na.action = na.action, treatment = w, 
						xlev = attr(object, "xlevels"))
	#print (data)

	if (!is.null(cl <- attr(Terms, "dataClasses")))
		.checkMFClasses(cl, data, TRUE)

	treatment <- data$`(treatment)`
	n <- nrow(data)
	Y <- model.response(data)
	where <- est.causalTree(object, causalTree.matrix(data))
	#return (where)
	#print (where)
	#return (where)


	#check the treatment condition:
	if (missing(treatment)) stop("You should import the treatment status for data.")
	if (length(treatment) != n) 
		stop("The length of treatment status vector should be same as number
			 of observations.")
			 if (length(which(treatment == 0)) == 0 || length(which(treatment == 1)) == 0)
				 stop("Can't make estimation since no treated cases or no control cases.")


			 unique_leaf <- unique(where)
			 causal_estimation <- rep(0, n)
			 treat <- which(treatment == 1)
			 control<- which(treatment == 0)
			 for (i in 1:length(unique_leaf)) {
				 leaf_id <- unique_leaf[i]
				 index <- which(where == leaf_id)
				 index1 <- intersect(index, treat)
				 index0 <- intersect(index, control)    
				 effect <- mean(Y[index1]) - mean(Y[index0])
				 parent_node_id <- leaf_id
				 while(is.na(effect)){
					 ## go back to parent node who can compute the value:
					 parent_node_id <- floor(parent_node_id / 2)
					 sons_node_id <- sons(parent_node_id, leaf)
					 obs_in_parent<- which(where %in% sons_node_id)
					 index1 <- intersect(obs_in_parent, treat)
					 index0 <- intersect(obs_in_parent, control)
					 effect <- mean(Y[index1]) - mean(Y[index0])
				 }
				 causal_estimation[as.numeric(index)] = effect
			 }
			 return(causal_estimation)  
}
