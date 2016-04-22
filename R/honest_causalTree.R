# call honest.causalTree.R to caculate the/ honest estimate
# value of causal effect in each node.
# When the honest estimate is NA in one node, we use its parent's value to replace it.

# find out all sons of parent in the leaf vector:
# only suitable for full binary tree:
sons <- function(parent, leaf) {
    max <- 20
    result <- sonshelper(parent, leaf, 0, max)
    return (result)
}

sonshelper <- function(parent, leaf, count, maxdepth) {
    if (count > maxdepth) {
        # rarely happen
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

# give the id of the node, calculate the causal_effect of the node
calculateEffect <- function(idx, leaf, where, treat, control, estY, estW){
    parent_node_id <- idx
    effect = NA
    while(is.na(effect)){
        sons_node_id <- sons(parent_node_id, leaf)
        obs_in_parent<- which(where %in% sons_node_id)
        index1 <- intersect(obs_in_parent, treat)
        index0 <- intersect(obs_in_parent, control)
        effect <- weighted.mean(estY[index1], estW[index1]) - weighted.mean(estY[index0], estW[index0])
        parent_node_id <- floor(parent_node_id / 2)
    }
    return (effect)
}

honest_causalTree <- function(object, estX, estY, estNobs, estW, estTreatment)
{
    leaf <- as.numeric(row.names(object$frame)[which(object$frame$var == "<leaf>")])
    treat <- which(estTreatment == 1)
    control<- which(estTreatment == 0)
    where <- est.causalTree(object, causalTree.matrix(estX))
    where_honest = NULL
    
    ## to be edited:
    object$where_honest <- where
    
    nodes <- leaf
    while (length(nodes) > 1) {
        for (idx in nodes) {
            effect <- calculateEffect(idx, leaf, where, treat, control, estY, estW)
            tmp <- which(as.numeric(as.character(rownames(object$frame))) == idx)
            object$frame$yval[tmp] <- effect
        }
        nodes <- floor(nodes / 2)
        nodes <- setdiff(unique(nodes), c(0))
    }
    # calculate for the root node and replace it with yval:
    effect <- calculateEffect(1, leaf, where, treat, control, estY, estW)
    tmp <- which(as.numeric(rownames(object$frame)) == 1)
    object$frame$yval[tmp] <- effect
    return(object)
}

