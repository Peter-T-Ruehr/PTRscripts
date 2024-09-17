#' Adds taxon as sister to a tip.
#' @param tree A phylogenetic tree of class `phylo`.
#' @param tip A single character variable containing the exact tip name of the
#' sister froup of the new taxon. Must be present as a tree$tip.label.
#' @param new_tips A single character variable or a string of characters
#' containing the tip names to add.
#' @param position A numeric value between`0` and `1` to set where on the edge
#' leading to the sister tip the new splitting node should be placed.
#' @return A Phylogenetic tree of class `phylo`.
#' @export
#' @examples
#' set.seed(1)
#' # library(tibble)
#' # library(phytools)
#' # library(ape)
#' library(dplyr)
#'
#' tree <- phytools::pbtree(n=10,scale=1)
#'
#' tree$tip.label <- LETTERS[ape::Ntip(tree):1]
#' # plot(tree)
#'
#' par(mfrow=c(1,4))
#'
#' tip_colors_orig <- tibble::tibble(tip = tree$tip.label, color = NA) %>%
#'   mutate(color = case_when(grepl("^N", tip) ~ "red",
#'                            tip == "C" ~ "green",
#'                            TRUE ~ "darkgreen"))
#'
#' ape::plot.phylo(tree,
#'            mar=c(0.1,0.1,1.1,0.1),
#'            cex=1.5,
#'            tip.color = tip_colors_orig$color)
#' title(main = paste0("original"), cex.main=1.5)
#'
#' for(i in c(.2,.5,.8)){ # c(.1,.25,.5,.75,.9)
#'   tree_new1 <- add_cherry_to_tip(tree = tree,
#'                                  tip = "C",
#'                                  new_tips = "N1",
#'                                  position = i)
#'
#'
#'   # get colors
#'   tip_colors <- tibble::tibble(tip = tree_new1$tip.label, color = NA) %>%
#'     mutate(color = case_when(grepl("^N", tip) ~ "red",
#'                              tip == "C" ~ "green",
#'                              TRUE ~ "darkgreen"))
#'
#'   ape::plot.phylo(tree_new1,
#'              mar=c(0.1,0.1,1.1,0.1),
#'              cex=1.5,
#'              tip.color = tip_colors$color)
#'
#'   title(main = paste0("pos.= ", i), cex.main=1.5)
#' }
#'
#' ape::plot.phylo(tree,
#'            mar=c(0.1,0.1,1.1,0.1),
#'            cex=1.5,
#'            tip.color = tip_colors_orig$color)
#' title(main = paste0("original"), cex.main=1.5)
#'
#' for(i in c(.2,.5,.8)){ # c(.1,.25,.5,.75,.9)
#' tree_new1 <- add_cherry_to_tip(tree = tree,
#'                                  tip = "C",
#'                                  new_tips = c("N1", "N2", "N3", "N4"),
#'                                  position = i)
#'
#'   # get colors
#'   tip_colors <- tibble::tibble(tip = tree_new1$tip.label, color = NA) %>%
#'     mutate(color = case_when(grepl("^N", tip) ~ "red",
#'                              tip == "C" ~ "green",
#'                              TRUE ~ "darkgreen"))
#'
#'   ape::plot.phylo(tree_new1,
#'              mar=c(0.1,0.1,1.1,0.1),
#'              cex=1.5,
#'              tip.color = tip_colors$color)
#'
#'   title(main = paste0("pos.= ", i), cex.main=1.5)
#' }
#'
#' par(mar=c(5.1,4.1,4.1,2.1))
#'
add_cherry_to_tip <- function(tree,
                              tip,
                              new_tips,
                              position = 0.5) {

  require(tidytree)
  require(phytools)
  require(tidytree)
  require(dplyr)

  # # testing
  # tree = tree
  # tip = "t3"
  # new_tips = "new1"
  # position <- 0.5
  # plot(tree)

  tip.labels.init <- tree$tip.label



  # find original edge length
  edge.length <- tree$edge.length[which(tree$edge[,2]==
                                          which(tree$tip.label==tip))]

  # Find the edge leading to the tip
  tip_id <- match(tip, tree$tip.label)

  # Create the new cherry
  correct.topo = F
  while(correct.topo == F){
    tree_to_add <- phytools::pbtree(n=length(new_tips)+1, nsim = 1)

    # Naming the tips
    tree_to_add$tip.label <- c(tip, new_tips)

    # plot(tree_to_add)
    # nodelabels()

    if(length(new_tips) == 1){
      MRCA.desired = 2
    } else {
      MRCA.desired = length(new_tips)+3
    }
    # create new trees until new tips are sister taxa and keep that tree
    if(tidytree::MRCA(tree_to_add, new_tips) == MRCA.desired){
      correct.topo = T

      # change branch lengths of tree
      if(length(new_tips) > 1){
        tree_to_add$edge.length <- rep(1, length(tree_to_add$edge.length))
        tree_to_add <- ape::chronopl(tree_to_add, lambda = 1, age.min = edge.length*(1-position), age.max = NULL, node = "root")
      } else if(length(new_tips) == 1){
        tree_to_add$edge.length <- rep((edge.length*(1-position)), length(tree_to_add$edge.length))
      }
    }
  }
  # ggtree(tree_to_add) +
  #   theme_tree2()

  # Binding both trees
  tree.bound <- ape::bind.tree(tree, tree_to_add, where = tip_id)

  # ggtree(tree.bound) +
  #   theme_tree2()

  # find MRCA of all tips (tip, new_tips)
  mrca.tips <- tidytree::MRCA(tree.bound, c(tip, new_tips))

  # get edge length of all tips (tip, new_tips)
  edge.number.all <- tree.bound$edge %>%
    tibble::as_tibble() %>%
    mutate(rn = row_number()) %>%
    filter(V2 == mrca.tips) %>%
    pull(rn)

  curr.tip.edge.length <- tree.bound$edge.length[edge.number.all]

  tree.bound$edge.length[edge.number.all] <- curr.tip.edge.length*position

  # ggtree(tree.bound) +
  #   theme_tree2()

  tip.labels.end <- tree.bound$tip.label
  diff.tip.labels <- setdiff(tip.labels.end, tip.labels.init)

  if(length(diff.tip.labels) > 0){
    print(paste0("Added the following ", length(diff.tip.labels), " tip label(s):"))
    print(diff.tip.labels)
  } else {
    message("!!! no tips were added !!!")
  }
  return(tree.bound)
}
