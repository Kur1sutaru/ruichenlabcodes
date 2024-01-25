library(eulerr)

# this_comb should be a vector of logicals indicating possible combinations
intersectAll <- function(...) {
    args <- as.list(...)
    nargs <- length(args)
    intersection <- args[[1]]
    # only loop if more than 1 arg
    if (nargs > 1) {
        for (i in seq(1, nargs-1)) {
            intersection <- intersect(intersection, args[[i+1]])
        }
    }
    return(intersection)
}

unionAll <- function(...) {
    args <- as.list(...)
    nargs <- length(args)
   
    # if there are no args (i.e. all args are in include, none in exclude), exit early
    if (nargs == 0) {
        return(NULL)
    }
   
    union <- args[[1]]
    if (nargs > 1) {
        for (i in seq(1, nargs-1)) {
            union <- union(union, args[[i+1]])
        }
    }
    return(union)
}


getSetComb <- function(set_list, this_comb) {
    include <- intersectAll(set_list[this_comb])
    exclude <- unionAll(set_list[!this_comb])
   
    set_comb <- setdiff(include, exclude)
    return(set_comb)

}
# pass in a named list, where the name is the category and the elements are vectors of things to count
getDisjointSets <- function(set_list) {
   
    nsets <- length(set_list)

    # to get all possible indices -- permutation with replacement, P(2, nsets) but one perm will be all F's - subtract 1
    # you can do this with expand.grid too but it's less intuitive imo
    combs_index <- gtools::permutations(n = 2, r = nsets, v = c(TRUE, FALSE), repeats.allowed = TRUE)
    combs_index <- combs_index[-1, ]
    ncombs <- nrow(combs_index)
   
    # rearrange order - necessary so combs stays ordered the same way as set_list
    # otherwise set_list may have names "a, b, c" but combs, based on T/F ordering, is "c, b, a". This won't mess up the
    # actual euler numbers, but will make labeling confusing
    combs_index <- as.data.frame(combs_index) %>%
        arrange(desc(across())) %>%
        as.matrix
   
    combs <- vector(mode = "list", length = ncombs)
    names <- vector(mode = "character", length = ncombs)
    for (i in 1:nrow(combs_index)) {
        this_comb <- unlist(combs_index[i, ])
        this_set <- getSetComb(set_list, this_comb)
        names[i] <- paste(names(set_list)[this_comb], collapse = "&")
        combs[[i]] <- this_set
    }
   
    names(combs) <- names
    return(combs)
}


set_list <- list("a" = df1_filt$ncbi_gene,
                     "b" = df2_filt$Entrez.Gene.ID.for.Human,
                     "c" = df3_filt$entrezgene_id)
disjoint_sets <- getDisjointSets(set_list)

disjoint_num <- sapply(disjoint_sets, length)

fit <- euler(disjoint_num, input = "disjoint")

plot(fit, quantities = TRUE, labels = c("GSEA immune", "IPA fibrosis", "Reactome ECM"))