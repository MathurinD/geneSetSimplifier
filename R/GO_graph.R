#' @import GO.db
#' @import tidyverse
#' @import Rgraphviz

    library(Rgraphviz)
    library(GO.db)
    library(tidyverse)
    #enrichGO(c(''), 'org.Hs.eg.db', ont="BP",universe=unique(phosphosites$entrezgene_id), pvalueCutoff=0.01) -> enrichment
    enrichment = go_e1
    threshold = 0.05
    max_depth = Inf

#' Fit text to specified width by introducing linebreaks
fitToWidth <- function(text, width=20) {
    words = unlist(str_split(text, ' '))
    final = c()
    cline = ""
    for (ww in words) {
        if (nchar(cline) + 1 + nchar(ww) <= width) {
            cline = paste(cline, ww, sep=' ')
        } else if (nchar(ww) >= width) {
            final = paste(final, cline, sep='\n')
            cline = ww
        } else {
            final = paste(final, cline, sep='\n')
            cline = ww
        }
    }
    final = paste(final, cline, sep='\n')
    return(substr(final, 3, nchar(final)))
}


#' @param enrichment Object of class enrichResult
#' @param threshold adjusted p-value threshold
GOGraph <- function(enrichment, threshold=0.05, max_depth = Inf) {
    enrich = enrichment@result %>% as_tibble
    leaf_terms = enrich %>% filter(p.adjust < threshold) %>% pull(ID)
    #GOBPPARENTS[leaf_terms] %>% as.list %>% unlist %>% {.[grepl('is_a',names(.))]} %>% {names(.)=gsub('.is_a', '', names(.)); .} %>% {graphNEL(unique(c(unlist(.),names(.))),.,'direct')} %>% layoutGraph %>% plot

    to_follow = leaf_terms
    edge_list = list()
    depth = 0
    while (length(to_follow) > 0 && depth < max_depth) {
        depth = depth + 1
        this_loop = unique(to_follow)
        to_follow = c()
        for (term in this_loop) {
            adjacency = tryCatch(GOBPPARENTS[term], error=function(e){list()}) %>% as.list
            if (length(adjacency) > 0 && ! term %in% c('all', 'GO:0008150')) {
                adjacency = adjacency %>% unlist %>% {.[grepl('is_a',names(.))]} %>% {names(.)=gsub('.is_a', '', names(.)); .}
                adjacency = adjacency %>% unlist %>% {tmp=names(.); names(tmp)=.; tmp} %>% as.list
                parents = names(adjacency)
                to_follow = c(to_follow, parents[!parents %in% names(edge_list)])
                edge_list = c(edge_list, adjacency)
            }
        }
    }
    goplot = graphNEL(nodes=unique(c(names(edge_list), unlist(edge_list))), edgeL=edge_list, 'directed') %>% as('graphAM')
    nodeRenderInfo(goplot)$shape = 'circle'
    nodeRenderInfo(goplot)$shape = enrich %>% filter(p.adjust < threshold, ID %in% nodes(goplot)) %>% {tmp=rep('rectangle', nrow(.)); names(tmp)=.$ID; tmp} # Rectangle for significant nodes
    gradient_size = 10
    nodeRenderInfo(goplot)$fillcolor = 'white'
    nodeRenderInfo(goplot)$fillcolor = enrich %>% filter(p.adjust < threshold, ID %in% nodes(goplot)) %>% {tmp=colorRampPalette(c('yellow', 'red'))(gradient_size)[sapply(-log10(.$p.adjust), max, gradient_size)]; names(tmp)=.$ID; tmp} # Color corresponding to p.adjust
    nodeRenderInfo(goplot)$label={tmp=paste(colnames(goplot@adjMat), GOTERM[colnames(goplot@adjMat)] %>% as.list %>% sapply(function(gt){gt@Term}) %>% sapply(fitToWidth), sep='\n'); names(tmp)=colnames(goplot@adjMat); tmp}
    #nodeRenderInfo(goplot)$
    nodeRenderInfo(goplot)$fontsize=30
    goplot %>% layoutGraph(layoutType='dot') %>% renderGraph() # Does not display the colors...
#    plot(goplot, 'dot', nodeAttrs=nodeRenderInfo(goplot)) # Does not display the labels...
}
