library(dplyr)
library(visNetwork)
library(RColorBrewer)
library(htmlwidgets)

# load("GSEA_results.Rdata")
# gsea_result = results$cp$GSEA_H

plot_network = function(gsea_result) {
  # if hallmark GSEA, then remove "HALLMARK" text from start
  gsea_result$pathway = gsub("^HALLMARK_", "", gsea_result$pathway)
  
  # nodes dataframe
  nodes_df = data.frame(
    id = as.character(gsea_result$pathway), 
    label = gsea_result$pathway,
    p_value = gsea_result$padj,
    nes = gsea_result$NES,
    shape = "dot",
    title = paste0(
      "Gene Set: ", gsea_result$pathway, "<br>",
      "adj.p-value: ", sprintf("%.3e", gsea_result$padj), "<br>",
      "NES: ", sprintf("%.2f", gsea_result$NES), "<br>"
    ),
    stringsAsFactors = FALSE
  )
  
  # assign colours to nodes
  fixed_breaks_nes = seq(-4, 4, length.out = 101)
  nodes_df$nes_clipped = pmax(-4, pmin(4, nodes_df$nes))
  
  color_bins = cut(nodes_df$nes_clipped, breaks = fixed_breaks_nes,
                   include.lowest = TRUE, right = TRUE,
                   labels = FALSE)
  
  nes_palette = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
  nodes_df$color.background = nes_palette[color_bins]
  
  # assign border and highlight colours
  nodes_df$color.border = nodes_df$color.background
  nodes_df$color.highlight.background = "white"
  nodes_df$color.highlight.border = nodes_df$color.background
  
  # Opacity calculation
  nodes_df$color.opacity = pmax(0.1, pmin(1, 1 - nodes_df$p_value))  # Min opacity of 0.1
  
  # edges dataframe
  edges_df = data.frame(from = character(), to = character(), value = numeric(), stringsAsFactors = FALSE)
  
  gene_sets = as.character(gsea_result$pathway) 
  leading_edges = gsea_result$leadingEdge
  
  # Iterate through all pairs of gene sets to calculate overlap
  edge_count = 0
  for (i in 1:(length(gene_sets) - 1)) {
    for (j in (i + 1):length(gene_sets)) {
      le1 = leading_edges[[i]]
      le2 = leading_edges[[j]]
      
      # Skip if either leading edge is NULL or empty
      if(is.null(le1) || is.null(le2) || length(le1) == 0 || length(le2) == 0) {
        next
      }
      
      # Calculate Jaccard Index
      intersection_size = length(intersect(le1, le2))
      union_size = length(union(le1, le2))
      jaccard_index = if (union_size > 0) intersection_size / union_size else 0
      
      # Add edges with threshold
      if (jaccard_index > 0) {
        edge_count = edge_count + 1
        edges_df = rbind(edges_df, data.frame(
          from = gene_sets[i],
          to = gene_sets[j],
          value = jaccard_index * 10,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Filter for top 100 edges
  edges_df = edges_df[order(edges_df$value, decreasing=T), ]
  edges_df = edges_df[1:100,]
  
  
  # Create network visualisation
  network_plot = visNetwork(nodes_df, edges_df, height = "700px", width = "100%") %>%
    visNodes(
      size = 20,
      borderWidth = 1,
      font = list(size = 15, color = "black"),
      shadow = list(enabled = F)
    ) %>%
    visEdges(
      smooth = list(enabled = TRUE, type = "continuous"),  # Changed from cubicBezier
      color = list(color = "lightgray", highlight = "darkgray"),
      arrows = list(to = list(enabled = FALSE))
    ) %>%
    visLayout(
      improvedLayout = TRUE,
      randomSeed = 0
    ) %>%
    visOptions(
      highlightNearest = TRUE,
      nodesIdSelection = TRUE
    ) %>%
    visPhysics(
      enabled = TRUE,
      solver = "barnesHut",
      barnesHut = list(
        gravitationalConstant = -2000,
        centralGravity = 0.3,  
        springLength = 135,    
        springConstant = 0.01,
        avoidOverlap = 1.5
      )
    ) %>%
    visInteraction(
      tooltipDelay = 0
    )
  
  return(network_plot)
}


