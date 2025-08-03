# miRNA-TF-DEG Network Analysis
# Author: Your Name
# Date: August 2025
# Description: Constructs a validated miRNA-target interaction network for a set of DEGs,
# identifies highly interacting miRNAs, and explores shared TF-miRNA regulatory targets.

# ======================= SETUP ============================
# Load necessary libraries
suppressPackageStartupMessages({
  library(multiMiR)
  library(dplyr)
  library(igraph)
})

# Set working directory to project root (user should modify this)
# setwd("path_to_your_project_directory")

# =================== INPUT: Differentially Expressed Genes ===================
DEG_genes <- c("CCL5", "CD244", "CD27", "CTLA4", "CXCR6", "FASLG", "GZMB",
               "IL2RB", "KLRD1", "LCK", "NCR3", "PRF1", "TBX21", "ZAP70")

# =============== STEP 1: Get Validated miRNA-Target Interactions ===============
miRNA_results <- get_multimir(target = DEG_genes, table = "validated", summary = TRUE)
miRNA_summary <- miRNA_results@summary

# Replace NA with 0
na_cols <- c("mirtarbase", "tarbase", "validated.sum")
miRNA_summary[na_cols] <- lapply(miRNA_summary[na_cols], function(x) replace(x, is.na(x), 0))

# Add interaction score
miRNA_summary$interaction_score <- rowSums(miRNA_summary[na_cols])

# Save all results
write.csv(miRNA_summary, "output/miRNA_summary.csv", row.names = FALSE)

# =============== STEP 2: Filter by Score Threshold ===============
threshold <- 4
miRNA_filtered <- miRNA_summary %>%
  filter(interaction_score >= threshold)

# =============== STEP 3: Network Preparation ===============
# Nodes
node_ids <- unique(c(miRNA_filtered$mature_mirna_id, miRNA_filtered$target_symbol))
nodes <- data.frame(
  id = node_ids,
  type = ifelse(node_ids %in% miRNA_filtered$mature_mirna_id, "miRNA", "DEG")
)

# Edges
edges <- miRNA_filtered %>%
  select(mature_mirna_id, target_symbol, interaction_score) %>%
  rename(from = mature_mirna_id, to = target_symbol)

# Add ranking
edges <- edges %>%
  arrange(desc(interaction_score)) %>%
  mutate(miRNA_rank = rank(-interaction_score, ties.method = "first"))

# =============== STEP 4: Build & Visualize Graph ===============
miRNA_graph <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# Set attributes
V(miRNA_graph)$label <- V(miRNA_graph)$name
V(miRNA_graph)$color <- ifelse(V(miRNA_graph)$type == "miRNA", "skyblue", "orange")

# Plot graph
plot(miRNA_graph,
     vertex.size = 5,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     vertex.label.dist = 1.5,
     main = paste("miRNA-Target Network (Threshold ≥", threshold, ")"))

# Export files
write.csv(edges, "output/miRNA_target_network_filtered.csv", row.names = FALSE)
write_graph(miRNA_graph, "output/miRNA_target_network.graphml", format = "graphml")

# =============== STEP 5: Network Summary ===============
cat("Network Summary:\n")
cat("Nodes:", vcount(miRNA_graph), "\n")
cat("Edges:", ecount(miRNA_graph), "\n")

# =============== STEP 6: TF-miRNA Regulatory Crosstalk ===============
# Load TF and miRNA target files
mirna_targets <- read.csv("data/miRNA_target_for_TF_miRNA.csv")
tf_targets <- read.csv("data/tf_targets_for_TF_miRNA.csv")

# Merge to find shared targets
tf_mirna <- merge(tf_targets, mirna_targets, by = "Target")[, c("TF", "miRNA")]
tf_mirna <- unique(tf_mirna)

# Save results
write.csv(tf_mirna, "output/tf_mirna.csv", row.names = FALSE)

cat("TF-miRNA pairs found:", nrow(tf_mirna), "\n")
