# Packages
library(tidyverse)
library(igraph)
library(ggraph)

# ---- 1) Your table (rename to `links`) ----
# If you already have it in R, skip this and just ensure it’s called `links`.
links <- tribble(
  ~from, ~to, ~time, ~support, ~generations,
  6,     1,   188, 1.0000000,   5,
  4,     2,   -12, 0.9261477,   1,
  16,     3,   -18, 0.5608782,   1,
  14,     4,   -31, 0.7445110,   1,
  16,     5,   -11, 0.5888224,   1,
  4,     6,   -10, 0.3892216,   1,
  14,     7,   -17, 0.9381238,   1,
  NA,     8,   -69,        NA,  NA,
  12,     9,   133, 0.9161677,   5,
  4,    10,   184, 0.9201597,   5,
  6,    11,   315, 1.0000000,   5,
  8,    12,   -51, 0.9680639,   1,
  4,    13,   -14, 0.8882236,   1,
  12,    14,   -39, 0.5409182,   1,
  14,    15,   -13, 0.5588822,   1,
  14,    16,   -24, 0.5828343,   1,
  4,    17,   -11, 0.9281437,   1,
  4,    18,    -9, 0.9281437,   1,
  4,    19,   -14, 0.8502994,   1,
  12,    20,   -23, 0.9680639,   1
)

# ---- 2) Build graph ----
# Identify imports (roots): rows with from = NA
imports <- links %>% filter(is.na(from)) %>% pull(to)

# Edges: drop NA-from for graph construction
edges <- links %>% filter(!is.na(from)) %>%
  transmute(from = as.character(from),
            to   = as.character(to),
            support)

# Node list
nodes <- tibble(id = unique(c(edges$from, edges$to))) %>%
  mutate(id = as.character(id),
         is_import = id %in% as.character(imports))

# igraph
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# ---- 3) Compute node depth (generations from import) ----
# Distance from (one or more) import roots
if (length(imports) > 0) {
  d <- distances(g, v = V(g)[as.character(imports)], to = V(g), mode = "out")
  # If multiple imports, take the minimum distance from any root
  depth <- apply(d, 2, min)
} else {
  depth <- rep(0, vcount(g))
}
V(g)$depth <- as.numeric(depth)
V(g)$label <- V(g)$name

# ---- 4) Plot (left-to-right tree layout) ----
# Use a grayscale palette for journals; edge width encodes support.
p <- ggraph(g, layout = "tree") +
  geom_edge_link(aes(edge_alpha = support,
                     edge_width = support),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm'),
                 arrow = arrow(type = "closed", length = unit(3, "mm")),
                 colour = "grey20") +
  geom_node_point(aes(shape = is_import, size = is_import),
                  stroke = 0.7, colour = "black", fill = "white") +
  geom_node_text(aes(label = label), vjust = -1.0, size = 3) +
  scale_edge_alpha(range = c(0.2, 1)) +
  scale_edge_width(range = c(0.2, 1.4)) +
  scale_shape_manual(values = c(`TRUE` = 22, `FALSE` = 21)) +  # squares for imports
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 3)) +
  guides(edge_width = "none", edge_alpha = "none", shape = "none", size = "none") +
  theme_void() 
#ggtitle("Inferred transmission network (arrows: infector → infectee)",
#subtitle = "Edge thickness and opacity ∝ posterior support; squares mark import roots")

print(p)

# Save figure
ggsave("transmission_network.png", p, width = 8, height = 5, dpi = 300)
ggsave("transmission_network.pdf", p, width = 8, height = 5)


