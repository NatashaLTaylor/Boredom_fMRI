# Using Chord Diagram
# install.packages("circlize")
# install.packages("scales")
library(circlize)
library(scales)
library(colorspace)

# load in the functional connectivity matrix
m <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\sig_permute_coefs_dmn_within_ordered.csv", header = FALSE, sep = ",")
m <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\fc_dmn_within_nsqip_interact.csv", header = FALSE, sep = ",")

data <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\plot_fc_chord_data.csv", header = TRUE, sep = ",")

node_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/dmn_new_order_labels.csv", header = FALSE, stringsAsFactors = FALSE)
rownames(m) <- node_names$V1
colnames(m) <- node_names$V1

chordDiagram(m) # plots the diagram as is

# Restart circular layout parameters
circos.clear()


# attempt the plot without the 0 weights
# Filter out zero-weighted connections
non_zero_connections <- which(m != 0, arr.ind = TRUE)

# Create a data frame of connections with non-zero weights
connections <- data.frame(
  from = non_zero_connections[, 1],
  to = non_zero_connections[, 2],
  weight = m[non_zero_connections]
)


# Plot the diagram so it's grouped into the sub-networks of DMN
# load in sub-nework IDs
sub_net_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/sub_networks_DMN.csv", header = TRUE, stringsAsFactors = FALSE)


non_zero_connections <- which(connections != 0, arr.ind = TRUE)

region_names <- sub_net_names$ROI_subs
# unique_region_names <- make.unique(region_names)
## Creates the plot
connectivity_matrix <- m
connectivity_matrix[lower.tri(connectivity_matrix)] <- 0
# Filter out zero-weighted connections
non_zero_connections <- which(connectivity_matrix != 0, arr.ind = TRUE)


# Create a data frame of connections with non-zero weights
connections <- data.frame(
  from = region_names[non_zero_connections[, 1]],
  to = region_names[non_zero_connections[, 2]],
  weight = connectivity_matrix[non_zero_connections]
)


# group list

order <- c(
  grep("^IPL", region_names, value = TRUE),
  grep("^LatVent", region_names, value = TRUE),
  grep("^Rsp", region_names, value = TRUE),
  grep("^Para_Hippo", region_names, value = TRUE),
  grep("^Dorsal_pFC", region_names, value = TRUE),
  grep("^Medial_pFC", region_names, value = TRUE),
  grep("^Temp", region_names, value = TRUE),
  grep("^pCun_PPC", region_names, value = TRUE)
)


# plot specific colours for each sub-network
colors <- c(
  IPL = "#00a79cdc",
  LatVent_pFC = "#279e96cc",
  Rsp = "#00a79c67",
  Para_Hippo = "#5ac0b9d3",
  Dorsal_pFC = "#5b5da9",
  Medial_pFC = "#9799d8",
  Temp = "#5b5ea9a4",
  pCun_PCC = "#0681ba"
)

# Create a named vector that maps each region to a color based on its group
region_colors <- sapply(region_names, function(region) {
  if (grepl("^IPL", region)) {
    return(colors["IPL"])
  } else if (grepl("^LatVent_pFC", region)) {
    return(colors["LatVent_pFC"])
  } else if (grepl("^Rsp", region)) {
    return(colors["Rsp"])
  } else if (grepl("^Para_Hippo", region)) {
    return(colors["Para_Hippo"])
  } else if (grepl("^Dorsal_pFC", region)) {
    return(colors["Dorsal_pFC"])
  } else if (grepl("^Medial_pFC", region)) {
    return(colors["Medial_pFC"])
  } else if (grepl("^Temp", region)) {
    return(colors["Temp"])
  } else if (grepl("^pCun_PCC", region)) {
    return(colors["pCun_PCC"])
  }
})

connection_colors <- ifelse(connections$weight < 0,
  sapply(connections$from, function(region) {
    if (grepl("^IPL", region)) {
      return(darken(colors["IPL"], amount = 0.6))
    } else if (grepl("^LatVent_pFC", region)) {
      return(darken(colors["LatVent_pFC"], amount = 0.6))
    } else if (grepl("^Rsp", region)) {
      return(darken(colors["Rsp"], amount = 0.6))
    } else if (grepl("^Temp", region)) {
      return(darken(colors["Temp"], amount = 0.6))
    } else if (grepl("^Para_Hippo", region)) {
      return(darken(colors["Para_Hippo"], amount = 0.6))
    } else if (grepl("^Dorsal_pFC", region)) {
      return(darken(colors["Dorsal_pFC"], amount = 0.6))
    } else if (grepl("^Medial_pFC", region)) {
      return(darken(colors["Medial_pFC"], amount = 0.6))
    } else if (grepl("^pCun_PCC", region)) {
      return(darken(colors["pCun_PCC"], amount = 0.6))
    }
  }),
  sapply(connections$from, function(region) {
    if (grepl("^IPL", region)) {
      return(colors["IPL"])
    } else if (grepl("^LatVent_pFC", region)) {
      return(colors["LatVent_pFC"])
    } else if (grepl("^Rsp", region)) {
      return(colors["Rsp"])
    } else if (grepl("^Temp", region)) {
      return(colors["Temp"])
    } else if (grepl("^Para_Hippo", region)) {
      return(colors["Para_Hippo"])
    } else if (grepl("^Dorsal_pFC", region)) {
      return(colors["Dorsal_pFC"])
    } else if (grepl("^Medial_pFC", region)) {
      return(colors["Medial_pFC"])
    } else if (grepl("^pCun_PCC", region)) {
      return(colors["pCun_PCC"])
    }
  })
)


order <- c(
  "IPL",
  "LatVent_pFC",
  "Rsp",
  "Para_Hippo",
  "Dorsal_pFC",
  "Medial_pFC",
  "Temp",
  "pCun_PCC"
)

circos.clear()
# Open an SVG graphics device
svg("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/Figures/DMN_Figures/Within_Dmn_FC/chord_diagram_nsqdip_interact.svg", width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  grid.col = colors,
  col = connection_colors,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
)
# Add labels for the sectors
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

dev.off()
#only plot positive connections
# Assuming 'connections' has numerical values and you want to keep rows where all values are positive
positive_connections <- connections[apply(connections > 0, 1, all), ]
svg("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/Figures/DMN_Figures/Within_Dmn_FC/chord_diagram_nsqdip_interact_positive.svg", width = 20, height = 20)

chordDiagram(
  positive_connections,
  grid.col = colors,
  order = order,
  #col = connection_colors,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
  )
  dev.off()
  # Assuming 'connections' has numerical values and you want to keep rows where all values are positive
  negative_connections <- connections[connections$weight < 0, ]
#make the connections darker version of existing colour
colors_negative <- 

chordDiagram(
  negative_connections,
  grid.col = colors,
  order = order,
  #col = connection_colors,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
  )







## Create Chord Plot for Whole brain to FC
fc_dmn_whole_brain <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\ordered_fc_dmn_lm.csv", header = FALSE, sep = ",")
avg_fc_dmn_whole_brain <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\avg_ordered_fc_dmn_lm.csv", header = FALSE, sep = ",")
# avg. beta coefs 7nets, dmn sub-regions in there too
avg_fc_dmn_7_nets <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\avg_7nets_fc_dmn_betas.csv", header = FALSE, sep = ",")
# load in sub-nework IDs
sub_net_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/sub_networks_DMN.csv", header = TRUE, stringsAsFactors = FALSE)
schaef_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/schaef_names.csv", header = TRUE, stringsAsFactors = FALSE)
schaef_17nets <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/schaef_17_networks.csv", header = TRUE, stringsAsFactors = FALSE)
schaef_7nets <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/schaef_7_dmn_networks.csv", header = TRUE, stringsAsFactors = FALSE)


avg_fc_dmn_whole_brain <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\Chord_Diagrams\\avg_across_dmn_subnets_wholebrain.csv", header = FALSE, sep = ",")


region_names <- sub_net_names$ROI_subs
region_names2 <- c("IPL", "Dorsal_pFC", "Medial_pFC", "LatVent_pFC", "pCun_PCC", "Temp", "Rsp", "Para_Hippo")
schaef_nets <- schaef_names$Networks
schaef_nets <- schaef_7nets$Network
# label for the matrix
rownames(fc_dmn_whole_brain) <- make.unique(schaef_nets)
colnames(fc_dmn_whole_brain) <- region_names
rownames(avg_fc_dmn_whole_brain) <- make.unique(schaef_nets)
colnames(avg_fc_dmn_whole_brain) <- region_names2

rownames(avg_fc_dmn_7_nets) <- schaef_names$Networks
colnames(avg_fc_dmn_7_nets) <- region_names2
# makes unique row names for the schaef rows
# unique_schaef_names <- make.unique(schaef_nets)
# connections <- fc_dmn_whole_brain
connections <- avg_fc_dmn_whole_brain
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = schaef_nets[non_zero_connections[, 1]],
  to = region_names2[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)
# Assuming 'connections' has numerical values and you want to keep rows where all values are positive
positive_connections <- connections[apply(connections > 0, 1, all), ]
# Assuming 'connections' has numerical values and you want to keep rows where all values are positive
negative_connections <- connections[connections$weight < 0, ]
# reduce the dimensions of the data;

# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  transparency = 0.5,
)
# include labels and colour etc..

circos.clear()
# plot specific colours for each sub-network
colors <- c(
  IPL = "#00a79cdc",
  Dorsal_pFC = "#5b5da9",
  Medial_pFC = "#9799d8",
  LatVent_pFC = "#279e96cc",
  pCun_PCC = "#0681ba",
  Temp = "#5b5ea9a4",
  Rsp = "#00a79c67",
  Para_Hippo = "#5ac0b9d3"
)
order <- c(
  grep("^IPL", region_names2, value = TRUE),
  grep("^Dorsal_pFC", region_names2, value = TRUE),
  grep("^Medial_pFC", region_names2, value = TRUE),
  grep("^LatVent", region_names2, value = TRUE),
  grep("^pCun_PCC", region_names2, value = TRUE),
  grep("^Temp", region_names2, value = TRUE),
  grep("^Rsp", region_names2, value = TRUE),
  grep("^Para_Hippo", region_names2, value = TRUE)
)
# Open an SVG graphics device
svg("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/Figures/DMN_Figures/DMN_FC_whole_brain/chord_diagram.svg", width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  grid.col = colors,
  # col = connection_colors,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  transparency = 0.2
)
# Add labels for the sectors
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

dev.off()

sector_order2 <- c(region_names2, schaef_nets)

svg("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/neg_avg_DMN_subnets_wholebrain_chord_diagram.svg", width = 20, height = 20)

chordDiagram(
  negative_connections,
  order = sector_order2,
  grid.col = colors,
  annotationTrack = c("name", "grid"),
  transparency = 0.2,
)

dev.off()


## Create chord plot with subcortex grouped together
avg_fc_dmn_7_subcort <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\lm_model\\avg_7nets_fc_dmn_subcort.csv", header = FALSE, sep = ",")

# load in sub-nework IDsschaef_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/schaef_names.csv", header = TRUE, stringsAsFactors = FALSE)

schaef_7nets <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/schaef_7_dmn_networks.csv", header = TRUE, stringsAsFactors = FALSE)

# remove any brainstem nuclei -
avg_fc_dmn_subcort <- avg_fc_dmn_7_subcort[-c(18:27), ]
# load in sub-net names
schaef_7nets2 <- data.frame(schaef_7nets[-c(24:33), ])
schaef_7nets2 <- data.frame(schaef_7nets2[-c(17:22), ])
colnames(schaef_7nets2) <- "Network"
# Change row 16 to "Subcortex"
schaef_7nets2[16, "Network"] <- "Subcortex"
# rename rows and cols
rownames(avg_fc_dmn_subcort) <- schaef_7nets2$Network
colnames(avg_fc_dmn_subcort) <- region_names2
# makes unique row names for the schaef rows
# unique_schaef_names <- make.unique(schaef_nets)
# connections <- fc_dmn_whole_brain
connections <- avg_fc_dmn_subcort
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = schaef_7nets2$Network[non_zero_connections[, 1]],
  to = region_names2[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)

colors <- c(
  IPL = "#00a79cdc",
  Dorsal_pFC = "#5b5da9",
  Medial_pFC = "#9799d8",
  LatVent_pFC = "#279e96cc",
  pCun_PCC = "#0681ba",
  Temp = "#5b5ea9a4",
  Rsp = "#00a79c67",
  Para_Hippo = "#5ac0b9d3",
  Visual = "#71a95ba4",
  Somato-Motor ="#71a95ba4",
  Dors-Atten = "#71a95ba4",
  SalVentAtten = "#71a95ba4",
  Limbic = "#71a95ba4",
  Cognitive = "#71a95ba4",
  TempPar = "#71a95ba4",
  Subcortex = "#71a95ba4",
  Cerebellum = "#71a95ba4"
  )

color2 <- c(
  Visual =  "#00a79cdc",
  Somato-Motor ="#71a95ba4",
  Dors-Atten = "#71a95ba4",
  SalVentAtten = "#71a95ba4",
  Limbic = "#71a95ba4",
  Cognitive = "#71a95ba4",
  TempPar = "#71a95ba4",
  Subcortex = "#71a95ba4",
  Cerebellum = "#71a95ba4"
)
# Define the order of the sectors
sector_order2 <- c(region_names2, schaef_7nets2$Network)

svg("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/Figures/DMN_Figures/DMN_FC_whole_brain/avg_nets_chord_diagram.svg", width = 20, height = 20)

chordDiagram(
  connections,
  order = sector_order2,
  grid.col = colors,
  annotationTrack = c("name", "grid"),
  transparency = 0.2,
)

dev.off()