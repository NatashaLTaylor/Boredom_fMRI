## Chord Diagrams for within avg. FC
library(circlize)
library(scales)
library(colorspace)
library(grDevices)

# load data
data <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\order_sig_diff_del_health.csv", header = FALSE, sep = ",") # nolint

order_8_nets <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\grouping_11_nets_fc.csv", header = FALSE, sep = ",") # nolint

schaef_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/schaef_names.csv", header = TRUE, stringsAsFactors = FALSE)
sub_regions_networks <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/sub_regions_networks.csv", header = TRUE, stringsAsFactors = FALSE)
# order schaefer_names
order_vector <- order_8_nets[, 2]
# Reorder schaef_names based on the extracted order
schaef_names_ordered <- schaef_names[order_vector, ]
sub_regions_networks_ordered <- sub_regions_networks[order_vector, ]

# column names for within networks
colnames(avg_within_net_delirium) <- c("Visual", "Somato-motor", "Dorsal Attention", "Sal-Vent Attention", "Limbic", "Control", "Default", "Temp-parietal", "Subcortex", "Cerebellum", "AAS Nuclei")

rownames(data) <- sub_regions_networks_ordered$ROI.Name
colnames(data) <- sub_regions_networks_ordered$Name_nets
# specify connections of interest

# network name changes based on which network specify
network_name <- "Limbic"
# Filter columns where column names contain "Visual" (case-insensitive)
net_columns <- grep(network_name, colnames(data), ignore.case = TRUE)
net_columns_names <- grep(network_name, colnames(data), value = TRUE, ignore.case = TRUE)
# Filter rows where row names contain "Vis" (case-insensitive)
net_rows <- grep(network_name, rownames(data), ignore.case = TRUE)
net_rows_names <- grep(network_name, colnames(data), value = TRUE, ignore.case = TRUE)

# get the network locations of "Vis" from the data
# Subset the data frame
#filtered_data <- data[net_columns, net_columns]
# for default
filtered_data <- data[net_rows, net_columns]
connections <- filtered_data
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = net_rows_names[non_zero_connections[, 1]], # row and column names will be the sub-network names
  to = net_columns_names[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)

# Extract unique names

order <- c("SalVentAttn_ParOper", "SalVentAttn_FrOpe","SalVentAttn_OFC","SalVentAttn_Ins", "SalVentAttn_ParMed","SalVentAttn_FrMed", "SalVentAttn_PrC","SalVentAttn_PFCl", "SalVentAttn_IPL")
#order <- c("Cont Temp", "Cont IPS", "Cont IPL","Cont PFCd", "Cont PFCm","Cont PFCl/v", "Cont Cing","Cont pCun")

order <- c("Cont Temp", "Cont IPS", "Cont pCun","Cont IPL","Cont PFCd", "Cont PFCm","Cont PFCl/v", "Cont Cing")

order <- c("DMN IPL", "DMN PFCd", "DMN PFCm", "DMN PFCl/v", "DMN PCC", "DMN Temp", "DMN Rsp", "DMN Para-hippo")
order <- c("VisCent_ExStr", "VisPeri_ExtrSup", "VisPeri_StriCal","VisPeri_ExtrInf")
order <-c("SomMotA","SomMotB_Ins","SomMotB_S2","SomMotB_Aud","SomMotB_Cent")
order <-c("DorsAttn_TempOcc","DorsAttn_ParOcc","DorsAttnA_SPL","DorsAttnB_FEF","DorsAttnB_PostC")


# Define a vector of colors (you can customize this as needed)
#color_palette <- c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f07882", "#f95d6a", "#ff7c43", "#ffa600")
# pastel colours
color_palette <- c("#F7A792", "#9ED9DC", "#ED7D87", "#BDDDA4", "#f1b6c5", "#99D2F2", "#BD93B8", "#ACB7DD", "#FFD368")

# color_palette <- c("#003f5c", "#7a5195", "#ff6361", "#ffa600")
# Create a named vector of colors for each unique name
unique_names <- unique(net_columns_names)
# Create a vector of colors for the connections based on their weights
sector_colors <- setNames(color_palette[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
     "#A9A9A9" #positive light grey
  } else {
    "#282828" #negative dark grey
  }
})

circos.clear()
# Open an SVG graphics device - this creates figure with bolder outer colors
filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/4chord_diagram_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  annotationTrackHeight = c(0.05, 0.05), 
  preAllocateTracks = 1,
  transparency = 0.5,
  grid.col = sector_colors,
  col = connection_colors,
  #order = order
)

dev.off()


## Colors the sections based on the connections
filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/3chord_diagram_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.5,
  grid.col = sector_colors,
  #col = connection_colors,
  #order = order
)

dev.off()



# network name changes based on which network specify
network_name <- "DMN"
# Filter columns where column names contain "Visual" (case-insensitive)
net_columns_names <- grep(network_name, sub_regions_networks_ordered$Name_nets, value = TRUE, ignore.case = TRUE)
# Filter rows where row names contain "Vis" (case-insensitive)

connections <- filtered_data
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = net_columns_names[non_zero_connections[, 1]], # row and column names will be the sub-network names
  to = net_columns_names[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)
# Replace DMNpCun_PCC with DMN PCC
connections$to <- gsub("DMN pCun_PCC", "DMN PCC", connections$to, ignore.case = TRUE)
connections$from <- gsub("DMN pCun_PCC", "DMN PCC", connections$from, ignore.case = TRUE)
# Extract unique names
unique_names <- unique(connections$from)
order <- c("DMN IPL", "DMN PFCd", "DMN PFCm", "DMN PFCl/v", "DMN pCun_PCC", "DMN PCC", "DMN Temp", "DMN Rsp", "DMN Para-hippo")
order <- c("DMN IPL", "DMN PFCd", "DMN PFCm", "DMN PFCl/v", "DMN PCC", "DMN Temp", "DMN Rsp", "DMN Para-hippo")
# Define a vector of colors (you can customize this as needed)
# pastel colours
color_palette <- c("#F7A792", "#9ED9DC", "#ED7D87", "#BDDDA4", "#f1b6c5", "#99D2F2", "#BD93B8", "#ACB7DD", "#FFD368")

# color_palette <- c("#003f5c", "#7a5195", "#ff6361", "#ffa600")
# Create a named vector of colors for each unique name
sector_colors <- setNames(color_palette[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
    sector_colors[connections$to[i]]
  } else {
    darken(sector_colors[connections$to[i]], amount = 0.4)
  }
})
circos.clear()
# Open an SVG graphics device
filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/chord_diagram_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.5,
  grid.col = sector_colors,
  col = connection_colors,
  order = order
)

dev.off()


# temp parietal is a bit different to the others

rownames(data) <- schaef_names_ordered$ROI.Name
colnames(data) <- schaef_names_ordered$Networks
# specify connections of interest

# network name changes based on which network specify
network_name <- "Temp_Parietal"
network_name2 <- "TempPar"
# Filter columns where column names contain "Visual" (case-insensitive)
net_columns <- grep(network_name, colnames(data), ignore.case = TRUE)
net_columns_names <- grep(network_name2, colnames(data), value = TRUE, ignore.case = TRUE)
# Filter rows where row names contain "Vis" (case-insensitive)
net_rows <- grep(network_name2, rownames(data), ignore.case = TRUE)
net_rows_names <- grep(network_name2, rownames(data), value = TRUE, ignore.case = TRUE)

# get the network locations of "Vis" from the data
# Subset the data frame
filtered_data <- data[net_rows, net_rows]
# for default
# filtered_data <- data[net_columns, net_columns]
connections <- filtered_data
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = net_rows_names[non_zero_connections[, 1]], # row and column names will be the sub-network names
  to = net_rows_names[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)

color_palette <- c("#F7A792", "#9ED9DC", "#ED7D87", "#BDDDA4", "#f1b6c5", "#99D2F2", "#BD93B8", "#ACB7DD", "#f95d6a")

# color_palette <- c("#003f5c", "#7a5195", "#ff6361", "#ffa600")
# Create a named vector of colors for each unique name
unique_names <- unique(connections$from)
sector_colors <- setNames(color_palette[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(connections$weight, function(weight) {
  if (weight > 0) {
    sector_colors[connections$to]
  } else {
    darken(sector_colors[connections$to], amount = 0.4)
  }
})

circos.clear()
# Open an SVG graphics device
filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/2chord_diagram_%s.svg", network_name2)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.5,
  grid.col = sector_colors,
  #col = connection_colors,
  # order = order
)

dev.off()




## Chord Diagram from network to whole brain regions

# load data in
data <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\order_diff_salvent_wholebrain.csv", header = FALSE, sep = ",") # nolint

# network name changes based on which network specify
sub_regions_networks <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\Chord_Diagrams\\sub_regions_networks.csv", header = TRUE, sep = ",") # nolint
sub_regions_networks_ordered <- sub_regions_networks[order_vector, ]

# specify connections of interest
data_network <- subset(sub_regions_networks, Order_8_nets == 4)

net_names <- data_network$Name_nets

# needs to be name of networks
rownames(data) <- schaef_names_ordered$ROI.Name[1:400]
# network names for cortex
colnames(data) <- net_names
cortex_nets <- sub_regions_networks_ordered$Cort_nets[1:400]

# filtered_data <- data[net_columns, net_columns]
connections <- data
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = cortex_nets[non_zero_connections[, 1]], # row and column names will be the sub-network names
  to = net_names[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)

# order of DMN
order <- c(
  "SalVentAttn_ParOper", "SalVentAttn_FrOpe", "SalVentAttn_OFC", "SalVentAttn_Ins",
  "SalVentAttn_ParMed", "SalVentAttn_FrMed", "SalVentAttn_PrC", "SalVentAttn_PFCl", "SalVentAttn_IPL", "Visual", "Somato-motor", "Dors-Atten", "Control", "Limbic", "DMN", "TempPar"
)


# specify colours cortex networks
color_palette2 <- c(
  Visual = "#F7A792",
  `Somato-motor` = "#9ED9DC",
  `Dors-Atten` = "#ED7D87",
  Control = "#BDDDA4",
  Limbic = "#f1b6c5",
  DMN = "#BD93B8",
  TempPar = "#FFD368"
)
# specify colour for salvent
colors_palette1 <- c(
  SalVentAttn_ParMed = "#00a79cdc",
  SalVentAttn_ParOper = "#5b5da9",
  SalVentAttn_FrOpe = "#9799d8",
  SalVentAttn_FrMed = "#279e96cc",
  SalVentAttn_PFCl = "#0681ba",
  SalVentAttn_IPL = "#5b5ea9a4",
  SalVentAttn_Ins = "#00a79c67",
  SalVentAttn_PrC = "#5ac0b9d3",
  SalVentAttn_OFC = "#00a79c"
)
# combine color together
color_combine <- c(colors_palette1, color_palette2)

unique_names <- unique(connections$from)

sector_colors <- setNames(color_combine[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
    color_combine[connections$from[i]]
  } else {
    darken(color_combine[connections$from[i]], amount = 0.4)
  }
})
# make it dark grey lines
sector_colors <- setNames(color_palette[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
     "#A9A9A9" #positive light grey
  } else {
    "#282828" #negative dark grey
  }
})

circos.clear()


filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/pos_neg_whole_brain_chord_diagram2_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  #annotationTrackHeight = c(0.05, 0.05), 
  transparency = 0.5,
  grid.col = color_combine,
  col = connection_colors,
  order = order
)

dev.off()

filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/pos_neg_whole_brain_chord_diagram4_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  annotationTrackHeight = c(0.05, 0.05), 
  transparency = 0.5,
  grid.col = color_combine,
  col = connection_colors,
  order = order
)

dev.off()


### DMN Chord Plots
# load data in

data <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\order_diff_dmn_wholebrain.csv", header = FALSE, sep = ",") # nolint
sub_net_names <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/sub_networks_DMN.csv", header = TRUE, stringsAsFactors = FALSE)
network_name <- "DMN"
# network name changes based on which network specify
colnames(data) <- sub_net_names$ROI_subs
names_columns <- sub_net_names$ROI_subs
# specify connections of interest
# needs to be name of networks
rownames(data) <- make.unique(schaef_names_ordered$Networks[1:400])
# network names for cortex
cortex_nets <- schaef_names_ordered$Networks[1:400]


# filtered_data <- data[net_columns, net_columns]
connections <- data
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = cortex_nets[non_zero_connections[, 1]], # row and column names will be the sub-network names
  to = names_columns[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)

connections$from <- gsub("Rsp", "RSP", connections$from, ignore.case = TRUE)
connections$to <- gsub("Rsp", "RSP", connections$to, ignore.case = TRUE)

# order of DMN
# order <- c("DMN IPL", "DMN PFCd", "DMN PFCm", "DMN PFCl/v", "DMN pCun_PCC", "DMN PCC", "DMN Temp", "DMN Rsp", "DMN Para-hippo")

order <- c(
  "DMN IPL", "DMN PFCl/v", "DMN PFCd", "DMN PFCm", "DMN Temp", "DMN PCC",
  "DMN RSP", "DMN Para-hippo", "Visual", "Somato-motor", "SalVentAttn",
  "Dors-Atten", "Control", "Limbic", "TempPar"
)


# specify colour for salvent
colors_palette1 <- c(
  "DMN IPL" = "#00a79cdc",
  "DMN PFCd" = "#5b5da9",
  "DMN PFCm" = "#9799d8",
  "DMN PFCl/v" = "#279e96cc",
  "DMN PCC" = "#0681ba",
  "DMN Temp" = "#5b5ea9a4",
  "DMN RSP" = "#00a79c67",
  "DMN Para-hippo" = "#5ac0b9d3"
)
# specify colours cortex networks
color_palette2 <- c(
  "Visual" = "#F7A792",
  "Somato-motor" = "#9ED9DC",
  "Dors-Atten" = "#ED7D87",
  "Control" = "#BDDDA4",
  "Limbic" = "#f1b6c5",
  "SalVentAttn" = "#BD93B8",
  "TempPar" = "#FFD368"
)
color_combine <- c(colors_palette1, color_palette2)

# combine color together
unique_sectors <- unique(c(connections$from, connections$to))
unique_names <- unique(connections$from)

sector_colors <- setNames(color_combine[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
    color_combine[connections$from[i]]
  } else {
    darken(color_combine[connections$from[i]], amount = 0.4)
  }
})



circos.clear()
sector_colors <- setNames(color_palette[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
     "#A9A9A9" #positive light grey
  } else {
    "#282828" #negative dark grey
  }
})

circos.clear()


filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/pos_neg_whole_brain_chord_diagram2_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  #annotationTrackHeight = c(0.05, 0.05), 
  transparency = 0.5,
  grid.col = color_combine,
  col = connection_colors,
  order = order
)

dev.off()

filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/pos_neg_whole_brain_chord_diagram4_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  annotationTrackHeight = c(0.05, 0.05), 
  transparency = 0.5,
  grid.col = color_combine,
  col = connection_colors,
  order = order
)

dev.off()


### Visual Chord Plots
# load data in
data <- read.csv("C:\\Users\\natas\\OneDrive - The University of Sydney (Staff)\\Postdoc_Rob\\Analysis\\Graph_Theory\\schaef_400\\FC\\order_diff_visual_wholebrain.csv", header = FALSE, sep = ",") # nolint
network_name <- "Vis"
schaef_names$Networks[149:194] <- "Default"
schaef_names$Networks[358:390] <- "Default"
schaef_names_ordered <- schaef_names[order_vector, ]
# network name changes based on which network specify
colnames(data) <- sub_regions_networks_ordered$Name_nets[1:47]

rownames(data) <- sub_regions_networks_ordered$Cort_nets
cortex_nets <- sub_regions_networks_ordered$Cort_nets
cortex_nets[1:47] <- sub_regions_networks_ordered$Name_nets[1:47]
names_columns <- colnames(data)
# specify connections of interest
# needs to be name of networks
rownames(data) <- make.unique(cortex_nets)

# filtered_data <- data[net_columns, net_columns]
connections <- data
# remove zero connections
non_zero_connections <- which(connections != 0, arr.ind = TRUE)
connections <- data.frame(
  from = cortex_nets[non_zero_connections[, 1]], # row and column names will be the sub-network names
  to = names_columns[non_zero_connections[, 2]],
  weight = connections[non_zero_connections]
)

# order of visual network
order <- c(
  "VisCent_ExStr", "VisPeri_ExtrSup", "VisPeri_ExtrInf", "VisPeri_StriCal",
  "Somato-motor", "SalVentAttn", "Dors-Atten", "Control", "Limbic", "DMN", "TempPar"
)

# specify colour for salvent
colors_palette1 <- c(
  VisCent_ExStr = "#00a79cdc",
  VisPeri_ExtrSup = "#5b5da9",
  VisPeri_ExtrInf = "#9799d8",
  VisPeri_StriCal = "#0681ba"
)
# specify colours cortex networks
color_palette2 <- c(
  "Somato-motor" = "#9ED9DC",
  "Dors-Atten" = "#ED7D87",
  "Control" = "#BDDDA4",
  "Limbic" = "#f1b6c5",
  "SalVentAttn" = "#BD93B8",
  "DMN" = "#F7A792",
  "TempPar" = "#FFD368"
)
color_combine <- c(colors_palette1, color_palette2)

# combine color together
unique_sectors <- unique(c(connections$from, connections$to))
unique_names <- unique(connections$from)

connection_names <- unique(c(connections$from, connections$to))
missing_in_order <- setdiff(connection_names, order)
missing_in_connections <- setdiff(order, connection_names)


sector_colors <- setNames(color_combine[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
    color_combine[connections$from[i]]
  } else {
    darken(color_combine[connections$from[i]], amount = 0.4)
  }
})

# make it dark grey lines
sector_colors <- setNames(color_palette[1:length(unique_names)], unique_names)
# Create a vector of colors for the connections based on their weights
connection_colors <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
     "#A9A9A9" #positive light grey
  } else {
    "#282828" #negative dark grey
  }
})

circos.clear()


filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/pos_neg_whole_brain_chord_diagram2_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  #annotationTrackHeight = c(0.05, 0.05), 
  transparency = 0.5,
  grid.col = color_combine,
  col = connection_colors,
  order = order
)

dev.off()

filename <- sprintf("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/Chord_Diagrams/pos_neg_whole_brain_chord_diagram4_%s.svg", network_name)
svg(filename, width = 20, height = 20)
# Plot the chord diagram
chordDiagram(
  connections,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  annotationTrackHeight = c(0.05, 0.05), 
  transparency = 0.5,
  grid.col = color_combine,
  col = connection_colors,
  order = order
)

dev.off()