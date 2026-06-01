# Using Chord Diagram
# install.packages("circlize")
# install.packages("scales")
library(circlize)
library(scales)
library(colorspace)

## Unique to your data
#set filepath for outputs
filepath = "/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/Analysis_2/"

# Set the working directory
setwd("/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/Analysis_2")
# Load the CSV file of fc matrix ordered
data <- read.csv("sig_diff_fc_grp2_grp0.csv", header = FALSE, sep = ",")

schaef_names <- read.csv("/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/Code/figure_code/schaef_17_nets.csv", header = TRUE, stringsAsFactors = FALSE)


colnames(data) <- schaef_names$Networks
rownames(data) <- schaef_names$ROI.Name



#chordDiagram(data) # plots the diagram as is

# Restart circular layout parameters
circos.clear()

# attempt the plot without the 0 weights
# Filter out zero-weighted connections
non_zero_connections <- which(data != 0, arr.ind = TRUE)

# Create a data frame of connections with non-zero weights
connections <- data.frame(
  from = non_zero_connections[, 1],
  to = non_zero_connections[, 2],
  weight = data[non_zero_connections]
)



non_zero_connections <- which(connections != 0, arr.ind = TRUE)

region_names <- schaef_names$Networks
# unique_region_names <- make.unique(region_names)
## Creates the plot
connectivity_matrix <- data
connectivity_matrix[lower.tri(connectivity_matrix)] <- 0
# Filter out zero-weighted connections
non_zero_connections <- which(connectivity_matrix != 0, arr.ind = TRUE)


# Create a data frame of connections with non-zero weights
connections <- data.frame(
  from = region_names[non_zero_connections[, 1]],
  to = region_names[non_zero_connections[, 2]],
  weight = connectivity_matrix[non_zero_connections]
)


# group list - 17 Networks
order <- c(
  grep("^Visual", region_names, value = TRUE),
  grep("^Somato-Motor", region_names, value = TRUE),
  grep("^Dors-Atten", region_names, value = TRUE),
  grep("^SalVentAtten", region_names, value = TRUE),
  grep("^Limbic", region_names, value = TRUE),
  grep("^Cont", region_names, value = TRUE),
  grep("^Default", region_names, value = TRUE),
  grep("^Temp_Parietal", region_names, value = TRUE),
  grep("^Subcortex", region_names, value = TRUE),
  grep("^Cerebellum", region_names, value = TRUE),
  grep("^AAS", region_names, value = TRUE)
)


# pastel colours
color_palette <- c("#F7A792", "#9ED9DC", "#ED7D87", "#BDDDA4", "#F1B6C5", 
                   "#99D2F2", "#BD93B8", "#ACB7DD", "#FFD368", "#86C5A6", "#E48F98")
#pastel colours for all 17 networks
color_palette <- c("#FAD2E1", "#F4B6C2",  # Soft pink  - Visual
                   "#D4E2D4", "#B8D8BE",  # Soft green  - Somato-motor
                   "#FBE4C6", "#FAD89C",  # Soft peach  - Dors-Atten
                    "#E2CFEA", "#CBA6DD",  # Soft purple - SalVentAtten
                   "#C6DEF1", "#A9C6E4",  # Soft blue  - Limbic
                   "#F7D9C4", "#EFC4A6", "#F4B89A",  # Soft coral  - Cont
                   "#F4A8A5", "#EFA2A1", "#E68C8A",  # Soft Red  - Default
                   "#C8E7E3",   # Extra soft teal - Temp_Parietal
                   "#8A98C8", #soft navy - subcortex
                   "#F6EAC2", #soft yellow - cerebellum
                    "#6B7AAE") #soft navy- AAS
       


# Create a named vector that maps each region to a color based on its group
# Create a named vector of colors for each unique name
unique_names <- unique(connections$from)
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
svg(file.path(filepath, "chord_diagram_17nets_all.svg"), width = 20, height = 20) # nolint

svg_file_path <- "/Users/ntaylor/Desktop/chord_diagram_17nets_all.svg"
svg("/Users/ntaylor/Desktop/chord_diagram_17nets_all.svg", width = 20, height = 20)
png("/Users/ntaylor/Desktop/chord_diagram_17nets_all.png", width = 800, height = 800)

# Plot the chord diagram
# Set cell padding to avoid summation error
#circos.par(cell.padding = c(0.02, 0, 0.02, 0))

chordDiagram(
  connections,
  grid.col = sector_colors,
  col = connection_colors,
  order = order,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
)
dev.off()

#positive light grey and negative dark
connection_colors2 <- sapply(seq_along(connections$weight), function(i) {
  weight <- connections$weight[i]
  if (weight > 0) {
     "#A9A9A9" #positive light grey
  } else {
    "#282828" #negative dark grey
  }
})
png("/Users/ntaylor/Desktop/chord_diagram_17nets_ligth_darkgrey.png", width = 800, height = 800)

chordDiagram(
  connections,
  grid.col = sector_colors,
  col = connection_colors2,
  order = order,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
)
dev.off()


## only plot positive connections


# Assuming 'connections' has numerical values and you want to keep rows where all values are positive
positive_connections <- connections[apply(connections > 0, 1, all), ]
png("/Users/ntaylor/Desktop/chord_diagram_17nets_pos_only.png", width = 800, height = 800)

chordDiagram(
  positive_connections,
  grid.col = sector_colors,
  order = order,
  col = connection_colors2,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
  )
  dev.off()
positive_connections <- connections[apply(connections > 0, 1, all), ]
# Set the color for positive connections to light grey
positive_connection_colors <- rep("#A9A9A9", nrow(positive_connections))

png("/Users/ntaylor/Desktop/chord_diagram_17nets_pos_only2.png", width = 800, height = 800)

chordDiagram(
  positive_connections,
  grid.col = sector_colors,
  order = order,
  col = positive_connection_colors,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
)
dev.off()



  # Assuming 'connections' has numerical values and you want to keep rows where all values are positive
  negative_connections <- connections[connections$weight < 0, ]
#make the connections darker version of existing colour

png("/Users/ntaylor/Desktop/chord_diagram_17nets_neg_only.png", width = 800, height = 800)
chordDiagram(
  negative_connections,
  grid.col = sector_colors,
  order = order,
  #col = connection_colors,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
  )
  dev.off()

negative_connections_colors <- rep("#282828", nrow(negative_connections))
png("/Users/ntaylor/Desktop/chord_diagram_17nets_neg_only2.png", width = 800, height = 800)
chordDiagram(
  negative_connections,
  grid.col = sector_colors,
  order = order,
  col = negative_connections_colors,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = 1,
  transparency = 0.2
  )
  dev.off()