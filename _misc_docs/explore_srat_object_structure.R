


# Reassign for simplicity
srat = small_srat

# Get Gene names
rownames(srat)

# Get Cell Names
colnames(srat)

# Access genes 1:100 and cells 1:5
srat[1:100, 1:5]

# Access genes 1:100 and cells 1:5
srat@meta.data

# Quick access a specific metadata column
srat$orig.ident
srat$nCount_RNA

# Access raw Counts
srat@assays$RNA$counts

# Access normalized data
srat@assays$RNA$data

# Access scaled data
srat@assays$RNA$scale.data

# Which cells are filtered
srat@assays$RNA@cells

# Which genes are filtered
srat@assays$RNA@features

# Metadata for finding variable features
srat@assays$RNA@meta.data

# Variable feature rank
srat@assays$RNA@meta.data$var.features.rank
