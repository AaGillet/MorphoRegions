# Load dataset; vertebral index in first column ('Vertebra' column)
data("dolphin")

# Compute PC scores with prcomp:
PCA <- prcomp(dolphin[-1], scale=T)

# Extract PC scores and eigenvalues:
PCA_scores <- as.data.frame(PCA$x)
PCA_eigenval <- PCA$sdev^2


# Process PC scores and morphological data:
pco_dolphin <- process_PC(data=dolphin, pos='Vertebra',
                          pcscores=PCA_scores, eigenvals=PCA_eigenval,
                          posPC = dolphin[,1])