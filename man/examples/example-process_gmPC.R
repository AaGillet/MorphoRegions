# Load coordinates of aligned landmarks, PC scores, and eigenvalues; vertebra index as names of array elements of aligned landmarks
data("seal_gmm")

# Process PC scores and morphological data:
seal_pca <- process_gmPC(data=seal_gmm$coord_gpa,
                         pcscores=seal_gmm$scores,
                         eigenvals = seal_gmm$eigenvals)

# Process multiple datasets; specimen names in 'specimens' vector, vertebrae index in 'vertebrae' vector
data("wolves_gmm")

# Process PC scores and geometric morphometric dataset:
wolves_pca <- process_gmPC(data=wolves_gmm$coord_gpa,
                           pos=wolves_gmm$vertebrae,
                           pcscores=wolves_gmm$scores,
                           eigenvals = wolves_gmm$eigenvals,
                           specimens=wolves_gmm$specimens)