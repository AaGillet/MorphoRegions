#' Measurements from the vertebral column of an alligator
#'
#' Linear and angular measurements from Alligator mississipiensis MCZ 81457
#' @docType data
#' @format A matrix with 22 vertebrae and 19 measurements. Column 1, vertebra, is the positional
#' information.
#'
#' @keywords datasets
"alligator"

#' Measurements from the vertebral column of a mouse
#'
#' Linear and angular measurements from Mus musculus MCZ 59560
#' @docType data
#' @format A matrix with 23 vertebrae and 19 measurements. Column 1, vertebra, is the positional
#' information.
#'
#' @keywords datasets
"musm"

#' Measurements from the vertebral column of a dolphin
#'
#' Linear and angular measurements from Platanista gangetica SMNS 45653
#' @docType data
#' @format A matrix with 40 vertebrae and 16 measurements. Column 1, vertebra, is the positional
#' information.
#'
#' @keywords datasets
"dolphin"

#' Measurements from the vertebral column of three porpoises
#'
#' Linear and angular measurements from Phocoena phocoena NRM 815072 (`porpoise1`), NRM 835011 (`porpoise2`), and NRM 855083 (`porpoise3`).
#' @docType data
#' @usage data("porpoise")
#' @format Each is a data frame with 58, 56, and 59 vertebrae (respectively) and 16 measurements. Column 1, Vertebra, contains the positional information.
#'
#' @keywords datasets
#' @references Gillet, A., Frederich, B., Pierce, S. E., & Parmentier, E. (2022). Iterative habitat transitions are associated with morphological convergence of the backbone in delphinoids. Journal of Mammalian Evolution, 29(4), 931-946.
#' @name porpoise
"porpoise1"

#' @rdname porpoise
#' @usage NULL
#' @format NULL
"porpoise2"

#' @rdname porpoise
#' @usage NULL
#' @format NULL
"porpoise3"

#' Geometric morphometric data and PC scores of the vertebral column a seal
#'
#' Procrusted aligned 3d landmarks (`coord_gpa`) and PC scores (`scores`) of the backbone of Phoca vitulina (ZMUC 160).
#' @docType data
#' @usage data("seal_gmm")
#' @format A list containing: a 3D array of Procrustes aligned 3D landmarks (`coord_gpa`) of 25 presacral vertebrae with 40 landmarks per vertebra, a matrix of PC scores from a PCA performed on Procrustes aligned landmarks (`scores`),
#' and a numeric vector with the eigenvalues from the PCA (`eigenvals`).
#'
#' @keywords datasets
#' @references Esteban, J. M., Martín-Serra, A., Pérez-Ramos, A., Mulot, B., Jones K. E., & Figueirido, B.  (2023). The impact of the land-to-sea transition on evolutionary integration and modularity of the pinniped backbone. Nature Communications, 6, 1141.
#' @name seal_gmm
"seal_gmm"

#' Geometric morphometric data and PC scores of the vertebral column two wolves
#'
#' Procrusted aligned 3d landmarks (`coord_gpa`) and PC scores (`scores`) from  10 subsampled vertebrae along the backbone of two specimens of Canis lupus (MU 071 and MU 073).
#' @docType data
#' @usage data("wolves_gmm")
#' @format A list containing: a 3D array of Procrustes aligned 3D landmarks (`coord_gpa`) of 10 subsampled vertebrae per specimen, and 30 landmarks per vertebra, a matrix of PC scores from a PCA performed on Procrustes aligned landmarks (`scores`),
#' a numeric vector with the eigenvalues from the PCA (`eigenvals`), a numeric vector with the position information of vertebrae (`vertebrae`), and a vector with the names of the specimens (`specimens`).
#'
#' @keywords datasets
#' @references Schwab, J., Figueirido, B., & Jones, K. E. (2026). Ecological Inference From Isolated Vertebrae: Evaluating Functional Signal Across the Carnivoran Spine. Journal of Morphology, 287(1), e70109.
#' @name wolves_gmm
"wolves_gmm"

