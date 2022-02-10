# compare sas geneig results at
# 
# https://support.sas.com/documentation/cdl/en/imlug/67502/HTML/default/viewer.htm#imlug_langref_sect159.htm
#
# to results from R's geigen
library(geigen)

A <- c(10,   2,   3,   1,   1,
        2,  12,   1,   2,   1,
        3,   1,  11,   1,  -1,
        1,   2,   1,   9,   1,
        1,   1,  -1,   1,  15)

B <- c(12,   1,  -1,   2,    1,
        1,  14,   1,  -1,    1,
       -1,   1,  16,  -1,    1,
        2,  -1,  -1,  12,   -1,
        1,   1,   1,  -1,   11)

A <- matrix(A, 5, 5)
B <- matrix(B, 5, 5)

# R's geigen is in reverse of SAS's geneig
g <- geigen(A, B)
rev(g$values)
g$vectors[,ncol(g$vectors):1]
# note that the first column matches output at the sas.com link above, but the
# remaining columns match with reversed sign


# try using the geneig function 
source("../multivarious-master/R/geneig.R")
source("../multivarious-master/R/projector.R")
source("../multivarious-master/R/all_generic.R")
source("../multivarious-master/R/pre_process.R")

library(Matrix)
g_robust <- geneig(A, B, ncomp=5, method="robust")
g_robust$values # same as line 25 above
matrix(g_robust$v@x, 5, 5)

g_sdiag <- geneig(A, B, ncomp=5, method="sdiag")
g_sdiag$values # same as line 25 above
g_sdiag$v # different set of wrong -1*

g_lapack <- geneig(A, B, ncomp=5, method="lapack")
g_lapack$values # same as line 25 above
g_lapack$v # different set of wrong -1*
