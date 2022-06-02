library(OpenMx)
library(metaSEM)

# input sample size
n <- c(424,215,531,394,255,240,779,282)

# change data format
study1<-as.matrix(study1)
study2<-as.matrix(study2)
study3<-as.matrix(study3)
study4<-as.matrix(study4)
study5<-as.matrix(study5)
study6<-as.matrix(study6)
study7<-as.matrix(study7)
study8<-as.matrix(study8)

# change rownames
rownames(study1) <- c('EC','TC','A','C','R','IM')
rownames(study2) <- c('EC','TC','A','C','R','IM')
rownames(study3) <- c('EC','TC','A','C','R','IM')
rownames(study4) <- c('EC','TC','A','C','R','IM')
rownames(study5) <- c('EC','TC','A','C','R','IM')
rownames(study6) <- c('EC','TC','A','C','R','IM')
rownames(study7) <- c('EC','TC','A','C','R','IM')
rownames(study8) <- c('EC','TC','A','C','R','IM')

study<-list(study1,study2,study3,study4,study5,study6,study7,study8)
names(study) <- c('study1','study2','study3','study4','study5','study6','study7','study8')

data <- study
study <-list(study,n)
names(study) <- c('data','n')

# fixed effect model
fixed1 <- tssem1(study$data, study$n, method = "FEM")
summary(fixed1)

# pooled correlation matrix
coef(fixed1)

# random effect model
random1 <- tssem1(study$data, study$n, method="REM", RE.type="Diag")
summary(random1)

# pooled correlation matrix
vec2symMat(coef(random1, select="fixed"), diag=FALSE )

# Observed factor matrix
F1 <- create.Fmatrix(c(1, 1, 1, 1, 1, 1), as.mxMatrix=FALSE)
dimnames(F1) <- list(c('EC','TC','A','C','R','IM'),
                     c('EC','TC','A','C','R','IM'))
F1

# Regression coefficient matrix
A1 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               '-.2*A_EC','.2*A_TC',0,0,0,0,
               '-.2*C_EC','.2*C_TC',0,0,0,0,
               '-.2*R_EC','.2*R_TC',0,0,0,0,
               '-.2*EC_IM','.2*TC_IM','.2*A_IM','.2*C_IM','.2*R_IM',0),
             ncol = 6,byrow = TRUE)

A2 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               '-.2*A_EC','.2*A_TC',0,0,0,0,
               '-.2*C_EC','.2*C_TC',0,0,0,0,
               '-.2*R_EC','.2*R_TC',0,0,0,0,
               0,0,'.2*A_IM','.2*C_IM','.2*R_IM',0),
             ncol = 6,byrow = TRUE)

dimnames(A2) <- list(c('EC','TC','A','C','R','IM'),
                     c('EC','TC','A','C','R','IM'))
A1

# Varience and covarience matrix
S1 <- matrix(c(1,'-.2*cor1',0,0,0,0,
               '-.2*cor1',1,0,0,0,0,
               0,0,'.2*e1','.2*cor2','.2*co3',0,
               0,0,'.2*cor2','.2*e2','.2*cor4',0,
               0,0,'.2*co3','.2*cor4','.2*e3',0,
               0,0,0,0,0,'.2*e4'),
               ncol = 6,byrow = TRUE)

S2 <- matrix(c(1,'.2*cor1',0,0,0,0,
               '.2*cor1',1,0,0,0,0,
               0,0,'.2*e1',0,0,0,
               0,0,0,'.2*e2',0,0,
               0,0,0,0,'.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)

S3 <- matrix(c(1,0,0,0,0,0,
               0,1,0,0,0,0,
               0,0,'.2*e1','.2*cor2','.2*co3',0,
               0,0,'.2*cor2','.2*e2','.2*cor4',0,
               0,0,'.2*co3','.2*cor4','.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)
S4 <- matrix(c(1,'-.2*cor1',0,0,0,0,
               '-.2*cor1',1,0,0,0,0,
               0,0,'.2*e1',0,'.2*co3',0,
               0,0,0,'.2*e2','.2*cor4',0,
               0,0,'.2*co3','.2*cor4','.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)

S5 <- matrix(c(1,'-.2*cor1',0,0,0,0,
               '-.2*cor1',1,0,0,0,0,
               0,0,'.2*e1','.2*cor2',0,0,
               0,0,'.2*cor2','.2*e2','.2*cor4',0,
               0,0,0,'.2*cor4','.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)

S6 <- matrix(c(1,'-.2*cor1',0,0,0,0,
               '-.2*cor1',1,0,0,0,0,
               0,0,'.2*e1','.2*cor2','.2*co3',0,
               0,0,'.2*cor2','.2*e2',0,0,
               0,0,'.2*co3',0,'.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)

S7 <- matrix(c(1,'-.2*cor1',0,0,0,0,
               '-.2*cor1',1,0,0,0,0,
               0,0,'.2*e1','.2*cor2',0,0,
               0,0,'.2*cor2','.2*e2',0,0,
               0,0,0,0,'.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)

S8 <- matrix(c(1,'-.2*cor1',0,0,0,0,
               '-.2*cor1',1,0,0,0,0,
               0,0,'.2*e1',0,0,0,
               0,0,0,'.2*e2','.2*cor4',0,
               0,0,0,'.2*cor4','.2*e3',0,
               0,0,0,0,0,'.2*e4'),
             ncol = 6,byrow = TRUE)

dimnames(S8) <- list(c('EC','TC','A','C','R','IM'),
                     c('EC','TC','A','C','R','IM'))
S2

# Stage 2
random2 <- tssem2(random1, Amatrix=A2, Smatrix=S3, Fmatrix=F1,
                 model.name="Li")
summary(random2)

# Plot
plot(random2, whatLabels="path", edge.label.cex=0.8)
plot(random2, color="yellow")

