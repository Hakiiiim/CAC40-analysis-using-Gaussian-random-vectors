install.packages("RVAideMemoire")

library(RVAideMemoire)
library(datasets) 
library(ggplot2)

# Script TP2 
# Corrigé partiel exercice 1
# ----------
# EXERCISE 1 
# ----------


# Partie 1 - mu = (1 2), sigma1 = sigma2 = 1 et rho = 0.8
# On utilise la décomposition de Cholesky explicite vue en cours de la matrice
# de covariance du vecteur X - voir compte-rendu pour les détails

# Méthode 1 : pas de boucle for sur le numéro i de simulation, directement par
# calcul matriciel élémentaire

n <- 100 # taille de l'échantillon (les individus)

# initialisation de la matrice des données simulées : n = 100 individus repérés 
# par p = 2 deux variables en colonne --> c'est l'usage en Statistique 
X <- matrix(0,nrow= n,ncol=2) 

eps1 <- rnorm(n,mean=0,sd=1) # bruit blanc N(0,1) de taille n (composantes sur l'axe x1)
X[,1] <- 1 + eps1
eps2 <- rnorm(n,mean=0,sd=1) # bruit blanc N(0,1) de taille n (composantes sur l'axe x2)
X[,2] <- 2 + 0.8*eps1 + 0.6*eps2

# Visualisation de l'échantillon simulé de taille du vecteur X = (X1, X2) = nuage de n points

plot(X[,1],X[,2],asp=1,xlim=c(1-3,1+3),ylim=c(2-3,2+3),xlab="x1",ylab="x2")
abline(v=1,lty=2)
abline(h=2,lty=2)
abline(a=1,b=1,col="red")
abline(a=3,b=-1,col="red")
title("Echantillon de taille n = 100 d'un VG 2-dimensionnel (corrélation = 0.8)",cex.main=0.9)

# Méthode 2 : boucle for à éviter si n grand mais moins de risque d'erreur et programme plus lisible

n <- 100
Xbis <- matrix(0,nrow= n,ncol=2) 
for (i in 1:n){
  eps1 <- rnorm(1) # par défaut mean=0,sd=1
  Xbis[i,1] <- 1 + eps1
  Xbis[i,2] <- 2 + 0.8*eps1 + 0.6*rnorm(1)
}

points(Xbis[,1],Xbis[,2],pch="+",col="blue")

# Méthode 3 : on peut bien sûr ne pas utiliser l'expression analytique de la décomposition de Cholesky
# et l'obtenir numériquement avec la fonction chol() de R

Gamma <- matrix(0,nrow=2,ncol=2)
Gamma[1,1] <- 1
Gamma[2,2] <- 1
Gamma[1,2] <- 0.8
Gamma[2,1] <- 0.8
print(Gamma)

help(chol)
chol_Gamma <- chol(Gamma)
print(chol_Gamma)
S <- t(chol_Gamma) # attention à bien transposer si Gamma = S*t(S)
print(S) 

# ensuite méthode 1 ou 2 en utilisant les valeurs S[2,1] = 0.8 et S[2,2] = 0.6 

# Question 2 : décomposition spectrale de la matrice de covariance du vecteur X, lien avec la  
# décomposition de Mahalanobis et l'Analyse en Composantes Principales

# On utilise l'expression analytique de la décomposition spectrale, voir le CR pour les détails
lambda1 <- 1 + 0.8 # 1 + rho
lambda2 <- 1 - 0.8 # 1 - rho
U1 <- (1/sqrt(2))*matrix(c(1,1),n=2,ncol=1) # premier vecteur propre en matrice colonne
U2 <- (1/sqrt(2))*matrix(c(-1,1),n=2,ncol=1) # second vecteur propre en matrice colonne 

# mu + (U1, U2) définit une nouvelle base orthonormée orientée positivement
# Nouvelles coordonnées qui correspondent aux deux axes orthonormés précédents y = x + 1 et y = -x + 3

# Il faut d'abord centrer les données
Xc <- matrix(0,nrow= n,ncol=2) 
Xc[,1] <- X[,1] - 1
Xc[,2] <- X[,2] - 2

# Calcul des nouvelles coordonnées
C1 <- Xc%*%U1 # produit matriciel avec %*%, attention au fait que les individus sont en ligne
C2 <- Xc%*%U2

plot(C1,C2,xlab="u1",ylab="u2",asp=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
title("Nouvelles variables indépendantes, de variances lambda1 = 1.8 et lambda2 = 0.2",cex.main=0.9)

# Vérification empirique sur les données simulées

# variances estimées 
var1 <- var(C1)
print(var1) # à comparer avec lambda1 = 1.8 variance théorique
var2 <- var(C2)
print(var2) # à comparer avec lambda2 = 0.2

# normalité des composantes
qqnorm(C1,main="Q-Q Plot C1");qqline(C1,probs = c(0.1, 0.9),col="red")
qqnorm(C2,main="Q-Q Plot C2");qqline(C2,probs = c(0.1, 0.9),col="red")

# Indépendance, c'est le graphique déjà obtenu et qui correspondrait à une ACP normée ou non
# puisque les variables initiales sont de variance 1 : axe 1 = axe de plus grande inertie ou de
# variance projetée 

plot(C1,C2,xlab="C1 (90%)",ylab="C2 (10%)",asp=1) # 90% = 100*1.8/(1.8+0.2) %
abline(h=0,lty=2)
abline(v=0,lty=2)
title("Nuage des individus décrits par les nouvelles variables",cex.main=0.9)

# Lien avec la décomposition de Mahalanobis

U <- cbind(U1,U2)
print(U)
sv <- c(sqrt(1.8),sqrt(0.2)) # valeurs singulières = écart-types des nouvelles composantes
SIGMA <- diag(sv)
print(SIGMA)

# Vérification :
print(U%*%SIGMA%*%SIGMA%*%t(U))

# simulation de X par la décomposition de Mahalanobis <-> simulation des nouvelles coordonnées indépendantes
# à éviter si la dimension d est grande (ici d = 2)

n <- 10000
Xsimu <- U%*%SIGMA%*%matrix(rnorm(2*n),nrow=2,ncol=n) # attention aux dimensions
Xsimu <- t(Xsimu) # cf. format standard d'un jeu de données en Stat
# Ne pas oublier de décentrer
Xsimu[,1] <- Xsimu[,1] + 1  
Xsimu[,2] <- Xsimu[,2] + 2 
Xsimu <- data.frame(Xsimu) # encore mieux = structure de données R de type data frame 

plot(Xsimu,asp=1);abline(v=1,lty=2);abline(h=2,lty=2)
pairs(Xsimu)
cor(Xsimu) # matrice de corrélation empirique

# Partie 2
# Below: Create a function to simulate a 2D Gaussian Vector X = (X1, X2)
# ------
# Inputs
# ------
#  mu   : a vector of size 2 giving the mean of X
#  rho  : a real number between [-1, 1] giving the correlation cor(X1,X2)
#  sig    : a vector of size 2 containing the standard deviation of X. Default is (1,1).
#  n    : an integer giving the sample size. Default is 1000.
#  plot : should we plot the result? Default is TRUE.
#  ...  : optional arguments (color, labels, etc.) to be passed to plot
# ------
# Output
# ------
# X     : A matrix of size nx2 containing a sample of size n from the Gaussian distribution of X


simu_VG <- function(mu, rho, sig = c(1,1), n = 1000, plot = TRUE, ...){
  
  # construction of the covariance matrix, such that :
  # cor(X1, X2) = rho, var(X1) = sig[1]^2, var(X2) = sig[2]^2
  
  Gamma <- matrix(0,nrow=2,ncol=2)
  Gamma[1,1] <- sig[1]*sig[1]
  Gamma[2,2] <- sig[2]*sig[2]
  Gamma[1,2] <- sig[1]*sig[2]*rho
  Gamma[2,1] <- Gamma[1,2] 
  
  # here simulation of a sample of size n drawn from N(mu, Gamma) 

  # -- MY CODE --
  
  X <- matrix(NA, n, ncol=2)   # the matrix of size nx2 containing the simulations
  
  # warning : to be computed out of the loop 
  chol_Gamma <- chol(Gamma)
  S <- t(chol_Gamma) # to be a lower matrix
  
  for (i in 1:n){
    eps <- rnorm(n=2,mean=0,sd=1) # a 2-dimensional Gaussian noise
    X[i, ] <- mu + S%*%eps   # a vector of length 2 (row or column, as you like)
  }
  
  # plot the results if argument 'plot' is equal to TRUE
  if (plot){
    par(mfrow = c(1,1))
    plot(X, asp=1, ...)  # asp = 1 --> same scale for the x and y axis
    abline(v = mu[1], h = mu[2])     
  }
    
  return(X)
}
    
# Run the function
X <- simu_VG(mu = c(1,2), rho = -0.8, sig=c(3,1), n = 1000)

# A vous de jouer et compléter


# ----------
# EXERCISE 2 
# ----------
# consider the points such that x1-h <= X1 <= x1+h

n <- 5000 # taille de l'échantillon (les individus)

# initialisation de la matrice des données simulées : n = 100 individus repérés 
# par p = 2 deux variables en colonne --> c'est l'usage en Statistique 
X <- matrix(0,nrow= n,ncol=2) 

eps1 <- rnorm(n,mean=0,sd=1) # bruit blanc N(0,1) de taille n (composantes sur l'axe x1)
X[,1] <- 1 + eps1
eps2 <- rnorm(n,mean=0,sd=1) # bruit blanc N(0,1) de taille n (composantes sur l'axe x2)
X[,2] <- 2 + 0.8*eps1 + 0.6*eps2

x1 <- 1

h <- 0.01

indices <- which(X[,1]>=x1-h & X[,1]<=x1+h)  # the indices of X that fulfill the condition

X2 <- X[indices, 2]

plot(X[,1],X[,2],asp=1,xlim=c(1-3,1+3),ylim=c(2-3,2+3),xlab="x1",ylab="x2")
abline(v=1,lty=2)
abline(h=2,lty=2)
abline(a=1,b=1,col="red")
abline(a=3,b=-1,col="red")
points(X[indices, 1], X[indices, 2], col="blue")
abline(v = c(x1-h, x1+h))

# study the distribution of X2
# Is it a normal one ?
# What can you say of its mean ? Its variance ?

qqnorm(X2,main="Q-Q Plot x2");qqline(X2,probs = c(0.1, 0.9),col="red")

mu <- mean(X2)
print(mu)
vari <- var(X2)
print(vari)

#print(var(X[,2]))

# -- MY CODE --

N <- length(X[,1])

esperance <- rep(0,N)
variance <- rep(0,N)

for (i in 1:N) {
  x1 <- X[i,1]
  
  h <- 0.1
  
  indices <- which(X[,1]>=x1-h & X[,1]<=x1+h)  # the indices of X that fulfill the condition
  
  X2 <- X[indices, 2]
  
  esperance[i] <- mean(X2)
  
  variance[i] <- var(X2)
} 

graphes <- data.frame(X[,1],esperance,variance)

p = ggplot() +
  geom_line(data = graphes, aes(x =X[,1], y =esperance , color = "Esperance"))+
  geom_line(data = graphes, aes(x =X[,1], y =variance , color = "Variance"))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)
plot(p)

# ----------
# EXERCISE 3 
# ----------

RENAULT_data <- read.csv("data CAC40/RENAULT_2019-09-24.txt", header=T, sep="", dec=".")
RENAULT <- RENAULT_data$clot
plot(RENAULT, type='l', ylab = "cours de clôture (en euros)", xlab="jours boursiers (du 24/09/2018 au 24/09/2019)")
title("Evolution de l'action RENAULT sur un an")

# Taux de hausse ou de baisse logarithmique

tx_RENAULT <- diff(log(RENAULT)) # diff pour la différence discrète log(S(t)) - log(S(t-1))
plot(100*tx_RENAULT, type='l',ylab="taux en %",xlab="Temps")
abline(h=0, col="red")
title("Taux de rendement logarithmiques action RENAULT sur un an (en %)",cex.main=0.9)


#1

SOGE_data <- read.csv("data CAC40/SOCIETEGENERALE_2019-09-24.txt", header=T, sep="", dec=".")
SOGE <- SOGE_data$clot

BNP_data <- read.csv("data CAC40/BNPPARIBAS_2019-09-24.txt", header=T, sep="", dec=".")
BNP <- BNP_data$clot

tx_SOGE <- diff(log(SOGE))
tx_BNP <- diff(log(BNP))

V <- cbind(tx_SOGE,tx_BNP)

graphes <- data.frame(1:255,100*tx_SOGE,100*tx_BNP)

p = ggplot() +
  geom_line(data = graphes, aes(x =1:255, y =100*tx_SOGE , color = "SOGE"))+
  geom_line(data = graphes, aes(x =1:255, y =100*tx_BNP , color = "BNP"))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)
plot(p)

# plot(100*tx_SOGE, type='l',ylab="taux en %",xlab="Temps", col="blue")
# lines(100*tx_BNP, type='l',ylab="taux en %",xlab="Temps", col="red")

plot(tx_SOGE,tx_BNP)

qqnorm(tx_SOGE,main="Q-Q Plot x2");qqline(tx_SOGE,probs = c(0.1, 0.9),col="red")
qqnorm(tx_BNP,main="Q-Q Plot x2");qqline(tx_BNP,probs = c(0.1, 0.9),col="red")

mqqnorm(V,main = "Multi-normal Q-Q Plot");

mshapiro.test(V)

# x <- rnorm(30)
# y <- rnorm(30)
# mshapiro.test(cbind(x,y))

n <- 500

mu_SOGE <- mean(tx_SOGE)
mu_BNP <- mean(tx_BNP)

var_SOGE <- var(tx_SOGE)
var_BNP <- var(tx_BNP)

covariance <- cov(tx_BNP,tx_SOGE)

mu <- c(mu_SOGE,mu_BNP)

Gamma <- cbind(c(var_SOGE,covariance),c(covariance,var_BNP))
print(Gamma)
print(det(Gamma))

correlation <- cor(tx_BNP,tx_SOGE)
print(correlation)

chol_Gamma <- chol(Gamma)
print(chol_Gamma)

S <- t(chol_Gamma)

X <- matrix(0,nrow= n,ncol=2) 

eps1 <- rnorm(n,mean=0,sd=1) # bruit blanc N(0,1) de taille n (composantes sur l'axe x1)
eps2 <- rnorm(n,mean=0,sd=1)
X[,1] <- mu[1] + S[1,1]*eps1+S[1,2]*eps2
 # bruit blanc N(0,1) de taille n (composantes sur l'axe x2)
X[,2] <- mu[2] + S[2,1]*eps1+S[2,2]*eps2

# Visualisation de l'échantillon simulé de taille du vecteur X = (X1, X2) = nuage de n points

plot(X[,1],X[,2],xlab="x1",ylab="x2")
points(tx_SOGE,tx_BNP,col="red")
#mqqnorm(V,main = "Multi-normal Q-Q Plot"); #qqline(V,probs = c(0.1, 0.9),col="red")

#2

print(covariance)

DASSAULT_data <- read.csv("data CAC40/DASSAULTSYSTEMES_2019-09-24.txt", header=T, sep="", dec=".")
DASSAULT <- DASSAULT_data$clot
tx_DASSAULT <- diff(log(DASSAULT))

LOREAL_data <- read.csv("data CAC40/LOREAL_2019-09-24.txt", header=T, sep="", dec=".")
LOREAL <- LOREAL_data$clot
tx_LOREAL <- diff(log(LOREAL))

VEOLIA_data <- read.csv("data CAC40/VEOLIAENVIRONNEM_2019-09-24.txt", header=T, sep="", dec=".")
VEOLIA <- VEOLIA_data$clot
tx_VEOLIA <- diff(log(VEOLIA))

VUITTON_data <- read.csv("data CAC40/LVMHMOETVUITTON_2019-09-24.txt", header=T, sep="", dec=".")
VUITTON <- VUITTON_data$clot
tx_VUITTON <- diff(log(VUITTON))

AIRL_data <- read.csv("data CAC40/AIRLIQUIDE_2019-09-24-1.txt", header=T, sep="", dec=".")
AIRL <- AIRL_data$clot
tx_AIRL <- diff(log(AIRL))

AXA_data <- read.csv("data CAC40/AXA_2019-09-24.txt", header=T, sep="", dec=".")
AXA <- AXA_data$clot
tx_AXA <- diff(log(AXA))

CAPGEMINI_data <- read.csv("data CAC40/CAPGEMINI_2019-09-24.txt", header=T, sep="", dec=".")
CAPGEMINI <- CAPGEMINI_data$clot
tx_CAPGEMINI <- diff(log(CAPGEMINI))

DANONE_data <- read.csv("data CAC40/DANONE_2019-09-24.txt", header=T, sep="", dec=".")
DANONE <- DANONE_data$clot
tx_DANONE <- diff(log(DANONE))

ORANGE_data <- read.csv("data CAC40/ORANGE_2019-09-24.txt", header=T, sep="", dec=".")
ORANGE <- ORANGE_data$clot
tx_ORANGE <- diff(log(ORANGE))

SCHNEIDEREL_data <- read.csv("data CAC40/SCHNEIDEREL_2019-09-24.txt", header=T, sep="", dec=".")
SCHNEIDEREL <- SCHNEIDEREL_data$clot
tx_SCHNEIDEREL <- diff(log(SCHNEIDEREL))

TOTAL_data <- read.csv("data CAC40/TOTAL_2019-09-24.txt", header=T, sep="", dec=".")
TOTAL <- TOTAL_data$clot
tx_TOTAL <- diff(log(TOTAL))

covariance <- cor(SOGE,BNP)
print(covariance)

V2 <- cbind(tx_SOGE,tx_DASSAULT)

#mqqnorm(V2,main = "Multi-normal Q-Q Plot"); #qqline(V2,probs = c(0.1, 0.9),col="red")

# plot(SOGE,ylim=c(0,400), type='l',ylab="taux en %",xlab="Temps", col="blue")
# lines(BNP, type='l',ylab="taux en %",xlab="Temps", col="red")
# lines(RENAULT, type='l',ylab="taux en %",xlab="Temps", col="green")
# lines(DASSAULT, type='l',ylab="taux en %",xlab="Temps", col="black")
# lines(LOREAL, type='l',ylab="taux en %",xlab="Temps", col="yellow")
# lines(VEOLIA, type='l',ylab="taux en %",xlab="Temps", col="pink")
# lines(VUITTON, type='l',ylab="taux en %",xlab="Temps", col="purple")

graphes <- data.frame(1:256,SOGE,BNP,RENAULT,DASSAULT,LOREAL,VEOLIA,VUITTON,AIRL,AXA,CAPGEMINI,DANONE,ORANGE,SCHNEIDEREL,TOTAL)

p = ggplot() +
  geom_line(data = graphes, aes(x =1:256, y =SOGE , color = "SOGE"))+
  geom_line(data = graphes, aes(x =1:256, y =BNP , color = "BNP"))+
  geom_line(data = graphes, aes(x =1:256, y =RENAULT , color = "RENAULT"))+
  geom_line(data = graphes, aes(x =1:256, y =DASSAULT , color = "DASSAULT"))+
  geom_line(data = graphes, aes(x =1:256, y =LOREAL , color = "LOREAL"))+
  geom_line(data = graphes, aes(x =1:256, y =VEOLIA , color = "VEOLIA"))+
  geom_line(data = graphes, aes(x =1:256, y =VUITTON , color = "VUITTON"))+
  geom_line(data = graphes, aes(x =1:256, y =AIRL , color = "AIRL"))+
  geom_line(data = graphes, aes(x =1:256, y =AXA , color = "AXA"))+
  geom_line(data = graphes, aes(x =1:256, y =CAPGEMINI , color = "CAPGEMINI"))+
  geom_line(data = graphes, aes(x =1:256, y =DANONE , color = "DANONE"))+
  geom_line(data = graphes, aes(x =1:256, y =ORANGE , color = "ORANGE"))+
  geom_line(data = graphes, aes(x =1:256, y =SCHNEIDEREL , color = "SCHNEIDEREL"))+
  geom_line(data = graphes, aes(x =1:256, y =TOTAL , color = "TOTAL"))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)+
  coord_cartesian(ylim =c(0, 100))
plot(p)

#3

Vchoisi <- cbind(tx_VUITTON,tx_LOREAL,tx_RENAULT,tx_DANONE,tx_TOTAL,tx_SOGE,tx_VEOLIA,tx_SCHNEIDEREL)

graphes <- data.frame(1:256,SOGE,RENAULT,LOREAL,VEOLIA,VUITTON,DANONE,SCHNEIDEREL,TOTAL)

p = ggplot() +
  geom_line(data = graphes, aes(x =1:256, y =SOGE , color = "SOGE"))+
  geom_line(data = graphes, aes(x =1:256, y =RENAULT , color = "RENAULT"))+
  geom_line(data = graphes, aes(x =1:256, y =LOREAL , color = "LOREAL"))+
  geom_line(data = graphes, aes(x =1:256, y =VEOLIA , color = "VEOLIA"))+
  geom_line(data = graphes, aes(x =1:256, y =VUITTON , color = "VUITTON"))+
  geom_line(data = graphes, aes(x =1:256, y =DANONE , color = "DANONE"))+
  geom_line(data = graphes, aes(x =1:256, y =SCHNEIDEREL , color = "SCHNEIDEREL"))+
  geom_line(data = graphes, aes(x =1:256, y =TOTAL , color = "TOTAL"))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)
  #coord_cartesian(ylim =c(0, 100))
plot(p)

mqqnorm(Vchoisi,main = "Multi-normal Q-Q Plot");

mshapiro.test(Vchoisi)

M <- cov(Vchoisi)
print(det(M))

decomp <- eigen(M, symmetric=TRUE)
values <- decomp$values
vectors <- decomp$vectors

Vf <- matrix(0,nrow= 255,ncol=8) 
Vf[,1] <- Vchoisi[,1] - mean(Vchoisi[,1])
Vf[,2] <- Vchoisi[,2] - mean(Vchoisi[,2])
Vf[,3] <- Vchoisi[,3] - mean(Vchoisi[,3])
Vf[,4] <- Vchoisi[,4] - mean(Vchoisi[,4])
Vf[,5] <- Vchoisi[,5] - mean(Vchoisi[,5])
Vf[,6] <- Vchoisi[,6] - mean(Vchoisi[,6])
Vf[,7] <- Vchoisi[,7] - mean(Vchoisi[,7])
Vf[,8] <- Vchoisi[,8] - mean(Vchoisi[,8])


# Calcul des nouvelles coordonnées
C1 <- Vf%*%vectors[,1] 
C2 <- Vf%*%vectors[,2]
C3 <- Vf%*%vectors[,3]
C4 <- Vf%*%vectors[,4]
C5 <- Vf%*%vectors[,5]
C6 <- Vf%*%vectors[,6]
C7 <- Vf%*%vectors[,7]
C8 <- Vf%*%vectors[,8]

Vnew <- cbind(C1,C2,C3,C4,C5,C6,C7,C8)
print(Vnew)

graphes <- data.frame(1:255,C1,C2,C3,C4,C5,C6,C7,C8)

p = ggplot() +
  geom_line(data = graphes, aes(x =1:255, y =C1 , color = "C1"))+
  geom_line(data = graphes, aes(x =1:255, y =C2 , color = "C2"))+
  geom_line(data = graphes, aes(x =1:255, y =C3 , color = "C3"))+
  geom_line(data = graphes, aes(x =1:255, y =C4 , color = "C4"))+
  geom_line(data = graphes, aes(x =1:255, y =C5 , color = "C5"))+
  geom_line(data = graphes, aes(x =1:255, y =C6 , color = "C6"))+
  geom_line(data = graphes, aes(x =1:255, y =C7 , color = "C7"))+
  geom_line(data = graphes, aes(x =1:255, y =C8 , color = "C8"))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)
  #coord_cartesian(ylim =c(0, 100))
plot(p)

sing <- sqrt(values)
print(sing)

print(summary(sing))

