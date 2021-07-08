# je suis parti de ce post
# https://darrenjw.wordpress.com/2012/11/20/getting-started-with-bayesian-variable-selection-using-jags-and-rjags/

# je n'ai utilisé que le code de la section "Basic variable selection" c'est a dire que
# j ai laisse tombé les modifs sur le code que l'auteur propose dans les sections
# "Variable selection with random effects" et "Variable selection with random effects and a prior on the inclusion probability"

# j ai rajoute des lignes de code a la fin pour calculer les probabilites a posteriori des modeles 
# (voir par exemple la table 2 dans https://dl.dropboxusercontent.com/u/23160641/my-pubs/Gimenezetal2009Evolution.pdf)

library(rjags)

################################################
# ETAPE 1 : on simule des donnees 
# regression lineaire 
# 5 variables explicatives
# et 500 observations
################################################

# nb obs
n <- 200
# nb var explicatives
p <- 5
# genere les p variables dans une loi normale
X <- matrix(rnorm(n*p),
            ncol = p)

# genere les coefficients de regression associes aux 5 variables 
beta <- 4^(0:(1-p))
print(beta)
plot(beta) # on voit que la significativite des coefficients est decroissante ; 
# a partir du coef 2, les beta se rapprochent de 0

#  ordonnee a l origine
alpha <- 3

# precision sur la residuelle
tau <- 2
# le vecteur de residus
eps <- rnorm(n,0,1/sqrt(tau))

# les observations
y <- alpha + as.vector(X %*% beta + eps)
length(y)

################################################
########### ETAPE 2 : ajustement ##################
################################################

#------------- methode classique
mod <- lm(y~X)
print(summary(mod)) # on retrouve le fait que seuls les premiers beta sont significatifs



# methode bayesienne sans la selection de variable
require(rjags)
data <- list(y = y , X = X, n = n, p = 5)
init <- list(tau = 1, alpha = 0, beta = rep(0,5))
modelstring="
model{
for (i in 1:n){
mean[i]<-alpha + beta[1]*X[i,1] + beta[2]*X[i,2] + beta[3]*X[i,3] + beta[4]*X[i,4] + beta[5]*X[i,5]
y[i]~dnorm(mean[i],tau)
}
for(j in 1:p){
beta[j]~dnorm(0,0.001)
}
alpha~dnorm(0,0.0001)
tau~dgamma(1,0.001)
}
"
model <- jags.model(textConnection(modelstring), 
                    data = data,
                    inits = init)
update(model, n.iter = 100)
output <- coda.samples(model = model,
                       variable.names = c("alpha", "beta", "tau"), 
                       n.iter = 10000,
                       thin = 1)
print(summary(output))
plot(output)

# on compare les estimations via le max de vrais versus les moyennes a posteriori
cbind(summary(mod)$coef[1:p,1], summary(output)[[1]][1:p,1]) # les estimations sont tres proches

#------------- approche bayesienne avec selection de variables a la Kuo et Mallick
#------------- et calcul de la probabilite des modeles a posteriori
# voir le lien au debut de ce script pour des explications detaillees
# voir aussi la section 2.4 de l article http://ba.stat.cmu.edu/journal/2009/vol04/issue01/ohara.pdf et plus specifiquement 
# la sous-section consacree a la methode de Kuo et Mallick


data <- list(y = y, X = X, n = n, p = 5)
init <- list(tau = 1, alpha = 0, beta = rep(0,p), ind = rep(0,p))
modelstring="
model {
for (i in 1:n){
mean[i]<- alpha + ind[1] * beta[1]*X[i,1] + 
           ind[2] * beta[2]*X[i,2] + 
           ind[3] * beta[3]*X[i,3] + 
           ind[4] * beta[4]*X[i,4] + 
           ind[5] * beta[5]*X[i,5]
y[i]~dnorm(mean[i],tau)
}

for (j in 1:p){
# model label
ind[j]~dbern(.2)
# prior for regression coef (but intercept)
beta[j]~dnorm(0,0.0001)
# Create indicators for the possible variables in the model
	TempIndicator[j]<-ind[j]*pow(2, j-1) 
}

# prior for intercept
alpha~dnorm(0,0.0001)

# prior for precision of the residual term
tau~dgamma(1,0.001)

#Create a model number for each possible model
mdl<- 1+sum(TempIndicator[]) 

# calculate the percentage of time each model is selected (2^(nb of covariates) = 2^5 = 32)
for (j in 1 : 32){ 
pmdl[j]<-equals(mdl, j) 
}
}
"
model <- jags.model(textConnection(modelstring), data = data, inits = init)
update(model, n.iter = 1000)
output <- coda.samples(model = model, 
                       variable.names = c("alpha","beta","ind","tau","pmdl"), 
                       n.iter = 10000,
                       thin = 1)
#print(summary(output))
#plot(output)

# fonction qui genere toutes les combinaisons possibles de covariables 
# (1 dans la col j si la cov j est presente, 0 sinon)
ind <- function(p){
  if(p == 0) {return(t <- 0)}
  else if(p == 1) {return(t <- rbind(0, 1))}
  else if(p == 2) {return(t <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1)))}
  else {
    t <- rbind(cbind(ind(p - 1), rep(0, 2^(p - 1))), 
               cbind(ind(p - 1), rep(1, 2^(p - 1))))
    return(t)
  }
}

# cree toutes les combi possibles si p=5 covariables
# la premiere ligne de mat.modele correspond au modele avec intercept seulement
# la derniere ligne de mat.modele correspond au modele avec toutes les covariables
mat.modeles <- ind(5)
mat.modeles

output <- as.mcmc(output)

# on recupere les probas a posterior de tous les modeles (moyennes)
# on cherche dans output ou sont les pmdl
output[1,]
grep("pmdl", names(output[1,]))
pmdl <- output[,grep("pmdl", names(output[1,]))]
pmp <- apply(pmdl,2,mean) # probabilite a posteriori des modeles

# classe les proba a posteriori des modeles de la plus grande a la plus petite
ii <- order(pmp, decreasing = T)

# affiche les modeles
res <- cbind(mat.modeles[ii,], pmp[ii])
res

# on voit clairement que le premier modele qui inclut les variables X1 et X2 est dominant, 
# avec une prob a posteriori tres tres forte (regarder la derniere colonne)
# dans ce cas-la, il suffit de faire retourner le modele avec ces 2 variables seulement
# pour avoir l estimation de leurs effets.