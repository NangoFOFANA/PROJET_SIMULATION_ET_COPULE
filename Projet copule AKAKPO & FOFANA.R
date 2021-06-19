

## FOFANA Nango & AKAKPO Crépin 
#Projet de Simulation et Copule Master 2 Probabilité et Statistiques des nouvelles Données 

## Sujet: Problème de couveture par ensemble 


library(sets)  ## Package pour construire les ensembles
# Un ensemble de 13 éléments
N=13
U=set(1,2,3,4,5,6,7,8,9,10,11,12,13)

## definition des sous ensembble de l'ensemble s
n=13 # Nombre de sous ensemble

s1=set(1,2,3,4,5);s2=set(1,2,3,6,7)

s3=set(1,2,3,4,7,8);s4=set(1,3,4,5,8)

s5=set(1,4,5,8,9);s6=set(2,6,7,10,11)

s7=set(2,3,6,7,8,11,12);s8=set(3,4,5,7,8,9,12,13)

s9=set(5,8,9,13);s10=set(6,10,11)

s11=set(6,7,10,11,12,13);s12=set(7,8,11,12,13)

s13=set(8,9,11,12,13)

S=list(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13) ## list des sous ensemble 

#Fonction de constructions des sous ensemples 
D <- structure(list(sets = structure(c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L,
                3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 
                6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 
                8L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 11L, 11L, 11L, 11L, 11L, 11L, 12L, 
                12L, 12L, 12L, 12L, 13L, 13L, 13L, 13L, 13L), .Label = c("s1", "s2", 
                "s3", "s4","s5", "s6","s7","s8","s9","s10","s11","s12","s13"),
                class="factor"), n = c(1, 2, 3, 4, 5, 1, 2, 3, 6, 7, 1, 2, 3, 4, 7,
                8, 1, 3, 4, 5, 8, 1, 4, 5, 8, 9, 2, 6, 7, 10, 11, 2, 3, 6, 7, 8, 11,
                12, 3, 4, 5, 7, 8, 9, 12, 13, 5, 8, 9, 13, 6, 10, 11, 6, 7, 10, 11,
                12, 13, 7, 8, 11, 12, 13, 8, 9, 11, 12, 13)),.Names = c("sets", "n"),
               row.names = c(NA, 69L), class = "data.frame")

#D
## %Matrix des sous ensembles (tableaux à double entré)
(matrix_couverture=t(table(D$sets,D$n))) 

Y=rep(1,13) # Chaque sous ensemble répresenté par la valeur 1 
dimnames(matrix_couverture) = list( items=1:13, sets=paste0("s", 1:13) )
#dimnames(matrix_couverture)
## les poids associés au sous ensembles (aleatoire)
W=c(floor(runif(13,10,20)))
#W

M=1:13 ## la liste des sous ensemble 
##construction de la fonction H
nmax=400
H=rep(0,nmax) ##Stockage de la valeur H

# construction de la fonction de transition

fonction_transition=function(X,beta){
  Xnew=X  ## On copie l'ensemble du vecteur
  P=sample(M,2,replace = TRUE) ## Tirage aléatoire de deux sous ensemble 
  p1=P[1]
  p2=P[2]
  # permutation des sous ensembles 
  Xnew[p1]=M[p2]
  Xnew[p2]=M[p1]
  #N1=W[p1]*Y[p1]+W[[p2]]*Y[p2]
  #N2=W[p1]*Y[p1]+W[[p2]]*Y[p1]
  #
  cou=c(1,X) # le poid du sous ensemble avec X
  new_cou=c(1,Xnew) #  le poid du sous ensemble avec nouveau X
  
  N1=W[cou[p1]]+W[cou[p2]]
  N2=W[new_cou[p1]]+W[new_cou[p2]]
  
  DELTA=N1-N2# difference entre la fonction H(X) et H(Xnew)
  ALPH=exp(DELTA*beta)  ## calcul de la valeurs d alpha
  # définition de la fonction de réjet 
  v=runif(1)  
  if(ALPH<v){ 
    Xnew=X  ## On revient a X
    DELTA=0   ## H reste le meme
  }
  return(c(Xnew,DELTA))
}

cover_algorithme= function(nmax, bet, REMO){
  s<-sample(M,13,replace=T) ## tirage aléatoires de tous les sous ensembles 
  
  H=rep(0,nmax) ## Stockage des H
  
  ## fonction H définition 
  
  H[1]=0 # initialisation de la fonction H
  for (i in 1:(13)){
    H[1]=H[1]+W[s[i]]*Y[i]
  }
  meilleur_H=H[1]
  meilleur_H
  s_best=s
  for (j in 1:nmax)
  {
    if ( REMO == 2) {bet0=bet*j^(1/2)} # Recuit simulé 
    if ( REMO == 1) {bet0=bet} # Metropolis
    
    f=fonction_transition(s,bet0)
    s=f[1:N]
    DELTA=f[N+1]
    
    H[j+1]=H[j]-DELTA
    
    ## ameloration denotre fonction H afin d'obtenir le meilleur H (meilleur poids)
    if (H[j+1]<meilleur_H)
    {
      meilleur_H=H[j+1]
      meilleur_s=s
    }
  } ## Fin de boucle d'iteration
  
  ## Affichage de l'evolution de H au cours des 500 iterations
  plot(H, type="l",lwd=4 ,main=paste("valeur minimale de la fonction H est : ",
                                     meilleur_H),col="blue",xlab="Nombre d'iterations de l'algorithme",
       ylab="La Valeurs de notre fonction H")
  
  return(c(meilleur_s,meilleur_H))
}
R1 = cover_algorithme(500, 15, 2)## recuit simulé
#R2 = cover_algorithme(500, 15, 1) #metropolis

############################################################
