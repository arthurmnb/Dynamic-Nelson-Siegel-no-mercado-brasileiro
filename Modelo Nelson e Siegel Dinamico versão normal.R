############################################# Definição dos parâmetros do modelo ################################################
set.seed(232)

setwd("INSERIR DIRETORIO DE TRABALHO") #Inserir diretorio de trabalho

#Matriz com as taxas a serem utilizadas para modelagem já em %
y <- as.matrix(read.table("y.txt")) #Abrir arquivo dos estados observaveis utilizados

#Vetor de maturidades a serem analisadas, no caso: 3   6   9  12  15  18  21  24  30  36  48  60  72  84  96 108 120 meses
maturidade <- append(append(seq(3,24,3),c(30,36)),seq(48,120,12))  #maturidades dos títulos em meses


#########################################################################################################################################################
######################################### Definição dos valores utilizados a priori para o modelo #######################################################
#########################################################################################################################################################

#Priori para a media dos estados nao observaveis mi
mi.priori <- c(7.43, -1.37, 0.54)
S.priori.mi <- 10*diag(3)


#Priori para a matriz de covariancia dos estados observaveis sigma
alpha.priori.sigma <- 5
beta.priori.sigma <- 0.2


#Priori para o parametro lambda da equação de NS
alpha.priori.lambda <- 0.01
beta.priori.lambda <- 0.1


#Priori para a matriz de covariancia dos estados nao observaveis W
S.priori.W <- 5*diag(3)
v.priori <- 5


#Priori para a matriz de troca de estados do processo nao observavel A
mi.A.priori <- 0.9*diag(3)
C.A.priori <- 100*diag(3)


###############################################################################################################################
############################################### Funções auxiliares ############################################################
###############################################################################################################################

#função para gerar a matriz lambda 
#Entradas: parâmetro lambda (escalar) e as maturidades (vetor)

matriz_lambda <- function(lambda,maturidade){
  matriz.lambda <- matrix(0,length(maturidade),3) #Cria a matriz para ser preenchida
  matriz.lambda[,1] <- rep(1,length(maturidade)) #Preenche a primeira coluna da matriz
  matriz.lambda[,2] <- (1-exp(-lambda*maturidade))/(lambda*maturidade) #Preenche a segunda coluna da matriz
  matriz.lambda[,3] <- (1-exp(-lambda*maturidade))/(lambda*maturidade)-exp(-lambda*maturidade) #Preenche a terceira coluna da matriz
  return(matriz.lambda)
}

#Função para calcular as estatisticas utilizadas na geração da proposta de lambda no Metropolis Hastings
#Entradas: lambda anterior (escalar), alpha a priori de lambda (escalar), beta a priori de lambda (escalar), maturidade utilizada (vetor), estados não observáveis f (matriz (T+1)x3), estados observaveis y (matrix Txr), matriz de covariancia sigma (matriz rxr), variancia para a variancia sigma U (vetor)

dados.proposta.lambda <- function(lambda, alpha.priori.lambda, beta.priori.lambda,maturidade,f,y,sigma){
  #primeira derivada da matriz lambda
  matriz.lambda.d1 <- matrix(0,length(maturidade),3) #Cria a matriz para ser preenchida
  matriz.lambda.d1[,2] <- exp(-lambda*maturidade)/lambda + (exp(-lambda*maturidade)-1)/(lambda**2*maturidade) #Preenche a segunda coluna da matriz
  matriz.lambda.d1[,3] <- matriz.lambda.d1[,2] + maturidade*exp(-lambda*maturidade) #Preenche a terceira coluna da matriz
  #segunda derivada da matriz lambda
  matriz.lambda.d2 <- matrix(0,length(maturidade),3) #Cria a matriz para ser preenchida
  matriz.lambda.d2[,2] <- 2*(1 - exp(-lambda*maturidade))/(maturidade*lambda**3) - exp(-lambda*maturidade)*(lambda*maturidade + 2)/(lambda**2) #Preenche a segunda coluna da matriz
  matriz.lambda.d2[,3] <- matriz.lambda.d2[,2] - maturidade**2*exp(-lambda*maturidade) #Preenche a terceira coluna da matriz
  
  f <- f[-1,] #Removendo f em t=0
  
  q.d1 <- (alpha.priori.lambda-1)/lambda - beta.priori.lambda + sum((matrix(rep(diag(sigma),nrow(y)),byrow = T,ncol = length(maturidade))**-1)*f%*%t(matriz.lambda.d1)*(y - f%*%t(matriz_lambda(lambda,maturidade))))
  
  q.d2 <- -(alpha.priori.lambda -1)/lambda**2 + sum((matrix(rep(diag(sigma),nrow(y)),byrow = T,ncol = length(maturidade))**-1)*((y - f%*%t(matriz_lambda(lambda,maturidade)))*(f%*%t(matriz.lambda.d2)) - (f%*%t(matriz.lambda.d1))*(f%*%t(matriz.lambda.d1))))
  
  med.sd <- matrix(c(mean = lambda - (q.d1/q.d2), sd = sqrt(-q.d2**-1)),ncol=2)
  return(med.sd)
}

#Função para calcular o log da função de densidade alvo de lambda
#Entradas: lambda (escalar), alpha da priori de lambda (escalar), beta da pripro de lambda (escalar), maturidades utilizadas (vetor), matriz de estados nao observaveis f (matriz (T+1)x3), matriz de estados observaveis y (matriz Txr), matriz de covariancia sigma (matriz rxr), variancia para a variancia sigma U (vetor)

log.alvo.lambda <- function(lambda, alpha.priori.lambda, beta.priori.lambda, maturidade, f, y, sigma){
  f <- f[-1,]
  log.densidade <- log(lambda)*(alpha.priori.lambda-1) - lambda*beta.priori.lambda - 0.5*sum(colSums((matrix(rep(diag(sigma),nrow(y)),byrow = T,ncol = length(maturidade))**-1)*(y - f%*%t(matriz_lambda(lambda,maturidade)))*(y - f%*%t(matriz_lambda(lambda,maturidade)))))
  return(log.densidade)
}

#Função para cálculo da matriz de covariancia de f em t=0
#Entradas: matriz de troca de estados A (matriz 3x3) e matriz de covariancia W (matriz 3x3)

C0.t <- function(A,W){
  C0 <- matrix(solve(diag(9) - kronecker(A,A))%*%matrix(c(W),ncol=1),byrow = F,ncol=3) #Calcula C0
  if (isSymmetric.matrix(C0)==FALSE){ #Força a simetria caso tenha algum erro de arredondamento
    C0[2,1] <- C0[1,2]
    C0[3,1] <- C0[1,3]
    C0[3,2] <- C0[2,3]
  }
  return(C0)
}

#Função para log da função auxiliar g(A,W) da densidade de f0
#Entradas: matriz de troca de estados A (matriz 3x3), matriz de covariancia W (matriz 3x3), matriz de estados nao observaveis f (matriz (T+1)x3) e matriz da média dos estados nao observaveis mi (matriz 3x1)

l.g.AW <- function(A,W,f,mi){
  C0 <- C0.t(A,W) #Calculo da covariancia de f em t=0
  g.AW <- -.5*log(det(C0)) -.5*((f[1,] - matrix(mi,ncol=3))%*%solve(C0)%*%t(f[1,] - matrix(mi,ncol=3))) #Log da densidade de f em t=0
  return(g.AW)
}

#Função para gerar propostas para a matriz de troca de estados A para o Metropolis Hastings
#Entradas: media da priori de A (matriz 3x3), matriz de covariancia a priori de A (matriz 3x3), matriz de estados nao observaveis f (matriz (T+1)x3), media dos estados nao observaveis mi (matriz 3x1), matriz de covariancia W (matriz 3x3)

proposta.alvo.A <- function(mi.A.priori, C.A.priori, f, mi, W){
  
  s <- f[-nrow(f),] - matrix(rep(mi,nrow(f)-1),byrow = T,ncol = 3)
  s <- matrix(colSums(cbind(s[,1]*s[,1],s[,1]*s[,2],s[,1]*s[,3],s[,2]*s[,1],s[,2]*s[,2],s[,2]*s[,3],s[,3]*s[,1],s[,3]*s[,2],s[,3]*s[,3])),ncol=3)
  C.til <- solve(solve(C.A.priori) + s)
  
  s <- f[-nrow(f),] - matrix(rep(mi,nrow(f)-1),byrow = T,ncol = 3)
  r <- f[-1,] - matrix(rep(mi,nrow(f)-1),byrow = T,ncol = 3)
  s <- matrix(colSums(cbind(s[,1]*r[,1],s[,1]*r[,2],s[,1]*r[,3],s[,2]*r[,1],s[,2]*r[,2],s[,2]*r[,3],s[,3]*r[,1],s[,3]*r[,2],s[,3]*r[,3])),ncol=3)
  A.til <- (mi.A.priori%*%solve(C.A.priori) + s)%*%C.til
  
  autovalormax <- 1
  while (autovalormax >=1){
    proposta <- matrix(mvtnorm::rmvnorm(1,c(t(A.til)),kronecker(W,C.til)),byrow = T, ncol = 3)
    autovalormax <- max(Mod(eigen(proposta)$values))
  }
  
  return(proposta)
}


#########################################################################################################################################################
###################################### Algoritmos MCMC para geração da amostra ##########################################################################
#########################################################################################################################################################

#Função do algoritmo MCMC Metropolis Hastings para lambda
#Entradas: valor anterior da cadeia lambda.anterior (escalar), alpha da priori de lambda (escalar), beta da priori de lambda (escalar), maturidades utilizadas (vetor), matriz de estados nao observaveis f (matriz (T+1)x3), matriz de estados observaveis y (matriz Txr), matriz de covariancia sigma (matriz rxr), variancia da variancia U (vetor)

metropolis.hastings.lambda <- function(lambda.anterior, alpha.priori.lambda, beta.priori.lambda, maturidade, f, y, sigma){
  
  est.anterior <- dados.proposta.lambda(lambda.anterior, alpha.priori.lambda, beta.priori.lambda,maturidade,f,y,sigma)
  
  proposta <- est.anterior[1,1] + est.anterior[1,2]*rnorm(1,0.0,1.0)
  
  est.proposta <- dados.proposta.lambda(proposta, alpha.priori.lambda, beta.priori.lambda,maturidade,f,y,sigma)
  
  alpha <- min(1,exp(log.alvo.lambda(proposta, alpha.priori.lambda, beta.priori.lambda, maturidade, f, y, sigma) - log.alvo.lambda(lambda.anterior, alpha.priori.lambda, beta.priori.lambda, maturidade, f, y, sigma) + dnorm(lambda.anterior, mean = est.proposta[1,1], sd = est.proposta[1,2], log = T) - dnorm(proposta, mean = est.anterior[1,1], sd = est.anterior[1,2], log = T)))
  
  if (runif(1)<alpha){
    cadeia <- proposta
  } else {
    cadeia <- lambda.anterior
  }
  return(cadeia)
}

#Função do algoritmo MCMC Metropolis Hastings para matriz de covariancia W
#Entradas: valor anterior de W (matriz 3x3), v priori de W (escalar), S priori de W (matriz 3x3), média da priori de A mi.a.priori (matriz 3x3), variância a priori de A (matriz 3x3), matriz de troca de estados A (matriz 3x3), matriz de estados nao observaveis f ((T+1)x3), média dos estados nao observaveis mi (matriz 3x1)

metropolis.hastings.W <- function(W.anterior, v.priori, S.priori.W, mi.A.priori, C.A.priori, A, f, mi){
  
  S1 <- (A - mi.A.priori)%*%solve(C.A.priori)%*%t(A - mi.A.priori)
  
  S2 <- f[-1,]  - matrix(rep(mi,nrow(f)-1),byrow = T,ncol=ncol(f)) - as.matrix(f[-nrow(f),] - matrix(rep(mi,nrow(f)-1),byrow = T,ncol=ncol(f)))%*%t(A)
  
  S2 <- colSums(cbind(S2[,1]*S2[,1],S2[,1]*S2[,2],S2[,1]*S2[,3],S2[,2]*S2[,1],S2[,2]*S2[,2],S2[,2]*S2[,3],S2[,3]*S2[,1],S2[,3]*S2[,2],S2[,3]*S2[,3]))
  
  proposta <- MCMCpack::riwish((nrow(f)-1) + v.priori + 3, S.priori.W + S1 + matrix(S2,ncol=3))
  
  alpha <- min(1,exp(l.g.AW(A,proposta,f,mi) - l.g.AW(A,W.anterior,f,mi)))
  
  if(runif(1)<alpha){
    cadeia <- proposta
  } else {
    cadeia <- W.anterior
  }
  return(cadeia)
}

#Função do algoritmo MCMC Metropolis Hastings para matriz de de troca de estados A
#Entradas: valor anterior de A (matriz 3x3), média da priori de A mi.a.priori (matriz 3x3), variância a priori de A C.A.priori (matriz 3x3), matriz de covariancia W (matriz 3x3), matriz de estados nao observaveis f (matriz (T+1)x3), média dos estados nao observaveis mi (matriz 3x1)

metropolis.hastings.A <- function(A.anterior,mi.A.priori, C.A.priori, W, f, mi){
  
  proposta <- proposta.alvo.A(mi.A.priori, C.A.priori, f, mi, W)
  
  alpha <- min(1,exp(l.g.AW(proposta,W,f,mi) - l.g.AW(A.anterior,W,f,mi)))
  
  if(runif(1)<alpha){
    cadeia <- proposta
  } else {
    cadeia <- A.anterior
  }
  return(cadeia)
}

#Função para geração de valores para matriz de covariancia sigma
#Entradas: alpha priori de sigma alpha.priori.sigma (escalar), beta da priori de sigma beta.priori.sigma (escalar), valor do parametro lambda (escalar), maturidades utilizadas (vetor), matriz de estados nao observaveis f (matriz (T+1)x3), matriz de estados observaveis y (matriz Txr), variancia da variancia U (vetor)

sigma.t <- function(alpha.priori.sigma,beta.priori.sigma,lambda,maturidade,f,y){
  f <- f[-1,] #Removendo f em t=0
  matriz.lambda <- matriz_lambda(lambda,maturidade)
  alpha <- nrow(y) + alpha.priori.sigma/2
  soma <- colSums(((y - f%*%t(matriz.lambda))**2))
  beta <- beta.priori.sigma/2 + soma
  matriz.sigma <- diag(mapply( "rgamma", n=1, shape=alpha/2, rate=beta/2 )**-1)
  return(matriz.sigma)
}

#Função para geração de valores da media dos estados nao observaveis mi
#Entradas: media da priori de mi mi.priori (matriz 3x1), covariancia da priori de mi S.priori (matriz 3x3), matriz de troca de estados A (matriz 3x3), matriz de estados nao observaveis f (matriz (T+1)x3), matriz de covariancia W (matriz 3x3)

mi.t <- function(mi.priori,S.priori,A,f,W){
  dif.f <- f[2:nrow(f),] -  f[1:nrow(f)-1,]%*%t(A)
  soma_para_mi0 <- dif.f%*%t(t(diag(3) - A)%*%solve(W))
  soma_para_mi0 <- matrix(c(sum(soma_para_mi0[,1]),sum(soma_para_mi0[,2]),sum(soma_para_mi0[,3])),nrow=3)
  
  variancia.mi = solve(solve(S.priori) +(nrow(y))*t((diag(rep(1,3)) - A))%*%solve(W)%*%(diag(rep(1,3)) - A) + solve(C0.t(A,W))) #variância da condicional completa de mi
  media.mi <- variancia.mi%*%(solve(S.priori)%*%mi.priori + soma_para_mi0 + solve(C0.t(A,W))%*%matrix(f[1,],ncol=1)) #Média da condicional completa de mi
  return(matrix(mvtnorm::rmvnorm(1,mean = media.mi, sigma = variancia.mi),ncol=1)) #Retorna uma amostra de um valor para mi (matrix 3x1)
}

#função para calcular o valor do proximo estado nao observavel f na simulação
#Entradas: matriz de estados nao observaveis f (matrix (t-1)x3), media dos estados nao observaveis mi (matriz 3x1), matriz de troca de estados A (matrix 3x3) e matriz de covariância W (matrix 3x3)

f.t <- function(mi,A,f,W){
  return(mvtnorm::rmvnorm(1, mean = mi + A%*%(matrix(f[nrow(f),],ncol=1)-mi), sigma = W)) #Retorna uma amostra de um valor para f no tempo t (matrix 3x1)
}

#função para calcular y(t) na simulação
#Entradas: parâmetro lambda (escalar), as maturidades (vetor), estado nao observavel f no tempo t (matrix 3x1) e matriz de covariância sigma (matrix 3x3)

y.t <- function(lambda,maturidade,f,sigma){
  return(c(matriz_lambda(lambda,maturidade)%*%matrix(f,nrow=3)) + mvtnorm::rmvnorm(1, mean = rep(0,17),sigma = sigma))
}

####################################################################################################################################################################
################################################## FILTRO DE KALMAN E FFBS #########################################################################################
####################################################################################################################################################################

#Função para calcular a media e variancia das distribuições do filtro de kalman
#Entradas: media a priori dos estados nao observaveis mi0 (matriz 3x1), covariancia a priori dos estados nao observaveis C0 (matriz 3x3), matriz de covariancia dos estados observaveis sigma (matriz rxr), matriz de covariancia dos estados nao observaveis W (matriz 3x3), matriz dos estados nao observaveis Ft (matriz (T+1)x3), matriz dos estados observaveis y (matriz Txr), variancia adicionada a variancia de y U (vetor)

filtro.kalman <- function(mi0, C0, sigma, W, Ft, Gt, y){
  
  at <- as.list(rep(NA,nrow(y)))
  at[[1]] <- Gt%*%mi0
  Rt <- as.list(rep(NA,nrow(y)))
  Rt[[1]] <- Gt%*%C0%*%t(Gt) + W
  
  ft <- as.list(rep(NA,nrow(y)))
  ft[[1]] <- Ft%*%at[[1]]
  Qt <- as.list(rep(NA,nrow(y)))
  Qt[[1]] <- Ft%*%Rt[[1]]%*%t(Ft) + sigma
  
  erro <- as.list(rep(NA,nrow(y)))
  erro[[1]] <- y[1,] - c(ft[[1]])
  mt <- as.list(rep(NA,nrow(y)))
  mt[[1]] <- at[[1]] + Rt[[1]]%*%t(Ft)%*%solve(Qt[[1]])%*%erro[[1]]
  Ct <- as.list(rep(NA,nrow(y)))
  Ct[[1]] <- Rt[[1]] - Rt[[1]]%*%t(Ft)%*%solve(Qt[[1]])%*%Ft%*%Rt[[1]]
  
  #Forçando simetria na matriz por conta da precisão do arredondamento do R 
  Ct[[1]][2,1] <- Ct[[1]][1,2]
  Ct[[1]][3,1] <- Ct[[1]][1,3]
  Ct[[1]][3,2] <- Ct[[1]][2,3]
  
  for (ii in 2:nrow(y)){
    
    at[[ii]] <- Gt%*%mt[[ii-1]]
    Rt[[ii]] <- Gt%*%Ct[[ii-1]]%*%t(Gt) + W
    
    ft[[ii]] <- Ft%*%at[[ii]]
    Qt[[ii]] <- Ft%*%Rt[[ii]]%*%t(Ft) + sigma
    
    erro[[ii]] <- y[ii,] - c(ft[[ii]])
    mt[[ii]] <- at[[ii]] + Rt[[ii]]%*%t(Ft)%*%solve(Qt[[ii]])%*%erro[[ii]]
    Ct[[ii]] <- Rt[[ii]] - Rt[[ii]]%*%t(Ft)%*%solve(Qt[[ii]])%*%Ft%*%Rt[[ii]]
    
    #Forçando simetria na matriz por conta da precisão do arredondamento do R
    Ct[[ii]][2,1] <- Ct[[ii]][1,2]
    Ct[[ii]][3,1] <- Ct[[ii]][1,3]
    Ct[[ii]][3,2] <- Ct[[ii]][2,3]
    
  }
  return(list(at,Rt,ft,Qt,mt,Ct))
}

#Função para gerar valores para os estados nao observaveis f atraves do metodo de FFBS
#Entradas: media a priori dos estados nao observaveis mi0 (matriz 3x1), covariancia a priori dos estados nao observaveis C0 (matriz 3x3), matriz de covariancia dos estados observaveis sigma (matriz rxr), matriz de covariancia dos estados nao observaveis W (matriz 3x3), matriz dos estados nao observaveis Ft (matriz (T+1)x3), matriz dos estados observaveis y (matriz Txr), variancia adicionada a variancia de y U (vetor) 

FFBS <- function(mi0, C0, sigma, W, Ft, Gt, y){
  
  filtro <- filtro.kalman(mi0, C0, sigma, W, Ft, Gt, y)
  
  at <- filtro[[1]]
  Rt <- filtro[[2]]
  mt <- filtro[[5]]
  Ct <- filtro[[6]]
  
  t.max <- length(at)
  
  theta <- as.list(rep(NA,t.max+1))
  theta[[t.max+1]] <- mvtnorm::rmvnorm(1,mean = mt[[t.max]], sigma = Ct[[t.max]])
  
  ht <- as.list(rep(NA,t.max))
  Ht <- as.list(rep(NA,t.max))
  
  for (tt in t.max:2){
    
    ht[[tt]] <- t(mt[[tt-1]] + Ct[[tt-1]]%*%t(Gt)%*%solve(Rt[[tt]])%*%t(theta[[tt+1]] - t(at[[tt]])))
    
    Ht[[tt]] <- Ct[[tt-1]] - Ct[[tt-1]]%*%t(Gt)%*%solve(Rt[[tt]])%*%Gt%*%Ct[[tt-1]]
    
    #Forçando simetria na matriz por conta da precisão do arredondamento do R
    Ht[[tt]][2,1] <- Ht[[tt]][1,2]
    Ht[[tt]][3,1] <- Ht[[tt]][1,3]
    Ht[[tt]][3,2] <- Ht[[tt]][2,3]
    
    
    theta[[tt]] <- mvtnorm::rmvnorm(1,mean = ht[[tt]], sigma = Ht[[tt]])
  }
  
  ht[[1]] <- t(mi0 + C0%*%t(Gt)%*%solve(Rt[[1]])%*%t(theta[[2]] - t(at[[1]])))
  
  Ht[[1]] <- C0 - C0%*%t(Gt)%*%solve(Rt[[1]])%*%Gt%*%C0
  
  #Forçando simetria na matriz por conta da precisão do arredondamento do R
  Ht[[1]][2,1] <- Ht[[1]][1,2]
  Ht[[1]][3,1] <- Ht[[1]][1,3]
  Ht[[1]][3,2] <- Ht[[1]][2,3]
  
  
  theta[[1]] <- mvtnorm::rmvnorm(1,mean = ht[[1]], sigma = Ht[[1]])
  
  theta <- matrix(unlist(theta),byrow = T,ncol = 3)
  return(theta)
}


########################################################################################################################################################
####################################################### Aplicação do modelo ############################################################################
########################################################################################################################################################

#Valores para o inicio do algoritmo

lambda.gibbs <- 0.04 #Valor inicial do parametro lambda
mi.gibbs <- matrix(c(7,-3,4),ncol=3) #valor inicial da media dos estados nao observaveis mi
sigma.gibbs <- rep(1,length(maturidade)) #valor inicial da matriz de covariancia dos estados observaveis sigma
W.gibbs <- diag(3) #valor inicial da matriz de covariancia dos estados nao observaveis W
A.gibbs <- matrix(c(0.5, -0.1, 0.04, 0.005, 0.5, 0.0002, -0.03, 0.09, 0.5), byrow = F, ncol = 3) #valor inicial da matriz de troca de estados A

#For para rodar as iterações do algoritmo
for (ii in 1:(1*10**5)){
  
  ffbs <- FFBS(t(mi.gibbs), C0.t(A.gibbs,W.gibbs), diag(c(sigma.gibbs)), W.gibbs, matriz_lambda(lambda.gibbs, maturidade), A.gibbs, y)
  
  lambda.gibbs <- metropolis.hastings.lambda(lambda.gibbs,alpha.priori.lambda, beta.priori.lambda, maturidade, ffbs, y, diag(c(sigma.gibbs)))
  mi.gibbs <- t(mi.t(mi.priori,S.priori.mi,A.gibbs,ffbs,W.gibbs))
  sigma.gibbs <- matrix(diag(sigma.t(alpha.priori.sigma,beta.priori.sigma,lambda.gibbs,maturidade,ffbs,y)), nrow = 1)
  W.gibbs <- metropolis.hastings.W(W.gibbs,v.priori, S.priori.W, mi.A.priori, C.A.priori, A.gibbs, ffbs, mi.gibbs)
  A.gibbs <- metropolis.hastings.A(A.gibbs,mi.A.priori, C.A.priori, W.gibbs, ffbs, mi.gibbs)
 
  if (ii>20000 & ii%%16 == 0){
    write.table(matrix(c(ffbs[,1]),nrow = 1) ,file = "nivel.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(matrix(c(ffbs[,2]),nrow = 1) ,file = "inclinacao.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(matrix(c(ffbs[,3]),nrow = 1) ,file = "curvatura.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(lambda.gibbs ,file = "lambda.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(mi.gibbs ,file = "mi.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(sigma.gibbs ,file = "sigma.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(matrix(c(W.gibbs),nrow = 1) ,file = "W.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    write.table(matrix(c(A.gibbs),nrow = 1) ,file = "A.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
    print(ii)
  }
}
