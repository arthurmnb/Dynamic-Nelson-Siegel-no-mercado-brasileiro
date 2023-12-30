import numpy as np
from itertools import compress
from math import gamma
from scipy.stats import invwishart

diretorio = "DIRETORIO DE TRABALHO"

#LEITURA DAS TAXAS - ALTERAR NOME DE Y.TXT PARA O NOME DO ARQUIVO DE TAXAS
y = 100*np.loadtxt(diretorio + '\\y.txt', delimiter=';',usecols=np.arange(1,18,1), skiprows=1)
ultimo_na = np.max(list(compress(range(len(y[:,(y.shape[1]-1)])), np.isnan(y[:,(y.shape[1]-1)]))))
y = y[(ultimo_na+1):,:]

maturidade = np.array([3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120])

#####################################################################################################################################################################################
######################################################################### DADOS DA PRIORI ###########################################################################################
#####################################################################################################################################################################################

#Priori para a media dos estados nao observaveis mu
mu_priori = np.array([7.43, -1.37, 0.54]).reshape(3,1)
S_priori = 10*np.identity(3)

#Priori para a matriz de covariancia dos estados observaveis sigma
alpha_priori_sigma = 5
beta_priori_sigma = 0.2

#Priori para o parametro lambda da equação de NS
alpha_priori_lambda = 0.01
beta_priori_lambda = 0.1

#Priori para a matriz de covariancia dos estados nao observaveis W
S_priori_W = 5*np.identity(3)
v_priori = 5

#Priori para a matriz de troca de estados do processo nao observavel A
mu_A_priori = 0.9*np.identity(3)
C_A_priori = 100*np.identity(3)

######################################################################################################################################################################################
################################################################ CRIAÇÃO DOS ARQUIVOS PARA SALVAR ####################################################################################
######################################################################################################################################################################################

open(diretorio+'\\lambda.txt', 'a').close()
open(diretorio+'\\mu.txt', 'a').close()
open(diretorio+'\\sigma.txt', 'a').close()
open(diretorio+'\\W.txt', 'a').close()
open(diretorio+'\\A.txt', 'a').close()
open(diretorio+'\\Nivel.txt', 'a').close()
open(diretorio+'\\Inclinacao.txt', 'a').close()
open(diretorio+'\\Curvatura.txt', 'a').close()

#####################################################################################################################################################################################
###################################################################### FUNÇÕES DO MODELO ############################################################################################
#####################################################################################################################################################################################

def dnorm_log(x:float = 1,mu:float = 0,sd:float = 1) -> float:
    '''Função para calcular a densidade do log de uma distribuição normal

    x: float, default 1
        Valor a ser avaliado na densidade de probabilidade
    mu: float, default 0
        Valor da média da distribuição
    sd: float, default 1
        Valor do desvio-padrão da distribuição
    return = ln: float
        Valor da densidade de probabilidade do log da normal com média mu e desvio-padrão sd avaliado em x
    '''
    ln = -0.5*((x-mu)/sd)**2 - np.log(sd) - 0.5*np.log(2*np.pi)
    return ln

def matriz_lambda(lambdA: float=1, maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120])-> list:
    '''Função para calcular a matriz de carrregamentos de decaimento do modelo de Nelson-Siegel

    LambdA: float, default 1
        Fator de decaimento
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = matriz: list
        A matriz com os carregamentos do modelo
    '''
    maturidade = np.array(maturidade)
    matriz = np.array([
        [1]*len(maturidade),
        list((1-np.exp(-lambdA*maturidade))/(lambdA*maturidade)),
        list((1-np.exp(-lambdA*maturidade))/(lambdA*maturidade)-np.exp(-lambdA*maturidade))
    ])
    return matriz.T

def dados_proposta_lambda(lambdA:float,y:list,f:list,sigma:list,alpha_priori_lambda: float = 0.01, beta_priori_lambda:float = 0.1,maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]) -> tuple:
    '''Função para calcular a média e variância para geração da proposta para lambda

    LambdA: float
        Fator de decaimento no estado anterior
    y: list
        Array com os valores das taxas dimensão: M x tau
    f: list
        Array com os fatores latentes dimensão: (M+1) x 3
    sigma: list
        Lista com a diagonal da matriz de covariância das taxas dimensão: tau elementos
    alpha_priori_lambda: float, default 0.01
        Valor do parâmetro alpha da priori de lambda
    beta_priori_lambda:float, default 0.1
        Valor do parâmetro beta da priori de lambda
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = media_sd: tuple
        tupla com os valores de media e variancia para a proposta
    '''
    maturidade = np.array(maturidade)
    
    tau = y.shape[1]
    M = y.shape[0]

    matriz_lambda_d1 = np.array([
        [0]*tau,
        list(np.exp(-lambdA*maturidade)/lambdA + (np.exp(-lambdA*maturidade)-1)/(lambdA**2*maturidade)),
        list(np.exp(-lambdA*maturidade)/lambdA + (np.exp(-lambdA*maturidade)-1)/(lambdA**2*maturidade) + maturidade*np.exp(-lambdA*maturidade))
    ])

    matriz_lambda_d2 = np.array([
        [0]*tau,
        list(2*(1 - np.exp(-lambdA*maturidade))/(maturidade*lambdA**3) - np.exp(-lambdA*maturidade)*(lambdA*maturidade + 2)/(lambdA**2)),
        2*(1 - np.exp(-lambdA*maturidade))/(maturidade*lambdA**3) - np.exp(-lambdA*maturidade)*(lambdA*maturidade + 2)/(lambdA**2) - maturidade**2*np.exp(-lambdA*maturidade)
    ])

    f = np.delete(f,(0),axis=0)

    q_d1 = (alpha_priori_lambda-1)/lambdA - beta_priori_lambda + np.sum(np.multiply(np.array(sigma*M, dtype=float).reshape(M,tau)**-1),np.multiply(np.matmul(f,matriz_lambda_d1),(y-np.matmul(f,matriz_lambda(lambdA,maturidade).T))))

    q_d2 = -(alpha_priori_lambda -1)/lambdA**2 + np.sum(np.multiply(np.array(sigma*M, dtype=float).reshape(M,tau)**-1,np.multiply((y - np.matmul(f,matriz_lambda(lambdA,maturidade).T)),np.matmul(f,matriz_lambda_d2)) - np.multiply(np.matmul(f,matriz_lambda_d1),np.matmul(f,matriz_lambda_d1))))

    media_sd = (lambdA - (q_d1/q_d2),  np.sqrt(-q_d2**-1))

    return media_sd

def log_alvo_lambda(lambdA:float,y:list,f:list,sigma:list,alpha_priori_lambda: float = 0.01, beta_priori_lambda:float = 0.1,maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]) -> float :
    '''Função para calcular o log da função de densidade alvo de lambda

    LambdA: float
        Fator de decaimento no estado anterior
    y: list
        Array com os valores das taxas dimensão: M x tau
    f: list
        Array com os fatores latentes dimensão: (M+1)x3
    sigma: list
        Lista com a diagonal da matriz de covariância das taxas dimensão: tau elementos
    alpha_priori_lambda: float, default 0.01
        Valor do parâmetro alpha da priori de lambda
    beta_priori_lambda:float, default 0.1
        Valor do parâmetro beta da priori de lambda
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = log_densidade: float
        Valor do log da densidade da distribuição alvo de lambda
    '''
    f = np.delete(f,(0),axis=0)

    M = y.shape[0]
    tau = y.shape[1]

    log_densidade = np.multiply(np.log(lambdA),(alpha_priori_lambda-1)) - lambdA*beta_priori_lambda - 0.5*np.sum(np.multiply(np.array(sigma*M, dtype=float).reshape(M,tau)**-1),(y - np.matmul(f,matriz_lambda(lambdA,maturidade).T))**2)
    
    return log_densidade

def C0_t(A:list,W:list) -> list:
    '''Função para cálculo da matriz de covariancia do fator latente f em t=0

    A: list
        Array referente a matriz de persistência de estados dimensão 3x3
    W: list
        Array referente a matriz de covariância dos fatores latentes dimensão 3x3
    return = C0: list
        Array da matriz de covariância do fator latente no tempo 0
    '''
    C0 = np.array(np.matmul(np.linalg.inv(np.identity(9) - np.kron(A,A)),W.reshape(9,1))).reshape(3,3)
    return C0

def l_g_AW(A:list,W:list,mu:list,f:list) -> float:
    '''Função para log da função auxiliar g(A,W) da densidade de f0

    A: list
        Array referente a matriz de persistência de estados dimensão: 3x3
    W: list
        Array referente a matriz de covariância dos fatores latentes dimensão: 3x3
    mu: list
      Array da média dos fatores latentes dimensão: 3x1
    return = g_AW: float
        Valor do log da função auxiliar para a densidade de f0
    '''
    C0 = C0_t(A,W)
    
    g_AW = -0.5*np.log(np.linalg.det(C0)) - 0.5*np.matmul(f[0] - mu.reshape(1,3),np.matmul(np.linalg.inv(C0),(f[0] - mu.reshape(1,3)).T))
    
    return g_AW

def proposta_alvo_A(W:list,mu:list,f:list,mu_A_priori:list,C_A_priori:list) -> list:
    '''Função para gerar a proposta para a matriz de persistência A

    W: list
      Array referente a matriz de covariância dos fatores latentes dimensão: 3x3
    mu: list
      Array da média dos fatores latentes dimensão: 3x1
    f: list
      Array com os estados não observáveis dimensão: (M+1) x 3
    mu_A_priori: list
      Array com a média a priori de A dimensão: 3x3
    C_A_priori: list
      Array com a primeira matriz de covariância da priori de A dimensão: 3x3
    return = g_AW: float
      Valor do log da função auxiliar para a densidade de f0
    '''
    s = np.delete(f,(f.shape[0]-1),axis=0) - np.repeat(mu.reshape(1,3),f.shape[0]-1,axis=0)

    s1 = s[:,0].reshape(s.shape[0],1)
    s2 = s[:,1].reshape(s.shape[0],1)
    s3 = s[:,2].reshape(s.shape[0],1)

    C_til = np.linalg.inv(np.linalg.inv(C_A_priori) + np.einsum('ij->i',np.array([s1**2,s1*s2,s1*s3,s2*s1,s2**2,s2*s3,s3*s1,s3*s2,s3**2]).reshape(9,s.shape[0])).reshape(3,3))

    r = np.delete(f,(0),axis=0) - np.repeat(mu.reshape(1,3),f.shape[0]-1,axis=0)

    r1 = r[:,0].reshape(r.shape[0],1)
    r2 = r[:,1].reshape(r.shape[0],1)
    r3 = r[:,2].reshape(r.shape[0],1)

    A_til = np.matmul(np.matmul(mu_A_priori,np.linalg.inv(C_A_priori)) + np.einsum('ij->i',np.array([s1*r1,s1*r2,s1*r3,s2*r1,s2*r2,s2*r3,s3*r1,s3*r2,s3*r3]).reshape(9,s.shape[0])).reshape(3,3), C_til)

    autovalormax = 1

    while autovalormax>=1:
        proposta = np.random.multivariate_normal(np.concatenate(A_til.reshape(1,9)),np.kron(W,C_til),1).reshape(3,3)
        autovalormax = max(abs(np.linalg.eig(proposta)[0]))
    
    return proposta

def MH_lambda(lambda_anterior:float,y:list,f:list,sigma:list,alpha_priori_lambda:float = 0.01, beta_priori_lambda:float = 0.1, maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]) -> float:
  '''Função do algoritmo MCMC Metropolis Hastings para geração de lambda

    Lambda_anterior: float
        Fator de decaimento no estado anterior
    y: list
        Array com os valores das taxas dimensão: M x tau
    f: list
        Array com os fatores latentes dimensão: (M+1)x3
    sigma: list
        Lista com a diagonal da matriz de covariância das taxas dimensão: tau elementos
    alpha_priori_lambda: float, default 0.01
        Valor do parâmetro alpha da priori de lambda
    beta_priori_lambda:float, default 0.1
        Valor do parâmetro beta da priori de lambda
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = log_densidade: float
        Valor do log da densidade da distribuição alvo de lambda
    '''
  est_anterior = dados_proposta_lambda(lambda_anterior,y,f,sigma,alpha_priori_lambda, beta_priori_lambda, maturidade)
  proposta = np.random.normal(est_anterior[0],est_anterior[1],1)
  est_proposta = dados_proposta_lambda(proposta,y,f,sigma,alpha_priori_lambda, beta_priori_lambda, maturidade)

  alpha = min(1,np.exp( log_alvo_lambda(proposta,y,f,sigma,alpha_priori_lambda,beta_priori_lambda,maturidade) - log_alvo_lambda(lambda_anterior,y,f,sigma,alpha_priori_lambda,beta_priori_lambda,maturidade) + dnorm_log(x=lambda_anterior,mu=est_proposta[0], sd=est_proposta[1]) -  dnorm_log(x=proposta,mu=est_anterior[0], sd=est_anterior[1])))

  if ( np.random.uniform(0,1,1) < alpha):
    prox_lambda = proposta
  else:
    prox_lambda = lambda_anterior
  
  return prox_lambda

def MH_W(W_anterior:list,f:list,A:list,mu:list,v_priori:float,S_priori_W:list, mu_A_priori:list, C_A_priori:list) -> list:
    '''Função do algoritmo MCMC Metropolis Hastings para geração da matriz de covariância W

    W_anterior: float
        Matriz de covariância dos estados não observáveis W no estado anterior
    f: list
        Array com os fatores latentes dimensão: (M+1)x3
    A: list
        Array da matriz de persistência A dimensão: 3x3
    mu: list
        Array com a matriz de média mu dos fatores latentes dimensão: 3x1
    v_priori: float
        Valor do parâmetro v da priori de W
    S_priori_W: list
        Valor do parâmetro S da priori de W
    mu_A_priori: list
      Array com a média a priori de A dimensão: 3x3
    C_A_priori: list
      Array com a primeira matriz de covariância da priori de A dimensão: 3x3
    return = cadeia: list
        Próximo estado da cadeia da geração de W
    '''
    S1 = np.matmul(A - mu_A_priori,np.matmul(np.linalg.inv(C_A_priori),(A - mu_A_priori).T))

    S2 = np.delete(f,(0),axis=0) - np.repeat(mu.reshape(1,3),f.shape[0]-1,axis=0) - np.matmul(np.delete(f,(f.shape[0]-1),axis=0) - np.repeat(mu.reshape(1,3),f.shape[0]-1,axis=0),A.T)
    
    S21 = S2[:,0].reshape(S2.shape[0],1)
    S22 = S2[:,1].reshape(S2.shape[0],1)
    S23 = S2[:,2].reshape(S2.shape[0],1)
    
    S2 = np.einsum('ij->i',np.array([S21**2,S21*S22,S21*S23,S22*S21,S22**2,S22*S23,S23*S21,S23*S22,S23**2]).reshape(9,S2.shape[0])).reshape(3,3)

    proposta = invwishart.rvs(df = (f.shape[0]-1) + v_priori + 3, scale = S_priori_W + S1 + S2)

    alpha = min(1,np.exp(l_g_AW(A,proposta,mu,f) - l_g_AW(A,W_anterior,mu,f)))

    if (np.random.uniform(0,1,1) < alpha):
      cadeia = proposta
    else:
      cadeia = W_anterior
    
    return cadeia

def MH_A(A_anterior:list,f:list,W:list,mu:list,mu_A_priori:list, C_A_priori:list) -> list:
    '''Função do algoritmo MCMC Metropolis Hastings para geração da matriz de persistência A

    A_anterior: float
        Matriz de persistência dos estados não observáveis A no estado anterior
    f: list
        Array com os fatores latentes dimensão: (M+1)x3
    W: list
        Array da matriz de covariância dos estados não observáveis W dimensão: 3x3
    mu: list
        Array com a matriz de média mu dos fatores latentes dimensão: 3x1
    mu_A_priori: list
      Array com a média a priori de A dimensão: 3x3
    C_A_priori: list
      Array com a primeira matriz de covariância da priori de A dimensão: 3x3
    return = cadeia: list
        Próximo estado da cadeia da geração de A
    '''
    proposta = proposta_alvo_A(W,mu,f,mu_A_priori,C_A_priori)

    alpha = min(1,np.exp(l_g_AW(proposta,W,mu,f) - l_g_AW(A_anterior,W,mu,f)))
    
    if ( np.random.uniform(0,1,1) < alpha):
      cadeia = proposta
    else:
      cadeia = A_anterior
    
    return cadeia

def gen_sigma(y:list,f:list,lambdA:float,alpha_priori_sigma:float = 5,beta_priori_sigma:float = 0.2,maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]) -> list:
    '''Função para geração da diagonal da matriz de covariância sigma para as taxas

    y: float
        Array com as taxas utilizadas no modelo dimensão: M x tau
    f: list
        Array com os fatores latentes dimensão: (M+1) x 3
    lambdA: float
        Valor do fator de decaimento lambda
    alpha_priori_sigma: float, default 5
        Valor do alpha a priori da covariância sigma
    beta_priori_sigma: float, default 0.2
        Valor do beta a priori da covariância sigma
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = matriz_sigma: list
        Diagonal da matriz sigma gerada da distribuição a posteriori
    '''
    f = np.delete(f,(0),axis=0)
    matriz_lambdA = matriz_lambda(lambdA,maturidade)
    alpha = y.shape[0] + alpha_priori_sigma/2
    soma = np.einsum('ij->j',(y - np.matmul(f,matriz_lambdA.T))**2)
    beta = beta_priori_sigma/2 + soma

    matriz_sigma = (np.array([np.random.gamma(shape=alpha/2, scale=2/x, size=1) for x in beta])**-1).reshape(1,len(maturidade))
    matriz_sigma = list(np.concatenate(matriz_sigma))

    return matriz_sigma

def gen_mu(f:list,A:list,W:list,mu_priori:list,S_priori:list) -> list:
    '''Função para geração de valores da media dos estados nao observaveis mu

    f: list
        Array com os fatores latentes dimensão: (M+1) x 3
    A: list
        Array da matriz de persistência dos fatores latentes A dimensão: 3 x 3
    W: list
        Array da matriz de covariância dos fatores latentes W dimensão: 3 x 3
    mu_priori: list
        Valor da média a priori de mu dimensão: 3 x 1
    S_priori: list
        Valor da covariância a priori de mu dimensão: 3 x 3
    return = novo_mu: list
        Novo estado da cadeia da geração de mu a partir da posteriori da condicional completa
    '''
    diff = np.delete(f,(0),axis=0) - np.delete(f,(f.shape[0]-1),axis=0)
    soma_para_mi0 = np.matmul(diff, np.matmul((np.identity(3) - A).T,np.linalg.inv(W)).T)
    soma_para_mi0 = np.einsum('ij->j',soma_para_mi0).reshape(3,1)
    
    variancia_mu = np.linalg.inv(np.linalg.inv(S_priori) + (y.shape[0])*np.matmul((np.identity(3) - A).T, np.matmul(np.linalg.inv(W), np.identity(3) - A)) + np.linalg.inv(C0_t(A,W)))
    media_mu = np.matmul(variancia_mu, np.matmul(np.linalg.inv(S_priori), mu_priori) + soma_para_mi0 + np.matmul(np.linalg.inv(C0_t(A,W)), f[0,:].reshape(3,1)))

    novo_mu = np.random.multivariate_normal(np.concatenate(media_mu),variancia_mu,1).reshape(3,1)

    return novo_mu

def f_t(f_anterior:list,A:list,W:list,mu:list) -> list:
    '''Função para calcular o valor do proximo estado nao observavel f na simulação

    f_anterior: list
        Array com os fatores latentes no tempo t-1 dimensão: 1 x 3
    A: list
        Array com a matriz de persistência dos fatores latentes A dimensão: 3 x 3
    W: list
        Array com a matriz de covariância dos fatores latentes W dimensão: 3 x 3
    mu: list
        Array com a matriz de médias dos fatores latentes mu dimensão: 3 x 1
    return = novo_f: list
        Lista com os valores dos fatores latentes no tempo t
    '''
    novo_f = np.random.multivariate_normal(np.concatenate(mu + np.matmul(A, f_anterior.reshape(3,1) - mu)), W, 1).reshape(1,3)
    return novo_f

def y_t(ft:list,lambdA:float,sigma:list,maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]) -> list:
    '''Função para calcular o valor das taxas no momento t na simulação

    ft: list
        Array com os fatores latentes no tempo t dimensão: 1 x 3
    lambdA: float
        Valor do fator de decaimento lambda
    sigma: list
        Lista da diagonal da matriz de covariância das taxas sigma dimensão: tau elementos
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = novo_f: list
        Lista com os valores dos fatores latentes no tempo t
    '''
    novo_y = np.matmul(matriz_lambda(lambdA,maturidade), ft.T) + np.random.multivariate_normal([0]*len(maturidade), np.diag(sigma) , 1).reshape(1,3)
    return novo_y

def FiltroKalman(y:list,A:list,W:list,sigma:list,mu:list,lambdA:float,maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]) -> tuple:
    '''Função para executar o filtro de Kalman

    y: list
        Array com as taxas utilizadas no modelo dimensão: M x tau
    A: list
        Array com a matriz de persistência dos fatores latentes A dimensão: 3 x 3
    W: list
        Array com a matriz de covariância dos fatores latentes W dimensão: 3 x 3
    sigma: list
        Lista da diagonal da matriz de covariância das taxas sigma dimensão: tau elementos
    mu: list
        Array com a matriz de média mu dos fatores latentes dimensão: 3x1
    lambdA: float
        Valor do fator de decaimento lambda
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = (at,Rt,ft,Qt,mt,Ct): tuple
        Tupla com as médias e variâncias das distribuições preditivas e de filtragem do filtro de Kalman
    '''
    Ft = matriz_lambda(lambdA,maturidade)
    C0 = C0_t(A,W)

    at = np.matmul(A,mu).reshape(1,3,1)
    Rt = (np.matmul(A,np.matmul(C0,A.T)) + W).reshape(1,3,3)

    ft = np.matmul(Ft,at[0,:,:]).reshape(1,1,y.shape[1])
    Qt = (np.matmul(Ft,np.matmul(Rt[0,:,:],Ft.T)) + np.diag(sigma)).reshape(1,y.shape[1],y.shape[1])

    erro = (y[0,:] - ft[0,0,:]).reshape(1,1,y.shape[1])
    mt = (at[0,:,:] + Rt[0,:,:] @ Ft.T @ np.linalg.inv(Qt[0,:,:]) @ erro[0,:,:].T).reshape(1,3,1)
    Ct = (Rt[0,:,:] - np.matmul(Rt[0,:,:], np.matmul(Ft.T, np.matmul(np.linalg.inv(Qt[0,:,:]), np.matmul(Ft, Rt[0,:,:]))))).reshape(1,3,3)

    for it in range(1,y.shape[0]):
        
        at = np.append(at,np.matmul(A, mt[it-1,:,:]).reshape(1,3,1),axis=0)
        Rt = np.append(Rt,(np.matmul(A, np.matmul(Ct[it-1,:,:], A.T)) + W).reshape(1,3,3),axis=0)

        ft = np.append(ft, np.matmul(Ft, at[it,:,:]).reshape(1,1,y.shape[1]),axis=0)
        Qt = np.append(Qt,(np.matmul(Ft, np.matmul(Rt[it,:,:], Ft.T)) + np.diag(sigma)).reshape(1,y.shape[1],y.shape[1]),axis=0)

        erro = np.append(erro,(y[it,:] - ft[it,:,:]).reshape(1,1,y.shape[1]),axis=0)
        
        mt = np.append(mt,(at[it,:,:] + Rt[it,:,:] @ Ft.T @ np.linalg.inv(Qt[it,:,:]) @ erro[it,:,:].T).reshape(1,3,1),axis=0)
        Ct = np.append(Ct, (Rt[it,:,:] - np.matmul(Rt[it,:,:], np.matmul(Ft.T, np.matmul(np.linalg.inv(Qt[it,:,:]), np.matmul(Ft, Rt[it,:,:]))))).reshape(1,3,3),axis=0)

    return (at,Rt,ft,Qt,mt,Ct)

def FFBS(y:list,A:list,W:list,sigma:list,mu:list,lambdA:float,maturidade:list = [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]):
    '''Função para executar o Forward Filtering Backward Sampling FFBS para gerar os estados não observáveis f

    y: list
        Array com as taxas utilizadas no modelo dimensão: M x tau
    A: list
        Array com a matriz de persistência dos fatores latentes A dimensão: 3 x 3
    W: list
        Array com a matriz de covariância dos fatores latentes W dimensão: 3 x 3
    sigma: list
        Lista da diagonal da matriz de covariância das taxas sigma dimensão: tau elementos
    mu: list
        Array com a matriz de média mu dos fatores latentes dimensão: 3x1
    lambdA: float
        Valor do fator de decaimento lambda
    maturidade: list, default [3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120]
        Maturidades escolhidas para realização do cálculo dimensão: tau elementos
    return = (at,Rt,ft,Qt,mt,Ct): tuple
        Tupla com as médias e variâncias das distribuições preditivas e de filtragem do filtro de Kalman
    '''
    Ft = matriz_lambda(lambdA,maturidade)
    C0 = C0_t(A,W)

    filtro = FiltroKalman(y,A,W,sigma,mu,lambdA,maturidade)

    at = filtro[0]
    Rt = filtro[1]
    mt = filtro[4]
    Ct = filtro[5]

    t_max = at.shape[0]-1

    theta = np.random.multivariate_normal(np.concatenate(mt[t_max,:,:]),Ct[t_max,:,:],1).reshape(1,3)

    for it in range(t_max-1,-1,-1):
        ht = mt[it,:,:] + np.matmul(np.matmul(np.matmul(Ct[it,:,:],A.T), np.linalg.inv(Rt[it+1,:,:])), (theta[theta.shape[0]-1,:] - at[it+1,:,:].T).T)
        Ht = Ct[it,:,:] - np.matmul(np.matmul(np.matmul(np.matmul(Ct[it,:,:], A.T), np.linalg.inv(Rt[it+1,:,:])), A), Ct[it,:,:])

        theta = np.append(theta, np.random.multivariate_normal(np.concatenate(ht),Ht,1).reshape(1,3), axis=0)
    
    h0 = mu + np.matmul(np.matmul(np.matmul(C0,A.T), np.linalg.inv(Rt[0,:,:])), (theta[theta.shape[0]-1,:] - at[it+1,:,:].T).T)
    H0 = C0 - C0 @ A.T @ np.linalg.inv(Rt[0,:,:]) @ A @ C0

    theta = np.append(theta, np.random.multivariate_normal(np.concatenate(h0),H0,1).reshape(1,3), axis=0)

    return theta[::-1,:]

######################################################################################################################################################################################
######################################################### APLICAÇÃO DO MODELO ########################################################################################################
######################################################################################################################################################################################

#VALORES PARA COMEÇAR O ALGORITMO DE GIBBS

lambda_gibbs = 0.04
mu_gibbs = np.array([7,-3,4]).reshape(3,1)
sigma_gibbs = [1]*len(maturidade)
W_gibbs = 3*np.identity(3)
A_gibbs = np.array([0.5, -0.1, 0.04, 0.005, 0.5, 0.0002, -0.03, 0.09, 0.5]).reshape(3,3).T

#ALGORITMO DE GIBBS

np.random.seed(458)

for it in range(10**5+1):
    
    ffbs = FFBS(y,A_gibbs,W_gibbs,sigma_gibbs,mu_gibbs,lambda_gibbs,maturidade)
    lambda_gibbs = MH_lambda(lambda_gibbs,y,ffbs,sigma_gibbs,alpha_priori_lambda, beta_priori_lambda, maturidade)
    mu_gibbs = gen_mu(ffbs,A_gibbs,W_gibbs,mu_priori,S_priori)
    sigma_gibbs = gen_sigma(y,ffbs,lambda_gibbs,alpha_priori_sigma,beta_priori_sigma,maturidade)
    W_gibbs = MH_W(W_gibbs,ffbs,A_gibbs,mu_gibbs,v_priori,S_priori_W, mu_A_priori, C_A_priori)
    A_gibbs = MH_A(A_gibbs,ffbs,W_gibbs,mu_gibbs,mu_A_priori,C_A_priori)
    
    if (it > 20000) and (it%16 == 0):
        open(diretorio + '\\lambda.txt', 'a').write(" ".join(map(str, lambda_gibbs))+'\n')
        open(diretorio + '\\mu.txt', 'a').write(" ".join(map(str, list(np.concatenate(mu_gibbs))))+'\n')
        open(diretorio + '\\sigma.txt', 'a').write(" ".join(map(str, sigma_gibbs))+'\n')
        open(diretorio + '\\W.txt', 'a').write(" ".join(map(str, list(np.concatenate(W_gibbs))))+'\n')
        open(diretorio + '\\A.txt', 'a').write(" ".join(map(str, list(np.concatenate(A_gibbs))))+'\n')
        open(diretorio + '\\Nivel.txt', 'a').write(" ".join(map(str, list(ffbs[:,0])))+'\n')
        open(diretorio + '\\Inclinacao.txt', 'a').write(" ".join(map(str, list(ffbs[:,1])))+'\n')
        open(diretorio + '\\Curvatura.txt', 'a').write(" ".join(map(str, list(ffbs[:,2])))+'\n')
        print(it)