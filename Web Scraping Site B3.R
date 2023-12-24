library(rvest)

DIRETORIO <- "DIRETORIO" #diretorio para salvar a base de dados, incluir "\\" no fim

anos <- seq(2003,2023,1) #Anos do periodo em que a busca sera realizada

variaveis.interesse <- c("DDI   - Cupom cambial", "DI1   - DI de 1 dia")  #Escolha dos instrumentos financeiros para realizar web scrapping



#Função para retirada da tabela das variaveis selecionadas no dia Data
web.scraping <- function(Data,variaveis.interesse){

  pagina <- html_form(read_html("https://www2.bmf.com.br/pages/portal/bmfbovespa/lumis/lum-ajustes-do-pregao-ptBR.asp"))[[1]] #puxa o primeiro form do html para selecionar os inputs, lendo a variavel vemos que a data é chamada pelo input dData1
  filtro <- pagina %>% html_form_set(dData1 = Data) #Define o input dData1 como a data de interesse
  resp <- html_form_submit(filtro) #submete o input da data de interesse na pagina
  pagina <- read_html(resp) #recarrega a pagina com a data de interesse

  summaries_xpath <- pagina %>%
    html_elements(xpath = '/html/body/div[1]/div[2]/div/table') #puxa o XPath da tabela de dados

  consulta <- as.data.frame(html_table(summaries_xpath, header = TRUE)) #faz a leitura da tabela de dados
  consulta$Data <- rep(Data,nrow(consulta)) #Criando a coluna de data para incluir na base final
  consulta$Mercadoria[consulta$Mercadoria==""] <- NA #Transformando em NA as linhas faltantes em mercadoria pra preencher no proximo comando
  consulta$Mercadoria <- zoo::na.locf(consulta$Mercadoria) #preenchendo as mercadorias vazias com a mercadoria anterior
  consulta <- consulta[consulta$Mercadoria %in% variaveis.interesse,] #Filtrando as variaveis nas variaveis de interesse
  
  return(consulta)
}


#Busca ano a ano
for (ano in anos){
  #lista com as datas no ano
  lista.datas <- chartr("-","/",format(seq(as.Date(paste0(ano,"-01-01")), as.Date(paste0(ano,"-12-31")), by="days"), format="%d-%m-%Y"))
  #base de dados para salvar consultas
  base <- data.frame(Mercadoria = c(), Vencimento = c(), Preço.de.ajuste.anterior = c(), Preço.de.ajuste.Atual = c() ,Variação = c(), Valor.do.ajuste.por.contrato..R.. = c(), Data = c())
  #Busca em cada data
  for (Data in lista.datas){
    consulta <- web.scraping(Data,variaveis.interesse)
    #verifica se a busca retornou com base no criterio de numero de colunas, se no futuro alterarem as colunas no site pode ser necessario realizar mudanças aqui
    if (ncol(consulta) == 7){
      base <- rbind(base,consulta)
    }
  }
  #Salva a base de dados anual
  write.table(base,paste0(DIRETORIO,ano," - BASE.txt"), sep = ";", dec = ".")

}

