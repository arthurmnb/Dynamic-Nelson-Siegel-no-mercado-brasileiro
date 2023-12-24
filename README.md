PORTUGUES:

Neste repositório estão os códigos para a geração dos parâmetros e variáveis a posteriori do modelo dinâmico de Nelson-Siegel.

O arquivo "Modelo Nelson e Siegel Dinamico versão normal.R" é a rotina em R para a geração de 5.000 amostras a posteriori a partir de 100.000 iterações utilizando o modelo dinâmico padrão proposto por Diebold e Li.

O arquivo "Modelo Nelson e Siegel Dinamico versão t-student.R" é a rotina em R para a geração de 5.000 amostras a posteriori a partir de 100.000 iterações utilizando o modelo com cauda pesada t-student.

O arquivo "Web Scraping Site B3.R" é a rotina de obtenção de dados financeiros no site da B3, ele está programado para retirar os dados de contratos futuros DI1 e de cupom cambial.



ENGLISH

This repository contains code for the generation of the posteriori distribution of the variables and parameters of the dynamic Nelson-Siegel model.

The file "Modelo Nelson e Siegel Dinamico versão normal.R" generates 5,000 samples of the posteriori distribution of the gaussian dynamic Nelson-Siegel model originally proposed by Diebold and Li after 100,000 iteractions.

The file "Modelo Nelson e Siegel Dinamico versão t-student.R" generates 5,000 samples of the posteriori distribution of the heavy tailed t-student extension of dynamic Nelson-Siegel model after 100,000 iteractions.

The file "Web Scraping Site B3.R" is for web scrapping brazillian market data to create a database to apply the model.
