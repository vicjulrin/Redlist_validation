---
title: "Validación RLI - IAvH"
author: 
  - name: "Diaz Corzo Camila; Rincon-Parra VJ"
    email: "rincon-v@javeriana.edu.co"
output: 
  github_document:
    md_extension: +gfm_auto_identifiers
    preserve_yaml: true
    toc: true
    toc_depth: 6
---

Validación RLI - IAvH
================

Este documento detalla la validación de los resultados para la
estimación del Red List Index (RLI) realizada por el Instituto de
Investigación de Recursos Biológicos Alexander von Humboldt (IAvH) en R
software
([mbi-colombia/RLI](https://github.com/PEM-Humboldt/mbi-colombia/tree/main/MBI/BI_subindex/RLI),
y
[biab/RLI_pipeline](https://github.com/PEM-Humboldt/biab-2.0/tree/RLI_pipeline)).
El RLI es un indicador que mide el riesgo de extinción de las especies a
lo largo del tiempo, y su cálculo se basa en los cambios en las
categorías de la Lista Roja de la UICN [(IUCN,
2024)](https://www.iucnredlist.org/assessment/red-list-index). La
validación se realizó revisando cada parte del código desarrollado por
el IAvH, desde la organización de los datos hasta la ejecución de la
ecuación del indicador, contrastando tanto los resultados obtenidos a
través de la librería red como mediante un proceso manual.

El ejemplo documentado estima el RLI para grupos taxonómicos de especies
clasificados y evaluados por la IUCN, utilizando la API de la libreria
[‘rredlist IUCN’ Red List
Client](https://cran.r-project.org/web/packages/rredlist/index.html)
para consultar sus bases de datos, y estimando el indice a través de la
librería [red IUCN Redlisting
Tools](https://cran.r-project.org/web/packages/red/index.html).


- [Cargar librerias/paquetes necesarios para el
  análisis](#cargar-libreriaspaquetes-necesarios-para-el-análisis)
- [Definir inputs](#definir-inputs)
- [Check lista de datos por taxon](#check-lista-de-datos-por-taxon)
- [Obtener evaluacion historica de las especies listadas por
  taxon](#obtener-evaluacion-historica-de-las-especies-listadas-por-taxon)
- [Ajustar como matriz](#ajustar-como-matriz)
- [Corregir categorias](#corregir-categorias)
- [Ajustar matriz](#ajustar-matriz)
- [Validación de resultados](#validación-de-resultados)
  - [Resultados de la libreria red](#resultados-de-la-libreria-red)
  - [Resultados - ejecución manual (Butchart et al, 2004;
    2007)](#resultados---ejecución-manual-butchart-et-al-2004-2007)
  - [Comparar resultados](#comparar-resultados)

### Cargar librerias/paquetes necesarios para el análisis

``` r
# Load libraries ####
packages_list<- c("this.path", "dplyr","terra","red", "rredlist","ggplot2", "pbapply", "tibble", "plyr", "ggpubr")
packagesPrev<- .packages(all.available = TRUE)
lapply(packages_list, function(x) {   if ( ! x %in% packagesPrev ) { install.packages(x, force=T)}    })
lapply(packages_list, library, character.only = TRUE)
```

### Definir inputs

Los inputs del código son un token de acceso válido otorgado por la UICN
(ej. `token= v11xxx22`), el nombre de un grupo taxonómico listado por la
UICN (ej. `taxonomic_group= "crocodiles_and_alligators"`) y el pais de
interes para el calculo (ej. `country_name<- Colombia`).

``` r
# Define inputs ####
output<- file.path(dirname(this.path::this.path()), "output"); dir.create(output)
token <- "v11xxx22"  # Token IUCN
country_name<- "Colombia" # Pais de interes
taxonomic_group<- "crocodiles_and_alligators" # grupo taxonomico de interes
```

### Check lista de datos por taxon

``` r
## Load sp country ####
UICN_isocode <- rredlist::rl_countries(key = token)$results %>% dplyr::filter(country == country_name) %>% {.$isocode}
UICN_country <- rredlist::rl_sp_country(country= UICN_isocode, key = token)$result

## Load sp taxonomic group ####
UICN_taxon <- rredlist::rl_comp_groups(group = taxonomic_group, key = token)$result

## Filter country list by taxonomic group ####
IUCN_sp_Taxon<- UICN_taxon %>% dplyr::filter(taxonid %in% UICN_country$taxonid)
```

|     | taxonid | scientific_name         | subspecies | rank | subpopulation | category |
|:----|--------:|:------------------------|-----------:|-----:|--------------:|:---------|
| 1   |   46584 | Caiman crocodilus       |         NA |   NA |            NA | LC       |
| 2   |    5659 | Crocodylus acutus       |         NA |   NA |            NA | VU       |
| 3   |    5661 | Crocodylus intermedius  |         NA |   NA |            NA | CR       |
| 4   |   13053 | Melanosuchus niger      |         NA |   NA |            NA | LR/cd    |
| 5   |   46587 | Paleosuchus palpebrosus |         NA |   NA |            NA | LC       |
| 6   |   46588 | Paleosuchus trigonatus  |         NA |   NA |            NA | LC       |

### Obtener evaluacion historica de las especies listadas por taxon

``` r
historyAssesment_data <- iucn_history_assessment_data <- pbapply::pblapply(IUCN_sp_Taxon[, "scientific_name"], function(x) {
  tryCatch({
    rredlist::rl_history(name = x, key = token)$result %>% dplyr::mutate(scientific_name= x) 
    }, error = function(e) {NULL})
}) %>% plyr::rbind.fill() %>% list(IUCN_sp_Taxon) %>% plyr::join_all(match = "first")
print(historyAssesment_data)
```

|     | category                          | scientific_name         | year | assess_year | code  | taxonid | subspecies | rank | subpopulation |
|:----|:----------------------------------|:------------------------|:-----|:------------|:------|--------:|-----------:|-----:|--------------:|
| 1   | Least Concern                     | Caiman crocodilus       | 2019 | 2016        | LC    |      NA |         NA |   NA |            NA |
| 2   | Lower Risk/least concern          | Caiman crocodilus       | 1996 | 1996        | LR/lc |      NA |         NA |   NA |            NA |
| 3   | Threatened                        | Caiman crocodilus       | 1988 | 1988        | T     |      NA |         NA |   NA |            NA |
| 4   | Threatened                        | Caiman crocodilus       | 1986 | 1986        | T     |      NA |         NA |   NA |            NA |
| 5   | Vulnerable                        | Crocodylus acutus       | 2022 | 2020        | VU    |      NA |         NA |   NA |            NA |
| 6   | Vulnerable                        | Crocodylus acutus       | 2021 | 2020        | VU    |      NA |         NA |   NA |            NA |
| 7   | Vulnerable                        | Crocodylus acutus       | 2012 | 2009        | VU    |      NA |         NA |   NA |            NA |
| 8   | Vulnerable                        | Crocodylus acutus       | 1996 | 1996        | VU    |      NA |         NA |   NA |            NA |
| 9   | Vulnerable                        | Crocodylus acutus       | 1994 | 1994        | V     |      NA |         NA |   NA |            NA |
| 10  | Endangered                        | Crocodylus acutus       | 1990 | 1990        | E     |      NA |         NA |   NA |            NA |
| 11  | Endangered                        | Crocodylus acutus       | 1988 | 1988        | E     |      NA |         NA |   NA |            NA |
| 12  | Endangered                        | Crocodylus acutus       | 1986 | 1986        | E     |      NA |         NA |   NA |            NA |
| 13  | Endangered                        | Crocodylus acutus       | 1982 | 1982        | E     |      NA |         NA |   NA |            NA |
| 14  | Critically Endangered             | Crocodylus intermedius  | 2018 | 2017        | CR    |      NA |         NA |   NA |            NA |
| 15  | Critically Endangered             | Crocodylus intermedius  | 1996 | 1996        | CR    |      NA |         NA |   NA |            NA |
| 16  | Endangered                        | Crocodylus intermedius  | 1994 | 1994        | E     |      NA |         NA |   NA |            NA |
| 17  | Endangered                        | Crocodylus intermedius  | 1990 | 1990        | E     |      NA |         NA |   NA |            NA |
| 18  | Endangered                        | Crocodylus intermedius  | 1988 | 1988        | E     |      NA |         NA |   NA |            NA |
| 19  | Endangered                        | Crocodylus intermedius  | 1986 | 1986        | E     |      NA |         NA |   NA |            NA |
| 20  | Endangered                        | Crocodylus intermedius  | 1982 | 1982        | E     |      NA |         NA |   NA |            NA |
| 21  | Lower Risk/conservation dependent | Melanosuchus niger      | 2000 | 2000        | LR/cd |      NA |         NA |   NA |            NA |
| 22  | Endangered                        | Melanosuchus niger      | 1996 | 1996        | EN    |      NA |         NA |   NA |            NA |
| 23  | Vulnerable                        | Melanosuchus niger      | 1994 | 1994        | V     |      NA |         NA |   NA |            NA |
| 24  | Endangered                        | Melanosuchus niger      | 1990 | 1990        | E     |      NA |         NA |   NA |            NA |
| 25  | Endangered                        | Melanosuchus niger      | 1988 | 1988        | E     |      NA |         NA |   NA |            NA |
| 26  | Endangered                        | Melanosuchus niger      | 1986 | 1986        | E     |      NA |         NA |   NA |            NA |
| 27  | Endangered                        | Melanosuchus niger      | 1982 | 1982        | E     |      NA |         NA |   NA |            NA |
| 28  | Least Concern                     | Paleosuchus palpebrosus | 2019 | 2018        | LC    |      NA |         NA |   NA |            NA |
| 29  | Lower Risk/least concern          | Paleosuchus palpebrosus | 1996 | 1996        | LR/lc |      NA |         NA |   NA |            NA |
| 30  | Least Concern                     | Paleosuchus trigonatus  | 2019 | 2018        | LC    |      NA |         NA |   NA |            NA |
| 31  | Lower Risk/least concern          | Paleosuchus trigonatus  | 1996 | 1996        | LR/lc |      NA |         NA |   NA |            NA |

### Ajustar como matriz

``` r
adjust_categories<- data.frame(Cat_IUCN= c("CR", NA, "EN", "EN", NA, NA, "LC", "LC", "LC", "NT",  "NT", "RE", "VU", "VU"),
                               code= c("CR", "DD", "E", "EN", "I", "K", "LC", "LR/cd", "LR/lc", "LR/nt",  "NT", "R", "V", "VU"))
print(adjust_categories)
```

### Corregir categorias

|     | Cat_IUCN | code  |
|:----|:---------|:------|
| 1   | CR       | CR    |
| 2   | NA       | DD    |
| 3   | EN       | E     |
| 4   | EN       | EN    |
| 5   | NA       | I     |
| 6   | NA       | K     |
| 7   | LC       | LC    |
| 8   | LC       | LR/cd |
| 9   | LC       | LR/lc |
| 10  | NT       | LR/nt |
| 11  | NT       | NT    |
| 12  | RE       | R     |
| 13  | VU       | V     |
| 14  | VU       | VU    |

``` r
form_matrix <- as.formula(paste0("scientific_name", "~", "assess_year"))


historyAssesment_matrix <-   reshape2::dcast(historyAssesment_data, form_matrix,  value.var = "code",
                                            fun.aggregate = function(x) {unique(x)[1]}) %>% tibble::column_to_rownames("scientific_name") %>% as.data.frame.matrix()

RedList_matrix<- historyAssesment_matrix %>% as.matrix()

for(i in seq(nrow(adjust_categories))){
  RedList_matrix[ which(RedList_matrix== adjust_categories[i,]$code, arr.ind = TRUE) ]<- adjust_categories[i,]$Cat_IUCN 
}

for(j in unique(adjust_categories$Cat_IUCN)){
  key<- c(tolower(j), toupper(j), j) %>% paste0(collapse = "|")
  RedList_matrix[ which(grepl(key, RedList_matrix), arr.ind = T) ]    <- j
}

RedList_matrix[which( (!RedList_matrix %in% adjust_categories$Cat_IUCN)  & !is.na(RedList_matrix) , arr.ind = TRUE )]<-NA
print(RedList_matrix)
```

|                         | 1982 | 1986 | 1988 | 1990 | 1994 | 1996 | 2000 | 2009 | 2016 | 2017 | 2018 | 2020 |
|:------------------------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|
| Caiman crocodilus       | NA   | NA   | NA   | NA   | NA   | LC   | NA   | NA   | LC   | NA   | NA   | NA   |
| Crocodylus acutus       | EN   | EN   | EN   | EN   | VU   | VU   | NA   | VU   | NA   | NA   | NA   | VU   |
| Crocodylus intermedius  | EN   | EN   | EN   | EN   | EN   | CR   | NA   | NA   | NA   | CR   | NA   | NA   |
| Melanosuchus niger      | EN   | EN   | EN   | EN   | VU   | EN   | LC   | NA   | NA   | NA   | NA   | NA   |
| Paleosuchus palpebrosus | NA   | NA   | NA   | NA   | NA   | LC   | NA   | NA   | NA   | NA   | LC   | NA   |
| Paleosuchus trigonatus  | NA   | NA   | NA   | NA   | NA   | LC   | NA   | NA   | NA   | NA   | LC   | NA   |

### Ajustar matriz

``` r
### Eliminar las especies que no aportan a la estimación. Menos de dos años con dato de categorización  ####
RedList_matrix_2<- RedList_matrix[rowSums(!is.na(RedList_matrix))>=2,]

### Asignar los años no evaluados NA con datos del año de ultima evaluacion ####
replace_na_with_previous <- function(df, target_col) {
  for (col in (target_col-1):2) {
    df[[target_col]] <- ifelse(is.na(df[[target_col]]), df[[col]], df[[target_col]])
  }; return(df) }

df = RedList_matrix_2 %>% as.data.frame.matrix()
for(k in 2:ncol(RedList_matrix_2)){ df <- replace_na_with_previous(df, k) }

print(df)
```

|                         | 1982 | 1986 | 1988 | 1990 | 1994 | 1996 | 2000 | 2009 | 2016 | 2017 | 2018 | 2020 |
|:------------------------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|
| Caiman crocodilus       | NA   | NA   | NA   | NA   | NA   | LC   | LC   | LC   | LC   | LC   | LC   | LC   |
| Crocodylus acutus       | EN   | EN   | EN   | EN   | VU   | VU   | VU   | VU   | VU   | VU   | VU   | VU   |
| Crocodylus intermedius  | EN   | EN   | EN   | EN   | EN   | CR   | CR   | CR   | CR   | CR   | CR   | CR   |
| Melanosuchus niger      | EN   | EN   | EN   | EN   | VU   | EN   | LC   | LC   | LC   | LC   | LC   | LC   |
| Paleosuchus palpebrosus | NA   | NA   | NA   | NA   | NA   | LC   | LC   | LC   | LC   | LC   | LC   | LC   |
| Paleosuchus trigonatus  | NA   | NA   | NA   | NA   | NA   | LC   | LC   | LC   | LC   | LC   | LC   | LC   |

## Validación de resultados

Para validar los resultados obtenidos con el paquete red, primero
ejecutamos la función red::rli, que calcula el Red List Index (RLI) de
acuerdo con los datos proporcionados y los métodos implementados en
dicho paquete. Posteriormente, estimamos manualmente el RLI utilizando
la ecuación descrita por [Butchart et al. (2004,
2007)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000140),
para verificar la coherencia entre ambos enfoques.

### Resultados de la libreria red

``` r
red_resuLt<- red::rli(df_filtered) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("Year") %>%  setNames(c("Year", "RLI_red")) %>% 
  dplyr::filter(!Year %in% "Change/year")

print(red_resuLt)
```

|     | Year |   RLI_red |
|:----|:-----|----------:|
| 1   | 1982 | 0.4000000 |
| 2   | 1986 | 0.4000000 |
| 3   | 1988 | 0.4000000 |
| 4   | 1990 | 0.4000000 |
| 5   | 1994 | 0.5333333 |
| 6   | 1996 | 0.7000000 |
| 7   | 2000 | 0.8000000 |
| 8   | 2009 | 0.8000000 |
| 9   | 2016 | 0.8000000 |
| 10  | 2017 | 0.8000000 |
| 11  | 2018 | 0.8000000 |
| 12  | 2020 | 0.8000000 |

``` r
red_resuLt_plot<- ggplot(red_resuLt, aes(x = Year, y = RLI_red)) +
  geom_line(group = 1, col= "red") +
  geom_point() +
  coord_cartesian(ylim = c(0,1))+
  theme_classic()+
  theme(    panel.grid.major = element_line(color = "gray"),
  ) + theme(text = element_text(size = 4)) + ggtitle("red_results")

print(red_resuLt_plot)
```

![](README_figures/red_resuLt_plot.jpg)

### Resultados - ejecución manual (Butchart et al, 2004; 2007)

La formula propuesta por este autor es una modificación de la
originalmente propuesta y asume que las especies categorizadas como
“data deficient” (DD) no son tomados en cuenta para el cálculo.
Adicionalmente, cambia respecto a la ecuación original en que el
resultado del tiempo anterior no es tenido en cuenta. Matemáticamente la
nueva equación se expresa así:

$$
RLI_{t} = \frac{(M - T_{t})}{M}
$$

Dónde M representa el máximo valor de amenaza, el cual se expresa de la
siguiente manera:

$$
M = W_{EX}*N
$$

Para lo cual:

N = numero de especies evaluadas = nrow(df_weights) calcula el número
total de especies evaluadas, es decir, el número de filas en la matriz
df_weights, que representa las especies y los pesos de las categorías de
amenaza a lo largo del tiempo.

W<sub>EX</sub> = Representa el peso maximo (5) de la categoría más alta
de amenaza, que es “Extinto” (EX).

En la ecuación general T hace referencia a el valor de amenaza actual y
es definido de la siguiente forma:

$$
T_t=\sum_{s}{W_{c(t,s)}}
$$

Para lo cual:

W<sub>c</sub> = Hace referencia al peso que se le asigna a cada
categoria de amenaza el cual se obtiene usando: colSums(df_weights) al
sumar los valores de los pesos de las categorías de amenaza para todas
las especies en cada año.

s = Especies

t = Tiempo

``` r
## Definir los pesos para cada categoría ####
category_weights <- c("LC" = 0, "NT" = 1, "VU" = 2, "EN" = 3, "CR" = 4, "DD" = NA)  # DD es NA porque no se incluye

### Reemplazar las categorías por sus pesos correspondientes ####
df_weights <- apply(df_filtered, 2, function(column) category_weights[column]) %>% as.data.frame.matrix() %>% 
  dplyr::mutate(ID_Fila= rownames(df_filtered)) %>% tibble::column_to_rownames(var = "ID_Fila")

print(df_weights)
```

|                         | 1982 | 1986 | 1988 | 1990 | 1994 | 1996 | 2000 | 2009 | 2016 | 2017 | 2018 | 2020 |
|:------------------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| Caiman crocodilus       |   NA |   NA |   NA |   NA |   NA |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| Crocodylus acutus       |    3 |    3 |    3 |    3 |    2 |    2 |    2 |    2 |    2 |    2 |    2 |    2 |
| Crocodylus intermedius  |    3 |    3 |    3 |    3 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
| Melanosuchus niger      |    3 |    3 |    3 |    3 |    2 |    3 |    0 |    0 |    0 |    0 |    0 |    0 |
| Paleosuchus palpebrosus |   NA |   NA |   NA |   NA |   NA |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| Paleosuchus trigonatus  |   NA |   NA |   NA |   NA |   NA |    0 |    0 |    0 |    0 |    0 |    0 |    0 |

``` r
### Ecuacion RLI Butchardt ####
#### Definir las variables de las ecuaciones ####
vars_period<- sapply(names(df_weights), function(y) {
  N= sum(!is.na(df_weights[,y])) # Numero de especies con evaluación
  W= 5 # peso maximo entre las especies evaluadas 5 para EX
  M = W * N # Maximum threath score
  Wc_ts<- sum(df_weights[,y], na.rm=T) # Suma del peso de la categoría de amenaza c para la especie s en el tiempo t
data.frame(N=N, W=W, M= M, Wc_ts=Wc_ts)
}) %>% as.data.frame()
print(vars_period)
```

|       | 1982 | 1986 | 1988 | 1990 | 1994 | 1996 | 2000 | 2009 | 2016 | 2017 | 2018 | 2020 |
|:------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|
| N     | 3    | 3    | 3    | 3    | 3    | 6    | 6    | 6    | 6    | 6    | 6    | 6    |
| W     | 5    | 5    | 5    | 5    | 5    | 5    | 5    | 5    | 5    | 5    | 5    | 5    |
| M     | 15   | 15   | 15   | 15   | 15   | 30   | 30   | 30   | 30   | 30   | 30   | 30   |
| Wc_ts | 9    | 9    | 9    | 9    | 7    | 9    | 6    | 6    | 6    | 6    | 6    | 6    |

``` r
#### Ejecutar ecuacion ####
RLI_Butchart_2007<-pblapply(Wc_ts, function(t) { (M-t)/M }) %>% unlist()
Butchart_result<- data.frame(Year= names(RLI_Butchart_2007), RLI_Butchart = (RLI_Butchart_2007))
print(Butchart_result)
```

|      | Year | RLI_Butchart |
|:-----|:-----|-------------:|
| 1982 | 1982 |    0.4000000 |
| 1986 | 1986 |    0.4000000 |
| 1988 | 1988 |    0.4000000 |
| 1990 | 1990 |    0.4000000 |
| 1994 | 1994 |    0.5333333 |
| 1996 | 1996 |    0.7000000 |
| 2000 | 2000 |    0.8000000 |
| 2009 | 2009 |    0.8000000 |
| 2016 | 2016 |    0.8000000 |
| 2017 | 2017 |    0.8000000 |
| 2018 | 2018 |    0.8000000 |
| 2020 | 2020 |    0.8000000 |

``` r
Butchart_resuLt_plot<-  ggplot(Butchart_result, aes(x = Year, y = RLI_Butchart)) +
  geom_line(group = 1, col= "blue") +
  geom_point() +
  coord_cartesian(ylim = c(0,1))+
  theme_classic()+
  theme(    panel.grid.major = element_line(color = "gray"),
  ) + theme(text = element_text(size = 4)) + ggtitle("Butchart_result")

print(Butchart_resuLt_plot)
```

![](README_figures/Butchart_resuLt_plot.jpg)

### Comparar resultados

``` r
compare_data<- list(red_resuLt, Butchart_result) %>% plyr::join_all()
print(compare_data)
```

|     | Year |   RLI_red | RLI_Butchart |
|:----|:-----|----------:|-------------:|
| 1   | 1982 | 0.4000000 |    0.4000000 |
| 2   | 1986 | 0.4000000 |    0.4000000 |
| 3   | 1988 | 0.4000000 |    0.4000000 |
| 4   | 1990 | 0.4000000 |    0.4000000 |
| 5   | 1994 | 0.5333333 |    0.5333333 |
| 6   | 1996 | 0.7000000 |    0.7000000 |
| 7   | 2000 | 0.8000000 |    0.8000000 |
| 8   | 2009 | 0.8000000 |    0.8000000 |
| 9   | 2016 | 0.8000000 |    0.8000000 |
| 10  | 2017 | 0.8000000 |    0.8000000 |
| 11  | 2018 | 0.8000000 |    0.8000000 |
| 12  | 2020 | 0.8000000 |    0.8000000 |

``` r
compare_plot<- ggpubr::ggarrange(plotlist= list(red_resuLt_plot, Butchart_resuLt_plot), ncol= 2)
print(compare_plot)
```

![](README_figures/compare_plot.jpg)
