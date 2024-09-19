# Load libraries ####
packages_list<- c("this.path", "dplyr","terra","red", "rredlist","ggplot2", "pbapply", "tibble", "plyr", "ggpubr")
packagesPrev<- .packages(all.available = TRUE)
lapply(packages_list, function(x) {   if ( ! x %in% packagesPrev ) { install.packages(x, force=T)}    })
lapply(packages_list, library, character.only = TRUE)

# Define inputs ####
##  IUCN token ####
output<- file.path(dirname(this.path::this.path()), "output"); dir.create(output)
token <- "f33e69dfa9b06a6495aca3b049606f6e08ceb37083ff88a9c7c9dfbcd56a9121" # Token IUCN
country_name<- "Colombia"
taxonomic_group<- "crocodiles_and_alligators" # grupo taxonomico de interes


# Script ####
## Load sp country ####
UICN_isocode <- rredlist::rl_countries(key = token)$results %>% dplyr::filter(country == country_name) %>% {.$isocode}
UICN_country <- rredlist::rl_sp_country(country= UICN_isocode, key = token)$result

## Load sp taxonomic group ####
UICN_taxon <- rredlist::rl_comp_groups(group = taxonomic_group, key = token)$result

## Filter country list by taxonomic group ####
IUCN_sp_Taxon<- UICN_taxon %>% dplyr::filter(taxonid %in% UICN_country$taxonid)



## Obtener evaluacion historica de las especies listadas por taxon ####
historyAssesment_data <- iucn_history_assessment_data <- pbapply::pblapply(IUCN_sp_Taxon[, "scientific_name"], function(x) {
  tryCatch({
    rredlist::rl_history(name = x, key = token)$result %>% dplyr::mutate(scientific_name= x) 
    }, error = function(e) {NULL})
}) %>% plyr::rbind.fill() %>% list(IUCN_sp_Taxon) %>% plyr::join_all(match = "first")



## Ajustar como matriz ####
form_matrix <- as.formula(paste0("scientific_name", "~", "assess_year"))
historyAssesment_matrix <-   reshape2::dcast(historyAssesment_data, form_matrix,  value.var = "code",
                                             fun.aggregate = function(x) {unique(x)[1]}) %>% tibble::column_to_rownames("scientific_name") %>% as.data.frame.matrix()

## Corregir categorias ####
adjust_categories<- data.frame(Cat_IUCN= c("CR", NA, "EN", "EN", NA, NA, "LC", "LC", "LC", "NT",  "NT", "RE", "VU", "VU"),
                               code= c("CR", "DD", "E", "EN", "I", "K", "LC", "LR/cd", "LR/lc", "LR/nt",  "NT", "R", "V", "VU"))



RedList_matrix<- historyAssesment_matrix %>% as.matrix()

for(i in seq(nrow(adjust_categories))){
  RedList_matrix[ which(RedList_matrix== adjust_categories[i,]$code, arr.ind = TRUE) ]<- adjust_categories[i,]$Cat_IUCN 
}

for(j in unique(adjust_categories$Cat_IUCN)){
  key<- c(tolower(j), toupper(j), j) %>% paste0(collapse = "|")
  RedList_matrix[ which(grepl(key, RedList_matrix), arr.ind = T) ]    <- j
}

RedList_matrix[which( (!RedList_matrix %in% adjust_categories$Cat_IUCN)  & !is.na(RedList_matrix) , arr.ind = TRUE )]<-NA




## Ajustar matriz

### Eliminar las especies que no aportan a la estimación. Menos de dos años con dato de categorización  ####
RedList_matrix_2<- RedList_matrix[rowSums(!is.na(RedList_matrix))>=2,]

### Asignar los años no evaluados NA con datos del año de ultima evaluacion ####
replace_na_with_previous <- function(df, target_col) {
  for (col in (target_col-1):2) {
    df[[target_col]] <- ifelse(is.na(df[[target_col]]), df[[col]], df[[target_col]])
  }; return(df) }

df = RedList_matrix_2 %>% as.data.frame.matrix()
for(k in 2:ncol(RedList_matrix_2)){ df <- replace_na_with_previous(df, k) }


# Resultados de la libreria red ####
red_resuLt<- red::rli(df) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("Year") %>%  setNames(c("Year", "RLI_red")) %>% 
  dplyr::filter(!Year %in% "Change/year")

red_resuLt_plot<- ggplot(red_resuLt, aes(x = Year, y = RLI_red)) +
  geom_line(group = 1, col= "red") +
  geom_point() +
  coord_cartesian(ylim = c(0,1))+
  theme_classic()+
  theme(    panel.grid.major = element_line(color = "gray"),
  ) + theme(text = element_text(size = 4)) + ggtitle("red_results")


# Resultados - ejecución manual (Butchart et al, 2004; 2007) ####
## Definir los pesos para cada categoría ####
category_weights <- c("LC" = 0, "NT" = 1, "VU" = 2, "EN" = 3, "CR" = 4, "DD" = NA)  # DD es NA porque no se incluye

### Reemplazar las categorías por sus pesos correspondientes ####
df_weights <- apply(df, 2, function(column) category_weights[column]) %>% as.data.frame.matrix() %>% 
  dplyr::mutate(ID_Fila= rownames(df)) %>% tibble::column_to_rownames(var = "ID_Fila")

### Ecuacion RLI Butchardt ####
#### Definir las variables de las ecuaciones ####
vars_period<- sapply(names(df_weights), function(y) {
  N= sum(!is.na(df_weights[,y])) # Numero de especies con evaluación
  W= 5 # peso maximo entre las especies evaluadas 5 para EX
  M = W * N # Maximum threath score
  Wc_ts<- sum(df_weights[,y], na.rm=T) # Suma del peso de la categoría de amenaza c para la especie s en el tiempo t
data.frame(N=N, W=W, M= M, Wc_ts=Wc_ts)
}) %>% as.data.frame()



#### Ejecutar ecuacion ####
RLI_Butchart_2007<-pblapply(vars_period, function(x) { (x$M-x$Wc_ts)/x$M }) %>% unlist()
RLI_Butchart_2007

Butchart_result<- data.frame(Year= names(RLI_Butchart_2007), RLI_Butchart = (RLI_Butchart_2007))
Butchart_resuLt_plot<-  ggplot(Butchart_result, aes(x = Year, y = RLI_Butchart)) +
  geom_line(group = 1, col= "blue") +
  geom_point() +
  coord_cartesian(ylim = c(0,1))+
  theme_classic()+
  theme(    panel.grid.major = element_line(color = "gray"),
  ) + theme(text = element_text(size = 4)) + ggtitle("Butchart_result")

# Comparar resultados ####
compare_data<- list(red_resuLt, Butchart_result) %>% plyr::join_all()

compare_plot<- ggpubr::ggarrange(plotlist= list(red_resuLt_plot, Butchart_resuLt_plot), ncol= 2)


# export results
openxlsx::write.xlsx(IUCN_sp_Taxon, file.path(output, paste0("IUCN_sp_Taxon", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(historyAssesment_data, file.path(output, paste0("historyAssesment_data", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(adjust_categories, file.path(output, paste0("adjust_categories", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(RedList_matrix, file.path(output, paste0("RedList_matrix", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(df, file.path(output, paste0("df", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(red_resuLt, file.path(output, paste0("red_resuLt", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(df_weights, file.path(output, paste0("df_weights", ".xlsx")), rowNames = TRUE ) 
openxlsx::write.xlsx(vars_period, file.path(output, paste0("vars_period", ".xlsx")), rowNames = TRUE )


openxlsx::write.xlsx(Butchart_result, file.path(output, paste0("Butchart_result", ".xlsx")), rowNames = TRUE )
openxlsx::write.xlsx(compare_data, file.path(output, paste0("compare_data", ".xlsx")), rowNames = TRUE )

ggsave(file.path(dirname(output), "README_figures", paste0("red_resuLt_plot", ".jpg")), red_resuLt_plot, height = 2, width = 4)
ggsave(file.path(dirname(output), "README_figures", paste0("Butchart_resuLt_plot", ".jpg")), Butchart_resuLt_plot, height = 2, width = 4)
ggsave(file.path(dirname(output), "README_figures", paste0("compare_plot", ".jpg")), compare_plot, height = 2, width = 4)




