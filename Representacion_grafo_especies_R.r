# Modelizaci�n de datos de geolocalizaci�n de especies arb�reas mediante 
# redes complejas. Aplicaci�n al caso andaluz.

# SCRIPT PARA OBTENER INFORMACI�N SOBRE LAS COMUNIDADES Y LOS NODOS M�S 
# IMPORTANTES DE UN GRAFO DE ESPECIES Y PARA SU REPRESENTACI�N GR�FICA.

# Toma como entrada una lista de adyacencia de un grafo de especies y el 
# el registro de la poblaci�n de �rboles correspondiente.
# Ofrece informaci�n sobre la proporci�n de individuos y el grado de cada
# especie, as� como el nombre de sus especies vecinas.
# Aplica el algoritmo de Louvain para detecci�n de comunidades.
# Muestra las especies m�s importantes en base a su centralidad en cada 
# comunidad.
# Representa el grafo completo con n�meros que identifican a las especies.
# Cada nodo aparece coloreado en base a la comunidad a la que pertenece. 


rm(list=ls())
options(scipen=999)
options(digits=15)
memory.limit(size=10000)


list.of.packages <- c("dplyr","kableExtra","dplyr","Imap", "igraph", "ggraph","data.table","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# Antes de crear el grafo, se obtiene informaci�n sobre el n�mero de individuos
# de cada especie de la poblaci�n. Para ello, se lee el registro.
registro <- read_delim("Registro_nevada.csv")

# Se obtienen las frecuencias absolutas de cada especie en la poblaci�n.
especies <- select(registro,Especie)
esp_table <- as.data.frame(table(especies))
colnames(esp_table) <- c("Especie","N�")
esp_table
# Se ordenan de forma decreciente
esp_frec_abs <- esp_table[order(-esp_table$N�),]
esp_frec_abs
esp_order <- esp_frec_abs$Especie
esp_order

# Se obtiene el porcentaje y la frecuencia relativa de cada una
porcent <- (esp_frec_abs$N�)*100/sum(esp_frec_abs$N�)
frec <- (esp_frec_abs$N�)/sum(esp_frec_abs$N�)

# Se crea y se guarda un data frame con esta informaci�n sobre las especies,
# obtenida �nicamente del registro.
info_especies<- cbind(esp_frec_abs,Porcentaje=porcent,Frecuencia=frec)
file.create("infoEsp_SNevada.csv")
write.table(format(info_especies, digits=22), file="infoEsp_SNevada.csv",sep=',',row.names=F, col.names = T)




# Se crea el grafo a partir de una lista de adyacencia (de especies).
adj <- fread("adj_0.04.csv", header = FALSE, sep = ",",dec = ".", col.names = c("nodo1","nodo2"))

# Para ello debemos obtener tambi�n los nodos (especies de la poblaci�n).
species_list <- distinct(registro,Especie)
species_list
nod <- as.data.frame(species_list)


# Se crea el grafo
graf <- graph_from_data_frame(adj, directed=FALSE, vertices=nod)

# Se obtiene informaci�n sobre el grafo.
# Nombres de los nodos:
especies_listadas <- nod$Especie
nod
# Para cada nodo, se muestran sus vecinos.
for(i in especies_listadas){
  inter <- as_ids(neighbors(graf, i))
  cat("----VECINOS DE ",i,": ")
  if (length(inter)==0){
    cat("-\n")} else{
      cat("\n    -")
      cat(inter,sep = "\n    -")
    }
}

# Se muestran las especies con m�s de 10 vecinos.
for(i in especies_listadas){
  inter <- as_ids(neighbors(graf, i))
  if (length(inter)>10){
    cat(i,"\n")} 
}

# Se muestra la especie con un mayor n�mero de vecinos.
V(graf)[degree(graf)==max(degree(graf))]

# Se obtienen las frecuencias absolutas de cada grado.
tabla_gm <- table(degree(graf))
grado <- as.integer(names(tabla_gm))
cuantos_nodos <- unname(tabla_gm)
gm_dist <- as.data.table(cbind(grado,cuantos_nodos))
gm_dist

# Se muestran los nombres de las especies que comparten cada valor de grado.
for (i in 1:nrow(gm_dist)){
  grado <- as.integer(unname(gm_dist[,1][i]))
  cat("\n-Especies con grado",grado,": ")
  cat(as_ids(V(graf)[degree(graf)==grado]),sep=", ")
  cat("\n")
}



# B�squeda de comunidades con el algoritmo de Louvain.

# Se obtiene un objeto tipo "communities" con el comando de cluster_louvain
# (que aplica el algoritmo de Louvain para buscar una partici�n del grafo
# que optimice el valor de modularidad).
louv_mod <- cluster_louvain(graph = graf)
ship <- louv_mod$membership # Pertenencia de cada nodo a una comunidad
ship
# N�mero de nodos
length(ship) 

graf$community <- louv_mod$membership
ncomun <- length(unique(graf$community))
# N�mero de comunidades que ha obtenido el algoritmo de Louvain
ncomun 


# Para la representaci�n gr�fica se asigna a cada especie un n�mero
# (para representar los n�meros y no los nombres, por ser estos �ltimos
# en general bastante largos).

# N�mero de nodos
num_esp <- gorder(graf) 
nums <- c(1:num_esp)
# Copia del grafo original
graf_original <- graf
# Grafo sobre el que vamos a cambiar los nombres de los nodos por n�meros
graf_num <- graf

for (i in 1:length(V(graf_num))) { 
  cat("\n",i,"-",V(graf_num)$name[i])
  V(graf_num)$name[i] <- i
} 

graf <- graf_num


# Representaci�n gr�fica

# Tama�o de nodos
V(graf)$size <- 4
# Color del fondo
V(graf)$frame.color <- "white" 
# Color de cada nodo en base a su comunidad
V(graf)$color <- graf$community 
# Se convierten los nombres del grafo en las etiquetas a representar en la
# imagen
V(graf)$label <- V(graf)$name 
# Tama�o de la fuente de las etiquetas
V(graf)$label.cex <- 3.2
# Tipo de fuente de las etiquetas (bold)
V(graf)$label.font <- 2
# Color de las etiquetas
V(graf)$label.color <- "black"

# Si no quisi�ramos etiquetas:
#for (i in 1:length(V(graf))) { # Descomentar para eliminar etiquetas.
#V(graf)$label[i] <- "" } 

# Se colorean los enlaces en base al nodo saliente (como no es un grafo 
# dirigido, en nuestro caso se colorea en base a cualquiera de los dos
# nodos). 
edge.start <- ends(graf, es = E(graf), names = F)[,1] 
E(graf)$color <- V(graf)$color[edge.start] 
E(graf)$arrow.mode <- 0 
V(graf)$size <- 5
# Layout en c�rculo (nodos ubicados en un c�rculo, no dentro de este)
l_circle <- layout_in_circle(graf, order = V(graf))
# Layout fr
l_fr <- layout_with_fr(graf)
# (Se pueden elegir muchos m�s layouts)
# En la consola de R no suele aparecer bien el grafo, es mejor exportarlo.
plot(graf, rescale = T, layout = l_fr)

# Se exporta la imagen en formato tiff:
current_dir <- getwd()
ruta_imagen <- paste0(current_dir,"/grafo_representacion.tiff")
tiff(ruta_imagen, units="in", width=30, height=35, res=300)
# Edge.width, otro par�metro �til que cambia el grosor de los enlaces.
plot(graf, edge.width=3,rescale = T, layout = l_fr)
dev.off()


# Creaci�n de tablas que dan informaci�n sobre los individuos m�s importantes
# en base a su centralidad (betweeness).

# Se vuelve a trabajar con el grafo original.
graf <- graf_original



# Centralidad de intermediaci�n de cada v�rtice (absoluta, no normalizada)
btw <- betweenness(graf)

btw_df <- data.frame(Especie = names(btw), Raw_Betweenness= unname(btw))

# Se ordena de forma decreciente el marco de datos.
btw_df_order <- btw_df[order(-btw_df$Raw_Betweenness),]
btw_df_order

# Se exporta a un archivo CSV.
file.create("Raw_Betweennesss.csv")
write.table(format(btw_df_order, digits=6), file="Raw_Betweennesss.csv",sep=',',row.names=F, col.names = T)

# Creaci�n de una tabla que muestra el n�mero de nodos que pertenecen a la 
# comunidad y los cinco nodos m�s importantes (con mayor centralidad) dentro de 
# ella.
communities <- data.frame()
for (i in unique(graf$community)) { 
  # Se crea un subgrafo para cada comunidad 
  subgraph <- induced_subgraph(graf, v = which(graf$community == i)) 
  # Se guarda el tama�o de cada subgrafo 
  size <- igraph::gorder(subgraph) 
  # Se calcula la centralidad de cada subgrafo 
  btwn <- igraph::betweenness(subgraph) 
  communities <- communities %>% 
    dplyr::bind_rows(data.frame(
      community = i, 
      n_nodos = size,
      most_important = names(which(btwn == max(btwn))) 
    ) 
    ) 
} 
colnames(communities) <- c("Comunidad","N� nodos","M�s importante")
com_mostt <- knitr::kable(
  communities %>% 
    dplyr::select("Comunidad", "N� nodos","M�s importante")
)


kbl(communities, booktabs = T)
kbl(communities, booktabs = T) %>%
  kable_styling(latex_options = "striped")


# Creaci�n de una tabla que muestra para cada comunidad, los cinco
# individuos m�s importantes (en base a su centralidad) si los tiene. 
top_five <- data.frame() 
for (i in unique(graf$community)) { 
  # Se crea un subgrafo para cada comunidad 
  subgraph <- induced_subgraph(graf, v = which(graf$community == i)) 
  # for larger communities 
  if (igraph::gorder(subgraph) > 4) { 
    # get degree 
    degree <- igraph::degree(subgraph) 
    # get top five degrees 
    top <- names(head(sort(degree, decreasing = TRUE), 5)) 
    result <- data.frame(community = i, rank = c("1�","2�","3�","4�","5�"), Species = top) 
  } else { 
    result <- data.frame(community = NULL, rank = NULL, Species = NULL) 
  } 
  top_five <- top_five %>% 
    dplyr::bind_rows(result) 
} 
colnames(top_five)[1] <- c("Comunidad")
kablee <- knitr::kable(
  top_five %>% 
    tidyr::pivot_wider(names_from = rank, values_from = Species) 
)
kablee
readr::write_file(kablee, "Most_imp_commun_and_species.html")
kablee %>%kable_styling(latex_options = "striped")

