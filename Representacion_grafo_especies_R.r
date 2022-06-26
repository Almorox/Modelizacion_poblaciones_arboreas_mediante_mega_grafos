### CÓDIGO PARA LA REPRESENTACIÓN GRÁFICA DE UN GRAFO DE ESPECIES.
# Este código toma un grafo de especies creado sobre una población de
# árboles dada, bajo cierto valor de distancia umbral y muestra información
# sobre él (qué especies son más importantes en base a su centralidad, etc.).
# También aplica el algoritmo de Louvain sobre el grafo (búsqueda de 
# comunidades) y representa el grafo por colores indicando distintas comunidades.
# Se exporta la imagen en formato tiff.


rm(list=ls())
options(scipen=999)
options(digits=15)
memory.limit(size=10000)


list.of.packages <- c("dplyr","kableExtra","dplyr","Imap", "igraph", "ggraph","data.table","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# Se crea el grafo a partir de una lista de adyacencia (de especies).
adj <- fread("adj_0.04.csv", header = FALSE, sep = ",",dec = ".", col.names = c("nodo1","nodo2"))

# Para ello debemos obtener también los nodos (especies de la población).
registro <- read_delim("Registro_nevada.csv")
species_list <- distinct(registro,Especie)
species_list
nod <- as.data.frame(species_list)

# Se crea el grafo
graf <- graph_from_data_frame(adj, directed=FALSE, vertices=nod)

# Se obtiene información sobre el grafo.
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

# Se muestran las especies con más de 10 vecinos.
for(i in especies_listadas){
  inter <- as_ids(neighbors(graf, i))
  if (length(inter)>10){
    cat(i,"\n")} 
}

# Se muestra la especie con un mayor número de vecinos.
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
  cat(as_ids(V(graf)[degree(graf)==grado]),sep=",")
  cat("\n")
}



# Búsqueda de comunidades con el algoritmo de Louvain.

# Se obtiene un objeto tipo "communities" con el comando de cluster_louvain
# (que aplica el algoritmo de Louvain para buscar una partición del grafo
# que optimice el valor de modularidad).
louv_mod <- cluster_louvain(graph = graf)
ship <- louv_mod$membership # Pertenencia de cada nodo a una comunidad
ship
# Número de nodos
length(ship) 

graf$community <- louv_mod$membership
ncomun <- length(unique(graf$community))
# Número de comunidades que ha obtenido el algoritmo de Louvain
ncomun 


# Para la representación gráfica se asigna a cada especie un número
# (para representar los números y no los nombres, por ser estos últimos
# en general bastante largos).

# Número de nodos
num_esp <- gorder(graf) 
nums <- c(1:num_esp)
# Copia del grafo original
graf_original <- graf
# Grafo sobre el que vamos a cambiar los nombres de los nodos por números
graf_num <- graf

for (i in 1:length(V(graf_num))) { 
  cat("\n",i,"-",V(graf_num)$name[i])
  V(graf_num)$name[i] <- i
} 

graf <- graf_num


# Representación gráfica

# Tamaño de nodos
V(graf)$size <- 4
# Color del fondo
V(graf)$frame.color <- "white" 
# Color de cada nodo en base a su comunidad
V(graf)$color <- graf$community 
# Se convierten los nombres del grafo en las etiquetas a representar en la
# imagen
V(graf)$label <- V(graf)$name 
# Tamaño de la fuente de las etiquetas
V(graf)$label.cex <- 3.2
# Tipo de fuente de las etiquetas (bold)
V(graf)$label.font <- 2
# Color de las etiquetas
V(graf)$label.color <- "black"

# Si no quisiéramos etiquetas:
#for (i in 1:length(V(graf))) { #esto dependiendo de si quieres labels o no
#V(graf)$label[i] <- "" } 

# Se colorean los enlaces en base al nodo saliente (como no es un grafo 
# dirigido, en nuestro caso se colorea en base a cualquiera de los dos
# nodos). 
edge.start <- ends(graf, es = E(graf), names = F)[,1] 
E(graf)$color <- V(graf)$color[edge.start] 
E(graf)$arrow.mode <- 0 
V(graf)$size <- 5
# Layout en círculo (nodos ubicados en un círculo, no dentro de este)
l_circle <- layout_in_circle(graf, order = V(graf))
# Layout fr
l_fr <- layout_with_fr(graf)
# (Se pueden elegir muchos más layouts)
# En la consola de R no suele aparecer bien el grafo, es mejor exportarlo.
plot(graf, rescale = T, layout = l_fr)

# Se exporta la imagen en formato tiff:
current_dir <- getwd()
ruta_imagen <- paste0(current_dir,"/grafo_representacion.tiff")
tiff(ruta_imagen, units="in", width=30, height=35, res=300)
# Edge.width, otro parámetro útil que cambia el grosor de los enlaces.
plot(graf, edge.width=3,rescale = T, layout = l_fr)
dev.off()


# Creación de tablas que dan información sobre los individuos más importantes
# en base a su centralidad (betweeness).

# Se vuelve a trabajar con el grafo original.
graf <- graf_original

# Creación de una tabla que muestra el número de nodos que pertenecen a la 
# comunidad y los cinco nodos más importantes (con mayor centralidad) dentro de 
# ella.
communities <- data.frame()
for (i in unique(graf$community)) { 
  # Se crea un subgrafo para cada comunidad 
  subgraph <- induced_subgraph(graf, v = which(graf$community == i)) 
  # Se guarda el tamaño de cada subgrafo 
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
colnames(communities) <- c("Comunidad","Nº nodos","Más importante")
com_mostt <- knitr::kable(
  communities %>% 
    dplyr::select("Comunidad", "Nº nodos","Más importante")
)


kbl(communities, booktabs = T)
kbl(communities, booktabs = T) %>%
  kable_styling(latex_options = "striped")


# Creación de una tabla que muestra para cada comunidad, los cinco
# individuos más importantes (en base a su centralidad) si los tiene. 
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
    result <- data.frame(community = i, rank = c("1º","2º","3º","4º","5º"), Species = top) 
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

