# Modelización de datos de geolocalización de especies arbóreas mediante 
# redes complejas. Aplicación al caso andaluz.

# SCRIPT PARA LA CREACIÓN DE ARCHIVOS RELACIONADOS CON GRAFOS ESPACIALES 
# QUE PUEDA TOMAR QGIS PARA SU REPRESENTACIÓN GEOGRÁFICA. 

# En este proyecto, los grafos espaciales corresponden a grafos de individuos
# creados con el script "Grafos_en_funcion_de_r.r" .

rm(list=ls())
options(scipen=999)
options(digits=22)
memory.limit(size=10000)

# Instalación y carga de paquetes.

list.of.packages <- c("Imap","ggplot2","dplyr","sqldf", "igraph", "ggraph", "shinyBS", "shinydashboard",
                      "sf", "rnaturalearth", "rnaturalearthdata", "tidygraph",
                      "ggpubr","data.table","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)



# 1. Generación de un archivo que consista en la lista de adyacencia
# de un grafo espacial pero donde para cada nodo aparezcan todos
# atributos del registro (identificador, coordenadas y especie, aunque para 
# la representación en QGIS solo se necesitan las coordenadas).

# El procedimiento para obtener dicho objeto de adyacencia es el mismo
# que en el script "Analisi_tp.r" (allí luego se filtran las especies).


# Se lee el registro de árboles de la población, para obtener los nodos
# del grafo.
registro <- fread("Registro_nevada.csv",header=T,col.names = c("Nodo","X","Y","Especie"))


# Se lee la lista de adyacencia del grafo espacial. 
adj_arb <- fread("adj_0.06.csv", header = FALSE, sep = ",",dec = ".", col.names = c("nodo1","nodo2"))
adj_arb

# Esta función es llamada por la función "añadir_atributos_a_lista_adj"
# y se aplicará, con la función apply, a cada tupla de la lista de adyacencia
# de individuos para obtener de cada uno de los dos nodos que aparecen en la
# tupla, todos sus atributos del registro (es decir, con cada identificador 
# de nodo consulta el registro y devuelve todos los campos).
para_una_tupla <- function(fila){
  n1 <- fila["nodo1"]
  n2 <- fila["nodo2"]
  resultado_consulta <- cbind(registro[registro$Nodo == n1],registro[registro$Nodo == n2])
  return(resultado_consulta)
}

# Esta función toma como entrada el  registro de árboles de la población 
# sobre la que se ha creado el grafo espacial, la lista de adyacencia de 
# dicho grafo y el valor de distancia umbral sobre el que se creó el 
# grafo, para añadirlo al nombre del archivo de adyacencia que crea
# # (se puede introducir en este parámetro cualquier cosa que se quiera
# añadir al nombre). El archivo lo crea en el mismo directorio desde donde 
# se ejecuta la función.
añadir_atributos_a_lista_adj <- function(registro,adj_arb,tresh){
  n_adj_arb <- nrow(adj_arb)
  # Se lee la lista de adyacencia de árboles por intervalos de 100 tuplas para 
  # aplicar la función "para una tupla" a cada intervalo. Los resultados 
  # de cada intervalo se van volcando en un archivo, para evitar crear un objeto 
  # demasiado grande en R.
  nombre_archivo <- paste0("adj_complet",tresh,".csv")
  cortes <- seq(from=1, to=n_adj_arb, by=1000)
  length(cortes)
  cortes[length(cortes)]
  cat("Creando archivo con",n_adj_arb,"tuplas\n")
  a <- 1
  for(un_corte in cortes){
    b <-un_corte
    cat("...tupla",b,"\n")
    interv_tuplas <- adj_arb[a:b,]
    Rconsulta <- apply(interv_tuplas[,1:2],MARGIN=1,FUN=para_una_tupla)
    Rconsulta_dt<- data.table()
    f <- length(Rconsulta)
    for (i in 1:f){
      una_tupla <- as.data.table(Rconsulta[i])
      Rconsulta_dt <-rbind(Rconsulta_dt,una_tupla)
    }
    if (a==1){
      file.create(nombre_archivo)
      colnames(Rconsulta_dt) <- c("nodo1", "X1", "Y1","Especie1", "nodo2", "X2", "Y2", "Especie2")
      fwrite(Rconsulta_dt, file = nombre_archivo, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    }
    else{
      fwrite(Rconsulta_dt, file = nombre_archivo, sep = ",", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE, showProgress = getOption("datatable.showProgress", interactive()))
    }
    a <- b+1
  }
  
  interv_final <- adj_arb[a:n_adj_arb,]
  Rconsulta_final <- apply(interv_final[,1:2],MARGIN=1,FUN=para_una_tupla)
  Rconsulta_dt<- data.table()
  f <- length(Rconsulta_final)
  for (i in 1:f){
    una_tupla <- as.data.table(Rconsulta_final[i])
    Rconsulta_dt <-rbind(Rconsulta_dt,una_tupla)
  }
  fwrite(Rconsulta_dt, file = nombre_archivo, sep = ",", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE, showProgress = getOption("datatable.showProgress", interactive()))
  
  cat("Archivo",nombre_archivo, "creado\n ")
  # COMPROBACIÓN DE QUE EL ARCHIVO DE ADYACENCIA DE ENTRADA TIENE EL MISMO
  # NÚMERO DE FILAS QUE EL DE SALIDA.
  
  prueba <- fread(nombre_archivo, header = T, sep = ",",dec = ".")
  prueba
  n_prueba <- nrow(prueba)
  if (n_adj_arb == n_prueba){
    cat("El número de filas de ambos objetos de adyacencia coinciden, todo ha ido bien.")
  } else {"Algo no ha funcionado, no coincide el número de filas de ambos objetos"}
  ####

}

# Ejecutamos la función.
añadir_atributos_a_lista_adj(registro,adj_arb,0.06)


# 2. Cálculo del camino mínimo entre dos nodos dados de un grafo
# y exportación de la ruta a un archivo representable por QGIS.


# Para la creación y exportación del archivo, se define la misma función 
# que aparece en el script"Grafos_en_funcion_de_r.r":

# FUNCIÓN PARA QUE EXPORTA UNA SECUENCIA DE NODOS FORMANDO UN CAMINO
# A UN ARCHIVO REPRESENTABLE EN QGIS.
# Esta función toma una secuencia (lista) de identificadores de nodos de 
# árboles que aparecen formando un camino y crea un archivo csv en el que 
# cada nodo aparece con sus coordenadas (tomadas del registro) y con un 
# número natural identificador del orden de aparición en la secuencia. 
create_path_doc <- function(path_sec, registro_c, nombre_archivo_path){
  n_sec <- length(path_sec)
  recop <- data.frame()
  orden <- c(1:n_sec)
  for (i in orden){
    nod <- path_sec[i]
    recop <- rbind(recop,registro_c[registro_c$Nodo == nod,])
  }
  recop$ruta_d <- orden
  write.table(recop, file = nombre_archivo_path, sep=",",row.names=F, col.names=T)
}

# Para el cálculo del camino (y creación y exportación del archivo a través 
# de una llamada a la función anterior) se crea la siguiente función, que
# toma como entrada un grafo (objeto igraph), los dos identificadores de 
# nodos cuyo camino mínimo se desea calcular y el valor de distancia umbral 
# sobre el que se creó el grafo, para añadirlo al nombre del archivo
# (se puede introducir en este parámetro cualquier cosa que se quiera
# añadir al nombre). El archivo lo crea en el mismo directorio desde donde 
# se ejecuta la función.
calc_export_shortest_path <- function(grafo_igraph,nodo1,nodo2,tresh){
  shp <- shortest_paths(graf, nodo1,nodo2)
  secuencia <- as_ids(shp$vpath[[1]])
  nombre_archivo <- paste0("sh_",nodo1,"_",nodo2,"_",tresh,".csv")
  create_path_doc(secuencia,registro,nombre_archivo)
}

# Para la creación del grafo se necesitan los nodos
nodos_a_id <- as.data.frame(unname(as.integer(registro$Nodo)))

# Y la lista de adyacencia, ya cargada de antes
graf <- graph_from_data_frame(adj_arb, directed=FALSE, vertices=nodos_a_id)


# Ejecutamos la función
calc_export_shortest_path(graf,27779,39007,0.06)
