# Modelizaci�n de datos de geolocalizaci�n de especies arb�reas mediante 
# redes complejas. Aplicaci�n al caso andaluz.

# SCRIPT PARA LA CREACI�N Y AN�LISIS DE GRAFOS EN FUNCI�N DEL VALOR
# DE Tp.

# A partir de un registro de �rboles (con las especies anotadas) de una 
# poblaci�n dada y de la lista de adyacencia de un grafo espacial de �rboles
# (constru�do con cierto valor de distancia umbral para determinar las 
# relaciones espaciales), este c�digo crea los grafos de especies 
# correspondientes a esa misma poblaci�n y a esa misma distancia umbral 
# para distintos valores de Tp, los cuales toma como entrada en un vector.
# (Tp es el m�nimo n�mero de pares de individuos que est�n a una distancia igual o
# menor que la distancia umbral que han de existir entre dos especies para que 
# estas se consideren adyacentes en el grafo de especies). Para cada uno de los 
# grafos as� creados, calcula diversas propiedades topol�gicas. Por cada grafo
# devuelve una carpeta con sus resultados para todas las propiedades, adem�s
# de una matriz de adyacencia para crear el grafo otra vez en cualquier momento.
# Por cada propiedad, devuelve un archivo de los resultados de todos los grafos,
# indicando su valor de Tp para estudiar las propiedades en funci�n de este 
# valor.


rm(list=ls())
options(scipen=999)
options(digits=22)
memory.limit(size=10000)

# Instalaci�n y carga de paquetes.

list.of.packages <- c("Imap","dplyr","sqldf", "igraph","sf","data.table","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)


# Se lee la lista de adyacencia de individuos para una poblaci�n y un valor
# de distancia umbral dado. 
# En este ejemplo, el archivo que se lee corresponde a la lista del grafo
# de individuos de Sierra Nevada con una distancia umbral de 0.34 km.
adj_arb <- fread("adj_0.34.csv", header = FALSE, sep = ",",dec = ".", col.names = c("nodo1","nodo2"))

# N�mero de enlaces entre �rboles (no todos implicar�n un enlace entre especies,
# solo aquellos en los que los dos �rboles son de especies diferentes).
n_adj_arb <- nrow(adj_arb)
n_adj_arb

# Se lee el registro de �rboles de la poblaci�n deseada.
registro <- fread("Registro_nevada.csv",header=T,col.names = c("Nodo","X","Y","Especie"))
# N�mero de �rboles en la poblaci�n
nrow(registro)

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

# Esta funci�n se aplicar� con la funci�n apply a cada tupla de la 
# lista de adyacencia de individuos para obtener de cada uno de los dos 
# nodos que aparecen en la tupla, todos sus atributos del registro 
# (es decir, con cada identificador de nodo consulta el registro y devuelve
# todos los campos).
para_una_tupla <- function(fila){
  n1 <- fila["nodo1"]
  n2 <- fila["nodo2"]
  resultado_consulta <- cbind(registro[registro$Nodo == n1],registro[registro$Nodo == n2])
  return(resultado_consulta)
}

# Se lee la lista de adyacencia de �rboles por intervalos de 100 tuplas para 
# aplicar la funci�n "para una tupla" a cada intervalo. Los resultados 
# de cada intervalo se van volcando en un archivo, para evitar crear un objeto 
# demasiado grande en R.
cortes <- seq(from=1, to=n_adj_arb, by=1000)
length(cortes)
cortes[length(cortes)]

a <- 1
for(un_corte in cortes){
  b <-un_corte
  interv_tuplas <- adj_arb[a:b,]
  Rconsulta <- apply(interv_tuplas[,1:2],MARGIN=1,FUN=para_una_tupla)
  Rconsulta_dt<- data.table()
  f <- length(Rconsulta)
  for (i in 1:f){
    una_tupla <- as.data.table(Rconsulta[i])
    Rconsulta_dt <-rbind(Rconsulta_dt,una_tupla)
  }
  if (a==1){
    file.create("adj_complet.csv")
    colnames(Rconsulta_dt) <- c("nodo1", "X1", "Y1","Especie1", "nodo2", "X2", "Y2", "Especie2")
    fwrite(Rconsulta_dt, file = "adj_complet.csv", sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
  }
  else{
    fwrite(Rconsulta_dt, file = "adj_complet.csv", sep = ",", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE, showProgress = getOption("datatable.showProgress", interactive()))
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
fwrite(Rconsulta_dt, file = "adj_complet.csv", sep = ",", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE, showProgress = getOption("datatable.showProgress", interactive()))

# COMPROBACI�N DE QUE EL ARCHIVO DE ADYACENCIA DE ENTRADA TIENE EL MISMO
# N�MERO DE FILAS QUE EL DE SALIDA.

prueba <- fread("adj_complet.csv", header = T, sep = ",",dec = ".")
prueba
n_prueba <- nrow(prueba)
if (n_adj_arb == n_prueba){
  cat("El n�mero de filas de ambos objetos de adyacencia coinciden, todo ha ido bien.")
}
####

# Se obtiene del objeto de adyacencia creado un objeto de adyacencia 
# equivalente pero �nicamente con el campo de especie para cada nodo 
# (L_adj_ind_prima).
# (el objeto con todos los campos queda guardado por si se quisieran representar
# los enlaces espaciales en QGIS, para lo que se necesitan las coordenadas).
colnames(prueba) <- c("nodo1", "X1", "Y1","Species1", "nodo2", "X2", "Y2","Species2")
L_adj_ind_prima <- select(prueba,c('Species1','Species2'))

# Se guarda este objeto para poder leerlo en el futuro.
file.create("L_adj_ind_prima.csv")
write.table(L_adj_ind_prima, file = "L_adj_ind_prima.csv", sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)

###### Para leer la lista en el futuro:
# L_adj_ind_prima <- fread("L_adj_ind_prima.csv", header = T, sep = ",",dec = ".", col.names = c("Species1","Species2"))
# L_adj_ind_prima
######

# Se convierte el objeto con los nombres de las especies que aparecen en el
# registro en un vector.
class(esp_order)
esp_order_v <- as.vector(esp_order)
length(esp_order_v)

# El vector indica los nombres de las filas y columnas de la matriz de
# frecuencias absolutas. Se crea tal matriz con ceros.
Mat_F <- matrix(nrow=length(esp_order_v),ncol=length(esp_order_v))
colnames(Mat_F) <- esp_order_v
row.names(Mat_F) <- esp_order_v
Mat_F[is.na(Mat_F)] <- 0
#Mat_F

# Se rellena la matriz de frecuencias absolutas (sim�trica y con digonal 
# principal compuesta por ceros) construida en base a la lista L_adj_ind_prima de
# adyacencia entre especies (donde cada tupla es un enlace entre individuos 
# indicando su especie).

# No se guardan enlaces entre una especie y ella misma en la matriz de 
# frecuencias. 

for (fila in 1:nrow(L_adj_ind_prima)){
  especie1 <- as.character(unname(L_adj_ind_prima[fila,"Species1"]))
  especie1 <- especie1[1]
  especie2 <- as.character(unname(L_adj_ind_prima[fila,"Species2"]))
  especie2 <- especie2[1]
  if (especie1 != especie2){
    Mat_F[especie1,especie2]<-Mat_F[especie1,especie2]+1
    Mat_F[especie2,especie1]<-Mat_F[especie2,especie1]+1
  }
}

Mat_F

# Se guarda la matriz de frecuencias absolutas para poder realizar otro
# an�lisis de Tp con ella en cualquier momento. 
# Es el objeto (para una poblaci�n y una distancia umbral dada) sobre el que
# hay aplicar la condici�n del valor de Tp para crear distintos grafos
# y estudiar el comportamiento de cada propiedad en base al valor de Tp.
file.create("Mat_F.csv")
write.table(Mat_F, file = "Mat_F.csv", sep = ",", dec = ".", row.names = TRUE, col.names = NA)

##### Para leer la matriz F en el futuro
# Mat_F <- as.matrix(fread("Mat_F.csv",sep = ",",dec = "."),rownames=1)
# Mat_F
#####

# Funci�n que dada un registro, una matriz de frecuencias y un vector de valores
# Tp, crea un grafo para cada valor de Tp a partir de la matriz de frecuencias
# y calcula las propiedades de todos ellos. Devuelve los resultados en un
# directorio cuyo nombre toma tambi�n como input.
MatF_to_MatAdj <- function(registro,Mat_F,array_tp,nombre_proyect){
  dir.create(paste0("Proyecto_tp_",nombre_proyect))
  setwd(paste0("Proyecto_tp_",nombre_proyect))
  dir_out <- getwd()
  
  # Se crean los archivos de los par�metros calculados. Cada uno recoge los 
  # valores de un par�metro calculados sobre todos los grafos.
  file.create("coef_clustering.csv")
  file.create("coef_clustering_zeros.csv")
  file.create("densidad.csv")
  file.create("LCC.csv")
  file.create("modularidad.csv")
  file.create("num_comunidades.csv")
  file.create("num_conexiones.csv")
  file.create("pc_nodos_no_aislados.csv")
  file.create("grado_medio.csv")
  
  # Se escriben los cabeceros de todos ellos.
  
  cc <- as.data.table(cbind("Tp","coef.clustering"))
  fwrite(cc, file = "coef_clustering.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  cc2 <- as.data.table(cbind("Tp","coef.clustering.zeros"))
  fwrite(cc2, file = "coef_clustering_zeros.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  dens <- as.data.table(cbind("Tp","densidad conexiones"))
  fwrite(dens, file = "densidad.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  LCC <- as.data.table(cbind("Tp", "LCC"))
  fwrite(LCC, file = "LCC.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  modul <- as.data.table(cbind("Tp", "modularidad"))
  fwrite(modul, file = "modularidad.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  comps <- as.data.table(cbind("Tp", "n� componentes"))
  fwrite(comps, file = "num_componentes.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  comuns <- as.data.table(cbind("Tp", "n� comunidades"))
  fwrite(comuns, file = "num_comunidades.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  conexs <- as.data.table(cbind("Tp", "n� conexiones"))
  fwrite(conexs, file = "num_conexiones.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  aisl <- as.data.table(cbind("Tp", "nodos conectados"))
  fwrite(aisl, file = "pc_nodos_no_aislados.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  gm <- as.data.table(cbind("rTp", "grado medio"))
  fwrite(gm, file = "grado_medio.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  
  #lista_especies <- distinct(registro,Especie)
  #n_especies <- nrow(lista_especies)
  
  for (tp in array_tp){
    dir_tp <- paste0("Tp",tp)
    dir.create(dir_tp)
    setwd(dir_tp)
    dir_tp <- getwd()
    n <- nrow(Mat_F)
    adj_matrix <- matrix(0, n, n)
    
    #Determinaci�n de la adyacencia entre dos especies.
    for (i in 1:n){
      for (j in 1:n){
        if (Mat_F[i,j] >= tp){
          adj_matrix[i,j] = 1
        }
      }
    }
    colnames(adj_matrix) <- colnames(Mat_F)
    row.names(adj_matrix) <- row.names(Mat_F)
    nomb_adj <- paste0("mat_adj_",tp,".csv")
    # Se guarda la matriz de adyacencia para poder crear el grafo en cualquier
    # momento.
    file.create(nomb_adj)
    write.table(as.matrix(Mat_F), file = nomb_adj, sep = ",", dec = ".", row.names = TRUE, col.names = NA)
 
  # C�LCULO DE PROPIEDADES Y PAR�METROS DEL GRAFO DE ESPECIES.
  # Cada propiedad calculada se escribe sobre un archivo com�n a todos los 
  # grafos de especies (creados bajo distinto valor de Tp) de la poblaci�n.
  cat("...EMPIEZA LA CREACI�N DEL GRAFO DE ESPECIES...")
 
  if (nrow(adj_matrix) != 0){
    
    #GRAFOS
    tramit <- as.matrix(adj_matrix)
    #nodos_s_id <- as.data.frame(lista_especies)
    cat("Creando grafo...")
    graf <- graph_from_adjacency_matrix(tramit, mode= "undirected", diag=F)
    cat("DONE...")
    
    cat("Calculando num v�rtices y ejes...")
    num_ver <- gorder(graf)
    num_edges <- gsize(graf)
    tresh_num_edges <- as.data.table(cbind(tp, num_edges))
    d_edges <- paste0(dir_out,"/","num_conexiones.csv")
    fwrite(tresh_num_edges, file = d_edges, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando modularidad...")
    louv_mod <- cluster_louvain(graph = graf)
    mod <- round(modularity(graf, membership = louv_mod$membership),3)
    tresh_mod <- as.data.table(cbind(tp, mod))
    d_mod <- paste0(dir_out,"/","modularidad.csv")
    fwrite(tresh_mod, file = d_mod, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("calculando coef clustering...")
    cc <- round(transitivity(graf, type = "average"),3)
    tresh_cc <- as.data.table(cbind(tp,cc))
    d_cc <- paste0(dir_out,"/","coef_clustering.csv")
    fwrite(tresh_cc, file = d_cc, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("calculando coef clustering ZEROS...")
    cc_zeros <- round(transitivity(graf, type = "average", isolates = "zero"),3)
    tresh_cc_zeros <- as.data.table(cbind(tp,cc_zeros))
    d_cc <- paste0(dir_out,"/","coef_clustering_zeros.csv")
    fwrite(tresh_cc_zeros, file = d_cc, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando densidad...")
    dens <- edge_density(graf)
    tresh_dens <- as.data.table(cbind(tp,dens))
    d_dens <- paste0(dir_out,"/","densidad.csv")
    fwrite(tresh_dens, file = d_dens, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando num comunidades...")
    graf$community <- louv_mod$membership
    ncomun <- length(unique(graf$community))
    tresh_ncomun <- as.data.table(cbind(tp, ncomun))
    d_comuns <- paste0(dir_out,"/","num_comunidades.csv")
    fwrite(tresh_ncomun, file = d_comuns, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando LCC...")
    LCC <- length(component_distribution(graf,mode="weak"))-1 
    tresh_LCC <- as.data.table(cbind(tp, LCC))
    d_LCC <- paste0(dir_out,"/","LCC.csv")
    fwrite(tresh_LCC, file = d_LCC, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando num de nodos aislados...")
    pc_parti <- length(V(graf)[degree(graf)!=0])*100/num_ver
    tresh_pc_parti <- as.data.table(cbind(tp, pc_parti))
    d_parti <- paste0(dir_out,"/","pc_nodos_no_aislados.csv")
    fwrite(tresh_pc_parti, file = d_parti, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando N� componentes conexas...")
    ncompons <- count_components(graf)
    tresh_ncompons <- as.data.table(cbind(tp, ncompons))
    d_compons <- paste0(dir_out,"/","num_componentes.csv")
    fwrite(tresh_ncompons, file = d_compons, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("Calculando grado medio...")
    gm <- mean(degree(graf))
    tresh_gm <- as.data.table(cbind (tp,gm))
    d_gm <- paste0(dir_out,"/","grado_medio.csv")
    fwrite(tresh_gm, file = d_gm, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    cat("DONE...")
    
    cat("calculando distribuci�n de grados...")
    tabla_g <- table(degree(graf))
    grado <- as.integer(names(tabla_g))
    cuantos_nodos <- unname(tabla_g)
    g_distribucion <- as.data.table(cbind(grado,cuantos_nodos))
    colnames(g_distribucion) <- c("Grado","N� nodos")
    file.create("distribucion_grados.csv")
    fwrite(g_distribucion, file = "distribucion_grados.csv", sep=",",row.names=F, col.names=T)
    cat("DONE - TODOS LOS PAR�METROS CALCULADOS...")
    
    
    # Una vez calculados todos los par�metros del grafo, se procede a recoger
    # todos estos valores y volcarlos en un documento correspondiente al grafo
    # (hasta ahora solo se han volcado en documentos separados para cada 
    # par�metro, comunes a todos los grafos creados).
    cat("Escribiendo carpeta del tp",tp,"...\n")
  
    params <- c(num_ver,tp,num_edges,mod,ncomun,LCC,ncompons,cc,cc_zeros,dens,pc_parti,gm)
    matt <- matrix(params,nrow=12,ncol=1)
    colnames(matt) <- c("VALOR")
    row.names(matt) <- c("N� nodos","Radio de vecindad (km)","N� conexiones","Modularidad","N� comunidades","LCC","N� componentes","CC","CC_zeros","Densidad de enlaces","% nodos conectados","Grado medio")
    matt <- as.data.table(matt, keep.rownames = TRUE)
    nomb_arxi_prop <- paste0("propiedades",tp,".csv")
    fwrite(matt, file = nomb_arxi_prop, sep=",",row.names=T, col.names=T)
    
    
    cat("Tp",tp,"COMPLETADO")
    # Si no se ha podido formar el grafo de especies (por no haber ning�n enlace)
    # se avisa en pantalla.
  }else{cat("\nCero conexiones, no se forma grafo de especies para tp",tp,"\n")}
  
  
  #return(graf)
  setwd(dir_out)
  
  }
}


# Ejemplo de an�lisis llevado a cabo sobre esta poblaci�n (Sierra Nevada con
# distancia umbral de 0.34km).
array_tp1 <- c(1,35,50)
array_tp2 <- seq(100,30000,by=1200)
array_tp <- c(array_tp1,array_tp2,30000)
array_tp
MatF_to_MatAdj(registro,Mat_F,array_tp,"1-30000")

