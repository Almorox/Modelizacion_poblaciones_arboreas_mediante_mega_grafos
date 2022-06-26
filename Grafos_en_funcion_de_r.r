# Modelizaci�n de datos de geolocalizaci�n de especies arb�reas mediante 
# redes complejas. Aplicaci�n al caso andaluz.

# SCRIPT PARA LA CREACI�N Y EL AN�LISIS DE GRAFOS DE INDIVIDUOS Y DE 
# ESPECIES BAJO DISTINTOS VALORES DE DISTANCIA UMBRAL A A PARTIR DE UN 
# REGISTRO DE �RBOLES GEORREFERENCIADOS.

# El an�lisis completo se lanza con la funci�n "Ejecutar_analisis"
# la cual toma como entrada:

# - un registro de �rboles en el que tienen que aparecer cuatro campos
# (identificador de �rbol, longitud, altitud (coordenadas decimales) y 
# especie), como la salida que produce el script 
# "Procesamiento_de_registros".

# - un vector cuyos elementos corresponden a distintos valores de distancia
# umbral con los que se van a crear los grafos.

# - un valor booleano (T o F) para el par�metro "CALC_DIAM", indicando
# si se desea calcular el mayor di�metro de una componente conexa 
# (no es el di�metro del grafo, sino la mayor distancia geod�sica no 
# infinita del grafo). En caso de que este par�metro sea T, se incluye
# el c�lculo de esta distancia como una propiedad m�s (para la red
# de �rboles) y para cada grafo se devuelve un archivo con la ruta
# compuesta por la secuencia de nodos que aparecen en ella. Cada nodo,
# adem�s, aparece con sus coordenadas y su indicador de orden de aparici�n 
# en la ruta, lo cual permite su posterior representaci�n en QGIS.

# Para cada valor de distancia umbral se crear� un grafo espacial de 
# �rboles (un enlace entre dos nodos indica que los dos �rboles 
# se encuentran a una distancia igual o menor que la distancia umbral) 
# y un grafo no espacial de especies (un enlace entre dos especies
# indica que hay al menos un par de individuos que, perteneciendo cada
# uno a una de las dos especies, se encuentran a una distancia igual o 
# menor que la distancia umbral).

# La salida consiste en una carpeta correspondiente al proyecto,
# que contiene dos carpetas: "Arboles" y "Especies". En cada una de ellas
# hay una carpeta por cada valor de distancia umbral, que contiene 
# la lista de adyacencia del grafo correspondiente, un archivo "propiedades"
# con los c�lculos de las propiedades de dicho grafo y un archivo 
# "distribucion_grados" con las frecuencias absolutas de cada grado en el 
# grafo. En caso de que se haya calculado el di�metro, aqu� aparece 
# tambi�n el archivo con la ruta (en el caso de los grafos de �rboles).
# En cada una de las carpetas "Arboles" y "Especies" aparece tambi�n
# un archivo para cada propiedad con los resultados de todos los grafos 
# de individuos o de especies respectivamente (se indica el valor de
# distancia umbral del grafo y el valor de la propiedad). Esto permite
# estudiar el comportamiento de cada propiedad en funci�n de la distancia
# umbral, para ambos tipos de grafos (individuos y especies).

# Este script est� preparado para modelizar grandes poblaciones de �rboles.
# Para ello aplica un procedimiento de divisi�n del territorio en 
# en parcelas, que permite disminuir considerablemente el n�mero de 
# distancias por pares a calcular. 

# Por ahora se ha probado su eficacia con problaciones de hasta mill�n y 
# medio de individuos (con la opci�n CALC_DIAM = F). Seguramente funcione 
# con muchos m�s, pero habr�a que ajustar los p�rametros que controlan la 
# parcelaci�n (crear parcelas m�s peque�as que las que crean los valores 
# por defecto, pero sin perder la conectividad entre parcelas).

# Para que la funci�n "ejecutar_analisis" pueda lanzar el an�lisis
# han de estar definidas el resto de funciones que aparecen en este
# script.

# La parcelaci�n se lleva a cabo sobre el territorio definido por 
# cuatro valores de coordenadas decimales: MAX_X, MIN_X, MAX_Y, MIN_Y 
# (X son valores de longitud e Y de latitud). En este script, estos 
# valores definen el territorio andaluz, por lo que para la modelizaci�n
# de otro territorio han de cambiarse. 

rm(list=ls())
options(scipen=999)
options(digits=15)
memory.limit(size=10000)

# Instalaci�n y carga de paquetes necesarios.

list.of.packages <- c("Imap","dplyr","sqldf", "igraph","sf", "tidygraph",
                      "data.table","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# Lectura del archivo procesado.
#(archivo de GBIF procesado a trav�s del script "procesamiento_de_registros.r")
registro <- read_delim("Registro_nevada.csv")
# N�mero de instancias que componen el registro.
nrow(registro)

# El territorio a considerar es un rect�ngulo formado por los valores m�ximos y 
# m�nimos de longitud y latitud de las ocurrencias de la REDIAM para los g�neros
# de inter�s. El territorio as� definido sirve tanto para la modelizaci�n
# de registros de �rboles correspondientes al territorio andaluz completo o a 
# cualquier regi�n andaluza, en este caso: Sierra Morena andaluza y Sierra 
# Nevada (se lleva a cabo la modelizaci�n de cada poblaci�n de forma 
# independiente, es decir, con ejecuciones diferentes del c�digo).
MAX_X <- -1.637031
MIN_X <- -7.521846
MAX_Y <- 38.72789
MIN_Y <- 36.01393

# FUNCI�N 4 (en base al orden de llamada de funciones).
# Determinaci�n de enlaces en la red individuos y de especies en una
# parcela perteneciente al territorio considerado.
# La funci�n calcula y eval�a las distancias por pares entre �rboles 
# para crear las listas de adyacencia de �rboles y de especies de la parcela
# de forma simult�nea.
# Toma como entrada el registro correspondiente a una parcela y el valor de
# distancia umbral.

adj_list_P <- function(registro_parcela, tresh) {
  # Si en la parcela hay menos de dos �rboles, no se da ning�n enlace en ella, 
  # por lo que se para la ejecuci�n de la funci�n.
  if (nrow(registro_parcela)<2){
    return()
  }
  
  # N�mero de �rboles en la parcela.
  nInstancias <- nrow(registro_parcela)
  
  # Primera columna de la lista de adyacencia de �rboles de la parcela.
  adj_a_1 <- vector()
  # Segunda columna de la lista de adyacencia de �rboles de la parcela.
  adj_a_2 <- vector()
  # Primera columna de la lista de adyacencia de especies de la parcela.
  adj_s_1 <- vector()
  # Segunda columna de la lista de adyacencia de especies de la parcela.
  adj_s_2 <- vector()
  # Contador de filas en la lista de adyacencia de �rboles.
  v_a <- 1
  # Contador de filas en la lista de adyacencia de especies.
  v_s <- 1
  
  for (i in 1:(nInstancias-1)){
    # De cada instancia de la parcela, se guarda:
    # La especie.
    Esp_i <- registro_parcela$Especie[i]
    # La longitud
    X_i <- registro_parcela$X[i]
    # La latitud
    Y_i <- registro_parcela$Y[i]
  
    # Contador que muestra en pantalla el n�mero de instancia i cada 1000 
    # instancias recorridas.
    #if(i %% 1000 == 0){
      #cat(i)
      #cat("..")
    #}
    
    # Para cada instancia i se recorren las instancias j por delante de ella en
    # el registro (es equivalente a recorrer la segunda mitad de la matriz de 
    # adyacencia cortada por la diagonal principal y sin incluir esta).
    for (j in (i+1):nInstancias){ 
      # De cada instancia j se guarda:
      # La especie.
      Esp_j <- registro_parcela$Especie[j]
      # La longitud
      X_j <- registro_parcela$X[j]
      # La latitud
      Y_j <- registro_parcela$Y[j]
      
      # Se calcula la distancia entre el �rbol i y el �rbol j
      dist_ij <- 6371*acos(cos((90-Y_i)*pi/180)*cos((90-Y_j)*pi/180)+sin((90-Y_i)*pi/180)*sin((90-Y_j)*pi/180)*cos((X_i-X_j)*pi/180))
      
      # Si no se ha podido calcular la distanica, se cont�nua con la siguiente 
      # instancia j
      if (is.null(dist_ij) || is.na(dist_ij)){
        next
      }
      
      # Si la distancia es menor o igual que la distancia umbral, tenemos un 
      # enlace entre los dos �rboles en la red de individuos y un enlace
      # entre las dos especies en la red de especies (siempre que las dos 
      # especies no sean la misma). 
      if (dist_ij <= tresh){
        # LISTA DE ADYACENCIA DE �RBOLES.
        # Para luego poder filtrar ocurrencias �nicas en la lista de la red
        # completa, se guarda el nodo con n�mero identificador m�s bajo 
        # en la primera columna.
        if (registro_parcela$Nodo[i] < registro_parcela$Nodo[j]){
          adj_a_1[v_a] = registro_parcela$Nodo[i]
          adj_a_2[v_a] = registro_parcela$Nodo[j]
        } else {
          adj_a_1[v_a] = registro_parcela$Nodo[j]
          adj_a_1[v_a] = registro_parcela$Nodo[i]
        }
        # Se actualiza el contador de enlaces entre �rboles de la parcela
        v_a <- v_a + 1
        
        #LISTA DE ADYACENCIA DE ESPECIES.
        # Para luego poder filtrar ocurrencias �nicas en la lista de la red
        # completa, se guarda en la primera columna el nodo que aparece antes al
        # ordenar ambas especies alfab�ticamente.
        if (Esp_i != Esp_j){
          par_s <- sort(c(Esp_i,Esp_j))
          adj_s_1[v_s] <- par_s[1]
          adj_s_2[v_s] <- par_s[2]
          # Se actualiza el contador de enlaces entre �rboles de la parcela.
          # En este caso, es muy probable que muchos enlaces est�n repetidos.
          v_s <- v_s + 1
        }
      }
    }
  }
  adj_arb <- unname(cbind(adj_a_1, adj_a_2))
  adj_esp <- unname(cbind(adj_s_1,adj_s_2))
  list_adj <- list("arboles" = adj_arb, "especies" = adj_esp)
  return(list_adj)  
}

# FUNCI�N 3.
# Esta funci�n toma un "bloque horizontal" (conjunto de ocurrencias filtradas
# por encontrarse entre ciertos valores de latitud), generado con la funci�n 
# "analisis_bajo_r".
# Realiza cortes verticales sobre el bloque de forma que genera parcelas que
# se encuentran limitadas por distintos valores de longitud.
# Recoge las ocurrencias que se encuentran dentro de cada parcela y llama a la 
# funci�n adj_list_P, que calcula los enlaces entre �rboles y entre especies
# que se dan en la parcela. 
horizontal <- function(bloque_horizontal,tresh){
  # Se inicializa un vector que guarda cada corte vertical (valor de longitud)
  # realizado sobre el bloque.
  cortes_V <- vector()
  
  # Extremos en la longitud del territorio.
  ini <- MIN_X
  fin <- MAX_X
  
  # Diferencia de longitud entre los extremos.
  long <- fin + (-ini)
  
  # Se realizan 80 cortes verticales sobre el territorio, cada uno con una
  # diferencia de "incremento_long" respecto al anterior.
  incremento_long <- long/80 
  
  # Se define un margen que rodee cada corte realizado a ambos lados
  # La longitud de este margen es el valor de distancia umbral m�s
  # un valor de seguridad que en este caso es de 200 m (bastante conservativo).
  margenkm <- tresh+0.2 
  # El valor anterior est� en kilom�tros. Para convertirlo a un valor de 
  # distancia geod�sica en grados decimales, hemos comprobado que en Andaluc�a
  # un grado de longitud equivale aproximadamente a 88 km, as� que se usa 
  # este n�mero como factor de conversi�n. Debido al error que se pueda estar
  # cometiendo en cada caso y a que hay que evitar que el margen suponga una 
  # longitud menor que la de la distancia umbral (para no perder enlaces entre
  # parcelas), se ha elegido la alta distancia de seguridad (200 m).
  margenG <- margenkm/88 
  
  adj_def_a <- data.table()
  adj_def_s <- data.table()
  
  # Si el margen que rodea cada corte es superior al incremento en latitud,
  # evidentemente, la parcelaci�n no va a ser efectiva computacionalmente
  # (va a producir un aumento del tiempo de computaci�n), por lo que se 
  # avisa en pantalla.
  if (margenG > incremento_long){
    cat("AVISO: MARGEN > INCREMENTO DE LONGITUD")
  }
  # Si la extensi�n corre el riesgo de ser m�s peque�a que el valor de distancia 
  # umbral, tambi�n se avisa en pantalla (se estar�an perdiendo enlaces).
  if(incremento_long*88<tresh){
    cat("ERROR: INCREMENTO DE LONGITUD PARA LA PARCELACI�N DEMASIADO PEQUE�O, ABORTA EL PROGRAMA")
  }
  
  # PRIMERA PARCELA DEL BLOQUE HORIZONTAL.
  # La primera parcela est� definida por un solo corte de longitud
  # (el primero), de forma que la parcela recoge todas las ocurrencias del
  # bloque horizontal que tienen un valor de longitud menor que el del corte.
  #cat("...P1...")
  f <- ini + incremento_long
  P1 <- subset(bloque_horizontal, X < f)
  
  # Con estas ocurrencias se llama a la funcion "adj_list_P" que devuelve
  # la lista de adyacencia de �rboles y de especies en la parcela.
  adj_P1 <- adj_list_P(P1,tresh)
  
  # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
  # definitiva de �rboles del bloque (de todas las parcelas del bloque 
  # horizontal).
  adj_P1_a <- adj_P1$arboles
  adj_P1_a <- as.data.table(adj_P1_a)
  adj_def_a <- rbind(adj_def_a,adj_P1_a)
  
  # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
  # definitiva de especies del bloque (de todas las parcelas del bloque 
  # horizontal).
  adj_P1_s <- adj_P1$especies
  adj_P1_s <- as.data.table(adj_P1_s)
  adj_def_s <- rbind(adj_def_s,adj_P1_s)

  # Se guarda la longitud del primer corte vertical.
  cortes_V[1] <- f
  
  # Para los cortes verticales del 2 al 80:
  for (i in 2:80) {
    # Las parcelas que no est�n en los extremos est�n definidas por
    # un valor inferior y superior de  longitud, siendo el valor inferior el 
    # valor superior de la divisi�n anterior...
    a <- f 
    # y el valor superior corresponde a su propio valor inferior m�s el valor de
    # incremento de longitud.
    b <- a+incremento_long 
    
    # La parcela en cuesti�n recoge las ocurrencias del bloque horizontal que 
    # quedan entre ambos valores de longitud.
    Pi <- subset(bloque_horizontal, X > a & X < b)
    
    # Con estas ocurrencias se llama a la funcion "adj_list_P" que devuelve
    # la lista de adyacencia de �rboles y de especies para la parcela.
    adj_Pi <- adj_list_P(Pi,tresh)
    
    # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
    # definitiva de �rboles (de todos las parcelas del bloque horizontal).
    adj_Pi_a <- adj_Pi$arboles
    adj_Pi_a <- as.data.table(adj_Pi_a)
    adj_def_a <- rbind(adj_def_a,adj_Pi_a)
    
    # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
    # definitiva de especies (de todos las parcelas del bloque horizontal).
    adj_Pi_s <- adj_Pi$especies
    adj_Pi_s <- as.data.table(adj_Pi_s)
    adj_def_s <- rbind(adj_def_s,adj_Pi_s)
  
    # Cuando se acaba con una parcela se guarda su valor superior de longitud
    # para que el siguiente bloque lo coja de referencia.
    f <- b 
    # Se guarda cada valor superior como un nuevo corte.
    cortes_V[i] <- b
  }
  
  # �LTIMA PARCELA DEL BLOQUE
  # Con el valor superior de longitud de la �ltima parcela  se crea la �ltima
  # parcela, que recoge todas las ocurrencias con un valor de longitud superior 
  # a ese valor.
  a <- f
  P81 <- subset(bloque_horizontal, X > a)
  
  # Con estas ocurrencias se llama a la funcion "adj_list_P" que devuelve
  # la lista de adyacencia de �rboles y de especies para la parcela.
  adj_P81 <- adj_list_P(P81,tresh)
  
  # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
  # definitiva de �rboles (de todos las parcelas del bloque horizontal).
  adj_P81_a <- adj_P81$arboles
  adj_P81_a <- as.data.table(adj_P81_a)
  adj_def_a <- rbind(adj_def_a,adj_P81_a)
  
  # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
  # definitiva de especies (de todos las parcelas del bloque horizontal).
  adj_P81_s <- adj_P81$especies
  adj_P81_s <- as.data.table(adj_P81_s)
  adj_def_s <- rbind(adj_def_s,adj_P81_s)
  
  # Eliminaci�n de duplicados en las listas.
  adj_def_a <- distinct(adj_def_a)
  adj_def_s <- distinct(adj_def_s)
  
  # Ahora para cada corte vertical realizado es necesario crear otro "parcela 
  # conectora" con l�mite de longitud inferior igual al valor de lonitud
  # del corte menos el margen y con l�mite superior igual al valor de latitud 
  # del corte m�s el margen.
  n_cortes <- length(cortes_V)
  for (k in 1:n_cortes){
    este <- cortes_V[k]+margenG
    oeste <- cortes_V[k]-margenG
    
    # La parcela conectora recoge las ocurrencias con valores de altitud 
    # comprendidos entre esos dos extremos. 
    Pconect <- subset(bloque_horizontal, X > oeste & X < este)
    
    # Con estas ocurrencias se llama a la funcion "adj_P_list" que devuelve
    # la lista de adyacencia de �rboles y de especies para la parcela
    # conectora.
    adj_Pconect <- adj_list_P(Pconect,tresh)
    
    # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
    # definitiva de �rboles (de todas las parcelas del bloque horizontal).
    adj_Pconect_a <- adj_Pconect$arboles
    adj_Pconect_a <- as.data.table(adj_Pconect_a)
    adj_def_a <- rbind(adj_def_a,adj_Pconect_a)
    # Se van eliminando repeticiones
    adj_def_a <- distinct(adj_def_a)
    
    # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
    # definitiva de especies (de todas las parcelas del bloque horizontal).
    adj_Pconect_s <- adj_Pconect$especies
    adj_Pconect_s <- as.data.table(adj_Pconect_s)
    adj_def_s <- rbind(adj_def_s,adj_Pconect_s)
    # Se van eliminando repeticiones
    adj_def_s <- distinct(adj_def_s)
  }
  
  # La funci�n devuelve el objeto adj_def que consiste en una lista que
  # recoge las listas de adyacencia de �rboles y de especies del bloque 
  # horizontal entero. 
  adj_def <- list("arboles"=adj_def_a,"especies"=adj_def_s)
  return(adj_def)
}

# FUNCI�N PARA QUE EXPORTA UNA SECUENCIA DE NODOS FORMANDO UN CAMINO
# A UN ARCHIVO REPRESENTABLE EN QGIS
# Esta funci�n toma una secuencia (lista) de identificadores de nodos de 
# �rboles que aparecen formando un camino y crea un archivo csv en el que 
# cada nodo aparece con sus coordenadas (tomadas del registro) y con un 
# n�mero natural identificador del orden de aparici�n en la secuencia.
# Esta salida puede ser tomada e interpretada por QGIS para dibujar el 
# camino sobre el mapa. 
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


# FUNCI�N 2. 
# Esta funci�n es la toma un �nico valor de distancia umbral y crea ambos
# grafos (de individuos y de especies) bajo dicho valor. Para ello realiza
# cortes horizontales sobre el territorio para crear "bloques horizontales"
# de ocurrencias y con cada uno de ellos llama a la funci�n "horizontal".
# Integrando las salidas de cada llamada a dicha funci�n, crea ambos grafos.
# Para los dos, calcula todas las propiedades y adem�s de incluirlas en la
# carpeta del grafo, vuelca los resultados sobre los archivos comunes a todos
# los grafos creados bajo distintos valores de distancia umbral en este
# proyecto (cada proyecto corresponde a una ejecuci�n de la funci�n
# "Ejecutar_analisis"). 
analisis_bajo_r <- function(registro,tresh,dir_arboles,dir_especies,dir_proj,CALC_diam){
 
  # N�mero de �rboles en el registro (n� nodos en el grafo de individuos)
  n_arboles <- nrow(registro)
  
  # Se inicializa el vector que almacena los cortes horizontales.
  cortes_H <- vector()
  
  # N�mero de especies diferentes en el registro (n� nodos en el grafo de especies).
  lista_especies <- distinct(registro,Especie)
  n_especies <- nrow(lista_especies)
  
  # Conjunto de nodos en el grafo de especies.
  nodos_a_id <- unname(as.integer(registro$Nodo))
  
  # Se crea la carpeta correspondiente al valor de distancia umbral.
  # Va a estar tanto en la carpeta de �rboles como en la de especies (con
  # distintos contenidos).
  dir_r <- paste0("tresh",tresh)
  setwd(dir_especies)
  dir.create(dir_r)
  setwd(dir_r)
  dir_esp_tresh <- getwd()
  
  setwd(dir_arboles)
  dir.create(dir_r)
  setwd(dir_r)
  dir_arb_tresh <- getwd()
  
  # Se define un margen que rodee cada corte realizado a ambos lados
  # La longitud de este margen es el valor de distancia umbral m�s
  # un valor de seguridad que en este caso es de 200 m (bastante conservativo).
  margenkm <- tresh+0.2 
  # El valor anterior est� en kilom�tros. Para convertirlo a un valor de 
  # distancia geod�sica en grados decimales, hemos comprobado que en Andaluc�a
  # un grado de altitud equivale aproximadamente a 107 km, as� que se usa 
  # este n�mero como factor de conversi�n. Debido al error que se pueda estar
  # cometiendo en cada caso y a que hay que evitar que el margen suponga una 
  # longitud menor que la de la distancia umbral (para no perder enlaces entre
  # parcelas), se ha elegido la alta distancia de seguridad (200 m).
  margenG <- margenkm/107 

  # Se crean las listas de adyacencia (vac�as) correspondientes al grafo de 
  # individuos y al de especies. En ellas se ir�n guardando los enlaces
  # calculados sobre cada parcela. Se ir�n filtrando los enlaces repetidos
  # por haber sido determinados en m�s de una parcela diferente.
  adj_DEF_a <- data.table()
  adj_DEF_s <- data.table()
  
  # Extremos en la latitud del territorio.
  bot <- MIN_Y
  top <- MAX_Y
  # Diferencia de latitud entre los extremos.
  tot_lat <- top-bot
  # Se realizan 40 cortes horizontales sobre el territorio, cada uno con una
  # diferencia de "incremento_lat" respecto al anterior.
  incremento_lat <- tot_lat/40
  
  # Si el margen que rodea cada corte es superior al incremento en latitud,
  # evidentemente, la parcelaci�n no va a ser efectiva computacionalmente
  # (va a producir un aumento del tiempo de computaci�n), por lo que se 
  # avisa en pantalla.
  if (margenG > incremento_lat){
    cat("AVISO: MARGEN > INCREMENTO DE LATITUD")
  }
  
  # Si la extensi�n corre el riesgo de ser m�s peque�a que el valor de distancia 
  # umbral, tambi�n se avisa en pantalla (se estar�an perdiendo enlaces).
  if(incremento_lat*107<tresh){
    cat("ERROR: INCREMENTO DE LATITUD PARA LA PARCELACI�N DEMASIADO PEQUE�O, ABORTA EL PROGRAMA")
  }
  
  cat("Empezando con el primer bloque horizontal H1...")
  
  # El primer bloque horizontal est� definido por un solo corte de altitud
  # (el primero), de forma que el bloque recoge todas las ocurrencias
  # por debajo de dicho valor de altitud.
  f <- bot + incremento_lat
  bloqueH1 <- subset(registro, Y < f)
  
  # Con estas ocurrencias se llama a la funcion "horizontal" que devuelve
  # la lista de adyacencia de �rboles y de especies para el bloque horizontal.
  adj_bloqueH1 <- horizontal(bloqueH1,tresh)
  
  # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
  # definitiva de �rboles (de todos los bloques).
  adj_bloqueH1_a <- adj_bloqueH1$arboles
  adj_DEF_a <- rbind(adj_DEF_a,adj_bloqueH1_a)
  
  # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
  # definitiva de especies (de todos los bloques).
  adj_bloqueH1_s <- adj_bloqueH1$especies
  adj_DEF_s <- rbind(adj_DEF_s,adj_bloqueH1_s)
  
  # Se guarda la altitud del primer corte horizontal.
  cortes_H[1] <- f 
  
  # Para los cortes horizontales del 2 al 39:
  for(i in 2:39) {
    # Los bloques horizontales que no est�n en los extremos est�n definidos por
    # un valor inferior y superior de  altitud, siendo el valor inferior el 
    # valor superior de la divisi�n anterior...
    a <- f 
    # y el valor superior corresponde a su propio valor inferior m�s el valor de
    # incremento de latitud.
    b <- a+incremento_lat
    # El bloque en cuesti�n recoge las ocurrencias del registro que quedan 
    # entre ambos valores de altitud.
    bloqueHi <- subset(registro, Y > a & Y < b)
    
    # Con estas ocurrencias se llama a la funcion "horizontal" que devuelve
    # la lista de adyacencia de �rboles y de especies para el bloque horizontal.
    adj_bloqueHi <- horizontal(bloqueHi,tresh)
    
    # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
    # definitiva de �rboles (de todos los bloques).
    adj_bloqueHi_a <- adj_bloqueHi$arboles
    adj_DEF_a <- rbind(adj_DEF_a,adj_bloqueHi_a)
    # Se van eliminando enlaces repetidos en la lista
    adj_DEF_a <- as.data.table(adj_DEF_a)
    adj_DEF_a <- distinct(adj_DEF_a) 
    
    # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
    # definitiva de especies (de todos los bloques).
    adj_bloqueHi_s <- adj_bloqueHi$especies
    adj_DEF_s <- rbind(adj_DEF_s,adj_bloqueHi_s)
    # Se van eliminando enlaces repetidos en la lista.
    adj_DEF_s <- as.data.table(adj_DEF_s)
    adj_DEF_s <- distinct(adj_DEF_s)
    
    # Se imprime en pantalla el n�mero del bloque horizontal de cada iteraci�n. 
    cat("H",i, "...")
    # Cuando se acaba con un bloque se guarda su valor superior de altitud
    # para que el siguiente bloque lo coja de referencia.
    f <- b 
    # Se guarda cada valor superior como un nuevo corte.
    cortes_H[i] <- f
  }
  
  # Con el valor superior de latitud del �ltimo bloque horizontal creado
  # se crea el �ltimo bloque horizontal, que recoge todas las ocurrencias
  # con un valor de altitud superior a ese valor.
  cat("�LTIMO BLOQUE HORIZONTAL H41...")
  b <- f
  bloqueH41 <- subset(registro, Y > b)
  
  # Con estas ocurrencias se llama a la funcion "horizontal" que devuelve
  # la lista de adyacencia de �rboles y de especies para el bloque horizontal.
  adj_bloqueH41 <- horizontal(bloqueH41,tresh)
  
  # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
  # definitiva de �rboles (de todos los bloques).
  adj_bloqueH41_a <- adj_bloqueH41$arboles
  adj_DEF_a <- rbind(adj_DEF_a,adj_bloqueH41_a)
  adj_DEF_a <- distinct(as.data.table(adj_DEF_a))
  
  # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
  # definitiva de especies (de todos los bloques).
  adj_bloqueH41_s <- adj_bloqueH41$especies
  adj_DEF_s <- rbind(adj_DEF_s,adj_bloqueH41_s)
  adj_DEF_s <- distinct(as.data.table(adj_DEF_s))
  
  cat("TERMINA EL �LTIMO BLOQUE HORIZONTAL (H41)...")
  
  # Ahora para cada corte horizontal realizado es necesario crear otro "bloque
  # horizontal conector" con l�mite inferior igual al valor de latitud del corte
  # menos el margen y con l�mite superior igual al valor de latitud del corte
  # m�s el margen. 
  
  cat("EMPIEZA EL PRIMER BLOQUE HORIZONTAL CONECTOR (HC1)")
  n_cortesH <- length(cortes_H)
  for (k in 1:n_cortesH){
    superior <- cortes_H[k]+margenG
    inferior <- cortes_H[k]-margenG
    
    # El bloque conector recoge las ocurrencias con valores de altitud 
    # comprendidos entre esos dos extremos. 
    BHCi <- subset(registro, Y > inferior & Y < superior)
    
    # Con estas ocurrencias se llama a la funcion "horizontal" que devuelve
    # la lista de adyacencia de �rboles y de especies para el bloque horizontal
    # conector.
    adj_BHCi <- horizontal(BHCi,tresh)
    
    # Los enlaces de la lista de �rboles se a�aden a la lista de adyacencia
    # definitiva de �rboles (de todos los bloques).
    adj_BHCi_a <- adj_BHCi$arboles
    adj_DEF_a <- rbind(adj_DEF_a,adj_BHCi_a)
    adj_DEF_a <- distinct(as.data.table(adj_DEF_a))
    
    # Los enlaces de la lista de especies se a�aden a la lista de adyacencia
    # definitiva de especies (de todos los bloques).
    adj_BHCi_s <- adj_BHCi$especies
    adj_DEF_s <- rbind(adj_DEF_s,adj_BHCi_s)
    adj_DEF_s <- distinct(as.data.table(adj_DEF_s))
    cat("HC",k)
  }
  cat("SE HA COMPLETADO EL �LTIMO BLOQUE HORIZONTAL CONECTOR HC40")
  cat("...TODOS LOS ENLACES RECOGIDOS Y FILTRADOS..")
  cat("...EMPIEZA LA CREACI�N DEL GRAFO DE ARBOLES...")
  
  # Ruta de la carpeta de la distancia umbral en cuesti�n dentro de la carpeta
  # de �rboles.
  setwd(dir_arb_tresh) 
  
  # Lo primero es guardar la lista de adyacencia definitiva para poder
  # volver a crear el grafo en cualquier momento.
  
  # Nombre de la lista de adyacencia de �rboles bajo ese valor de r.
  nomb_adj <- paste0("adj_",tresh,".csv")
  file.create(nomb_adj)
  cat("Escribiendo lista de adyacencia...")
  fwrite(adj_DEF_a, file = nomb_adj, sep=",",row.names=F, col.names=F)
  cat("DONE")
  
  # CREACI�N DEL GRAFO DE �RBOLES.
  
  # Data frame de la lista de adyacencia.
  tramit <- as.data.frame(adj_DEF_a)
  # Data frame del conjunto de identificadores de nodos (�rboles).
  nodos_a_id <- as.data.frame(nodos_a_id)
  cat("Creando grafo...")
  # Grafo no dirigido.
  graf <- graph_from_data_frame(tramit, directed=FALSE, vertices=nodos_a_id)
  cat("DONE...")
  
  # C�LCULO DE PROPIEDADES Y PAR�METROS DEL GRAFO DE �RBOLES.
  # Capa propiedad calculada se escribe sobre un archivo com�n a todos los 
  # grafos de �rboles (con distinto valor de distancia umbral) de la poblaci�n.
  
  cat("Calculando n� v�rtices y ejes...")
  num_ver <- gorder(graf)
  num_edges <- gsize(graf)
  tresh_num_edges <- as.data.table(cbind(tresh, num_edges))
  d_edges <- paste0(dir_arboles,"/","num_conexiones.csv")
  fwrite(tresh_num_edges, file = d_edges, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando modularidad...")
  louv_mod <- cluster_louvain(graph = graf)
  mod <- round(modularity(graf, membership = louv_mod$membership),3)
  tresh_mod <- as.data.table(cbind(tresh, mod))
  d_modu <- paste0(dir_arboles,"/","modularidad.csv")
  fwrite(tresh_mod, file = d_modu, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("calculando coef clustering...")
  cc <- round(transitivity(graf, type = "average"),3)
  tresh_cc <- as.data.table(cbind(tresh,cc))
  d_cc <- paste0(dir_arboles,"/","coef_clustering.csv")
  fwrite(tresh_cc, file = d_cc, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("calculando coef clustering ZEROS...")
  cc_zeros <- round(transitivity(graf, type = "average", isolates = "zero"),3)
  tresh_cc_zeros <- as.data.table(cbind(tresh,cc_zeros))
  d_cc2 <- paste0(dir_arboles,"/","coef_clustering_zeros.csv")
  fwrite(tresh_cc_zeros, file = d_cc2, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  # El di�metro se calcula solo si lo ha solicitado el usuario (CALC_diam = T)
  if (CALC_diam){
    cat("Calculando di�metro...")
    path_d <-get_diameter(graf, directed= FALSE, unconnected = TRUE, weights = NULL)[[]]
    path_d <- as_ids(path_d)
    diam <- length(path_d)-1 
    tresh_diam <- as.data.table(cbind(tresh, diam))
    diam_dir <- paste0(dir_arboles,"/","diametro.csv")
    fwrite(tresh_diam, file = diam_dir, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    setwd(dir_arb_tresh)
    
    # Se crea, tambi�n en la carpeta del tresh en la carpeta de �rboles, un 
    # archivo que contiene la secuencia de nodos que componenen el di�metro.
    # Para ello se llama a la funci�n create_path_doc que toma como entrada
    # el objeto con la secuencia (path_d) y el nombre del archivo (que contiene
    # el valor de distancia umbral correspondiente para identificar el grafo).
    nombre_diam_doc <- paste0("diametro",tresh,".csv")
    create_path_doc(path_d,registro,nombre_diam_doc)
  }
   
  
  cat("Calculando densidad de enlaces...")
  dens <- edge_density(graf)
  tresh_dens <- as.data.table(cbind(tresh,dens))
  d_dens <- paste0(dir_arboles,"/","densidad.csv")
  fwrite(tresh_dens, file = d_dens, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando N� de comunidades...")
  graf$community <- louv_mod$membership
  ncomun <- length(unique(graf$community))
  tresh_ncomun <- as.data.table(cbind(tresh, ncomun))
  d_comuns <- paste0(dir_arboles,"/","num_comunidades.csv")
  fwrite(tresh_ncomun, file = d_comuns, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando LCC...")
  nrow_adj <- nrow(adj_DEF_a)
  LCC <- length(component_distribution(graf,mode="weak"))-1 
  tresh_LCC <- as.data.table(cbind(tresh, LCC))
  d_LCC <- paste0(dir_arboles,"/","LCC.csv")
  fwrite(tresh_LCC, file = d_LCC, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando % de nodos no aislados...")
  pc_parti <- length(V(graf)[degree(graf)!=0])*100/num_ver
  tresh_pc_parti <- as.data.table(cbind(tresh, pc_parti)) 
  d_parti <- paste0(dir_arboles,"/","pc_nodos_no_aislados.csv")
  fwrite(tresh_pc_parti, file = d_parti, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  
  cat("DONE...")
  
  cat("Calculando num componentes conexas...")
  #num de componentes conexas
  ncompons <- count_components(graf)
  tresh_ncompons <- as.data.table(cbind(tresh, ncompons))
  d_compons <- paste0(dir_arboles,"/","num_componentes.csv")
  fwrite(tresh_ncompons, file = d_compons, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando grado medio...")
  #grado medio
  gm <- mean(degree(graf))
  tresh_gm <- as.data.table(cbind (tresh,gm))
  d_gm <- paste0(dir_arboles,"/","grado_medio.csv")
  fwrite(tresh_gm, file = d_gm, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("calculando distribuci�n de grados...")
  #distribuci�n de grados - ESTO NO ES EL DE LA CARPETA ARBOLES SINO DE LA DEL TRESH (dentro de �rboles)
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
  cat("Escribiendo carpeta del treshold ",tresh, " para el grafo de �rboles...")
  options(scipen=999)
  options(digits=10)
  
  dist_tr <- tresh
  # Se recogen todos los par�metros.
  if (CALC_diam){
    params <- c(num_ver,dist_tr,num_edges,mod,ncomun,LCC,ncompons,cc,cc_zeros,dens,pc_parti,gm,diam)
  } else {
    params <- c(num_ver,dist_tr,num_edges,mod,ncomun,LCC,ncompons,cc,cc_zeros,dens,pc_parti,gm)
  }
  
  nparams = length(params)
  # Se guardan en una matriz
  matt <- matrix(params,nrow=nparams,ncol=1)
  colnames(matt) <- c("VALOR")
  if (CALC_diam){
    row.names(matt) <- c("N� nodos","Radio de vecindad (km)","N� conexiones","Modularidad","N� comunidades","LCC","N� componentes","CC","CC_zeros","Densidad de enlaces","% nodos conectados","Grado medio","Di�metro")
  } else {
    row.names(matt) <- c("N� nodos","Radio de vecindad (km)","N� conexiones","Modularidad","N� comunidades","LCC","N� componentes","CC","CC_zeros","Densidad de enlaces","% nodos conectados","Grado medio")
  }
  matt <- as.data.table(matt, keep.rownames = TRUE)
  nomb_archivo_prop <- paste0("propiedades",tresh,".csv")
  fwrite(matt, file = nomb_archivo_prop, sep=",",row.names=T, col.names=T)
  
  # ------------------------------------------------------------
  # C�LCULO DE PROPIEDADES Y PAR�METROS DEL GRAFO DE ESPECIES.
  # Capa propiedad calculada se escribe sobre un archivo com�n a todos los 
  # grafos de especies (con distinto valor de distancia umbral) de la poblaci�n.
  cat("...EMPIEZA LA CREACI�N DEL GRAFO DE ESPECIES...")
  setwd(dir_esp_tresh) 
  
  # Lo primero es guardar la lista de adyacencia definitiva para poder
  # volver a crear el grafo en cualquier momento.
  file.create(nomb_adj)
  if (nrow(adj_DEF_s) != 0){
  cat("Escribiendo lista de adyacencia...")
  fwrite(adj_DEF_s, file = nomb_adj, sep=",",row.names=F, col.names=F)
  cat("DONE")
  #GRAFOS
  tramit <- as.data.frame(adj_DEF_s)
  nodos_s_id <- as.data.frame(lista_especies)
  cat("Creando grafo...")
  graf <- graph_from_data_frame(tramit, directed=FALSE, vertices=nodos_s_id)
  cat("DONE...")
  
  cat("Calculando num v�rtices y ejes...")
  num_ver <- gorder(graf)
  num_edges <- gsize(graf)
  tresh_num_edges <- as.data.table(cbind(tresh, num_edges))
  d_edges <- paste0(dir_especies,"/","num_conexiones.csv")
  fwrite(tresh_num_edges, file = d_edges, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando modularidad...")
  louv_mod <- cluster_louvain(graph = graf)
  mod <- round(modularity(graf, membership = louv_mod$membership),3)
  tresh_mod <- as.data.table(cbind(tresh, mod))
  d_mod <- paste0(dir_especies,"/","modularidad.csv")
  fwrite(tresh_mod, file = d_mod, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("calculando coef clustering...")
  cc <- round(transitivity(graf, type = "average"),3)
  tresh_cc <- as.data.table(cbind(tresh,cc))
  d_cc <- paste0(dir_especies,"/","coef_clustering.csv")
  fwrite(tresh_cc, file = d_cc, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("calculando coef clustering ZEROS...")
  cc_zeros <- round(transitivity(graf, type = "average", isolates = "zero"),3)
  tresh_cc_zeros <- as.data.table(cbind(tresh,cc_zeros))
  d_cc <- paste0(dir_especies,"/","coef_clustering_zeros.csv")
  fwrite(tresh_cc_zeros, file = d_cc, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando densidad...")
  dens <- edge_density(graf)
  tresh_dens <- as.data.table(cbind(tresh,dens))
  d_dens <- paste0(dir_especies,"/","densidad.csv")
  fwrite(tresh_dens, file = d_dens, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando num comunidades...")
  graf$community <- louv_mod$membership
  ncomun <- length(unique(graf$community))
  tresh_ncomun <- as.data.table(cbind(tresh, ncomun))
  d_comuns <- paste0(dir_especies,"/","num_comunidades.csv")
  fwrite(tresh_ncomun, file = d_comuns, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando LCC...")
  nrow_adj <- nrow(adj_DEF_s)
  LCC <- length(component_distribution(graf,mode="weak"))-1 
  tresh_LCC <- as.data.table(cbind(tresh, LCC))
  d_LCC <- paste0(dir_especies,"/","LCC.csv")
  fwrite(tresh_LCC, file = d_LCC, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando num de nodos aislados...")
  pc_parti <- length(V(graf)[degree(graf)!=0])*100/num_ver
  tresh_pc_parti <- as.data.table(cbind(tresh, pc_parti))
  d_parti <- paste0(dir_especies,"/","pc_nodos_no_aislados.csv")
  fwrite(tresh_pc_parti, file = d_parti, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando N� componentes conexas...")
  ncompons <- count_components(graf)
  tresh_ncompons <- as.data.table(cbind(tresh, ncompons))
  d_compons <- paste0(dir_especies,"/","num_componentes.csv")
  fwrite(tresh_ncompons, file = d_compons, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  cat("DONE...")
  
  cat("Calculando grado medio...")
  gm <- mean(degree(graf))
  tresh_gm <- as.data.table(cbind (tresh,gm))
  d_gm <- paste0(dir_especies,"/","grado_medio.csv")
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
  cat("Escribiendo carpeta del tresh",tresh)
  dist_tr <- tresh
  params <- c(num_ver,dist_tr,num_edges,mod,ncomun,LCC,ncompons,cc,cc_zeros,dens,pc_parti,gm)
  matt <- matrix(params,nrow=12,ncol=1)
  colnames(matt) <- c("VALOR")
  row.names(matt) <- c("N� nodos","Radio de vecindad (km)","N� conexiones","Modularidad","N� comunidades","LCC","N� componentes","CC","CC_zeros","Densidad de enlaces","% nodos conectados","Grado medio")
  matt <- as.data.table(matt, keep.rownames = TRUE)
  nomb_arxi_prop <- paste0("propiedades",tresh,".csv")
  fwrite(matt, file = nomb_arxi_prop, sep=",",row.names=T, col.names=T)
  
  
  cat("tresh",tresh,"COMPLETADO")
  # Si no se ha podido formar el grafo de especies (por no haber ning�n enlace)
  # se avisa en pantalla.
  }else{cat("\nCero conexiones, no se forma grafo de especies para tresh",tresh,"\n")}
  
  # Se vuelve a la carpeta del proyecto.
  setwd(dir_proj)
  #return(graf)
}


# FUNCI�N 1.
# La �nica funci�n que el usuario necesita llamar para lanzar el an�lisis.
# Toma los tres par�metros ya comentados como entrada. Crea las carpetas
# necesarias y llama con cada valor de distancia umbral (del vector
# vector_tresh) a la funci�n "analisis_bajo_r" que es la que crea y 
# analiza los grafos.
Ejecutar_analisis <- function(vector_tresh,registro,CALC_diam){
  
  # Se crea en el directorio desde el que se ha llamado a la funci�n,
  # un directorio con la fecha del proyecto.
  old_wd <- getwd()
  date <- Sys.Date()
  dir_proj <- paste0("Proyecto",date)
  dir.create(dir_proj)
  
  # Se trabaja dentro de ese directorio a partir de ahora.
  setwd(dir_proj)
  # Se guarda la ruta.
  dir_proj <- getwd()
  
  # All� se crean dos carpetas, la de �rboles y la de especies, en las que
  # se guardar�n los resultados de cada una de las redes.
  carpetas <- c("Arboles","Especies")
  
  for (i in 1:2){
    dir.create(carpetas[i])
    setwd(carpetas[i])
    if (i==1){
      # Se guarda la ruta de la carpeta de �rboles.
      dir_arboles <- getwd()
      # En el caso de que se haya marcado la opci�n de calcular el di�metro (para
      # la red de �rboles) se crea en este momento el archivo correspondiente a 
      # dicha propiedad dentro de la carpeta de �rboles.
      if(CALC_diam){
        file.create("diametro.csv")
        # Se escribe el encabezado de los campos del archivo del di�metro.
        dimm <- as.data.table(cbind("radio (km)", "di�metro"))
        fwrite(dimm, file = "diametro.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
      }
      
    } else{
      # Se guarda la ruta de la carpeta de especies.
      dir_especies <- getwd()
    }
    
    # Y se crean todos estos archivos en ambas carpetas
    
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
    
    cc <- as.data.table(cbind("radio(km)","coef.clustering"))
    fwrite(cc, file = "coef_clustering.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    cc2 <- as.data.table(cbind("radio(km)","coef.clustering.zeros"))
    fwrite(cc2, file = "coef_clustering_zeros.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    dens <- as.data.table(cbind("radio (km)","densidad conexiones"))
    fwrite(dens, file = "densidad.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    LCC <- as.data.table(cbind("radio (km)", "LCC"))
    fwrite(LCC, file = "LCC.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    modul <- as.data.table(cbind("radio (km)", "modularidad"))
    fwrite(modul, file = "modularidad.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    comps <- as.data.table(cbind("radio (km)", "n� componentes"))
    fwrite(comps, file = "num_componentes.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    comuns <- as.data.table(cbind("radio (km)", "n� comunidades"))
    fwrite(comuns, file = "num_comunidades.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    conexs <- as.data.table(cbind("radio (km)", "n� conexiones"))
    fwrite(conexs, file = "num_conexiones.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    aisl <- as.data.table(cbind("radio (km)", "nodos conectados"))
    fwrite(aisl, file = "pc_nodos_no_aislados.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    gm <- as.data.table(cbind("radio (km)", "grado medio"))
    fwrite(gm, file = "grado_medio.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    
    
    #Se vuelve al directorio del proyecto.
    setwd(dir_proj)
  }
  
  # Para cada valor de distancia umbral (r) del vector "vector_tresh", se llama a 
  # una segunda funci�n para el an�lisis de la red de individuos y de especies 
  # de la poblaci�n con este valor de radio de vecindad.
  n_tresh <- length(vector_tresh)
  for (t in 1:n_tresh){
    actual_tresh <- vector_tresh[t]
    cat("\n EMPEZANDO CON RADIO DE VECINDAD",actual_tresh,"...")
    act_graf <- analisis_bajo_r(registro,actual_tresh,dir_arboles,dir_especies,dir_proj,CALC_diam)
  }
}




# Ejemplo de ejecuci�n de an�lisis
vect_tresh <- c(0.002,0.005,0.01,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38)

Ejecutar_analisis(vect_tresh,registro,CALC_diam = T)





