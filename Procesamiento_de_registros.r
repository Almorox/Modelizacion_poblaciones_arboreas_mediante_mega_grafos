# Modelizaci�n de datos de geolocalizaci�n de especies arb�reas mediante 
# redes complejas. Aplicaci�n al caso andaluz.

# SCRIPT PARA EL PROCESAMIENTO DEL ARCHIVOS DE OCURRENCIAS DE GBIF.

rm(list=ls())
options(scipen=999) 


# Instalaci�n y carga de paquetes.
list.of.packages <- c("dplyr","data.table","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# Funci�n que toma como entrada el archivo CSV de GBIF y devuelve otro archivo con solo 3
# atributos: longitud, latitud y especie. Adem�s, elimina instancias con coordeandas 
# repetidas. 

# A la hora de guardar, con write.table, un conjunto de coordenadas �nicas en un archivo CSV 
# y luego leerlo, un peque�o porcentaje de instancias aparecen como repetidas
# (seguramente por alg�n tipo de tratamiento de los decimales), por ello se decidi� primero
# exportar todas las ocurrencias a un archivo temporal controlando el n�mero de decimales (6)
# y luego leer este archivo para proceder a la eliminaci�n de instancias repetidas.

# Se consideran 6 coordenadas decimales debido a que son m�s que suficientes para
# para distinguir dos �rboles cercanos (aunque depende de d�nde nos encontremos,
# un grado se suele aproximar como 111,1 km, por lo que con seis decimales seremos capaces
# de distinguir distancias que difieren de 11,11 cm - no vamos a encontrar dos �rboles a esta
# distancia -)

Procesamiento <- function(archivo_GBIF,nombre_final_csv){
  # Lee el archivo de GBIF.
  long_checklist <- fread(archivo_GBIF)
  # Selecciona �nicamente estos tres atributos.
  XYS <- long_checklist[ , c("decimalLongitude", "decimalLatitude","species")]
  # N�mero de instancias en el archivo sin procesar.
  principio <- nrow(XYS)
  # Si hubiera alguna instancia con uno de los 3 atributos de inter�s vac�os, se
  # elimina la instancia.
  XYS <- na.omit(XYS)
  colnames(XYS) <- c("x","y","s")
  file.create("XYS_TEMP.csv")
  # Se vuelca el objetos con los 3 atributos de inter�s sobre un archivo temporal CSV
  # con seis decimales en las coordenadas (contiene instancias repetidas).
  write.table(format(XYS,nsmall=6), file = "XYS_TEMP.csv", sep=",",row.names=F, col.names=T)
  # Se lee dicho archivo y se eliminan las instancias repetidas.
  XYS_read <- read_delim("XYS_TEMP.csv")
  XYS_distinct <- na.omit(distinct(XYS_read, x, y, .keep_all = TRUE))
  num_nodos <- nrow(XYS_distinct)
  # Se asocia un n�mero entero identificador de nodo a cada instancia (�rbol).
  nodos_id <- c(1:num_nodos)
  # El registro final tiene los 3 campos anteriores m�s el campo de identificador de nodo.
  # Se crea y se vuelca el registro sobre otro archivo y se borra el archivo temporal anterior.
  registro <- data.table(Nodo=nodos_id,X=XYS_distinct$x,Y=XYS_distinct$y,Especie=XYS_distinct$s)
  file.create(nombre_final_csv)
  write.table(registro, file = nombre_final_csv, sep=",",row.names=F, col.names=T)
  unlink("XYS_TEMP.csv")
  
  #Se comprueba que efectivamente todas las ocurrencias del archivo sean �nicas.
  prueba <- read_delim(nombre_final_csv)
  cat("---------------------------------------------------- \n")
  cat("N�mero de ocurrencias en el archivo sin procesar:", principio,"\n")
  cat("Comprobando que en el archivo procesado no haya ocurrencias repetidas...\n")
  cat("N� de ocurrencias tras leer el archivo:",nrow(prueba),"\n")
  XYS_def <- na.omit(distinct(prueba, X, Y, .keep_all = TRUE))
  cat("N� de ocurrencias tras eliminar ocurrencias repetidas:",nrow(XYS_def),"\n")
  
  if (nrow(prueba) == nrow(XYS_def)){
    cat("COMPROBADO.")
  } else {cat("Algo no funciona")}
  
}


