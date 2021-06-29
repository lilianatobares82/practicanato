#Cargo la función que me automatiza la busqueda de pqs en un conjunto de secuencias
source("C:/Users/BANGHO/Dropbox/practicas_Liliana/funciones/getPQSfromFasta.R")
#Cargo la función que me guarda en forma segura los resultados obtenidos
source("C:/Users/BANGHO/Dropbox/practicas_Liliana/funciones/savePQSfromFasta.R")

library(pqsfinder)
#Ejemplo con los datos de T.evansis
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.evansis")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_TevansiSTIB805_Genome.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "T.evansis")
list.files()
########################################################################
#Ejemplo con los datos de L.majorFriedlin
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/L.major")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_LmajorFriedlin.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "L.major")
########################################################################
#Ejemplo con los datos de T.cruzi(TcruziCLBrenerNon-Esmeraldo-like)
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.cruzi")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_TcruziCLBrenerNon-Esmeraldo-like.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "T.cruzi")

########################################################################
#Ejemplo con los datos de T.cruzi(TcruziCLBrenerEsmeraldo-like)
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.cruzi")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_TriTrypDB-48_TcruziCLBrenerEsmeraldo-like.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "T.cruzi")
########################################################################
#Ejemplo con los datos de L.donovani
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/L.donovani")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_LdonovaniBPK282A1.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "L.donovani")
########################################################################
#Ejemplo con los datos de T.brucei
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.brucei")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_TbruceigambienseDAL972.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "T.brucei")
########################################################################
#Ejemplo con los datos de L.infantum
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/L.infantum")
inicio = Sys.time()
datosPQSs = getPQSfromFasta("IGIs_LinfantumJPCM5.fasta")
final = Sys.time()
# Guardo los resultados...solo si no existen previamente!!!
savePQSfromFasta(lista = datosPQSs, nombreBase = "L.infantum")

