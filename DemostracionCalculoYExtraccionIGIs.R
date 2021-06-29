#CARGO MIS FUNCIONES ARTESANALES
source("/home/diego/Dropbox/practicas_General/CodigoR/funciones/igiCalculator.R")
source("/home/diego/Dropbox/practicas_General/CodigoR/funciones/igiExtractor.R")
##################################################################################
# Demostración con L.major
setwd("/home/diego/Dropbox/practicas_General/datosTritryp/LmajorFriedlin")
#Calculo y guardo las IGIs del primer organismos
igisDF = igiCalculator(gffFile="TriTrypDB-48_LmajorFriedlin.gff")
#Puedo guardar la tabla en formato csv (o el que considere más adecuado)
write.csv(igisDF, "IGIs_LmajorFriedlin.csv", row.names = FALSE)
#xtraigo y guardo las secuencias de las IGIS calculadas previamente
igisSequences = igiExtractor(fastaFile="TriTrypDB-48_LmajorFriedlin_Genome.fasta", positions=igisDF)
#Puedo guardar las secuencias en formato fasta
Biostrings::writeXStringSet(igisSequences, filepath="IGIs_LmajorFriedlin.fasta")
#####################################################################################