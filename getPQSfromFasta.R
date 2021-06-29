#Función para obtener las posiciones de los pqs a partir de un set de secuencias (.fasta)
#input: archivo fasta con las secuencias IGIs a evaluar
#output: lista conteniendo dos data frames:
#-tabla con el número de pqs detectado en cada IGI
#-tabla con los detalles de posición, hebra y score de cada pqs
library("pqsfinder")
getPQSfromFasta = function (fastaFile=NULL){
  IGIs = Biostrings::readDNAStringSet(fastaFile)
  pqsPositions = data.frame()
  pqsByIGIs = data.frame()
  #Bucleo cada secuencia, obtengo sus pqs y almaceno en dos tablas (número de PQSs por secuencia IGI,
  #Detalles de PQSs)
  for (i in 1:length(IGIs)){
    print(paste0("Procesando IGI ",i," de ",length(IGIs)))
    #Realizo la búsqueda de pqs
    pqsIGI = pqsfinder(IGIs[[i]])
    seq_id = names(IGIs[i])
    #Obtengo el número de pqs obtenidos anteriormente
    pqsNumber = length(pqsIGI)
    
    pqsByIGIs=rbind(pqsByIGIs,data.frame(seq_id=seq_id, pqsNumber = pqsNumber))
    if(pqsNumber > 0){
      temp = data.frame(seq_id=seq_id, inicio=start(pqsIGI), fin=end(pqsIGI),
                        strand=strand(pqsIGI), score=score(pqsIGI))
      pqsPositions = rbind(pqsPositions, temp)
    }
  }
  #Cargo tablas en una lista
  salida = list(pqsByIGIs = pqsByIGIs, pqsPositions = pqsPositions)
  return(salida)
}

#Ejecuto la función (alrededor de 10 min si no hay problemas)
# lista = getPQSfromFasta("IGIs_LmajorFriedlin.fasta")
# #Guardo los resultados
# write.csv(lista$pqsByIGIs, "pqsByIGIs_LmajorFriedlin.csv", row.names = FALSE)
# write.csv(lista$pqsPositions, "pqsPositions_LmajorFriedlin.csv", row.names = FALSE)


