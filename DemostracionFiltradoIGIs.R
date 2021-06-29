setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/L.major")
pqsByIGIs = read.csv("pqsByIGIs_Lmajor.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_Lmajor.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_Lmajor.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)

#############################################################
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.cruzi")
pqsByIGIs = read.csv("pqsByIGIs_T.cruziCLBrener-Esmeraldo-like.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_T.cruziCLBrener-Esmeraldo-like.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_T.cruziCLBrener-Esmeraldo-like.csv", row.names = F)
#############################################################
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.cruzi")
pqsByIGIs = read.csv("pqsByIGIs_T.cruziCLBrenerNon-Esmeraldo-like.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_T.cruziCLBrenerNon-Esmeraldo-like.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_T.cruziCLBrenerNon-Esmeraldo-like.csv", row.names = F)
#############################################################
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/L.donovani")
pqsByIGIs = read.csv("pqsByIGIs_L.donovani.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_L.donovani.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_L.donovani.csv", row.names = F)
#############################################################
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.evansis")
pqsByIGIs = read.csv("pqsByIGIs_T.evansis.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_T.evansis.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_T.evansis.csv", row.names = F)
#############################################################
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/T.brucei")
pqsByIGIs = read.csv("pqsByIGIs_T.brucei.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_T.brucei.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_T.brucei.csv", row.names = F)
#############################################################
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana/L.infantum")
pqsByIGIs = read.csv("pqsByIGIs_L.infantum.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_L.infantum.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 90#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_L.infantum.csv", row.names = F)

