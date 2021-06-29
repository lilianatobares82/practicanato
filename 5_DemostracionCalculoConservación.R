#El análisis se realiza en 4 pasos
#1) Carga de datos de IGIs
#2) ALINEAMIENTOS
#3) BUSQUEDA DE PQSs en las secuencias consenso
#4) Generación de LOGOS de los PQSs de interés
################################################

#RETOQUES PRELIMINARES:####
#Creo un vector con los nombres resumidos de los genomas, que me ayudará a trabajar
#ATENCION: Estos nombres deben coincidir entre los archivos a utilizar (de la carpeta outputs_1)
# y las columnas de la tabla de ortólogos
genomesNames = c("Ldonovani", "Linfantum", "Lmajor", "Tbrucei", "TcruziE", "TcruziNE", "Tevansi")
#Abro la tabla de salida de Roary con los nómbres de las proteínas ortólogas
setwd("C:/Users/BANGHO/Dropbox/practicas_Liliana")
ortologos = read.csv("C:/Users/BANGHO/Dropbox/practicas_Liliana/script5/Galaxy45-Roary on data 7,data 6and others Gene Presence Absence.csv",stringsAsFactors = F)
#Reviso la distribución de secuencias por grupo de ortólogos

table(ortologos$No..sequences)
#Grafico la distribución de secuencias por grupo de ortólogos
plot(table(ortologos$No..sequences), main="ortologos, blast70, MCL0.5")
####
#Renombro las columnas con el vector de nombres creado al principio
colnames(ortologos)[15:21]
#Me aseguro que el vector de nombres esté en el orden correcto (mismo que las columnas)
#(Esto debe ser verificado manualmente)
genomesNamesReordered=genomesNames
#Renombro las columnas
colnames(ortologos)[15:21]=genomesNamesReordered
#Obtengo los ortólogos con representantes en los x genomas y reduzco la tabla
x = 7
ortologos = ortologos[ortologos$No..sequences==x,c(15:21)]

#1)Carga de datos de IGIs ############################################################################
#Cargo la lista de archivos de la carpeta 1 (Primero los csv y luego los fasta)
IGIsCsvs=list.files("../script1", pattern = "*.csv")
IGIsFasta=list.files("../script1", pattern = "*.fasta")
#Leo las secuencias IGIs completas y sus anotaciones (me interesa conocer las direcciones de las hebras)list.files("./Outputs_1")
#Cargo todo en una lista
IGIsData=list()
IGIsData$csv=list()
IGIsData$fasta=list()
setwd("./script1")
for (i in 1:7){
  #Tomo los nombres individuales de los genomas
  genomeName=genomesNamesReordered[i]
  #Obtengo los paths de cada archivo (de a un genoma por vez)
  pathCsv=IGIsCsvs[grep(genomeName,IGIsCsvs)]
  pathFasta=IGIsFasta[grep(genomeName,IGIsFasta)]
  #Abro y guardo en la lista los archivos csv y fasta de cada genoma
  print(paste0("Cargando datos de ", genomeName, "..."))
  IGIsData$csv[[genomeName]]=read.csv(pathCsv, stringsAsFactors = FALSE)
  IGIsData$fasta[[genomeName]]=Biostrings::readDNAStringSet(pathFasta)
  
}

#2)ALINEAMIENTOS###########################################################################
# Para cada grupo de ortólogos, tomo las secuencias de cada genoma, 
# las agrupo en un mismo objeto DNAStringSet y aplico los diferentes
# algoritmos de alineamiento
library(msa)
alignments=list()
alignments$clustalOmega=list()
alignments$clustalW=list()
alignments$muscle=list()
for(i in 1:nrow(ortologos)){#bucle 1 para ir analizando cada grupod e ortólogos
  multiSeq=Biostrings::DNAStringSet()
  k=0;#para controlar cuantas secuencias se invierten
  for (j in 1:7){#bucle 2, toma los IGIs IDs de cada columna de la tabla de ortólogos
    #Tomo los nombres individuales de los genomas
    genomeName=genomesNamesReordered[j]
    #Tomo el ID del IGI para el grupo de ortólogos (j) y genóma actual (i)
    IGIname=ortologos[i,colnames(ortologos)==genomeName]
    if(IGIname!=""){#Extraigo la secuencia
      seq=IGIsData$fasta[[genomeName]][IGIname]
      #IMPORTANTE: hago corrección de hebra (reverseComplement para hebras negativas)
      strand=as.character(IGIsData$csv[[genomeName]][IGIsData$csv[[genomeName]]$ID==IGIname,]$strand)
      if(strand=="-"){
        seq=Biostrings::reverseComplement(seq)
        k=k+1
      }
      multiSeq=c(multiSeq, seq)
   }
      
  }
  print(paste0("Alineando ortologos numero ",i,". (Se han invertido ",k," secuencias antes del alineamiento)..."))
  print(multiSeq)
  # alignments$clustalOmega[[i]] = msa(multiSeq,"ClustalOmega")
  # alignments$clustalW[[i]] = msa(multiSeq,"ClustalW")
  ti = Sys.time()
  alignments$muscle[[i]] = msa(multiSeq,"Muscle")
  te=Sys.time()
  print(paste0("Alineamiento numero ",i,"exitoso. Se ha demorado: ", te-ti, " min"))
}
#Guardo los alineamientos como objeto R###
# save(alignments,file="../script6/alignments.Rdata")
#save(alignments,file="../script6/alignments_6sec.Rdata")
# load("../Outputs_6/alignments.Rdata")
###

#3) BUSQUEDA DE PQSs en las secuencias consenso##########################################################
#Obtengo las secuencias consenso de los alineamientos y hago
#una búsqueda de PQSs en dichas secuencias.
#Adicionalmente obtengo los valores de conservación de cada nt
#(Variarán según como se construya la matriz de sustitución "mat")
mat <- nucleotideSubstitutionMatrix(4, -1)
mat <- cbind(rbind(mat, "-"=-8), "-"=-8)
library(pqsfinder)
#Creo df donde almacenaré toda la info de los pqs
PQSfullData=data.frame(ortologo="",pqsID="",start="",width="",score="",
                       strand="", pqs="", conservationMeanIn="",
                       conservationMeanOut="", conservationDiff="",
                       stringsAsFactors = F)
#Bucleo para cada alineamiento y calculo los PQSs
for(i in 1:length(alignments$muscle)){
  print(paste0("Calculando pqs de alineamiento ",i))
  alineamiento = alignments$muscle[[i]]
  consenso = msaConsensusSequence(alineamiento)
  #Obtengo los valores de conservaciónd e cada nucleótido
  conservation = msaConservationScore(alineamiento, mat, gapVsGap=0)
  #Calculo pqs a partir de consenso
  consenso=gsub("\\?","N",consenso)
  consenso=DNAString(consenso)
  consensusPQSs=pqsfinder(consenso)
  print(paste0("Se han encontrado ", length(consensusPQSs), " cuadruples en la secuencia consenso ",i))
  if(length(consensusPQSs)>0){
    pqsData=data.frame()
    for(j in 1:length(consensusPQSs)){
      start=start(consensusPQSs[j])
      end=start(consensusPQSs[j])+width(consensusPQSs[j])-1
      #Obtengo las conservaciones medias para el pqs y su entorno inmediato (+-25nt)
      conservationMeanIn=mean(conservation[start:end])
      conservationMeanOut=mean(c(conservation[(start-25):start],conservation[end:(end+25)]))
      conservationDiff=conservationMeanIn-conservationMeanOut
      #Cargo todo en un vector y lo agrego a la tabla general
      x=c(ortologo=i, pqsID=j, start=start, width=width(consensusPQSs[j]),
          score=score(consensusPQSs[j]), strand=strand(consensusPQSs[j]),
          pqs=as.character(subseq(consensusPQSs[j]@subject, start, end)),
          conservationMeanIn=conservationMeanIn, conservationMeanOut=conservationMeanOut,
          conservationDiff=conservationDiff)
      PQSfullData=rbind(PQSfullData,x)
    }
  }
}
PQSfullData=PQSfullData[-(PQSfullData$ortologo==""),]
#Guardo en formato CSV
write.csv(PQSfullData,"../Outputs_6/PQSs_fullData_6.csv", row.names = FALSE)
list.files("..")
#4) Generación de LOGOS de los PQSs de interés############################################################
#install.packages("ggseqlogo")
library(ggseqlogo)
#Grafico el entorno los PQSs consenso
for(i in 1:nrow(PQSfullData)){
  alignNumber=as.integer(PQSfullData[i,]$ortologo)
  pqsID=as.integer(PQSfullData[i,]$pqsID)
  start=as.integer(PQSfullData[i,]$start)
  end=start+as.integer(PQSfullData[i,]$width)-1
  alineamiento = alignments$muscle[[alignNumber]]
  subseqs=c()
  for(j in 1:length(alineamiento@unmasked)){
    temp=as.character(subseq(alineamiento@unmasked[j], start-10,end+10))
    subseqs=c(subseqs,temp)
  }
  png(filename=paste0("../Outputs_6/PQSlogo_ortoGroup_",alignNumber,"_pqsID_",pqsID,".png"),width = 1000, height = 400)
  print(ggseqlogo(subseqs, method = 'bits', seq_type='DNA'))
  dev.off()
}

#OPCIONAL
#Grafico el entorno los PQSs consenso (solo en genomas de Leishmanias)
# for(i in 1:nrow(PQSfullData)){
#   alignNumber=as.integer(PQSfullData[i,]$ortologo)
#   start=as.integer(PQSfullData[i,]$start)
#   end=start+as.integer(PQSfullData[i,]$width)-1
#   alineamiento = alignments$muscle[[alignNumber]]
#   subseqs=c()
#   for(j in 5:7){#Las 3 últimas secuencias del alineamiento corresponden a los Leishmania
#     temp=as.character(subseq(alineamiento@unmasked[j], start-10,end+10))
#     subseqs=c(subseqs,temp)
#   }
#   png(filename=paste0("../Outputs_6/PQSlogo_ortoGroup_",alignNumber,"_pqs_",i,".png"),width = 1000, height = 400)
#   print(ggseqlogo(subseqs, method = 'bits', seq_type='DNA'))
#   dev.off()
# }
