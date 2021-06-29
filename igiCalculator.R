#Función simple para abrir gff de un organismo y calcular las IGIs
#Nota: No se tienen en cuenta posibles solapamientos de genes
#Input: anotaciones en formato GFF
#Output: dataframe con las posiciones start-end de todas las IGIs
igiCalculator = function (gffFile=NULL){
  #Controlo que se hayan ingresado el único argumento (gff)
  if(is.null(gffFile)){
    #Si falta algun argumento, salgo de la función y aviso del error
    stop("Debe ingresar anotaciones en formato GFF.")
  }
  #Abro el GFF
  print("Cargando archivo gff...")
  annotations = rtracklayer::readGFF(gffFile)
  #Proceso las anotaciones para obtener una tabla############################################################
  #Paso a data.frame para operar más rápido (no es lo óptimo)
  df = as.data.frame(annotations)
  #Mantengo unicamente las filas del tipo gene (secuencias génicas...proteínas y otras)
  cdsRows=df$type=="gene"
  print(paste0("Se reducen las anotaciones de ",length(cdsRows), " registros a ",sum(cdsRows)," registros"))
  df = df[cdsRows,]
  #Dejo las columnas que me interesan: 
  columns = colnames(df)%in%c("seqid", "type", "start", "end", "strand", "ID", "description")
  df = df[,columns]
  
  #creo copia de la tabla de anotaciones para ir llenando con las posiciones de las IGIs
  igisTable =df
  igisTable$type=as.character(igisTable$type)
  
  #Identifico todos los cromosomas
  sequencesIDs = levels(df$seqid)
  for (j in 1:length(sequencesIDs)){
    chromosome = sequencesIDs[j]
    print(paste0("Obteniendo IGIs del cromosoma ",chromosome))
    temporalSeq = df[df$seqid==chromosome,]
    #Trabajo de a un cromosoma por vez######
    #Ordeno las secuencias génicas según su posición
    temporalSeq = temporalSeq[order(temporalSeq$start),]
    #Obtengo la igi número 1 (start, end e ID)
    igiStart=1
    igiEnd = temporalSeq[2,]$start-1
    igiId = temporalSeq[1,]$ID
    #Reemplazo los nuevos valores start y end en la tabla Igis
    igisTable[igisTable$ID==igiId,c(2,3,4)] = c("IGI",igiStart,igiEnd)
    #Obtengo las igis segunda a ante-última
    for (i in 2:(nrow(temporalSeq)-1)){
      print(i)
      #Calculo los start y end de cada igi y extraigo su id
      igiStart = temporalSeq[i-1,]$end+1
      igiEnd = temporalSeq[i+1,]$start-1
      igiId = temporalSeq[i,]$ID
      #Reemplazo los nuevos valores start y end en la tabla Igis
      igisTable[igisTable$ID==igiId,c(2,3,4)] = c("IGI",igiStart,igiEnd)
    }
    #Obtengo la igi final
    igiStart=temporalSeq[nrow(temporalSeq)-1,]$end+1
    igiEnd = 0#Desconozco la longitud del cromosoma en este momento...lo seteo a cero y después actualizo el valor real
    igiId = temporalSeq[nrow(temporalSeq),]$ID
    #Reemplazo los nuevos valores start y end en la tabla Igis
    igisTable[igisTable$ID==igiId,c(2,3,4)] = c("IGI",igiStart,igiEnd)
    
  }
  ###########################################################################################
  #devuelvo el resultado fuera de la función
  igisTable$start = as.numeric(igisTable$start)
  igisTable$end = as.numeric(igisTable$end)
  return(igisTable)
}



#Prueba local de la función (luego se elimina o comenta)
#setwd("/home/diego/MEGA/UNLaR/Practicanato/datos Tritryp/Leishmania/")
#igisDF = igiCalculator(gffFile="TriTrypDB-48_LmajorFriedlin.gff")
#Puedo guardar la tabla en formato csv (o el que considere más adecuado)
#...



