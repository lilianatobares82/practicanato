#Función simple para extraer las secuencias de ADN a partir de un .fasta y una tabla con posiciones
#Input: genoma en formato fasta, tabla de posiciones de IGIs calculada con la función igiCalculator()
#Output: DNAStringSet con las IGIs individuales
igiExtractor = function (fastaFile=NULL, positions=NULL){
  #Controlo que se hayan ingresado los dos argumentos (fasta y tabla con posiciones)
  if(is.null(fastaFile) || is.null(positions)){
    #Si falta algun argumento, salgo de la función y aviso del error
    stop("Debe ingresar las secuencias genómicas en formato fasta y la tabla con posiciones de las IGIs")
  }
  #Abro el genoma
  print("Cargando genoma...")
  genome = Biostrings::readDNAStringSet(fastaFile)
  #Creo un contenedor para las secuencias que iré recortando
  igisExtracted = Biostrings::DNAStringSet()
  #Recorto de a un cromosoma por vez
  seqNumber = length(genome)
  for (i in 1:seqNumber){
    chromosome = genome[i]
    #Extraigo el nombre del cromosoma
    chrName = strsplit(names(chromosome),"|", TRUE)[[1]][1]
    chrName = trimws(chrName)#Remuevo espacioes en blanco
    print(paste0("Extrayendo secuencias IGIs del cromosoma ",chrName))
    #Obtengo un subset de la tabla igis
    chrAnnotations = igisDF[igisDF$seqid==chrName,]
    #Recuerdo que debo corregir el "end" de la última secuencia, que figura como cero.
    #Reemplazo por la longitud del cromosoma
    chrLong = Biostrings::width(chromosome)
    chrAnnotations[chrAnnotations$end==0,]$end=chrLong
    #Obtengo todas las IGIs del cromosoma
    if(sum(chrAnnotations$end-chrAnnotations$start<1)>0){#Si hay algun rango negativo, lo elimino
      chrAnnotations = chrAnnotations[chrAnnotations$end-chrAnnotations$start>=1,]
      print("Detectadas secuencias con rangos negativos. No se consideran para la extracción de IGIs")
    }
    view = Biostrings::Views(chromosome[[1]], start=chrAnnotations$start, end=chrAnnotations$end, width=NULL, names=chrAnnotations$ID)
    #Paso a Objeto DNAStringSet
    temporalDSS = Biostrings::DNAStringSet(view)
    #Cargo en el objeto contenedor principal
    igisExtracted = c(igisExtracted,temporalDSS)
  }
  #Devuelvo el resultado
  return(igisExtracted)
}


#Prueba local de la función (luego se elimina o comenta)
#setwd("/home/diego/MEGA/UNLaR/Practicanato/datos Tritryp/Leishmania/")
#igisSequences = igiExtractor(fastaFile="TriTrypDB-48_LmajorFriedlin_Genome.fasta", positions=igisDF)


#Puedo guardar las secuencias en formato fasta
#Biostrings::writeXStringSet(igisSequences, filepath="TriTrypDB-48_LmajorFriedlin_IGIs.fasta")
