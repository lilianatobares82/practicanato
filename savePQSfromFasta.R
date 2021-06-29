#Función para obtener las posiciones de los pqs a partir de un set de secuencias (.fasta)
#input: archivo fasta con las secuencias IGIs a evaluar
#output: lista conteniendo dos data frames:
#-tabla con el número de pqs detectado en cada IGI
#-tabla con los detalles de posición, hebra y score de cada pqs
savePQSfromFasta = function (lista=NULL, nombreBase=NULL){
  name1 = paste0("pqsByIGIs_",nombreBase,".csv")
  name2 = paste0("pqsPositions_",nombreBase,".csv")
  a = list.files(pattern = name1)
  b = list.files(pattern = name2)
  if(length(a)==0 & length(b)==0){
    print(paste0("Guardando '",name1, "' y '",name2,"'"))
    write.csv(lista$pqsByIGIs, name1, row.names = FALSE)
    write.csv(lista$pqsPositions, name2, row.names = FALSE)
  }else{
    i = 1
    while(length(a)>0 || length(b)>0){
      #cambio el nombre y vuelvo a chequear ocurrencia
      name1 = paste0("pqsByIGIs_",nombreBase,"_",i,".csv")
      name2 = paste0("pqsPositions_",nombreBase,"_",i,".csv")
      a = list.files(pattern = name1)
      b = list.files(pattern = name2)
      i=i+1
    }
    #Cuando salgo del bucle, guardo
    warning(paste0("Ya hay un archivo con el mismo nombre.\n  Los datos actuales se guardan con nombres nuevos:\n '",name1, "' y '",name2,"'"))
    
    print(paste0("Guardando '",name1, "' y '",name2,"'"))
    write.csv(lista$pqsByIGIs, name1, row.names = FALSE)
    write.csv(lista$pqsPositions, name2, row.names = FALSE)
  }
}


#Ejecuto la función 
# savePQSfromFasta(lista=salida, nombreBase="LmajorFriedlin")
# #Guardo los resultados
# write.csv(lista$pqsByIGIs, "pqsByIGIs_LmajorFriedlin.csv", row.names = FALSE)
# write.csv(lista$pqsPositions, "pqsPositions_LmajorFriedlin.csv", row.names = FALSE)


