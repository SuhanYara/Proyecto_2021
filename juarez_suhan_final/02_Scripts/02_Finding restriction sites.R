###################################
###### Locating Restriction Sites
####################################

#Ejercicio a resolver:

#A partir de una secuencia de al menos 1000 pb ubicar la posición y longitud de cada palindromo inverso en la cadena
##Palindromo inverso: serán secuencias que tienen la misma secuencia en su complemento(de dercha a izquierda)

#Primero cargamos una secuencia de 1000 pb
library(Biostrings)
Bacillus<-readDNAStringSet("01_Raw_Data/Bacillus_pseudomycoides.fasta")
Bacillus

Bac_matrix<-as.matrix(Bacillus)
Bac_matrix
Sec_Bacillus<-as.vector(Bac_matrix)
Sec_Bacillus
length(Sec_Bacillus)


##Primero necesitamos la cadena complementaria (C-G y A-T), yo usé lo siguiente:
#Busqué la posición de las bases y después un reemplazo por el complementario
Adenin<-which(Sec_Bacillus =="A")
Adenin
Timin<-which(Sec_Bacillus =="T")
Timin
Citocin<-which(Sec_Bacillus =="C")
Citocin
Guanin<-which(Sec_Bacillus =="G")
Guanin

Complement<-replace(Sec_Bacillus,Adenin,"T")
Complement<-replace(Complement,Timin,"A")
Complement<-replace(Complement,Citocin,"G")
Complement<-replace(Complement,Guanin,"C")

#Guardé en un objeto esto:
Complement

#Para hacer una comparación con la cadena principal, reacomodar la cadena complementaria de final a inicio
ReversoC<-Complement[1000:1]
ReversoC

#Ahora hay que buscar los palindromos
Posicion<-1
Palindromo4<-c()

for (Posicion in 1:length(Sec_Bacillus)) {
  if (all(Sec_Bacillus[Posicion:(Posicion+3)]==ReversoC[Posicion:(Posicion+3)])) {
    Palindromo4<-c(Palindromo4, Posicion)
  }
  Posicion<-Posicion +1
}

##Aquí vemos el número de palindromos
Palindromo4_Número<-length(Palindromo4)
Palindromo4_Número

##Aquí las posiciones
Posiciónes_Palindromos4<-Palindromo4
Posiciónes_Palindromos4

##Con esto veo la secuencia de cada palindromo
Sec_Bacillus[117:120]
Sec_Bacillus[881:884]


###Si sólo existen un límitado número de palindromos de tamaño 4, como el ejemplo, sólo podríamos buscar
#en esas posiciones algunos con tamaño mayor. O aplicar el comando anterior y adicionar para 5 y hasta 12

##Aquí únicamente haré una comparación simple, porque sólo tengo 2

Posiciónes_Palindromos4

##Marca 117 y 881, y esas buscaré con una comparación simple para ver si son iguales

##Para palindromos de 5
all(Sec_Bacillus[117:121]==ReversoC[117:121])
all(Sec_Bacillus[881:885]==ReversoC[881:885])

##Ambas son falsas, así que no podremos continuar si sólo tenemos estos 

#Aquí cargamos algunas de las secuencias de de enzimas más comunes

Pac1<-c("TTAATTAA")
SgFI<-c("GCGATCGC")
Sse83871<-c("CCTGCAGG")
Notl<-c("GCGGCCGC")
RpaI<-c("GTYGGAG")
RpaBI<-c("CCCGCAG")
RpaB51<-c("CGRGGAC")
HindIII<-c("AAGCTT")
EcoRV<-c("GATATC")
XbaI<-c("TCTAGA")
FokI<-c("GGATG")
Hind4II<-c("CCTTC")
DpnI<-c("GATC")
TaqI<-c("TCGA")
HaeIII<-c("GGCC")
AluI<-c("AGCT")

##Aquí el programita para determinar si alguna secuencia se puede utilizar para cortar nuestra secuencia

W<-0
while(W == 0) {
  Extension<-readline(prompt = "Ingresa le extensión de tu palindromo (4:8): ")
  Extension<-as.numeric(Extension)
  Insertar<-readline(prompt = "Ingresa la secuencia,(usa mayúsculas, no insertes espacios o caracteres especiales):")
  Insertar<-as.character(Insertar)
  
if (Extension == 8) {
  if (Insertar == Pac1) {
    print("Necesitas la enzima Pac1 para cortar tu secuencia")
  } else if (Insertar == SgFI) {
    print("Necesitas la enzima SgFI para cortar tu secuencia")
  } else if (Insertar == Sse83871) {
    print("Necesitas la enzima Sse83871 para cortar tu secuencia")
  } else if (Insertar == Notl) {
    print("Necesitas la enzima Notl para cortar tu secuencia")
  } else {"No se encuentra en el catálogo"}
  
} else if (Extension == 7) {
  if (Insertar == RpaI) {
    print("Necesitas la enzima RpaI para cortar tu secuencia")
  } else if (Insertar == RpaBI) {
    print("Necesitas la enzima RpaBI para cortar tu secuencia")
  } else if (Insertar == RpaB51) {
    print("Necesitas la enzima RpaB51 para cortar tu secuencia")
  } else {"No se encuentra en el catálogo"}
  
} else if (Extension == 6) {
  if (Insertar == HindIII) {
    print("Necesitas la enzima HindIII para cortar tu secuencia")
  } else if (Insertar == EcoRV) {
    print("Necesitas la enzima RpaBI para cortar tu secuencia")
  } else if (Insertar == XbaI) {
    print("Necesitas la enzima XbaI para cortar tu secuencia")
  } else {"No se encuentra en el catálogo"}
  
} else if (Extension == 5) {
  if (Insertar == FokI) {
    print("Necesitas la enzima FokI para cortar tu secuencia")
  } else if (Insertar == Hind4II) {
    print("Necesitas la enzima Hind4II para cortar tu secuencia")
  } else {"No se encuentra en el catálogo"}
  
} else if (Extension == 4) {
  if (Insertar == DpnI) {
    print("Necesitas la enzima DpnI para cortar tu secuencia")
  } else if (Insertar == TaqI) {
    print("Necesitas la enzima TaqI para cortar tu secuencia")
  } else if (Insertar == HaeIII) {
    print("Necesitas la enzima HaeIII para cortar tu secuencia")
  } else if (Insertar == AluI) {
    print("Necesitas la enzima AluI para cortar tu secuencia")
  } else {"No se encuentra en el catálogo"}
}
  
  W<-readline(prompt = "¿Continuar comparando? (Si=0, No=1): ")
  W<-as.numeric(W)
  
}
