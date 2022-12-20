install.packages("ggplot2")
install.packages("ggrepel")
install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "DESeq2", "org.Hs.eg.db", "biomaRt", "edgeR", "tweeDEseq", "GOstats", "tweeDEseqCountData", "annotate", "biomaRt"))
#abrimos el archivo SRE
load(file = "path_to_data_file/rse_gene_pleura.Rdata")
stage <- rse_gene$gdc_cases.diagnoses.tumor_stage
ids.early <-grep(paste("stage i$", "stage ia$","stage ib$","stage ic$", "stage ii$", "stage iia$","stage iib$", "stage iic$",sep="|"), stage)
ids.late <-grep(paste("stage iii$", "stage iiia$", "stage iiib$","stage iiic$", "stage iv$", "stage iva$","stage ivb$", "stage ivc$",sep="|"), stage)
# create an empty column named GROUP
colData(rse_gene)$GROUP <-rep(NA, ncol(rse_gene))
# add early for those patients with tumours at stages i-ii
colData(rse_gene)$GROUP[ids.early] <- "early"
# add late for those patients with tumours at stages iii-iv
colData(rse_gene)$GROUP[ids.late] <- "late"
#creamos un data frame para representarlo en un grafico
dataGroups <- data.frame(Stage = rse_gene$gdc_cases.diagnoses.tumor_stage, Group = rse_gene$GROUP)
library(ggplot2)
ggplot(data=dataGroups, mapping = aes(x=Stage, y=Group))+geom.point()
#no entiendo como he hecho este grafico. Geom_col es para cuando tienes la tabla precomuntada (ya tengo una tabla en la que indica que en Stage1 hay 2 personas, en el 2 hay 3,...)
ggplot(data=dataGroups, mapping = aes(x=Stage, fill=Group))+
  geom_bar()
#ahora voy a hacer lo mismo pero poniendo en el eje y que se lea Number of patients en lugar de count
ggplot(data=dataGroups, mapping = aes(x=Stage, fill=Group))+
  geom_bar()+
  labs(y="Number of patients")

#comprobamos que tenemos la información del stage del cancer de todos los pacientes
#en naData guardan cuales de los pacientes del dataset rseo_gene que tienen valor na en la variable GROUP
naData <- is.na(rse_gene$GROUP)#q me escanee la columba GROUP em busca de "na" (ver dataGroups para entenderlo)
#visualizo naData
table(naData)
#Me da FALSE 87 --> hay 87 que es FALSO que son na --> hay 87 que tienen datos
dim(rse_gene)
#hay 58037 genes analizados (filas) de 87 personas (columnas)
#no eliminamos nada porque no tenemos ningún "na"

table(rse_gene$GROUP)
#hay 26 con un tumor en early stage y 61 en late stage

counts <- assay(rse_gene, "counts")#con counts veo para cada gen cuántos reads han mapeado para cada muestra
counts[1:5, 1:2]#quiero ver las 5 primeras filas y las 2 primeras columnas

phenotype <- colData(rse_gene)
phenotype[1:5,1:2]

#compruebo si los individuos están en lso dos datasets (quiero comprobar que tengo el genotipo (counts) y fenotipo(phenotype) de todo el mundo)
identical(colnames(counts),rownames(phenotype))
#como pone TRUE es q tengo las dos infos de todos los indvs

#guardo en un vector el id de las personas que coinciden en ambos (q son todos)
individualsCommon <- intersect(colnames(counts),rownames(phenotype))
#filtramos el dataset de count para quedarme solo con los que están en ambos (q son todos)
count <- counts[,individualsCommon]
#filtramos el dataset de fenotipos para quedarme solo con los que están en ambos (q son todos)
phenotype <- phenotype[individualsCommon,]


annotation <- rowData(rse_gene)
annotation
#ahora he guardado los datos de rse_gene, q es un archivo random en un data frame q puedo visualizar

#Can you explore which information is available for each gene
head(annotation)#aqui vemos q información hay para cada gen (id, length y symbol)
head(dataGroups)#aquí vemos cada individuo q tumor tiene (en q estadio) y el grupo al q pertenece (early o late), entonces el número de individuos será el número de filas
dim(dataGroups)#así veo cuántos indv tengo
table(dataGroups$Stage)#ahora hago una tabla de frecuencias de individuos q me agrupe en función de las posibilidades que haya en la variable stage


#2 NORMALIZATION

#necesitamos normalizar?? 
#para eso haré un MA-plot
maPlot(counts[,1],
       counts[,2],
       pch=19,
       cex=.5,
       ylim=c(-8,8),
       allCol="darkgray",
       lowess=TRUE,
       xlab=expression(A == log[2] (sqrt(S1/N %.% S2/N))),
       ylab=expression(M == log[2](S1/N)-log[2](S2/N)))
grid(col="black")
title("Raw data")
#espero q los genes estén distribuidos simétricamente, es decir, que la línea roja vaya justo por el 0. Podemos ver que hay una región en la que no está exactamente sobre 0, entonces sí que hay q normalizar (por esta mini desviación del a línea roja)
#primero nos guardamos la longitud del gen
geneLength <- annotation$bp_length
#en annotation teníamos una columna bp_length, entonces en gene_Length me guardo esa columna y las filas (es decir, se ha quedado con la iformación de la longuitud de cada gen)

#ahora normalizamos por el método RPKM
counts.rpkm <- t(t(counts/geneLength*1000)/colSums(counts)*1e6)
#he guardado en counts.rpkm las longitudes normalizadas según la fórmula
counts.rpkm[1:5,1:2]

maPlot(counts.rpkm[,1],
       counts.rpkm[,2],
       pch=19,
       cex=.5,
       ylim=c(-8,8),
       allCol="darkgray",
       lowess=TRUE,
       xlab=expression( A==log[2] (sqrt(S1/N%.%S2/N))),
       ylab=expression(M==log[2](S1/N)-log[2](S2/N)))
grid(col="black")
title("RPKM")
#la línea sigue sin estar en 0 exactamente

#ahora normalizo por el método TMM
counts.tmm <- normalizeCounts(counts, method = "TMM")
counts.tmm[1:5,1:2]
maPlot(counts.tmm[,1],
       counts.tmm[,2],
       pch=19, cex=.5,
       ylim=c(-8,8),
       allCol="darkgray",
       lowess=TRUE,
       xlab=expression( A==log[2] (sqrt(S1/N%.%S2/N))),
       ylab=expression(M==log[2](S1/N)-log[2](S2/N)))
grid(col="black")
title("TMM")


#entre las 3 opciones: sin normaliza, RPKM y TMM queremos ver cuándo los datos siguen una distribución simétrica y eso lo vemos cuando la líne roja va justo por 0
#la que tiene la línea roja sobre el 0 durante más tiempo es la TMM


#3DIFFERENTIAL EXPRESSION ANALYSIS  

# stage of each patient
pheno.stage <- subset(phenotype, select=GROUP)

# recreate the counts in a new matrix (una fila por cada gen y una columna por cada muestra)
counts.adj <- matrix((as.vector(as.integer(counts))), nrow=nrow(counts), ncol=ncol(counts))
#la unica diferencia q hay entre count y count.adj es q, si me fijo, en counts los datos están guardados como caracter (chr) y en counts.adj como número (int)

rownames(counts.adj) = rownames(counts)#quiro q las filas se llamen igual
colnames(counts.adj) = colnames(counts)#quiero q las columnas se llamen igual

# check information
identical(colnames(counts.adj), rownames(pheno.stage))
#efectivamente los genotipos y  los fenotipos siguen perteneciendo a cada persona
counts[1:5,1:2]
counts.adj[1:5,1:2]

