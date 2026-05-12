#Instalamos las librerías necesarias si no las tenemos ya: 
if (!require("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")

if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")  

if (!require("edgeR", quietly = TRUE)) 
  BiocManager::install("edgeR")

if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

if (!require("cmapR", quietly = TRUE))
  BiocManager::install("cmapR")

if (!require("DT", quietly = TRUE))
  install.packages("DT")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!require("ggrepel", quietly = TRUE))
  install.packages("ggrepel")

if (!require("stringr", quietly = TRUE))
  install.packages("stringr")


#Cargamos las librerías que vamos a usar:
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(dplyr)
library(edgeR)
library(limma)
library(ggplot2)
library(biomaRt)
library(pheatmap)
library(ggrepel)
library(cmapR)
library(stringr)

#Definimos la ruta al directorio en el que se descargarán los datos y lo creamos:
dir_gdc <- "C:/Users/anaal/OneDrive - UNIVERSIDAD DE GRANADA/TCGA/GDCdata" #CAMBIAR EN EL FUTURO
dir.create(dir_gdc, recursive = TRUE, showWarnings = FALSE) #Con recursive true nos aseguramos de que se creen las carpetas intermediarias que puedan faltar

#Descargamos datos de expresión de LUAD y LUSC
#y generamos datatables de ambos que podremos guardar para consulta:

query_expr_LUAD<-GDCquery(
  project="TCGA-LUAD",
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  workflow.type="STAR - Counts"
)
GDCdownload(query_expr_LUAD, directory = dir_gdc)
expr_LUAD_raw<-GDCprepare(query_expr_LUAD, directory = dir_gdc)

datatable(
  as.data.frame(colData(expr_LUAD_raw)),
  options = list(scrollX = TRUE, pageLength = 5),
  rownames = FALSE
)

query_expr_LUSC<-GDCquery(
  project="TCGA-LUSC",
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  workflow.type="STAR - Counts"
)
GDCdownload(query_expr_LUSC, directory = dir_gdc)
expr_LUSC_raw<-GDCprepare(query_expr_LUSC, directory = dir_gdc)

datatable(
  as.data.frame(colData(expr_LUSC_raw)),
  options = list(scrollX = TRUE, pageLength = 5),
  rownames = FALSE
)


#Descargamos datos de mutación de LUAD y LUSC:

query_mut_LUAD<-GDCquery(
  project="TCGA-LUAD",
  data.category="Simple Nucleotide Variation",
  data.type="Masked Somatic Mutation",
  workflow.type="Aliquot Ensemble Somatic Variant Merging and Masking",
  access="open"
)
GDCdownload(query_mut_LUAD, directory = dir_gdc)
mut_LUAD<-GDCprepare(query_mut_LUAD, directory = dir_gdc)

query_mut_LUSC<-GDCquery(
  project="TCGA-LUSC",
  data.category="Simple Nucleotide Variation",
  data.type="Masked Somatic Mutation",
  workflow.type="Aliquot Ensemble Somatic Variant Merging and Masking",
  access="open"
)
GDCdownload(query_mut_LUSC, directory = dir_gdc)
mut_LUSC<-GDCprepare(query_mut_LUSC, directory = dir_gdc)


#Descargamos los datos clínicos de LUAD y LUSC:

query_clinical_LUAD<-GDCquery(
  project="TCGA-LUAD",
  data.category="Clinical",
  data.type="Clinical Supplement",
  data.format="BCR Biotab"
)
GDCdownload(query_clinical_LUAD, directory = dir_gdc)
clinical_LUAD_BCRtab<-GDCprepare(query_clinical_LUAD, directory = dir_gdc)

query_clinical_LUSC<-GDCquery(
  project="TCGA-LUSC",
  data.category="Clinical",
  data.type="Clinical Supplement",
  data.format="BCR Biotab"
)
GDCdownload(query_clinical_LUSC, directory = dir_gdc)
clinical_LUSC_BCRtab<-GDCprepare(query_clinical_LUSC, directory = dir_gdc)

#Guardo los objetos que he obtenido en un nuevo directorio al que llamaré "data":
dir.create("data")
saveRDS(expr_LUAD_raw, file="data/expr_LUAD_raw.rds")
saveRDS(expr_LUSC_raw, file="data/expr_LUSC_raw.rds")
saveRDS(mut_LUAD, file="data/mut_LUAD_raw.rds")
saveRDS(mut_LUSC, file="data/mut_LUSC_raw.rds")
saveRDS(clinical_LUAD_BCRtab, file="data/clinical_LUAD_BCRtab.rds")
saveRDS(clinical_LUSC_BCRtab, file="data/clinical_LUSC_BCRtab.rds")

#---------------------------------------
#DATOS DE EXPRESIÓN: eliminamos pacientes duplicados
#expr_LUAD es un SummarizedExperiment en el que las filas son los genes y las columnas las muestras
#de RNA-seq, cada columna corresponde a un archivo .tsv
#cuando TCGAbiolinks hace gdcprepare() pone el barcode de la muestra como nombre de la columna
#----------------------

procesar_expr<-function(expr, 
                        clinical_patients){
  #Extraer identificador de paciente:
  barcodes<-colnames(expr)
  patient_id<-substr(barcodes,1,12)
  colData(expr)$patient_id<-patient_id
  
  #Extraemos sample type:
  sample_type_code<-substr(barcodes,14,15)
  colData(expr)$sample_type_code<-sample_type_code
  
  sample_type<-case_when(
    sample_type_code=="01"~"Tumor",
    sample_type_code=="11"~"Normal"
  )
  colData(expr)$sample_type<-sample_type
  expr<-expr[,!is.na(colData(expr)$sample_type)] #como solo guardamos los sample type de las 01 y 11, quitamos en las que no haya nada, porque podrán ser por ejemplo, 02 (tumor solido recurrente) que no nos interesan
  
  #Actualizamos barcodes porque hemos quitado columnas de expr
  barcodes<-colnames(expr)
  
  #Calcular tamaño de librería:
  library_size<-colSums(assay(expr))
  colData(expr)$library_size<-library_size
  
  #Dataframe auxiliar:
  coldata_df<-as.data.frame(colData(expr))
  coldata_df$barcode<-barcodes
  
  #Seleccionamos una muestra por paciente y por tipo de tumor (la de mayor library size, en este caso):
  #De cada paciente se guardarán como máximo los datos de dos muestras: una de tumor y una sana:
  samples_to_keep<-pull(
    slice_max(
      group_by(coldata_df, patient_id, sample_type),
      order_by=library_size, n=1, with_ties=FALSE
    ),
    barcode
  )
  
  #Nuestros datos de expresión sin duplicados:
  expr_unique<-expr[,samples_to_keep]
  
  #Cruce con datos clínicos:
  patients_expr<-colData(expr_unique)$patient_id
  common_patients<-intersect(patients_expr, clinical_patients)
  
  expr_def<-expr_unique[,colData(expr_unique)$patient_id %in% common_patients]
  
  #Si un paciente tiene TUMOR y NORMAL, quedarnos solo con NORMAL, puesto que es la menos abundante en este dataset:
  
  df<-as.data.frame(colData(expr_def))
  
  #Agrupa las filas por paciente; se queda solo con los que tengan mas de un tipo de muestra, es decir, más de una fila (n>1);
  #extrae la columna de patient_id y elimina duplicados. Por tanto pareados será un vector con los patient_id de los pacientes que tengan 
  #más de una muestra:
  pareados<- df %>%
    group_by(patient_id) %>%
    filter(n()>1) %>%
    pull(patient_id) %>%
    unique()
  
  #se queda con todo menos con las muestras tumorales de pacientes pareados: 
  expr_def<-expr_def[,!(
    colData(expr_def)$patient_id %in% pareados & 
      colData(expr_def)$sample_type=="Tumor"  
  )]
  
  list(
    expr_def=expr_def,
    common_patients=common_patients
  )
  
}

#Aplicamos la función que acabamos de crear a LUAD para obtener los perfiles de expresión definitivos:
res_LUAD<-procesar_expr(
  expr=expr_LUAD_raw,
  clinical_patients=clinical_LUAD_BCRtab$clinical_patient_luad$bcr_patient_barcode
)

expr_LUAD_def<-res_LUAD$expr_def
common_patients_LUAD<-res_LUAD$common_patients

saveRDS(expr_LUAD_def,"data/expr_LUAD_def.rds")
saveRDS(common_patients_LUAD, "data/common_patients_LUAD.rds")

#Hacemos lo mismo para LUSC:
res_LUSC<-procesar_expr(
  expr=expr_LUSC_raw,
  clinical_patients=clinical_LUSC_BCRtab$clinical_patient_lusc$bcr_patient_barcode
)

expr_LUSC_def<-res_LUSC$expr_def
common_patients_LUSC<-res_LUSC$common_patients

saveRDS(expr_LUSC_def,"data/expr_LUSC_def.rds")
saveRDS(common_patients_LUSC, "data/common_patients_LUSC.rds")

preprocess_edgeR<-function(expr_data,                    #Objeto Summarized Experiment
                           group_variable="sample_type", #Nombre de la columna en colData
                           reference_level="Normal",
                           plot_prefix="Edge R",
                           plots=list(boxplot=TRUE, mds=TRUE, pca=TRUE)){
  #Extraemos matriz de counts:
  #si hay NA que se quede con 0 y se asegura de que todo sea numerico 
  counts<-assay(expr_data)  #Assay es una función del paquete SummarizedExperiment. Devuelve la matriz de expresión (filas=genes, columnas=muestras, valores=counts)
  counts<-as.matrix(counts) #para asegurar que sea matriz estandar
  mode(counts)<-"numeric"
  counts[is.na(counts)]<-0  #busca los NA y los convierte en 0
  
  #Definimos grupo:
  group<-factor(colData(expr_data)[[group_variable]]) #la variable "tumor" o "normal" se convierte en variable categórica al usar factor. Con el doble corchete se extrae el contenido de esa columna
  group<-relevel(group,ref=reference_level) #define el grupo de referencia en el modelo, en este caso, las muestras "normales"
  
  #Creamos DGElist,estructura base para el análisis de edgeR:
  dgl<-DGEList(counts=counts,group=group) 
  
  info_genes<-as.data.frame(rowData(expr_data))
  info_genes$id_cruce<-sub("\\..*", "", rownames(info_genes)) #el id_cruce es el Ensembl ID sin el nº de versión
  
  #Añadimos tipo de gen:
  dgl$genes <- info_genes[rownames(dgl), c("id_cruce", "gene_name", "gene_type")]
  
  #Filtrado de genes poco expresados:
  keep<-filterByExpr(dgl, group=group)  #elimina genes con muy baja expresión, reduce ruido, mejora potencia estadística etc. Es la que más se adapta a las condiciones experimentales
  #umbral_muestras<-ceiling(0.10*ncol(counts))
  #keep<-rowSums(counts>=10)>=umbral_muestras
  
  cat("Genes totales:", nrow(counts), "\n")
  cat("Genes eliminados:", sum(!keep), "\n")
  cat("Genes retenidos:", sum(keep), "\n")
  
  dgl<-dgl[keep,,keep.lib.sizes=FALSE] #nos quedamos solo con los genes que hayan pasado el filtro, hacemos que se recalculen los tamaños de biblioteca
  rm(keep)
  
  #Normlización TMM:
  dgl<-calcNormFactors(dgl, method="TMM") #Normalización Trimmed Mean of M-Values, corrige diferencias de composición, de profundidad de secuenciación...
  
  #Generamos los diferentes gráficos indicados por el usuario:
  if (plots$boxplot){   #el boxplot muestra que el rango de counts y las posiciones centrales son similares en muestras diferentes después de la normalizacion
    logCPM<-cpm(dgl, log=TRUE) #log-transformed counts per million
    boxplot(logCPM, 
            xaxt="n",
            main=paste(plot_prefix, "- Normalized logCPM"),
            xlab="Samples",
            ylab="Log (normalized counts)")
  }
  
  if (plots$mds) {
    plotMDS(dgl, col=as.numeric(group), pch=16,
            main=paste(plot_prefix, "- MDS Plot"))
    legend("topright", legend=levels(group),
           col=1:length(levels(group)),
           pch=16)
  }
  
  if (plots$pca) {
    logCPM<-cpm(dgl, log=TRUE)
    
    pca<-prcomp(t(logCPM), scale.=TRUE)
    
    pca_df<-data.frame(
      PC1=pca$x[,1],
      PC2=pca$x[,2],
      group=group
    )
    
    #Varianza explicada:
    var_expl<- round(100*pca$sdev^2/sum(pca$sdev^2),1) 
    
    
    p<-ggplot(pca_df, aes(PC1,PC2,color=group))+
      geom_point(size=2)+
      labs(title=paste(plot_prefix,"- PCA"),
           x=paste0("PC1 (", var_expl[1],"%)"),
           y=paste0("PC2 (", var_expl[2],"%)"))+
      theme_minimal()
    print(p)
  }
  
  return(list(
    dgl=dgl,
    group=group
  ))
}

analysis_edgeR<-function(dgl,
                         group,
                         lfc_threshold=1, 
                         fdr_threshold=0.05, 
                         plots=list(bcv=TRUE, ma=TRUE, volcano=TRUE, hm=TRUE),
                         base_folder="data",
                         plot_prefix="Edge R"){ 
  #Diseño:
  design<-model.matrix(~group) #convierte el modelo en matriz de diseño
  
  #Estimamos dispersión:
  dgl<-estimateDisp(dgl,design) #se calcula la dispersión, que mide la variabilidad biologica
  
  #Ajustamos modelo: 
  fit<-glmQLFit(dgl,design) #se ajusta al modelo QL Quasi-likelihood
  
  #Test de umbral de fold change:  
  treat<-glmTreat(fit,coef=2,lfc=lfc_threshold) 
  
  #Resultados completos:
  all_results<-topTags(treat, n=nrow(dgl), sort.by="none")$table #devuelve table con logFC, logCPM, F, PValue, FDR. n=nrow(gdl) significa que se hace para todos los genes
  
  if(!"id_cruce" %in% colnames(all_results)){
    all_results$id_cruce <- sub("\\..*","", rownames(all_results))
  } 
  
  
  #Resultados significativos (FDR menor que el umbral establecido):
  significant<-subset(all_results, FDR<fdr_threshold)
  
  #Separamos los genes en up y down-regulated:
  up<-subset(significant, logFC>0) #Sobreexpresado en tumor
  down<-subset(significant, logFC<0) #Infraexpresado en tumor
  
  cat("Número de DEGs (EdgeR):", nrow(significant), " (Up:", nrow(up), ", Down:", nrow(down), ")\n")
  
  #Generamos las gráficas indicadas por el usuario:
  if(plots$bcv) {
    plotBCV(dgl, main=paste(plot_prefix,"- BCV Plot"))
  }
  
  if (plots$ma){
    plotSmear(treat, de.tags=rownames(significant),
              panel.first=NULL, main=paste(plot_prefix, "- MA Plot at FDR <", fdr_threshold))
  }
  
  if (plots$volcano){
    
    volcano_df<-all_results
    volcano_df$DEG<-"Not significant"
    volcano_df$DEG[volcano_df$logFC > lfc_threshold &
                     volcano_df$FDR<fdr_threshold]<-"Up"
    
    volcano_df$DEG[volcano_df$logFC < -lfc_threshold &
                     volcano_df$FDR<fdr_threshold]<-"Down"
    top_genes<-volcano_df %>%
      dplyr::filter(DEG != "Not significant") %>%
      dplyr::arrange(FDR) %>%
      head(10)
    
    n_deg_volcano<-sum(volcano_df$DEG!="Not significant")
    
    p<-ggplot(volcano_df, aes(x=logFC, y=-log10(FDR), color=DEG))+
      geom_point(alpha=0.7,size=1.2)+
      scale_color_manual(values=c(
        "Not significant"="grey", "Up"="firebrick3", "Down"="dodgerblue3"))+
      geom_vline(xintercept=c(-lfc_threshold, lfc_threshold),
                 linetype="dashed", color="black")+
      geom_hline(yintercept=-log10(fdr_threshold),
                 linetype="dashed", color= "black")+
      geom_text_repel(
        data=top_genes,
        aes(label=gene_name),
        size=2.5,
        nudge_y = 0.5,
        direction="both",
        point.padding = 0.2,
        box.padding = 0.4,          
        min.segment.length = 0,
        segment.size = 0.2,
        max.overlaps = Inf,
        seed=42
      )+
      labs(title=paste(plot_prefix,"- Volcano Plot"),
           x="Log2 Fold Change",
           y="-Log10 FDR")+
      theme_minimal()
    
    print(p) #porque ggplot no muestra los gráficos si no lo indicas explícitamente
    
    volcano_deg<-subset(volcano_df, DEG!="Not significant")                         
  }
  
  if (plots$hm && nrow(significant) > 0) {
    #Seleccionamos los IDs de los top genes:
    top_genes_hm <- significant %>%
      dplyr::arrange(FDR) %>% 
      head(50) %>%
      dplyr::pull(id_cruce)
    
    #Extraemos los counts (para edgeR usamos logCPM):
    heatmap_data <- cpm(dgl, log=TRUE)
    rownames(heatmap_data) <- sub("\\..*", "", rownames(heatmap_data))
    heatmap_data <- heatmap_data[rownames(heatmap_data) %in% top_genes_hm, ]
    
    #Ordenamos las muestras por grupo:
    orden_muestras <- order(group) 
    heatmap_data <- heatmap_data[, orden_muestras] 
    
    #Preparamos la anotación:
    annotation_col <- data.frame(Group = group[orden_muestras]) 
    rownames(annotation_col) <- colnames(heatmap_data)
    
    #Dibujamos el heatmap:
    pheatmap::pheatmap(
      heatmap_data, 
      scale = "row", 
      annotation_col = annotation_col,
      show_rownames = (nrow(heatmap_data) < 50), #Mostrar nombres si son pocos
      show_colnames = FALSE,
      cluster_cols = FALSE,           #Mantiene el orden Normal vs Tumor
      cluster_rows = TRUE,            
      clustering_distance_rows = "correlation", 
      main = paste(plot_prefix, "- Heatmap")
    )
  }
  
  #Creamos carpeta de resultados:
  results_folder <- file.path(base_folder, "analysis", "results")
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive=TRUE)
  }
  
  write.table(significant,
              file=file.path(results_folder,
                             paste0(plot_prefix,"_DEGs_FDR_",fdr_threshold,".txt")),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
  
  write.table(all_results,
              file=file.path(results_folder,
                             paste0(plot_prefix,"_all_genes.txt")),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
  
  #Guardamos los genes significativos que codifiquen para proteína:
  significant_pc <- significant[significant$gene_type == "protein_coding" & !is.na(significant$gene_type), ]
  up_pc <- significant_pc[significant_pc$logFC > 0, ]
  down_pc <- significant_pc[significant_pc$logFC < 0, ]
  
  return(list(
    all_genes=all_results,
    significant=significant,
    up_genes=up,
    down_genes=down,
    significant_pc=significant_pc,
    up_genes_pc=up_pc,
    down_genes_pc=down_pc
  ))
}

complete_edgeR_analysis<-function(expr_data,
                                  group_variable="sample_type",
                                  reference_level="Normal",
                                  lfc_threshold=2, 
                                  fdr_threshold=0.01, 
                                  preprocess_plots=list(boxplot=TRUE, mds=TRUE, pca=TRUE), 
                                  analysis_plots=list(bcv=TRUE, ma=TRUE, volcano=TRUE, hm=TRUE),
                                  plot_prefix="Edge R"){
  #Preprocesamiento:
  pre<-preprocess_edgeR(
    expr_data=expr_data,
    group_variable=group_variable,
    reference_level=reference_level,
    plots=preprocess_plots,
    plot_prefix = plot_prefix
  )
  
  #Análisis diferencial:
  results<-analysis_edgeR(
    dgl=pre$dgl,
    group=pre$group,
    lfc_threshold=lfc_threshold,
    fdr_threshold=fdr_threshold,
    plots=analysis_plots,
    plot_prefix = plot_prefix
  )
  return(results)
}

#Función de análisis análoga a la de EdgeR:
analysis_voom<-function(expr_data,
                        group_variable="sample_type",
                        reference_level="Normal",
                        lfc_threshold=2,
                        fdr_threshold=0.01,
                        plots=list(volcano=TRUE, hm=TRUE), 
                        base_folder="data",
                        plot_prefix="Voom"){
  
  info_genes_raw <- as.data.frame(rowData(expr_data))
  info_genes_raw$id_cruce <- sub("\\..*","", rownames(info_genes_raw))
  
  #GENE_TYPE:
  info_genes <- info_genes_raw[, c("id_cruce", "gene_name", "gene_type")]

  counts<-assay(expr_data)
  counts<-as.matrix(counts) #aseguramos que sea matriz
  mode(counts)<-"numeric" #forzamos a numerico
  counts[is.na(counts)]<-0 #sustituimos NAs por 0
  
  
  group<-factor(colData(expr_data)[[group_variable]])
  group<-relevel(group, ref=reference_level)
  
  #al menos 10 counts en el 10% de las muestras: (filtrado)
  umbral_muestras<-ceiling(0.10*ncol(counts))
  keep<-rowSums(counts>=10)>=umbral_muestras
  
  cat("Genes totales:", nrow(counts), "\n")
  cat("Genes eliminados:", sum(!keep), "\n")
  cat("Genes retenidos:", sum(keep), "\n")
  
  counts<-counts[keep,]
  
  info_genes_filtered <- info_genes[rownames(counts), ]
  
  dge <- DGEList(counts=counts, group=group, genes=info_genes_filtered)
  dge<-calcNormFactors(dge, method="TMM")
  
  design<-model.matrix(~group)
  colnames(design)[2]<- paste0(levels(group)[2], "_vs_", levels(group)[1])
  coef_name<-colnames(design)[2]
  
  v<-voom(dge,design, plot=FALSE)
  fit<-lmFit(v, design)
  fit<-eBayes(fit)
  
  results <- topTable(fit, coef=2, number=Inf, adjust.method="fdr")
  
  results$DEG <- "Not significant"
  results$DEG[results$logFC > lfc_threshold & results$adj.P.Val < fdr_threshold] <- "Up"
  results$DEG[results$logFC < -lfc_threshold & results$adj.P.Val < fdr_threshold] <- "Down"
  
  #Para ponerlo en el mismo formato que la de edgeR:
  cols_principales <- c("id_cruce", "gene_name", "logFC", "AveExpr", "P.Value", "adj.P.Val")
  cols_restantes <- setdiff(colnames(results), c(cols_principales, "DEG")) 
  results <- results[, c(cols_principales, cols_restantes, "DEG")] 
  
  significant <- subset(results, DEG != "Not significant")
  up_genes <- subset(results, DEG == "Up")
  down_genes <- subset(results, DEG == "Down")
  
  cat("Número de DEGs (voom):", nrow(significant), " (Up:", nrow(up_genes), ", Down:", nrow(down_genes), ")\n")
  
  #Generamos las gráficas:
  if (plots$volcano){
    volcano_df<-results
    
    top_genes<-volcano_df %>%
      dplyr::filter(DEG != "Not significant") %>%
      dplyr::arrange(adj.P.Val) %>% 
      head(10)
    
    p<-ggplot(volcano_df, aes(x=logFC, y=-log10(adj.P.Val), color=DEG))+ 
      geom_point(alpha=0.6, size=1)+
      scale_color_manual(values=c("Not significant"="grey", "Up"="firebrick3", "Down"="dodgerblue3"))+
      geom_vline(xintercept=c(-lfc_threshold, lfc_threshold), linetype="dashed", color="black")+
      geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="black")+
      geom_text_repel(
        data=top_genes,
        aes(label=gene_name),
        size=2.5,
        nudge_y = 0.5,
        direction="both",
        point.padding = 0.2,
        box.padding = 0.4,          # Empuja hacia los bordes
        min.segment.length = 0,
        segment.size = 0.2,
        max.overlaps = Inf,
        seed=42
      )+
      labs(title=paste(plot_prefix,"- Volcano Plot"),
           x="Log2 Fold Change",
           y="-Log10 adj.P.Val")+ 
      theme_minimal()
    
    print(p) 
  }
  
  if (plots$hm && nrow(significant) > 0) {
    #Seleccionamos los IDs de los top genes:
    top_genes_hm <- significant %>%
      arrange(adj.P.Val) %>% # O FDR en edgeR
      head(50) %>%
      pull(id_cruce)
    
    #Extraer los counts (v$E para voom o logCPM para edgeR):
    heatmap_data <- v$E 
    rownames(heatmap_data) <- sub("\\..*", "", rownames(heatmap_data))
    heatmap_data <- heatmap_data[rownames(heatmap_data) %in% top_genes_hm, ]
    
    #Ordenamos las muestras por grupo:
    # Creamos un índice para ordenar: primero Normal, luego Tumor (según tus niveles)
    orden_muestras <- order(group)
    heatmap_data <- heatmap_data[, orden_muestras] 
    
    #Preparamos la anotación con el nuevo orden:
    annotation_col <- data.frame(Group = group[orden_muestras]) # <<-- CAMBIO
    rownames(annotation_col) <- colnames(heatmap_data)
    
    #Generamos el heatmap:
    pheatmap(
      heatmap_data, 
      scale = "row", 
      annotation_col = annotation_col,
      show_rownames = FALSE, 
      show_colnames = FALSE,
      cluster_cols = FALSE,          
      cluster_rows = TRUE,            # Los genes sí los dejamos agrupados por parecido
      clustering_distance_rows = "correlation", 
      main = paste(plot_prefix, "- Heatmap")
    )
  }
  
  results_folder <- file.path(base_folder, "analysis", "results") 
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive=TRUE)
  }
  
  #Guardamos DEGs significativos:
  write.table(significant, 
              file=file.path(results_folder, paste0(plot_prefix, "_DEGs_FDR_",fdr_threshold,"_LFC_",lfc_threshold,".txt")),
              quote=FALSE, sep="\t", row.names=FALSE) 
  
  # Guardamos tabla completa:
  write.table(results, 
              file=file.path(results_folder, paste0(plot_prefix, "_all_genes.txt")),
              quote=FALSE, sep="\t", row.names=FALSE)
  
  #Guardamos los significativos que codifiquen para proteína:
  significant_pc <- significant[significant$gene_type == "protein_coding" & !is.na(significant$gene_type), ]
  up_pc <- significant_pc[significant_pc$DEG == "Up", ]
  down_pc <- significant_pc[significant_pc$DEG == "Down", ]
  
  return(list(
    all_genes=results,
    significant=significant,
    up_genes=up_genes,     
    down_genes=down_genes,
    significant_pc = significant_pc,
    up_genes_pc = up_pc,
    down_genes_pc = down_pc
  ))
}

#Análisis EdgeR LUAD:
res_LUAD_edgeR<-complete_edgeR_analysis(expr_LUAD_def,
                                        lfc_threshold = 2,
                                        fdr_threshold = 0.01,
                                        preprocess_plots = list(boxplot=FALSE, mds=FALSE, pca=FALSE),
                                        analysis_plots = list(bcv=FALSE, ma=FALSE, volcano=FALSE, hm=FALSE),
                                        plot_prefix="LUAD - Edge R")
saveRDS(res_LUAD_edgeR,"data/res_LUAD_edgeR_pcg.rds")

#Análisis voom LUAD:
res_LUAD_voom<-analysis_voom(expr_LUAD_def,
                             lfc_threshold = 2,
                             fdr_threshold = 0.01,
                             plots = list(volcano=FALSE, hm=FALSE),
                             plot_prefix="LUAD - Voom")
saveRDS(res_LUAD_voom,"data/res_LUAD_voom_pcg.rds")

common_deg_LUAD <- intersect(
  res_LUAD_edgeR$significant$id_cruce,
  res_LUAD_voom$significant$id_cruce
)
saveRDS(common_deg_LUAD,"data/common_deg_LUAD.rds")

#Análisis EdgeR LUSC:
res_LUSC_edgeR<-complete_edgeR_analysis(expr_LUSC_def, 
                                        lfc_threshold = 2,
                                        fdr_threshold = 0.01,
                                        preprocess_plots = list(boxplot=FALSE, mds=FALSE, pca=FALSE),
                                        analysis_plots = list(bcv=FALSE, ma=FALSE, volcano=FALSE, hm=FALSE),
                                        plot_prefix="LUSC - Edge R")
saveRDS(res_LUSC_edgeR,"data/res_LUSC_edgeR_pcg.rds")

#Análisis Voom LUSC:
res_LUSC_voom<-analysis_voom(expr_LUSC_def, 
                             lfc_threshold = 2,
                             fdr_threshold = 0.01,
                             plots = list(volcano=FALSE, hm=FALSE),
                             plot_prefix="LUSC - Voom")
saveRDS(res_LUSC_voom,"data/res_LUSC_voom_pcg.rds")

common_deg_LUSC <- intersect(
  res_LUSC_edgeR$significant$id_cruce,
  res_LUSC_voom$significant$id_cruce
)
saveRDS(common_deg_LUSC,"data/common_deg_LUSC.rds")

common_lung_deg<-intersect(common_deg_LUAD, common_deg_LUSC)
saveRDS(common_lung_deg,"data/common_lung_deg.rds")


#Preparación de los Inputs para las herramientas de reposicionamiento de fármacos:
#Función que genera los inputs para ShinyDeepDR:
#Para shinydeepDR tienes que hacer como si solo fuera una muestra
#Serviría muy bien para medicina personalizada
#Extrae la mediana de los TPMs de las muestras tumorales:

inputs_shinyDeepDR<-function(expr_data, file_name, folder="data/analysis/results") {
  if(!dir.exists(folder)){
    dir.create(folder, recursive=TRUE)
  }  
  es_tumor <- colData(expr_data)$sample_type == "Tumor"
  expr_tumor <- expr_data[, es_tumor]
  
  #extraer la matriz de TPM:
  if("tpm_unstrand" %in% assayNames(expr_tumor)) {
    tpm_matrix <- assay(expr_tumor, "tpm_unstrand")
  } else {
    stop("No se ha encontrado el ensayo 'tpm_unstrand'.")
  }
  
  #calculamos la mediana por gen:
  mediana_tpm <- apply(tpm_matrix, 1, median)
  
  #nombres de los genes:
  info_genes <- as.data.frame(rowData(expr_tumor))
  
  df_shiny <- data.frame(
    Gene = info_genes$gene_name,
    Promedio_Tumoral = mediana_tpm
  )
  
  #Eliminamos NAs:
  df_shiny <- df_shiny[!is.na(df_shiny$Gene), ]
  
  #Eliminamos duplicados:
  #Agrupamos por símbolo genético y sumamos los TPMs de las variantes/isoformas (en el caso de que varios ensembl id apunten a un mismo genesymbol)
  df_shiny <- df_shiny %>%
    group_by(Gene) %>%
    summarise(Promedio_Tumoral = sum(Promedio_Tumoral)) %>%
    ungroup() %>%
    as.data.frame()
  
  #Guardamos:
  ruta_salida <- file.path(folder, file_name)
  write.table(df_shiny, file = ruta_salida, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Archivo shinyDeepDR guardado en:", ruta_salida, "\n")
  cat("Tumores promediados:", sum(es_tumor), "\n")
  
}

inputs_shinyDeepDR(
  expr_data = expr_LUAD_def,
  file_name="shinyDeepDR_LUAD.txt"
)

inputs_shinyDeepDR(
  expr_data = expr_LUSC_def,
  file_name="shinyDeepDR_LUSC.txt"
)

#Inputs para ShinyDeepDR de MUTACIONES:
maf_shinyDeepDR_LUAD <- mut_LUAD_raw %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode) %>%
  mutate(Tumor_Sample_Barcode = "LUAD_Global") %>% # Forzamos una única muestra
  distinct() # Eliminamos duplicados si el mismo gen muta igual en varios pacientes

# Guardamos el MAF "unificado":
write.table(maf_shinyDeepDR_LUAD, 
            file = "data/LUAD_mutations_ShinyDeepDR.maf", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

maf_shinyDeepDR_LUSC <- mut_LUSC_raw %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode) %>%
  mutate(Tumor_Sample_Barcode = "LUSC_Global") %>% # Forzamos una única muestra
  distinct() # Eliminamos duplicados si el mismo gen muta igual en varios pacientes

write.table(maf_shinyDeepDR_LUSC, 
            file = "data/LUSC_mutations_ShinyDeepDR.maf", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

#Preparamos los inputs para CDRPipe:
inputs_cdrpipe<-function(df_degs, file_name, folder="data/analysis/results"){
  
  #creamos directorio si no existe:
  if(!dir.exists(folder)){
    dir.create(folder, recursive=TRUE)
  }
  
  #detectar si es Voom o edgeR para elegir la columna del FDR:
  if ("adj.P.Val" %in% colnames(df_degs)) {
    col_pval <- "adj.P.Val"   # Es Voom
  } else if ("FDR" %in% colnames(df_degs)) {
    col_pval <- "FDR"         # Es edgeR
  } else {
    stop("No se encuentra adj.P.Val ni FDR en los datos.")
  }
  
  #Seleccionamos y renombramos las columnas según lo que pide CDRpipe:
  df_cdr <- df_degs %>%
    dplyr::select(
      SYMBOL = gene_name,
      log2FC_1 = logFC,
      p_val_adj = all_of(col_pval)
    )
  
  #Limpieza de NAs
  df_cdr <- df_cdr %>% filter(!is.na(SYMBOL) & SYMBOL != "")
  
  #Manejo de duplicados(Nos quedamos con la variante que tenga el p-valor más significativo (el menor)):  
  df_cdr <- df_cdr %>%
    group_by(SYMBOL) %>%
    slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  
  #Guardamos en formato CSV:
  ruta_salida <- file.path(folder, file_name)
  write.csv(df_cdr, file = ruta_salida, row.names = FALSE, quote = FALSE)
  
  cat("Archivo para CDRpipe guardado en:", ruta_salida, "\n")
  cat("Total de genes exportados:", nrow(df_cdr), "\n")
  
  return(df_cdr)
}

degs_voom_luad <- read.delim("C:/Users/anaal/OneDrive - UNIVERSIDAD DE GRANADA/TFG/TFG-AAV/data/analysis/results/FDR0.01LFC2/LUAD - Voom_DEGs_FDR_0.01_LFC_2.txt", sep = "\t", header = TRUE)

input_limpio_luad<-inputs_cdrpipe(
  df_degs = degs_voom_luad, 
  file_name = "CDRpipe_LUAD_Voom.csv"
)

degs_voom_lusc <- read.delim("C:/Users/anaal/OneDrive - UNIVERSIDAD DE GRANADA/TFG/TFG-AAV/data/analysis/results/FDR0.01LFC2/LUSC - Voom_DEGs_FDR_0.01_LFC_2.txt", sep = "\t", header = TRUE)

input_limpio_lusc<-inputs_cdrpipe(
  df_degs = degs_voom_lusc, 
  file_name = "CDRpipe_LUSC_Voom.csv"
)

#Generación de inputs para iLINCS basándonos en el objeto de CDRpipe:
inputs_ilincs <- function(df_cdr, file_name, folder = "data/analysis/results") {
  
  #Creamos directorio si no existe:
  if(!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }

  ruta_salida_ilincs <- file.path(folder, file_name)
  
  #Guardamos en formato TXT delimitado por tabuladores (\t)
  write.table(
    df_cdr, 
    file = ruta_salida_ilincs, 
    sep = "\t", 
    row.names = FALSE, 
    col.names = FALSE, #sin nombres de columna
    quote = FALSE
  )
  cat("Archivo para iLINCS guardado en:", ruta_salida_ilincs, "\n")
}

inputs_ilincs(
  df_cdr = input_limpio_luad, 
  file_name = "iLINCS_LUAD_voom.txt"
)

inputs_ilincs(
  df_cdr = input_limpio_lusc, 
  file_name = "iLINCS_LUSC_voom.txt"
)

#Creamos los objetos GCT para CMap (para su herramienta de Query):

#LUAD:
#Ordenamos UP por significancia (adj.P.Val de menor a mayor)
#El p-valor más pequeño es el más significativo.
voom_up_LUAD <- res_LUAD_voom$up_genes %>%
  arrange(adj.P.Val) %>% 
  head(150) %>%
  pull(gene_name)

#Ordenamos DOWN por significancia (adj.P.Val de menor a mayor)
voom_down_LUAD <- res_LUAD_voom$down_genes %>%
  arrange(adj.P.Val) %>% 
  head(150) %>%
  pull(gene_name)

#Guardamos los genes en archivos .txt para subir a CLUE.io
write.table(voom_up_LUAD, file = "data/CMap_LUAD_voom_up_top150.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(voom_down_LUAD, file = "data/CMap_LUAD_voom_down_top150.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Listas de inputs para CMap generadas")

#LUSC:
voom_up_LUSC <- res_LUSC_voom$up_genes %>%
  arrange(adj.P.Val) %>% 
  head(150) %>%
  pull(gene_name)

voom_down_LUSC <- res_LUSC_voom$down_genes %>%
  arrange(adj.P.Val) %>% 
  head(150) %>%
  pull(gene_name)

write.table(voom_up_LUSC, file = "data/CMap_LUSC_voom_up_top150.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(voom_down_LUSC, file = "data/CMap_LUSC_voom_down_top150.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

#------------
#TRATAMIENTO DE LOS RESULTADOS DE CLUE:
analizar_resultados_clue <- function(ruta_gct, nombre_archivo_csv) {
  
  gct_data <- cmapR::parse.gctx(ruta_gct)
  
  #Extraemos metadatos y matriz:
  metadatos_farmacos <- as.data.frame(gct_data@rdesc)
  scores_ncs <- as.data.frame(gct_data@mat)
  
  #Creamos la tabla final combinada:
  tabla_resultados <- data.frame(
    Nombre_Farmaco = metadatos_farmacos$pert_iname,
    Mecanismo_Accion = metadatos_farmacos$moa,
    Target = metadatos_farmacos$target_name,
    NCS = scores_ncs[,1] #Tomamos la primera columna de resultados
  )
  
  #Filtrar los TOP HITS (Los más negativos son los mejores)
  mis_candidatos <- tabla_resultados %>%
    filter(NCS < 0) %>% 
    arrange(NCS)
  
  #Guardamos los resultados en la carpeta de trabajo:
  write.csv(mis_candidatos, nombre_archivo_csv, row.names = FALSE)
  
  #Devolvemos el objeto de los fármacos candidatos:
  return(mis_candidatos)
}


#Ejecutamos la función para LUAD:
res_LUAD_CMap <- analizar_clue_resultados(
  ruta_gct = "C:/Users/anaal/OneDrive - UNIVERSIDAD DE GRANADA/TFG/TFG-AAV/data/analysis/results/Resultados reposicionamiento/Clue Query (CMap)/query_LUAD/my_analysis.sig_queryl1k_tool.69f74eb58ed24f0014280d76/ncs.gct",
  nombre_archivo_csv = "Resultados_Reposicionamiento_LUAD.csv"
)
saveRDS(mis_candidatos,"data/candidatos_LUAD.rds")

res_LUSC_CMap <- analizar_clue_resultados(
  ruta_gct = "data/analysis/results/Resultados reposicionamiento/Clue Query (CMap)/query_LUSC/my_analysis.sig_queryl1k_tool.69f7596e8ed24f0014280d7f/ncs.gct",
  nombre_archivo_csv = "Resultados_Reposicionamiento_LUAD.csv"
)
saveRDS(mis_candidatos,"data/candidatos_LUSC.rds")


#FUNCIÓN para obtener fármacos comunes a los 4 métodos:

obtener_farmacos_consenso <- function(ruta_cdr, ruta_cmap, ruta_ilincs, ruta_shiny, archivo_salida) {
  
  #Carga de resultados de reposicionamiento:
  df_cdr_LUAD <- read.csv(ruta_cdr)
  df_cmap_LUAD <- read.csv(ruta_cmap)
  df_ilincs_LUAD <- read.delim(ruta_ilincs, sep = "\t")
  df_shiny_LUAD <- read.csv(ruta_shiny)
  
  #Normalización y firltrado inicial:
  #CDRpipe: Columna 'name'
  df_cdr_clean <- df_cdr_LUAD %>%
    mutate(drug_lower = tolower(trimws(name))) %>%
    filter(cmap_score < 0) # Filtramos por reversión
  
  #CMap (Clue): Columna 'Nombre_Farmaco'
  df_cmap_clean <- df_cmap_LUAD %>%
    mutate(drug_lower = tolower(trimws(Nombre_Farmaco))) %>%
    filter(NCS < 0) # Filtramos por conectividad negativa
  
  #iLINCS: Columna 'Perturbagen'
  df_ilincs_clean <- df_ilincs_LUAD %>%
    mutate(drug_lower = tolower(trimws(Perturbagen))) %>%
    filter(Concordance < 0) # Solo los que revierten (negativos)
  
  #ShinyDeepDR: Columna 'Drug.Name' 
  df_shiny_clean <- df_shiny_LUAD %>%
    mutate(drug_lower = tolower(trimws(Drug.name))) %>%
    filter(Predicted.IC50..log.uM. < 0) 
  
  #Preparar datos con su métrica de significancia (limpiando nombres antes)
  res_cdr <- df_cdr_clean %>% 
    group_by(drug_lower) %>% 
    summarise(Score_CDR = min(cmap_score), .groups = 'drop')
  
  res_cmap <- df_cmap_clean %>% 
    group_by(drug_lower) %>% 
    summarise(Score_CMap = min(NCS), .groups = 'drop')
  
  res_ilincs <- df_ilincs_clean %>% 
    group_by(drug_lower) %>% 
    summarise(Score_iLINCS = min(Concordance), pVal_iLINCS = min(pValue), .groups = 'drop')
  
  res_shiny <- df_shiny_clean %>% 
    group_by(drug_lower) %>% 
    summarise(IC50_Shiny = min(Predicted.IC50..log.uM.), .groups = 'drop')
  
  lista_cdr    <- unique(df_cdr_clean$drug_lower)
  lista_cmap   <- unique(df_cmap_clean$drug_lower)
  lista_ilincs <- unique(df_ilincs_clean$drug_lower)
  lista_shiny  <- unique(df_shiny_clean$drug_lower)
  
  presencia_cdr    <- data.frame(Farmaco = lista_cdr,    CDRpipe = 1)
  presencia_cmap   <- data.frame(Farmaco = lista_cmap,   CMap = 1)
  presencia_ilincs <- data.frame(Farmaco = lista_ilincs, iLINCS = 1)
  presencia_shiny  <- data.frame(Farmaco = lista_shiny,  ShinyDeepDR = 1)
  
  #Intersección y significancia:
  tabla_final <- presencia_cdr %>% 
    full_join(presencia_cmap, by = "Farmaco") %>%
    full_join(presencia_ilincs, by = "Farmaco") %>%
    full_join(presencia_shiny, by = "Farmaco") %>%
    # Añadir los valores numéricos
    left_join(res_cdr, by = c("Farmaco" = "drug_lower")) %>%
    left_join(res_cmap, by = c("Farmaco" = "drug_lower")) %>%
    left_join(res_ilincs, by = c("Farmaco" = "drug_lower")) %>%
    left_join(res_shiny, by = c("Farmaco" = "drug_lower"))
  
  #Sustituir NAs en las columnas de presencia (del 2 al 5)
  tabla_final[, 2:5][is.na(tabla_final[, 2:5])] <- 0
  
  
  #Cálculo de totales y orden final:
  tabla_final <- tabla_final %>%
    mutate(Total_Metodos = CDRpipe + CMap + iLINCS + ShinyDeepDR) %>%
    # Prioridad: 1º Más métodos, 2º Mejor reversión en iLINCS, 3º Mejor IC50 en Shiny
    arrange(desc(Total_Metodos), Score_iLINCS, IC50_Shiny)
  
  #Mostrar los resultados top
  print(head(tabla_final, 20))
  
  write.csv(tabla_final, archivo_salida, row.names = FALSE)
  
  cat("\nAnálisis finalizado. Archivo guardado como:", archivo_salida, "\n")
  
  return(tabla_final)
}

#Ejecutamos función de fármacos consenso para LUAD:
farmacos_consenso_LUAD<- obtener_farmacos_consenso(
  ruta_cdr = "data/analysis/results/Resultados reposicionamiento/CDRPipe/LUAD_1_005.csv",
  ruta_cmap = "data/analysis/results/Resultados reposicionamiento/Clue Query (CMap)/Resultados_LUAD/Resultados_Reposicionamiento_LUAD.csv",
  ruta_ilincs = "data/analysis/results/Resultados reposicionamiento/iLINCS/iLINCS_Connectivity_Results_LUAD.xls",
  ruta_shiny = "data/analysis/results/Resultados reposicionamiento/ShinyDeepDR/shinyDeepDR_LUAD.csv",
  archivo_salida = "Farmacos_Consenso_LUAD.csv"
  )

#LUSC:
farmacos_consenso_LUSC<-obtener_farmacos_consenso(
  ruta_cdr = "data/analysis/results/Resultados reposicionamiento/CDRPipe/LUSC_CDRPipe_1_005.csv",
  ruta_cmap = "data/analysis/results/Resultados reposicionamiento/Clue Query (CMap)/Resultados_LUSC/Resultados_Reposicionamiento_LUSC.csv",
  ruta_ilincs = "data/analysis/results/Resultados reposicionamiento/iLINCS/iLINCS_resultados_LUSC.xls",
  ruta_shiny = "data/analysis/results/Resultados reposicionamiento/ShinyDeepDR/shinyDeepDR_LUSC.csv",
  archivo_salida = "Farmacos_Consenso_LUSC.csv"
)

##NUEVO APPROACH: generar los inputs para CDRPipe, iLINCS y Clue
# a partir de los genes COMUNES a EdgeR y Voom
#a partir de los datos de expresión con FDR 0.01 y LFC 2 (los más estrictos)


#Ya habíamos sacado un vector con los comunes, tomamos los genes de ese vector como criterio de inclusión:
df_comun_luad <- res_LUAD_voom$significant %>%
  filter(id_cruce %in% common_deg_LUAD)

df_comun_lusc <- res_LUSC_voom$significant %>%
  filter(id_cruce %in% common_deg_LUSC)

#CDRPipe LUAD y LUSC:
input_cdr_luad_comun <- inputs_cdrpipe(
  df_degs = df_comun_luad, 
  file_name = "CDRpipe_LUAD_Common_Voom_edgeR.csv"
)

input_cdr_lusc_comun <- inputs_cdrpipe(
  df_degs = df_comun_lusc, 
  file_name = "CDRpipe_LUSC_Common_Voom_edgeR.csv"
)

#iLINCS LUAD y LUSC:
inputs_ilincs(
  df_cdr = input_cdr_luad_comun, 
  file_name = "iLINCS_LUAD_Common_Voom_edgeR.txt"
)

inputs_ilincs(
  df_cdr = input_cdr_lusc_comun, 
  file_name = "iLINCS_LUSC_Common_Voom_edgeR.txt"
)

preparar_cmap_consenso <- function(df_comun, prefix) {
  
  up_genes <- df_comun %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(150) %>% pull(gene_name)
  down_genes <- df_comun %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% head(150) %>% pull(gene_name)
  
  write.table(up_genes, file = paste0("data/CMap_", prefix, "_common_up.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(down_genes, file = paste0("data/CMap_", prefix, "_common_down.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("Archivos CMap generados para", prefix, "con", length(up_genes), "up y", length(down_genes), "down genes.\n")
}

preparar_cmap_consenso(df_comun_luad, "LUAD")
preparar_cmap_consenso(df_comun_lusc, "LUSC")