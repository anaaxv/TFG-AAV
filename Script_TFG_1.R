#Instalamos las librerías necesarias si no las tenemos ya: 
if (!require("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")

if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")  #REVISAR

if (!require("edgeR", quietly = TRUE)) #REVISAR
  BiocManager::install("edgeR")

if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

if (!require("DT", quietly = TRUE))
  install.packages("DT")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!require("ggrepel", quietly = TRUE))
  install.packages("ggrepel")

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

#Guardo los objetos que he obtenido, para poder cargarlos en cualquier momento, sin necesidad de volver a correr el código:
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
#cuando tcgabiolinks hace gdcprepare() pone el barcode de la muestra como nombre de la columna
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
  
  #actualizamos barcodes porque hemos quitado columnas de expr
  barcodes<-colnames(expr)
  
  #Calcular tamaño de librería:
  library_size<-colSums(assay(expr))
  colData(expr)$library_size<-library_size
  
  #Dataframe auxiliar:
  coldata_df<-as.data.frame(colData(expr))
  coldata_df$barcode<-barcodes
  
  #Seleccionamos una muestra por paciente y por tipo de tumor (la de mayor library size, en este caso):
  #De cada paciente se guardarán como máximo los datos de dos muestras: una de tumor y una normal
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
  
  #Si un paciente tiene TUMOR y NORMAL, quedarnos solo con NORMAL:
  
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

#Aplicamos la función que acabamos de crear a LUAD:
res_LUAD<-procesar_expr(
  expr=expr_LUAD_raw,
  clinical_patients=clinical_LUAD_BCRtab$clinical_patient_luad$bcr_patient_barcode
  )

expr_LUAD_def<-res_LUAD$expr_def
common_patients_LUAD<-res_LUAD$common_patients

saveRDS(expr_LUAD_def,"data/expr_LUAD_def.rds")
saveRDS(common_patients_LUAD, "data/common_patients_LUAD.rds")

#Aplicamos la función a LUSC:
res_LUSC<-procesar_expr(
  expr=expr_LUSC_raw,
  clinical_patients=clinical_LUSC_BCRtab$clinical_patient_lusc$bcr_patient_barcode
)

expr_LUSC_def<-res_LUSC$expr_def
common_patients_LUSC<-res_LUSC$common_patients

saveRDS(expr_LUSC_def,"data/expr_LUSC_def.rds")
saveRDS(common_patients_LUSC, "data/common_patients_LUSC.rds")

#table(colData(expr_LUAD_def)$sample_type)
#table(colData(expr_LUAD_def)$sample_type_code)
#table(colData(expr_LUSC_def)$sample_type)
#table(colData(expr_LUSC_def)$sample_type_code)

preprocess_edgeR<-function(expr_data,                    #Objeto Summarized Experiment
                           group_variable="sample_type", #Nombre de la columna en colData
                           reference_level="Normal",
                           plot_prefix="Edge R",
                           plots=list(boxplot=TRUE, mds=TRUE, pca=TRUE)){
  #Extraer matriz de counts:
  #si hay NA que se quede con 0 y se asegura de que todo sea numerico 
  counts<-assay(expr_data)  #Assay es una función del paquete SummarizedExperiment. Devuelve la matriz de expresión (filas=genes, columnas=muestras, valores=counts)
  counts<-as.matrix(counts) #para asegurar que sea matriz estandar
  mode(counts)<-"numeric"
  counts[is.na(counts)]<-0  #busca los NA y los convierte en 0
  
  #counts[]<-lapply(as.data.frame(counts), function(x) as.numeric(as.character(x))) #convierte a data frame, recorre columna por columna, fuerza que sea numerico, y lo vuelve a meter en la matriz

  
  #Definir grupo:
  group<-factor(colData(expr_data)[[group_variable]]) #la variable "tumor" o "normal" se convierte en variable categórica al usar factor. Con el doble corchete se extrae de verdad lo que hay en esa columna
  group<-relevel(group,ref=reference_level) #define el grupo de referencia en el modelo, en este caso, las muestras "normales"
  
  #Crear DGElist:
  dgl<-DGEList(counts=counts,group=group) #creamos la dgelist, estructura base para edgeR
  
  info_genes<-as.data.frame(rowData(expr_data))
  info_genes$id_cruce<-sub("\\..*", "", rownames(info_genes))
  dgl$genes<-info_genes[, c("id_cruce", "gene_name")]
  
  #Filtrado de genes poco expresados:
  keep<-filterByExpr(dgl, group=group)  #elimina genes con muy baja expresión, reduce ruido, mejora potencia estadística etc
  #umbral_muestras<-ceiling(0.10*ncol(counts))
  #keep<-rowSums(counts>=10)>=umbral_muestras
  
  #cat: como print pero mejor 
  cat("Genes totales:", nrow(counts), "\n")
  cat("Genes eliminados:", sum(!keep), "\n")
  cat("Genes retenidos:", sum(keep), "\n")
  
  dgl<-dgl[keep,,keep.lib.sizes=FALSE] #nos quedamos solo con los genes que hayan pasado el filtro, haces que se recalculen los tamaños de biblioteca
  rm(keep)
  
  #Normlización TMM:
  dgl<-calcNormFactors(dgl, method="TMM") #Normalización Trimmed Mean of M-Values, corrige diferencias de composición, de profundidas de secuenciación...
 
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
    
    var_expl<- round(100*pca$sdev^2/sum(pca$sdev^2),1) #varianza explicada
    
    
    p<-ggplot(pca_df, aes(PC1,PC2,color=group))+
      geom_point(size=3)+
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
                   lfc_threshold=1, #cambiar a 1.5???
                   fdr_threshold=0.05, #cambiar a 0.05???
                   plots=list(bcv=TRUE, ma=TRUE, volcano=TRUE, hm=TRUE),
                   base_folder="data",
                   plot_prefix="Edge R"){ 
  #Diseño:
  design<-model.matrix(~group) #convierte el modelo en matriz de diseño
  
  #Estimar dispersión:
  dgl<-estimateDisp(dgl,design) #se calcula la dispersión, que mide la variabilidad biologica
  
  #Ajustar modelo: 
  fit<-glmQLFit(dgl,design) #se ajusta al modelo QL Quasi-likelihood
  
  #Test de umbral de fold change:  Incluir también glmQLFTest??? para ver cualquier cambio, es decir, logFC distinto de 0
  treat<-glmTreat(fit,coef=2,lfc=lfc_threshold) #pregunta si logFC>1, o sea foldchange>2, mas exigente que si buscasemos si logFc distinto de 0
  
  #Resultados completos:
  all_results<-topTags(treat, n=nrow(dgl), sort.by="none")$table #devuelve table con logFC, logCPM, F, PValue, FDR. n=nrow(gdl) significa que se hace para todos los genes
  
  all_results$id_cruce<-sub("\\..*","", rownames(all_results))
  
  #if (!is.null(dgl$genes)){
  #  all_results<-dplyr::left_join(all_results, 
   #                               dgl$genes,
  #                                by="id_cruce")
  #}

  
  #Resultados significativos:
  significant<-subset(all_results, FDR<fdr_threshold)
  
  #Separar up y down:
  up<-subset(significant, logFC>0) #Sobreexpresado en tumor
  down<-subset(significant, logFC<0) #Infraexpresado en tumor
  
  
  if(plots$bcv) {
    plotBCV(dgl, main=paste(plot_prefix,"- BCV Plot"))
  }
  
  if (plots$ma){
    plotSmear(treat, de.tags=rownames(significant),
              panel.first=NULL, main=paste(plot_prefix, "- MA Plot at FDR <", fdr_threshold))
  }
  #REVISARRR
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
    cat("Número de DEGs para Volcano (|log2FC| >", lfc_threshold,
        "y FDR <", fdr_threshold, "):", n_deg_volcano, "\n")
    
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
        size=3
      )+
      labs(title=paste(plot_prefix,"- Volcano Plot"),
           x="Log2 Fold Change",
           y="-Log10 FDR")+
      theme_minimal()
  
  print(p) #pq en ggplot no salen los graficos automaticamente
        
  volcano_deg<-subset(volcano_df, DEG!="Not significant")                         
  }
  
  if (plots$hm) {
    top_genes<- significant %>%
      arrange(FDR) %>%
      head(50) %>%
      pull(id_cruce)
    
    heatmap_counts<-cpm(dgl, log=TRUE)
    rownames(heatmap_counts)<-sub("\\..*","",rownames(heatmap_counts))
    heatmap_counts<-heatmap_counts[rownames(heatmap_counts)  %in% top_genes,]
    annotation_col<-data.frame(Group=group)
    rownames(annotation_col)<-colnames(heatmap_counts)
    
    pheatmap(
      heatmap_counts,
      scale="row",
      annotation_col=annotation_col,
      show_rownames= FALSE,
      show_colnames = FALSE,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      main= paste(plot_prefix, "- Heatmap Top DEGs")
    )
  }
  

  
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
    
    if (plots$volcano){
      write.table(volcano_deg,
                  file=file.path(results_folder,
                                 paste0(plot_prefix,"_Volcano_DEGs_.txt")),
                  quote=FALSE,
                  sep="\t",
                  row.names=FALSE)
    }
    
return(list(
  all_genes=all_results,
  significant=significant,
  up_genes=up,
  down_genes=down
))
}

complete_edgeR_analysis<-function(expr_data,
                                  group_variable="sample_type",
                                  reference_level="Normal",
                                  lfc_threshold=1, #añadir tb para el bvolcano otra que sea para logfc 2
                                  fdr_threshold=0.05, #mejor 5%??
                                  preprocess_plots=list(boxplot=FALSE, mds=FALSE, pca=FALSE), #CAMBIARRR
                                  analysis_plots=list(bcv=FALSE, ma=FALSE, volcano=FALSE, hm=TRUE)){
  #Preprocesamiento:
  pre<-preprocess_edgeR(
    expr_data=expr_data,
    group_variable=group_variable,
    reference_level=reference_level,
    plots=preprocess_plots
  )
  
  #Análisis diferencial:
  results<-analysis_edgeR(
    dgl=pre$dgl,
    group=pre$group,
    lfc_threshold=lfc_threshold,
    fdr_threshold=fdr_threshold,
    plots=analysis_plots
  )
  return(results)
}

analysis_voom<-function(expr_data,
                        group_variable="sample_type",
                        reference_level="Normal",
                        lfc_threshold=1,
                        fdr_threshold=0.05,
                        plots=list(volcano=FALSE, hm=TRUE), #CAMBIARRR
                        plot_prefix="Voom"){
  
  info_genes <- as.data.frame(rowData(expr_data))
  info_genes$id_cruce <- sub("\\..*","", rownames(info_genes))
  info_genes <- info_genes[, c("id_cruce", "gene_name")]
  
  counts<-assay(expr_data)
  counts<-as.matrix(counts) #aseguramos que sea matriz
  mode(counts)<-"numeric" #forzamos a numerico
  counts[is.na(counts)]<-0 #sustituimos NAs por 0
  
  
  group<-factor(colData(expr_data)[[group_variable]])
  group<-relevel(group, ref=reference_level)
  
  umbral_muestras<-ceiling(0.10*ncol(counts))
  keep<-rowSums(counts>=10)>=umbral_muestras
  
  #cat: como print pero mejor 
  cat("Genes totales:", nrow(counts), "\n")
  cat("Genes eliminados:", sum(!keep), "\n")
  cat("Genes retenidos:", sum(keep), "\n")
  
  counts<-counts[keep,]
  
  dge<-DGEList(counts=counts, group=group)
  dge<-calcNormFactors(dge, method="TMM")
  
  design<-model.matrix(~group)
  colnames(design)[2]<- paste0(levels(group)[2], "_vs_", levels(group)[1])
  coef_name<-colnames(design)[2]
  
  v<-voom(dge,design, plot=FALSE)
  fit<-lmFit(v, design)
  fit<-eBayes(fit)
  
  results<-topTable(fit,
                    coef=coef_name,
                    number=Inf,
                    adjust.method="fdr")
  
  results$id_cruce<-sub("\\..*","", rownames(results))
  
  results<-dplyr::left_join(results, 
                 info_genes,
                 by="id_cruce")
  
  significant<-subset(results, adj.P.Val<fdr_threshold & abs(logFC)>lfc_threshold)
  
  cat("Número de DEGs (voom):", nrow(significant), "\n")
  
  if (plots$volcano){
    
    volcano_df<-results
    volcano_df$DEG<-"Not significant"
    
    # en limma el p-valor ajustado se llama adj.P.Val
    volcano_df$DEG[volcano_df$logFC>lfc_threshold &
                     volcano_df$adj.P.Val<fdr_threshold]<-"Up"
    volcano_df$DEG[volcano_df$logFC< -lfc_threshold &
                     volcano_df$adj.P.Val<fdr_threshold]<-"Down"
    
    top_genes<-volcano_df %>%
      dplyr::filter(DEG != "Not significant") %>%
      dplyr::arrange(adj.P.Val) %>%
      head(10)
    
    #color azul para diferenciarlo del rojo de edgeR
    p<-ggplot(volcano_df, aes(x=logFC, y=-log10(adj.P.Val), color=DEG))+
      geom_point(alpha=0.6, size=1)+
      scale_color_manual(values=c("Not significant"="grey", "Up"="firebrick3", "Down"="dodgerblue3"))+
      geom_vline(xintercept=c(-lfc_threshold, lfc_threshold), linetype="dashed", color="black")+
      geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="black")+
      geom_text_repel(
        data=top_genes,
        aes(label=gene_name),
        size=3
      )+
      labs(title=paste(plot_prefix,"- Volcano Plot"),
           x="Log2 Fold Change",
           y="-Log10 adj.P.Val")+
      theme_minimal()
    
    print(p) #Usamos print() para asegurar que el gráfico se dibuja al ejecutar la función
  }
  
  if (plots$hm){
    top_genes<- significant %>%
      arrange(adj.P.Val) %>%
      head(50) %>%
      pull(id_cruce)
    
    heatmap_counts<-v$E
    rownames(heatmap_counts)<-sub("\\..*","",rownames(heatmap_counts))
    heatmap_counts<-heatmap_counts[rownames(heatmap_counts) %in% top_genes,]
    
    annotation_col<-data.frame(Group=group)
    rownames(annotation_col)<-colnames(heatmap_counts)
    
    pheatmap(
      heatmap_counts,
      scale="row",
      annotation_col=annotation_col,
      show_rownames=FALSE,
      show_colnames = FALSE,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      main=paste(plot_prefix,"- Heatmap Top DEGs")
    )
    
  }
  
  return(list(
    all_genes=results,
    significant=significant
  ))
}

#Análisis EdgeR y voom para LUAD:
res_LUAD_edgeR<-complete_edgeR_analysis(expr_LUAD_def)
res_LUAD_voom<-analysis_voom(expr_LUAD_def)
saveRDS(res_LUAD_edgeR,"data/res_LUAD_edgeR.rds")
saveRDS(res_LUAD_voom,"data/res_LUAD_voom.rds")

common_deg_LUAD <- intersect(
  res_LUAD_edgeR$significant$id_cruce,
  res_LUAD_voom$significant$id_cruce
)
saveRDS(common_deg_LUAD,"data/common_deg_LUAD.rds")

#Análisis EdgeR y voom para LUSC:
res_LUSC_edgeR<-complete_edgeR_analysis(expr_LUSC_def)
res_LUSC_voom<-analysis_voom(expr_LUSC_def)
saveRDS(res_LUSC_edgeR,"data/res_LUSC_edgeR.rds")
saveRDS(res_LUSC_voom,"data/res_LUSC_voom.rds")

common_deg_LUSC <- intersect(
  res_LUSC_edgeR$significant$id_cruce,
  res_LUSC_voom$significant$id_cruce
)
saveRDS(common_deg_LUSC,"data/common_deg_LUSC.rds")

common_lung_deg<-intersect(common_deg_LUAD, common_deg_LUSC)
saveRDS(common_lung_deg,"data/common_lung_deg.rds")
