
#Instalamos las librerías necesarias si no las tenemos ya: 

paquetes_cran <- c("DT", "dplyr", "pheatmap", "ggrepel", "stringr", "here", "UpSetR", "patchwork", "tidyr", "ggvenn", "ggplot2")
paquetes_bioc <- c("TCGAbiolinks", "edgeR", "limma", "biomaRt", "cmapR", "SummarizedExperiment")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (x in paquetes_cran) {
  if (!requireNamespace(x, quietly = TRUE)) install.packages(x)
}
for (x in paquetes_bioc) {
  if (!requireNamespace(x, quietly = TRUE)) BiocManager::install(x, update = FALSE)
}

#Los cargamos en la sesión:
paquetes_totales <- c(paquetes_cran, paquetes_bioc)
invisible(lapply(paquetes_totales, library, character.only = TRUE))


#expr_LUAD es un SummarizedExperiment en el que las filas son los genes y las columnas las muestras
#de RNA-seq, cada columna corresponde a un archivo .tsv
#cuando TCGAbiolinks hace gdcprepare() pone el barcode de la muestra como nombre de la columna

#FUNCIÓN PARA PROCESAR DATOS DE EXPRESIÓN:
#Eliminación de duplicados, cruce con clínicos, quedarnos solo con dos tipos de muestras...
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

#-----------------
#Funciones de análisis:

#EdgeR:

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
      geom_point(size=0.8)+
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
 
    volcano_df <- volcano_df %>% 
      dplyr::filter(gene_type == "protein_coding")
    
    top_up <- volcano_df %>%
      dplyr::filter(DEG == "Up") %>%
      dplyr::arrange(FDR) %>%
      head(5)
    
    top_down <- volcano_df %>%
      dplyr::filter(DEG == "Down") %>%
      dplyr::arrange(FDR) %>%
      head(5)
    
    top_genes <- rbind(top_up, top_down)
    
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
        seed=42,
        show.legend = FALSE
      )+
      labs(title=paste(plot_prefix,"- Volcano Plot (Protein-Coding)"),
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
  #Para cada uso de la función se crea una carpeta nueva llamada como su combinación de umbrales:
  threshold_dir<- paste0("FDR_", fdr_threshold, "_LFC_", lfc_threshold)
  
  results_folder <- here(base_folder, "analysis", "results", threshold_dir)
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive=TRUE)
  }
  
  write.table(significant,
              file=file.path(results_folder,
                             paste0(plot_prefix,"_DEGs.txt")),
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
                                  base_folder="data",
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
    base_folder=base_folder,
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
                        plots=list(pca=TRUE, volcano=TRUE, hm=TRUE), 
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
  cols_principales <- c("id_cruce", "gene_name", "gene_type", "logFC", "AveExpr", "P.Value", "adj.P.Val")
  cols_restantes <- setdiff(colnames(results), c(cols_principales, "DEG")) 
  results <- results[, c(cols_principales, cols_restantes, "DEG")] 
  
  significant <- subset(results, DEG != "Not significant")
  up_genes <- subset(results, DEG == "Up")
  down_genes <- subset(results, DEG == "Down")
  
  cat("Número de DEGs (voom):", nrow(significant), " (Up:", nrow(up_genes), ", Down:", nrow(down_genes), ")\n")
  
  if (plots$pca) {
    pca <- prcomp(t(v$E), scale.=TRUE)
    
    pca_df <- data.frame(
      PC1 = pca$x[,1],
      PC2 = pca$x[,2],
      group = group
    )
    
    # Varianza explicada:
    var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1) 
    
    p_pca <- ggplot(pca_df, aes(PC1, PC2, color=group)) +
      geom_point(size=0.8) +
      labs(title=paste(plot_prefix, "- PCA"),
           x=paste0("PC1 (", var_expl[1], "%)"),
           y=paste0("PC2 (", var_expl[2], "%)")) +
      theme_minimal()
    
    print(p_pca)
  }
  
  if (plots$volcano){
    volcano_df<-results
    
    volcano_df <- volcano_df %>% 
      dplyr::filter(gene_type == "protein_coding")
    
    top_up <- volcano_df %>%
      dplyr::filter(DEG == "Up") %>%
      dplyr::arrange(adj.P.Val) %>% 
      head(5)
    
    top_down <- volcano_df %>%
      dplyr::filter(DEG == "Down") %>%
      dplyr::arrange(adj.P.Val) %>% 
      head(5)
    
    top_genes <- rbind(top_up, top_down)
    
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
        seed=42,
        show.legend=FALSE
      )+
      labs(title=paste(plot_prefix,"- Volcano Plot (Protein-Coding)"),
           x="Log2 Fold Change",
           y="-Log10 adj.P.Val")+ 
      theme_minimal()
    
    print(p) 
  }
  
  if (plots$hm && nrow(significant) > 0) {
    #Seleccionamos los IDs de los top genes:
    top_genes_hm <- significant %>%
      arrange(adj.P.Val) %>% 
      head(50) %>%
      pull(id_cruce)
    
    #Extraer los counts (v$E para voom):
    heatmap_data <- v$E 
    rownames(heatmap_data) <- sub("\\..*", "", rownames(heatmap_data))
    heatmap_data <- heatmap_data[rownames(heatmap_data) %in% top_genes_hm, ]
    
    #Ordenamos las muestras por grupo:
    orden_muestras <- order(group)
    heatmap_data <- heatmap_data[, orden_muestras] 
    
    #Preparamos la anotación con el nuevo orden:
    annotation_col <- data.frame(Group = group[orden_muestras])
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
  
  threshold_dir<- paste0("FDR_", fdr_threshold, "_LFC_", lfc_threshold)
  results_folder <- file.path(base_folder, "analysis", "results", threshold_dir) 
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive=TRUE)
  }
  
  #Guardamos DEGs significativos:
  write.table(significant, 
              file=file.path(results_folder, paste0(plot_prefix, "_DEGs_.txt")),
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
#__________
#Gráfica de distribución de DEGS:
plot_deg_distribution <- function(cancers = c("LUAD", "LUSC"),
                                  metodos = c("edgeR", "voom"),
                                  fdr_values = c("0.05", "0.01"),
                                  lfc_values = c("1", "2"),
                                  max_y = 7500) {
  

  df_general <- data.frame()
  
  #Bucle de extracción automatizada buscando en el entorno
  for (can in cancers) {
    for (met in metodos) {
      for (curr_fdr in fdr_values) {
        for (curr_lfc in lfc_values) {
          
          nombre_objeto <- paste0("res_", met, "_", can, "_FDR", curr_fdr, "_LFC", curr_lfc)
          
          # Especificamos envir = .GlobalEnv para que encuentre tus objetos cargados
          if (exists(nombre_objeto, envir = .GlobalEnv)) {
            obj <- get(nombre_objeto, envir = .GlobalEnv)
            
            up_total   <- nrow(obj$up_genes)
            up_pc      <- nrow(obj$up_genes_pc)
            down_total <- nrow(obj$down_genes)
            down_pc    <- nrow(obj$down_genes_pc)
            
            temp_df <- rbind(
              # UP
              data.frame(Cancer = can, Metodo = ifelse(met == "edgeR", "EdgeR", "Voom"),
                         Umbral = paste0("LFC ", curr_lfc, " | FDR ", curr_fdr),
                         Regulacion = "Up", Tipo_Gene = "Protein Coding", Conteo = up_pc, Total_Label = up_total, PC_Label = up_pc),
              data.frame(Cancer = can, Metodo = ifelse(met == "edgeR", "EdgeR", "Voom"),
                         Umbral = paste0("LFC ", curr_lfc, " | FDR ", curr_fdr),
                         Regulacion = "Up", Tipo_Gene = "Otros DEGs", Conteo = up_total - up_pc, Total_Label = up_total, PC_Label = up_pc),
              # DOWN
              data.frame(Cancer = can, Metodo = ifelse(met == "edgeR", "EdgeR", "Voom"),
                         Umbral = paste0("LFC ", curr_lfc, " | FDR ", curr_fdr),
                         Regulacion = "Down", Tipo_Gene = "Protein Coding", Conteo = down_pc, Total_Label = down_total, PC_Label = down_pc),
              data.frame(Cancer = can, Metodo = ifelse(met == "edgeR", "EdgeR", "Voom"),
                         Umbral = paste0("LFC ", curr_lfc, " | FDR ", curr_fdr),
                         Regulacion = "Down", Tipo_Gene = "Otros DEGs", Conteo = down_total - down_pc, Total_Label = down_total, PC_Label = down_pc)
            )
            df_general <- rbind(df_general, temp_df)
          }
        }
      }
    }
  }
  
  #Control de seguridad por si acaso no encuentra ningún objeto
  if (nrow(df_general) == 0) {
    stop("No se encontraron objetos con el patrón 'res_metodo_cancer_FDR_LFC' en tu entorno.")
  }
  
  #Generamos el orden de los niveles del eje X de forma dinámica según tus inputs
  lista_umbrales <- c()
  for (l in lfc_values) {
    for (f in fdr_values) {
      lista_umbrales <- c(lista_umbrales, paste0("LFC ", l, " | FDR ", f))
    }
  }
  df_general$Umbral <- factor(df_general$Umbral, levels = unique(lista_umbrales))
  
  #Formateo de posiciones y etiquetas
  df_general <- df_general %>%
    mutate(
      Metodo_num = ifelse(Metodo == "EdgeR", 1, 2),
      X_pos = ifelse(Regulacion == "Up", Metodo_num - 0.2, Metodo_num + 0.2),
      Etiqueta = paste0(Total_Label, "\n(", PC_Label, ")")
    )
  
  df_general$Tipo_Gene <- factor(df_general$Tipo_Gene, levels = c("Protein Coding", "Otros DEGs"))
  
  #Construcción de la gráfica de ggplot2
  grafica_degs <- ggplot(df_general, aes(x = X_pos, y = Conteo, fill = Regulacion, alpha = Tipo_Gene)) +
    geom_col(position = "stack", width = 0.35, color = "black") +
    
    geom_text(data = df_general %>% filter(Tipo_Gene == "Protein Coding"), 
              aes(x = X_pos, y = Total_Label, label = Etiqueta),
              vjust = -0.15, lineheight = 0.85, size = 2.5, alpha = 1, inherit.aes = FALSE) +
    
    scale_fill_manual(values = c("Up" = "firebrick3", "Down" = "dodgerblue3")) +
    scale_alpha_manual(values = c("Protein Coding" = 1.0, "Otros DEGs" = 0.45)) +
    
    scale_x_continuous(breaks = seq_along(metodos), labels = tools::toTitleCase(metodos)) +
    scale_y_continuous(limits = c(0, max_y), expand = expansion(mult = c(0, 0.02))) +
    
    facet_grid(Cancer ~ Umbral) +
    labs(title = "Distribución de DEGs Totales y Codificantes para Proteína",
         x = "Método", y = "Número de DEGs", fill = "Regulación", alpha = "Fracción de la Barra") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"), 
          legend.position = "bottom")

  return(grafica_degs)
}


generar_venn_comparativo <- function(res_edgeR, res_voom, tipo_cancer="LUAD", base_folder="data") {
  #Extraemos los IDs de los genes (id_cruce) que son UP y DOWN en edgeR
  edger_up   <- res_edgeR$up_genes$id_cruce
  edger_down <- res_edgeR$down_genes$id_cruce
  
  #Extraemos los IDs de los genes que son UP y DOWN en voom
  voom_up    <- res_voom$up_genes$id_cruce
  voom_down  <- res_voom$down_genes$id_cruce
  
  lista_up <- list(
    "edgeR (Up)" = edger_up,
    "Voom (Up)"  = voom_up
  )
  
  lista_down <- list(
    "edgeR (Down)" = edger_down,
    "Voom (Down)"  = voom_down
  )
  
  colores <- c("firebrick3", "orange")
  colores_down <- c("dodgerblue3", "turquoise3")
  
  #Grafica UP
  p_up <- ggvenn(
    lista_up, 
    fill_color = colores,
    stroke_size = 0.5, 
    set_name_size = 4,  
    text_size = 3
  ) + 
    labs(title = paste0(tipo_cancer, " - Genes Sobreexpresados (UP)")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  
  #Grafica DOWN
  p_down <- ggvenn(
    lista_down, 
    fill_color = colores_down,
    stroke_size = 0.5, 
    set_name_size = 4,   
    text_size = 3
  ) + 
    labs(title = paste0(tipo_cancer, " - Genes Infraexpresados (DOWN)")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) 
  
  print(p_up)
  print(p_down)
  
  output_dir <- file.path(base_folder, "analysis", "results", "Venn_Diagrams")
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
  
  ggsave(file.path(output_dir, paste0(tipo_cancer, "_Venn_UP.png")), plot=p_up, width=6, height=5, dpi=300)
  ggsave(file.path(output_dir, paste0(tipo_cancer, "_Venn_DOWN.png")), plot=p_down, width=6, height=5, dpi=300)
  
  return(list(
    shared_up = intersect(edger_up, voom_up),
    shared_down = intersect(edger_down, voom_down)
  ))
}



#----------------------------------

#Funciones para generar los inputs de reposicionamiento:

#Función que genera los inputs para ShinyDeepDR:
#Para shinydeepDR tienes que hacer como si solo fuera una muestra
#Serviría muy bien para medicina personalizada
#Extrae la mediana de los TPMs de las muestras tumorales:

inputs_shinyDeepDR<-function(expr_data, file_name, folder=here("data","analysis","inputs_repo", "ShinyDeepDR")) {
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
  ruta_salida <- file.path(folder, paste0(file_name, ".txt"))
  write.table(df_shiny, file = ruta_salida, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Archivo shinyDeepDR guardado en:", ruta_salida, "\n")
  cat("Tumores promediados:", sum(es_tumor), "\n")
}

inputs_cdrpipe<-function(df_degs, file_name, folder=here("data","analysis","inputs_repo", "CDRPipe")){
  
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
  ruta_salida <- file.path(folder, paste0(file_name, ".csv"))
  write.csv(df_cdr, file = ruta_salida, row.names = FALSE, quote = FALSE)
  
  cat("Archivo para CDRpipe guardado en:", ruta_salida, "\n")
  cat("Total de genes exportados:", nrow(df_cdr), "\n")
  
  return(df_cdr)
}


#Generación de inputs para iLINCS basándonos en el objeto de CDRpipe:
inputs_ilincs <- function(df_cdr, file_name, folder = here("data","analysis","inputs_repo", "iLINCS")) {
  
  #Creamos directorio si no existe:
  if(!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }
  
  ruta_salida_ilincs <- file.path(folder, paste0(file_name, ".txt"))
  
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


#Clue Query (CMap):
inputs_clue <- function(df_comun, prefix, folder= here("data","analysis","inputs_repo", "Clue")) {
  
  if(!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }
  
  up_genes <- df_comun %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(150) %>% pull(gene_name)
  down_genes <- df_comun %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% head(150) %>% pull(gene_name)
  
  write.table(up_genes, file = file.path(folder, paste0(prefix, "_up.txt")), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(down_genes, file = file.path(folder, paste0(prefix, "_down.txt")), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("Archivos CMap generados para", prefix, "con", length(up_genes), "up y", length(down_genes), "down genes.\n")
}

#-----------------
##TRATAMIENTO DE LOS RESULTADOS DE CLUE:
analizar_resultados_clue <- function(ruta_gct, nombre_archivo, folder= here("data","analysis","results","resultados_repo", "CMap")) {
  
  if(!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }
  
  gct_data <- cmapR::parse_gctx(ruta_gct)
  
  metadatos_farmacos <- as.data.frame(gct_data@rdesc)
  scores_ncs <- as.data.frame(gct_data@mat)
  
  # Si el archivo tiene muchas columnas, promediamos por fila para tener un score global
  # O seleccionamos la columna 'summary' si existe.
  ncs_promedio <- rowMeans(scores_ncs) 
  
  tabla_resultados <- data.frame(
    Nombre_Farmaco = metadatos_farmacos$pert_iname,
    Mecanismo_Accion = metadatos_farmacos$moa,
    Target = metadatos_farmacos$target_name,
    NCS = ncs_promedio 
  )
  
  # Filtro de calidad: Solo fármacos con impacto real (NCS < -1.5 es un estándar común)
  mis_candidatos <- tabla_resultados %>%
    filter(NCS < -1.5) %>% 
    group_by(Nombre_Farmaco) %>% 
    summarise(Mecanismo_Accion = first(Mecanismo_Accion),
              Target = first(Target),
              NCS = min(NCS)) %>% # Si hay duplicados, nos quedamos con el más potente
    arrange(NCS)
  
  ruta_final<-file.path(folder, paste0(nombre_archivo, ".csv"))
  write.csv(mis_candidatos, ruta_final, row.names = FALSE)
  return(mis_candidatos)
}

#FUNCIÓN PARA OBTENER FÁRMACOS CONSENSO:

obtener_farmacos_consenso <- function(ruta_cdr, ruta_cmap, ruta_ilincs, ruta_shiny, archivo_salida,
                                      folder=here("data","analysis","results","resultados_repo", "Consenso")) {
  
  if(!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }
  
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
  
  #Sustituir NAs en las columnas de presencia (2 a 5)
  tabla_final[, 2:5][is.na(tabla_final[, 2:5])] <- 0
  
  
  #Cálculo de totales y orden final:
  tabla_final <- tabla_final %>%
    mutate(Total_Metodos = CDRpipe + CMap + iLINCS + ShinyDeepDR) %>%
    #Prioridad: 1º Más métodos, 2º Mejor reversión en iLINCS, 3º Mejor IC50 en Shiny
    arrange(desc(Total_Metodos), Score_iLINCS, IC50_Shiny)
  
  ruta_final<-file.path(folder, paste0(archivo_salida, ".csv"))
  write.csv(tabla_final, ruta_final, row.names = FALSE)
  
  nombre_rds <- paste0("data/", archivo_salida, ".rds" )
  saveRDS(tabla_final, nombre_rds)
  
  #GRÁFICOS:
  
  #UpSet
  listas_interseccion <- list(
    CDRpipe     = tabla_final$Farmaco[tabla_final$CDRpipe == 1],
    CMap        = tabla_final$Farmaco[tabla_final$CMap == 1],
    iLINCS      = tabla_final$Farmaco[tabla_final$iLINCS == 1],
    ShinyDeepDR = tabla_final$Farmaco[tabla_final$ShinyDeepDR == 1]
  )
  
  print(upset(fromList(listas_interseccion), 
              order.by = "freq", 
              main.bar.color = "#2c3e50", 
              sets.bar.color = "#e74c3c",
              matrix.color = "#2c3e50",
              text.scale = 1.2,
              set_size.show = FALSE))
  
  #Boxplots:
  p1 <- ggplot(tabla_final, aes(x = as.factor(Total_Metodos), y = Score_iLINCS, fill = as.factor(Total_Metodos))) +
    geom_boxplot(alpha = 0.7, outlier.color = "red") +
    theme_light() +
    labs(x = "Nivel de Consenso", y = "Score iLINCS", title = "Reversión por Consenso") +
    theme(legend.position = "none")
  
  p2 <- ggplot(tabla_final, aes(x = as.factor(Total_Metodos), y = IC50_Shiny, fill = as.factor(Total_Metodos))) +
    geom_boxplot(alpha = 0.7, outlier.color = "red") +
    theme_light() +
    labs(x = "Nivel de Consenso", y = "IC50 Shiny (log uM)", title = "IC50 por Consenso") +
    theme(legend.position = "none")
  
  print(p1 + p2)
  
  cat("Fármacos consenso obtenidos\n")
  
  return(tabla_final)
}


