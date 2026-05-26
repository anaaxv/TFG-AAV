if(!require("here", quietly = TRUE))
  install.packages("here")

library(here)

source(here("funciones", "Funciones_TFG.R"))

#________________________________________________________________________

#BLOQUE 1: DESCARGA DE DATOS TRANSCRIPTÓMICOS, CLÍNICOS Y DE MUTACIONES: 

#Definimos la ruta al directorio en el que se descargarán los datos y lo creamos:
dir_gdc <- here("TCGA","GDCdata") 

if (!dir.exists(dir_gdc)){
  dir.create(dir_gdc, recursive = TRUE) #Con recursive true nos aseguramos de que se creen las carpetas intermediarias que puedan faltar
}

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

dir.create("data")
saveRDS(expr_LUAD_raw, file="data/expr_LUAD_raw.rds")
saveRDS(expr_LUSC_raw, file="data/expr_LUSC_raw.rds")
saveRDS(mut_LUAD, file="data/mut_LUAD_raw.rds")
saveRDS(mut_LUSC, file="data/mut_LUSC_raw.rds")
saveRDS(clinical_LUAD_BCRtab, file="data/clinical_LUAD_BCRtab.rds")
saveRDS(clinical_LUSC_BCRtab, file="data/clinical_LUSC_BCRtab.rds")

#Dejamos el espacio de trabajo lo más limpio posible
#De esta manera trabajaremos solo con lo imprescindible en cada momento, para ahorrar RAM:
rm(query_expr_LUAD, query_expr_LUSC, 
   query_mut_LUAD, query_mut_LUSC,
   query_clinical_LUAD, query_clinical_LUSC,
   mut_LUAD, mut_LUSC)
gc()
#________________________________________________________________________
#BLOQUE 2: PROCESAMIENTO DE DATOS DE EXPRESIÓN

#Obtenemos los perfiles de expresión definitivos para LUAD y LUSC:
res_LUAD<-procesar_expr(
  expr=expr_LUAD_raw,
  clinical_patients=clinical_LUAD_BCRtab$clinical_patient_luad$bcr_patient_barcode
)

expr_LUAD_def<-res_LUAD$expr_def
common_patients_LUAD<-res_LUAD$common_patients

saveRDS(expr_LUAD_def,"data/expr_LUAD_def.rds")
saveRDS(common_patients_LUAD, "data/common_patients_LUAD.rds")

rm(expr_LUAD_raw, res_LUAD, common_patients_LUAD)
gc()

#Hacemos lo mismo para LUSC:
res_LUSC<-procesar_expr(
  expr=expr_LUSC_raw,
  clinical_patients=clinical_LUSC_BCRtab$clinical_patient_lusc$bcr_patient_barcode
)

expr_LUSC_def<-res_LUSC$expr_def
common_patients_LUSC<-res_LUSC$common_patients

saveRDS(expr_LUSC_def,"data/expr_LUSC_def.rds")
saveRDS(common_patients_LUSC, "data/common_patients_LUSC.rds")

rm(expr_LUSC_raw,res_LUSC,common_patients_LUSC)
gc()

#________________________________________________________________________
#BLOQUE 3: ANÁLISIS DE EXPRESIÓN (EdgeR y Voom) 

#Definimos los parámetros de prueba
cancers <- c("LUAD", "LUSC")
fdr_values <- c(0.01, 0.05)
lfc_values <- c(1, 2)

#Creamos todas las combinaciones posibles (8 en total: 2 cánceres * 2 FDR * 2 LFC)
params <- expand.grid(cancer = cancers, fdr = fdr_values, lfc = lfc_values, stringsAsFactors = FALSE)

#Lista para almacenar los DEGs comunes de LUAD y LUSC, por combinación de umbrales
common_deg_lung <- list()

#Iniciamos el bucle
for (i in 1:nrow(params)) {
  
  # Extraer parámetros actuales
  can <- params$cancer[i]
  curr_fdr <- params$fdr[i]
  curr_lfc <- params$lfc[i]
  
  cat("\n--- Procesando:", can, "| FDR:", curr_fdr, "| LFC:", curr_lfc, "---\n")
  
  # Seleccionar el objeto de expresión correspondiente:
  expr_data <- get(paste0("expr_", can, "_def"))
  
  # A.Análisis EdgeR
  res_edgeR <- complete_edgeR_analysis(
    expr_data,
    lfc_threshold = curr_lfc,
    fdr_threshold = curr_fdr,
    plot_prefix = paste(can, "- EdgeR")
  )
  
  # B.Análisis Voom 
  res_voom <- analysis_voom(
    expr_data,
    lfc_threshold = curr_lfc,
    fdr_threshold = curr_fdr,
    plot_prefix = paste(can, "- Voom")
  )
  
  # C.Intersección por cáncer (Common DEGs del cáncer actual)
  common_can <- intersect(
    res_edgeR$significant$id_cruce,
    res_voom$significant$id_cruce
  )
  
  #Guardamos:
  suffix <- paste0(can, "_FDR", curr_fdr, "_LFC", curr_lfc)
  saveRDS(res_edgeR, paste0("data/res_edgeR_", suffix, ".rds"))
  saveRDS(res_voom, paste0("data/res_voom_", suffix, ".rds"))
  saveRDS(common_can, paste0("data/common_deg_", suffix, ".rds"))
  
  #Guardamos en una lista temporal para el cruce final entre LUAD y LUSC
  # Usamos una clave para agrupar por umbrales
  threshold_key <- paste0("FDR", curr_fdr, "_LFC", curr_lfc)
  common_deg_lung[[threshold_key]][[can]] <- common_can
}

rm(res_edgeR, res_voom, expr_data, common_can)
gc()

#Calculamos la intersección FINAL (LUAD vs LUSC) para cada umbral
cat("\n--- Calculando Intersecciones Finales (Common Lung) ---\n")

for (key in names(common_deg_lung)) {
  # Intersectamos lo que hay en LUAD con lo que hay en LUSC para este umbral
  common_lung <- intersect(
    common_deg_lung[[key]][["LUAD"]],
    common_deg_lung[[key]][["LUSC"]]
  )
  
  cat("Umbral", key, "- Genes comunes totales:", length(common_lung), "\n")
  saveRDS(common_lung, paste0("data/common_lung_deg_", key, ".rds"))
}

rm(cancers, fdr_values, lfc_values, common_deg_lung, common_lung)
gc()

#________________________________________________________________________
#BLOQUE 3.5 (Opcional): Generamos gráfica de distribución de DEGs:

# Bucle para cargar TODOS los .rds resultantes del análisis de expresión (Edge R y Voom) en el environment:
for (i in 1:nrow(params)) {
  can      <- params$cancer[i]
  curr_fdr <- params$fdr[i]
  curr_lfc <- params$lfc[i]
  
  #Creamos el sufijo y los nombres de las variables en formato texto
  suffix <- paste0(can, "_FDR", curr_fdr, "_LFC", curr_lfc)
  var_edgeR <- paste0("res_edgeR_", suffix)
  var_voom  <- paste0("res_voom_", suffix)
  
  #Construimos la ruta hacia los archivos .rds
  ruta_edgeR <- paste0("data/res_edgeR_", suffix, ".rds")
  ruta_voom  <- paste0("data/res_voom_", suffix, ".rds")
  
  #Si los archivos existen, los leemos y asignamos su contenido a la variable dinámica
  if (file.exists(ruta_edgeR)) {
    assign(var_edgeR, readRDS(ruta_edgeR), envir = .GlobalEnv)
  }
  if (file.exists(ruta_voom)) {
    assign(var_voom, readRDS(ruta_voom), envir = .GlobalEnv)
  }
}

rm(params, can, curr_fdr, curr_lfc, i, ruta_edgeR, ruta_voom, suffix, var_edgeR, var_voom)
gc()

grafica_degs <- plot_deg_distribution()
print(grafica_degs)

luad_compartidos <- generar_venn_comparativo(
  res_edgeR = res_edgeR_LUAD_FDR0.01_LFC2, 
  res_voom = res_voom_LUAD_FDR0.01_LFC1, 
  tipo_cancer = "LUAD"
)

lusc_compartidos <- generar_venn_comparativo(
  res_edgeR = res_edgeR_LUSC_FDR0.01_LFC2, 
  res_voom = res_voom_LUSC_FDR0.01_LFC2, 
  tipo_cancer = "LUSC"
)

objetos_a_borrar <- ls(pattern = "^res_(edgeR|voom)_")

#Si encuentra objetos, los borra en lote y libera memoria RAM
if (length(objetos_a_borrar) > 0) {
  rm(list = objetos_a_borrar, envir = .GlobalEnv)
  cat("\n--- Se han eliminado", length(objetos_a_borrar), "objetos de DEGs del entorno. ---\n")
} else {
  cat("\n--- No se encontraron objetos de DEGs para borrar. ---\n")
}

rm(objetos_a_borrar)
gc()

#________________________________________________________________________
#BLOQUE 4: GENERAMOS LOS INPUTS PARA LAS HERRAMIENTAS DE REPOSICIONAMIENTO

#ShinyDeepDR, expresión:
inputs_shinyDeepDR(
  expr_data = expr_LUAD_def,
  file_name="input_ShinyDeepDR_LUAD"
)

inputs_shinyDeepDR(
  expr_data = expr_LUSC_def,
  file_name="input_ShinyDeepDR_LUSC"
)

rm(expr_LUAD_def,expr_LUSC_def)
gc()

#ShinyDeepDR, mutaciones:
mut_LUAD_raw<-readRDS("data/mut_LUAD_raw.rds")

maf_shinyDeepDR_LUAD <- mut_LUAD_raw %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode) %>%
  mutate(Tumor_Sample_Barcode = "LUAD_Global") %>% # Forzamos una única muestra
  distinct() # Eliminamos duplicados si el mismo gen muta igual en varios pacientes

# Guardamos el MAF "unificado":
write.table(maf_shinyDeepDR_LUAD, 
            file = "data/analysis/inputs_repo/ShinyDeepDR/LUAD_mutations.maf", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

rm(mut_LUAD_raw)
gc()

mut_LUSC_raw<-readRDS("data/mut_LUSC_raw.rds")

maf_shinyDeepDR_LUSC <- mut_LUSC_raw %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode) %>%
  mutate(Tumor_Sample_Barcode = "LUSC_Global") %>% # Forzamos una única muestra
  distinct() # Eliminamos duplicados si el mismo gen muta igual en varios pacientes

write.table(maf_shinyDeepDR_LUSC, 
            file = "data/analysis/inputs_repo/ShinyDeepDR/LUSC_mutations.maf", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

rm(mut_LUSC_raw)
gc()

#CDRPipe:

res_LUAD_voom<- readRDS("data/res_voom_LUAD_FDR0.01_LFC2.rds")
common_deg_LUAD<-readRDS("data/common_deg_LUAD_FDR0.01_LFC2.rds")
res_LUSC_voom<- readRDS("data/res_voom_LUSC_FDR0.01_LFC2.rds")
common_deg_LUSC<-readRDS("data/common_deg_LUSC_FDR0.01_LFC2.rds")

df_comun_luad <- res_LUAD_voom$significant %>%
  filter(id_cruce %in% common_deg_LUAD)

input_cdr_luad <- inputs_cdrpipe(
  df_degs = df_comun_luad, 
  file_name = "input_CDRpipe_LUAD"
)

df_comun_lusc <- res_LUSC_voom$significant %>%
  filter(id_cruce %in% common_deg_LUSC)

input_cdr_lusc <- inputs_cdrpipe(
  df_degs = df_comun_lusc, 
  file_name = "input_CDRpipe_LUSC"
)

rm(res_LUAD_voom, res_LUSC_voom, common_deg_LUAD, common_deg_LUSC)
gc()

#iLINCS:
inputs_ilincs(
  df_cdr = input_cdr_luad, 
  file_name = "input_iLINCS_LUAD"
)

inputs_ilincs(
  df_cdr = input_cdr_lusc, 
  file_name = "input_iLINCS_LUSC"
)

#Clue Query (CMap):
inputs_clue(df_comun_luad, "LUAD")
inputs_clue(df_comun_lusc, "LUSC")

#________________________________________________________________________
#BLOQUE 5: TRATAMIENTO DE RESULTADOS DE CLUE QUERY 
res_LUAD_CMap <- analizar_resultados_clue(
  ruta_gct = "data/analysis/results/resultados_repo/LUAD_Clue_comunes/my_analysis.sig_queryl1k_tool.6a03285d8ed9720013827997/ncs.gct",
  nombre_archivo = "CMap_Resultados_LUAD"
)
saveRDS(res_LUAD_CMap,here("data", "candidatos_LUAD_CMap.rds")) 

res_LUSC_CMap <- analizar_resultados_clue(
  ruta_gct = "data/analysis/results/resultados_repo/LUSC_Clue_comunes/my_analysis.sig_queryl1k_tool.6a0330648ed24f001428125f/ncs.gct",
  nombre_archivo = "CMap_Resultados_LUSC"
)
saveRDS(res_LUAD_CMap,here("data", "candidatos_LUSC_CMap.rds")) 


#________________________________________________________________________
#BLOQUE 6: OBTENER FÁRMACOS CONSENSO 

#Ejecutamos función de fármacos consenso para LUAD:
farmacos_consenso_LUAD<- obtener_farmacos_consenso(
  ruta_cdr = "data/analysis/results/resultados_repo/LUAD_CDRPipe_Resultados.csv",
  ruta_cmap = "data/analysis/results/resultados_repo/CMap/CMap_Resultados_LUAD.csv",
  ruta_ilincs = "data/analysis/results/resultados_repo/iLINCS/LUAD/LUAD_exemplar_top100_resultados.xls",
  ruta_shiny = "data/analysis/results/resultados_repo/ShinyDeepDR/shinyDeepDR_LUAD.csv",
  archivo_salida = "Farmacos_Consenso_LUAD"
)

#LUSC:
farmacos_consenso_LUSC<-obtener_farmacos_consenso(
  ruta_cdr = "data/analysis/results/resultados_repo/LUSC_CDRPipe_Resultados.csv",
  ruta_cmap = "data/analysis/results/resultados_repo/CMap/CMap_Resultados_LUSC.csv",
  ruta_ilincs = "data/analysis/results/resultados_repo/iLINCS/LUSC/LUSC_exemplar_top100_resultados.xls",
  ruta_shiny = "data/analysis/results/resultados_repo/ShinyDeepDR/shinyDeepDR_LUSC.csv",
  archivo_salida = "Farmacos_Consenso_LUSC"
)




