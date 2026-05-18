---
editor_options: 
  markdown: 
    wrap: 72
---

# Análisis de datos ómicos y reposicionamiento farmacológico en enfermedades complejas

> Trabajo de Fin de Grado — Ana Alameda Valenzuela — Biotecnología UGR\
> Análisis de expresión diferencial e identificación de candidatos
> terapéuticos en **LUAD** y **LUSC** usando datos del TCGA

------------------------------------------------------------------------

## Descripción

Este repositorio contiene el código desarrollado para mi TFG, cuyo
objetivo es identificar fármacos candidatos para el reposicionamiento en
cáncer de pulmón, aunque es fácilmente adaptable al análisis de
cualquier otra enfermedad del TCGA. Los dos subtipos de cáncer de pulmón
que se abordan son:

-   **LUAD** — Adenocarcinoma de pulmón (*Lung Adenocarcinoma*)
-   **LUSC** — Carcinoma de células escamosas de pulmón (*Lung Squamous
    Cell Carcinoma*)

La estrategia combina análisis de expresión diferencial con varias
herramientas de reposicionamiento farmacológico para obtener candidatos
consenso.

------------------------------------------------------------------------

## Estructura de carpetas:

```         
├── funciones/
│   └── Funciones_TFG.R        # Todas las funciones reutilizables del proyecto
├── main.R                     # Script de ejecución principal 
├── TCGA/
│   └── GDCdata/               # Datos descargados del GDC 
├── data/
│   ├── *.rds                  # Objetos R serializados (expresión, mutaciones, resultados)
│   └── analysis/
│       ├── results/# Resultados del análisis diferencial y reposicionamiento
│           └── resultados_repo/
│       └── inputs_repo/       # Inputs para las herramientas externas
│           ├── ShinyDeepDR/
│           ├── CDRPipe/
│           ├── iLINCS/
│           └── Clue/
├── README.md
└──...
```

Todas las carpetas aquí recogidas se generarán automáticamente al ir
ejecutando el script principal.

------------------------------------------------------------------------

## Flujo de trabajo

El script principal (`main.R`) está organizado en 6 bloques
secuenciales:

```         
BLOQUE 1 ─► Descarga de datos (TCGA)
    │
BLOQUE 2 ─► Procesamiento de expresión
    │
BLOQUE 3 ─► Análisis diferencial (EdgeR + Voom)
    │
BLOQUE 4 ─► Generación de inputs para herramientas de reposicionamiento
    │
BLOQUE 5 ─► Procesamiento de resultados de CMap (Clue Query)
    │
BLOQUE 6 ─► Integración y obtención de fármacos consenso
```

### Bloque 1 — Descarga de datos

Usa `TCGAbiolinks` para descargar desde el GDC:

| Tipo de dato               | Proyecto             | Formato                 |
|------------------------|------------------------|------------------------|
| Expresión génica (RNA-seq) | TCGA-LUAD, TCGA-LUSC | STAR - Counts           |
| Mutaciones somáticas       | TCGA-LUAD, TCGA-LUSC | Masked Somatic Mutation |
| Datos clínicos             | TCGA-LUAD, TCGA-LUSC | BCR Biotab              |

### Bloque 2 — Procesamiento de expresión

La función `procesar_expr()` realiza:

1.  Extracción del ID de paciente y tipo de muestra desde el barcode
    TCGA (`01` = Tumor, `11` = Normal).
2.  Eliminación de duplicados por paciente y tipo (se conserva la
    muestra con mayor tamaño de librería).
3.  Cruce con los datos clínicos disponibles.
4.  Para pacientes pareados (Tumor + Normal), se conserva únicamente la
    muestra Normal para el análisis.

### Bloque 3 — Análisis de expresión diferencial

Se aplican dos métodos complementarios con distintas combinaciones de
umbrales (FDR ∈ {0.01, 0.05} × LFC ∈ {1, 2}):

-   **EdgeR** (`complete_edgeR_analysis`): - Filtrado con
    `filterByExpr` - Normalización TMM - Modelo QL (*Quasi-Likelihood*)
    con `glmTreat` para test con umbral de fold change.

-   **Voom/Limma** (`analysis_voom`): - Filtrado: ≥10 counts en al menos
    el 10% de las muestras - Normalización TMM + transformación voom -
    Modelo lineal con `eBayes` y corrección FDR.

Los DEGs finales por cáncer son la **intersección** entre ambos métodos.
Los DEGs comunes a LUAD y LUSC se obtienen cruzando ambas listas.

**Gráficos generados:** Boxplot de logCPM normalizado, MDS, PCA, BCV, MA
plot, Volcano plot, Heatmap (top 50 genes).

### Bloque 4 — Inputs para reposicionamiento

| Herramienta | Función | Input generado |
|------------------------|------------------------|------------------------|
| **ShinyDeepDR** | `inputs_shinyDeepDR()` | Mediana de TPMs tumorales (expresión) + archivo MAF (mutaciones) |
| **CDRPipe** | `inputs_cdrpipe()` | CSV con SYMBOL, log2FC y p-valor ajustado |
| **iLINCS** | `inputs_ilincs()` | TXT tabulado (formato CDRPipe sin cabecera) |
| **CMap / Clue Query** | `inputs_clue()` | TXT con top 150 genes up y top 150 genes down |

Los inputs para CDRPipe, iLINCS y CMap se generan a partir de los DEGs
comunes entre EdgeR y Voom con los umbrales FDR = 0.01 y LFC = 2.

### Bloque 5 — Procesamiento de resultados CMap

La función `analizar_resultados_clue()` parsea los archivos `.gct`
devueltos por Clue Query y filtra los fármacos con **NCS \< −1.5**
(conectividad negativa fuerte = reversión del perfil tumoral).

### Bloque 6 — Fármacos consenso

`obtener_farmacos_consenso()` integra los resultados de las cuatro
herramientas:

-   Normaliza los nombres de fármacos (minúsculas, eliminación de
    espacios).
-   Filtra por dirección de reversión en cada herramienta:
    -   CDRPipe: `cmap_score < 0`
    -   CMap: `NCS < 0`
    -   iLINCS: `Concordance < 0`
    -   ShinyDeepDR: `IC50 (log µM) < 0`
-   Construye una tabla de presencia (0/1) para cada fármaco en cada
    herramienta.
-   Ordena por número de herramientas que validan el candidato, score
    iLINCS e IC50.

**Gráficos generados:** UpSet plot de intersecciones, boxplots de Score
iLINCS e IC50 por nivel de consenso.

------------------------------------------------------------------------

## Dependencias

El script instala automáticamente los paquetes ausentes al inicio. Las
dependencias son:

**Bioconductor:**

```         
TCGAbiolinks, edgeR, limma, biomaRt, cmapR
```

**CRAN:**

```         
DT, dplyr, ggplot2, pheatmap, ggrepel, stringr, here, UpSetR, patchwork
```

La instalación de `TCGAbiolinks` y otras dependencias de Bioconductor
requiere `BiocManager`:

``` r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
```

------------------------------------------------------------------------

## Uso

1.  Clona el repositorio y abre el proyecto en RStudio.
2.  Asegúrate de tener conexión a internet para la descarga de datos del
    GDC (Bloque 1).
3.  Ejecuta `main.R` bloque a bloque. Se recomienda no ejecutar el
    script completo de una vez, ya que:
    -   El **Bloque 1** puede tardar varias horas dependiendo de la
        conexión.
    -   El **Bloque 3** puede tardar aproximadamente un par de horas
        generando todos los gráficos.
    -   Los bloques 5 y 6 requieren que los resultados de las
        herramientas externas estén disponibles en las rutas indicadas.
4.  Los objetos intermedios se guardan como `.rds` en `data/` para poder
    retomar la ejecución sin repetir pasos costosos.

> **Nota:** Las herramientas externas (CMap/Clue Query, CDRPipe, iLINCS,
> ShinyDeepDR) son servicios web. Los inputs generados en el Bloque 4
> deben subirse manualmente a cada plataforma y los resultados
> descargarse antes de ejecutar los Bloques 5 y 6.

> **Advertencia:** Se recomienda no limpiar el entorno de R manualmente
> entre bloques. Algunos objetos se mantienen en memoria de forma
> deliberada porque son reutilizados en pasos posteriores. Si se
> necesita liberar memoria, reinicia desde el bloque correspondiente
> cargando los `.rds` de `data/`.

------------------------------------------------------------------------

## Herramientas externas de reposicionamiento

| Herramienta | URL | Descripción |
|------------------------|------------------------|------------------------|
| **Clue Query (CMap)** | [clue.io](https://clue.io) | Conectividad de firma génica con perfiles de fármacos (LINCS L1000) |
| **CDRPipe** | [cdrpipe.org](https://www.cdrpipe.org/) | Pipeline de reposicionamiento basado en expresión diferencial |
| **iLINCS** | [ilincs.org](http://ilincs.org) | Integración de firmas transcriptómicas con perturbágenos |
| **ShinyDeepDR** | [shinydeepdr](https://shiny.crc.pitt.edu/shinydeepdr/) | Predicción de sensibilidad a fármacos basada en deep learning |

------------------------------------------------------------------------

## Resultados esperados

Para cada tipo de cáncer (LUAD/LUSC) y combinación de umbrales se
generan:

-   Tabla completa de todos los genes analizados (`*_all_genes.txt`).
-   Tabla de DEGs significativos (`*_DEGs.txt`).
-   Inputs para todas las herramientas de reposicionamiento.
-   Tabla de fármacos consenso con scores de cada herramienta
    (`Farmacos_Consenso_*.csv`).
