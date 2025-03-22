## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
Sys.setenv(AZURE_AUTH_USE_FILE = "FALSE")

library(PloidyPeaks)
library(flowCore)
library(here)
library(AzureStor)
library(AzureRMR)
library(BiocFileCache)

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("PloidyPeaks")

## ----warning = FALSE, eval=TRUE, echo=TRUE------------------------------------
library(PloidyPeaks)

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# install.packages("devtools", dependencies = TRUE)

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# library(devtools)
# install_github("MicroStatsLab/PloidyPeaks", dependencies = TRUE,
#                build_vignettes = TRUE)

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# install.packages("AzureStor")
# install.packages("AzureRMR")
# 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("BiocFileCache")
# 
# library(AzureStor)
# library(AzureRMR)
# library(BiocFileCache)

## ----warning = FALSE, eval=TRUE, echo=TRUE------------------------------------
library(BiocFileCache)
original_wd <- getwd()
cachePath <- tempfile(pattern = "bioc_cache_")
dir.create(cachePath, showWarnings = FALSE)
bfc <- BiocFileCache(cachePath, ask = FALSE)
on.exit({
   if (normalizePath(getwd()) != normalizePath(original_wd)) {
     setwd(original_wd) 
   }
}, add = TRUE)
rawFolder <- file.path(cachePath, "fcs_data")
dir.create(rawFolder, showWarnings = FALSE, recursive = TRUE)
rawUrls <- c(
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_1.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_2.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_3.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_4.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_5.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_6.fcs"
)
downloadedFiles <- character(length(rawUrls))
for (i in seq_along(rawUrls)) {
  rid <- bfcquery(bfc, rawUrls[i])$rid
  if (length(rid) == 0){
    rid <- bfcadd(bfc, paste0("raw_data_", i), rawUrls[i])
  }
  fpath <- bfcrpath(bfc, rawUrls[i])
  newFPath <- file.path(rawFolder, basename(rawUrls[i]))  
  file.copy(fpath, newFPath, overwrite = TRUE)
  downloadedFiles[i] <- newFPath
}
gatedFolder <- file.path(cachePath, "fcs_data")
dir.create(gatedFolder, showWarnings = FALSE, recursive = TRUE)
gatedUrls <- c(
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_1.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_2.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_3.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_4.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_5.fcs",
    "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_6.fcs"
)
downloadedFiles <- character(length(gatedUrls))
for (i in seq_along(gatedUrls)) {
  fpath <- bfcrpath(bfc, gatedUrls[i])
  newFPath <- file.path(gatedFolder, basename(gatedUrls[i]))  
  file.copy(fpath, newFPath, overwrite = TRUE)
  downloadedFiles[i] <- newFPath
}

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# rectGateFlowFrame(
#   rawDir = NA,
#   flowName = "Sample1",
#   xVariable = "FL1-A",
#   yVariable = "SSC-A",
#   xMinValue = 30000,
#   xMaxValue = 700000,
#   yMinValue = 30000,
#   yMaxValue = 750000,
#   savePlot = TRUE
# )

## ----echo=FALSE, fig.pos = "H", out.width="80%",fig.align="center"------------
knitr::include_graphics(here("vignettes/images/V2_gated2.png"))

## ----echo=FALSE, fig.pos = "H", out.width="80%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/findPairs.png"))

## ----echo=FALSE, fig.pos = "H", out.width="80%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/doublets.jpeg"))

## ----echo=FALSE, fig.pos = "H", out.width="80%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/cellProp.png"))

## ----echo=FALSE, fig.pos = "H", out.width="80%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/wrapper_workflow.png"))

## ----echo=FALSE, fig.pos = "H", out.width="60%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/20220807_workflow.png"))

## ----echo=FALSE, fig.pos = "H", out.width="80%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/flowPeakDetectionOutput.png"))

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# peakCorrection(
#   xVariable = "FL1-A",
#   flowDir,
#   sampleName,
#   numSubPop,
#   savePlot
# )

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# flowLineGraph(
#   flowControl = NA,
#   flowSamples,
#   flowColours = NA,
#   xVariable = "FL1-A",
#   flowDir,
#   grid = FALSE,
#   annotations = c("RSE", 10.87, 6.43, 34.38, 29.87, 4.54),
#   vertLine1 = c(180, 178, 156, 197, 195),
#   vertLine2 = c(254, 245, 234, 286, 312),
#   fileName = "Gated_Grid",
#   samplePeaks = NA
# )

## ----echo=FALSE, fig.pos = "H", out.width="50%", out.height="20%", fig.show='hold', fig.align="center"----
knitr::include_graphics(here("vignettes/images/flowLineGraph.png"))

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# massFlowGraph(
#   xVariable = "FL1-A",
#   flowDir,
#   filePath,
#   fileName = "Gated_Grid"
# )

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# RSEOutlierDetection(
#     xVariable = "FL1-A",
#     flowDir,
#     filePath,
#     fileName = "RSEOutlier",
#     alpha = 0.05
# )

## ----echo=FALSE, fig.pos = "H", out.width="30%", fig.align="center"-----------
image1 <- here("vignettes/images/RSE_Vignette_plot.jpeg")
image2 <- here("vignettes/images/RSECol.png")
image3 <- here("vignettes/images/RSEGrid.png")

html_code <- sprintf('<div style="text-align: center; margin-bottom: 10px;">
    <img src="%s" style="width: 35%%; height: auto; margin-right: 10px;" />
    <img src="%s" style="width: 30%%; height: auto;" />
    </div>
    <div style="text-align: center;">
    <img src="%s" style="width: 90%%; height: auto;" />
    </div>', image1, image2, image3)

htmltools::HTML(html_code)

## ----echo=FALSE, results='asis', table.pos='H'--------------------------------
library(magrittr)
library(knitr)
library(kableExtra)
data <- data.frame(
  Function = c("rectGateFlowFrame", "rectGateFlowSet", "flowPeakDetection", "peakCorrection", "flowLineGraph", "massFlowGraph", "RSEOutlierDetection"),
  Outputs = c("Single .fcs file of gated data and plots of the gated data (before and after gating).", "Set of .fcs files of gated data, plots of the gated data (before and after gating), and .csv file (named experimentName_percentOfCellsGatedOut). The user can use this file as an indication of if a sample is messy (a significantly greater proportion of cells gated out than the other samples), if the chosen gating parameters were chosen properly (there were no cells gated out), or the proportion was significantly large for all.", "A .csv file (named experimentName_ploidyPeaksOutput). This file contains all of the information from the peak detection analysis. Each row is a different sample. Column information is peak means ($G_1$ and $G_2$ peaks), peak heights, and includes the number of identified subpopulations. The user can also specify if they want similar information on doublets, whether to include an indicator of whether the sample should be investigated further by the user, which peak algorithm was used to analyze each sample, and the residual standard error value from the models.", "A .csv file (named sample_ploidyPeaksOutput) and line graph of the sample with identified peaks. The information in the *.csv* file is the same as the peak detection analysis.", "A single line graph in the console or a PDF file of multiple samples in a grid.", "PDF file of samples as determined by the user (all samples, flagged samples, or doublet flagged samples). Samples are plotted in a grid using `flowLineGraph()`.", "Provides up to three things to the user. First is a histogram of the distribution of RSE values for a set of samples. Second is a new column in the *.csv* file identifying outliers based on provided threshold. Third is a *PDF* file of the plotted outlier samples.")
)

# Print the data frame using kableExtra
kable(data, format = "html") %>% 
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
  column_spec(1, width = "20%") %>%
  column_spec(2, width = "80%")

## ----echo=FALSE, fig.pos = "H", out.width="50%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/folder_structure.png"))

## ----warning = FALSE, eval=TRUE, echo=TRUE------------------------------------
rectGateFlowFrame(
  rawDir = rawFolder,
  flowName = "A1_1.fcs",
  xVariable = "FITC-A",
  yVariable = "SSC-A",
  xMinValue = 50,
  xMaxValue = 800,
  yMinValue = 50,
  yMaxValue = 800,
  savePlot = TRUE
)

## ----echo=FALSE, fig.pos = "H", out.width="50%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/A1_1-propergate.png"))

## ----warning = FALSE, eval=TRUE, echo=TRUE------------------------------------
rectGateFlowSet(
  rawDir = rawFolder,
  xVariable = "FITC-A",
  yVariable = "SSC-A",
  xMinValue = 50,
  xMaxValue = 800,
  yMinValue = 50,
  yMaxValue = 800,
  savePlot = TRUE
)

## ----warning=FALSE, eval=TRUE, out.width="75%", echo=FALSE--------------------
flowLineGraph(flowSamples = c("A1_1.fcs", "A1_2.fcs", "A1_3.fcs", "A1_4.fcs",
                              "A1_5.fcs", "A1_6.fcs"), flowColours =
                c("#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC",
                  "#F5C710"), xVariable = "FITC-A",
              flowDir = gatedFolder, grid = FALSE)

## ----echo=FALSE, fig.pos = "H", out.width="75%", fig.align="center"-----------
knitr::include_graphics(here("vignettes/images/prelimVis.jpeg"))

## ----warning = FALSE, eval=TRUE, echo=TRUE------------------------------------
flowPeakDetection(
  flowDir = gatedFolder,
  xVariable = "FITC-A",
  doublet = FALSE,
  singleThreshold = 8.5,
  usedCellsThreshold = 86
)  

## ----warning = FALSE, eval=TRUE, echo=TRUE------------------------------------
peakCorrection(xVariable = "FITC-A",
               flowDir = gatedFolder,
               sampleName = "A1_4.fcs", numSubPop = 2)

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# massFlowGraph(
#   xVariable = "FITC-A",
#   flowDir = gatedFolder,
#   filePath = paste0(dirname(gatedFolder), "/analysis", "/BiocFileCache_ploidyPeaksOutput.csv"),
#   fileName = "allSamples"
# )

## ----warning = FALSE, eval=FALSE, echo=TRUE-----------------------------------
# RSEOutlierDetection(xVariable = "FITC-A",
#                     flowDir = gatedFolder,
#                     filePath = paste0(dirname(gatedFolder), "/analysis", "/BiocFileCache_ploidyPeaksOutput.csv"),
#                     fileName = "RSEOutlier",
#                     alpha = 0.05)

## ----warning = FALSE, echo=TRUE-----------------------------------------------
sessionInfo()

