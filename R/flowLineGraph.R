#' flowLineGraph
#'
#' A function that permits the users to plot multiple flow frames, either all
#' on the same graph panel or separately on a grid. Each colourful line on the
#' same panel represents a different flow frame, and the black line represents
#' the control sample used for comparison (if indicated by the user). If plotted
#' on a grid, each flow frame is the black line, with the optional control
#' sample as the grey line. Annotations as well as vertical lines can be added
#' to each graph in the grid. If the user chooses to only plot a control sample,
#' they can also specify annotations and vertical lines, but the 'grid'
#' parameter will be ignored.
#'
#' @param flowDir The directory of the .fcs files. If NA, then a window
#' will prompt the user to select a folder (string file path)
#' @param flowControl The control sample (string name of the .fcs file)
#' @param flowSamples A vector of samples that the user wants to analyze if they
#'  do not want to analyze the full dataset; default value is NA (vector
#'  containing the names of each file in flowDir as a string)
#' @param xVariable The fluorescence channel on the x-axis (determined by the
#'  set-up of your flow cytometer; ex. FLA-1, FITC-A)
#' @param flowColours Vector of colours for the flowSample to be plotted (only
#'  used if grid = FALSE; should be a vector of strings of the colour names)
#' @param grid Whether you want to plot samples individual samples in a grid or
#'  all together on one graph panel (TRUE or FALSE)
#' @param annotations The vector of annotations to plot on each graph in the
#'  grid where the first entry is the label for the annotation (vector containing
#'  strings of the values you want to display)
#' @param vertLine1 The vector of values for a vertical line to be plotted on
#'  each graph in the grid (vector containing a numerical value for each sample)
#' @param vertLine2 The vector of values for a second vertical line to be
#'  plotted on each graph in the grid (vector containing a numerical value for
#'  each sample)
#' @param fileName To label the PDF file if plotting as a grid (string)
#' @param samplePeaks The matrix of G1 and G2 peak values for each sample (NOTE:
#'  this is only used in the `RSEOutlierDetection()` function and should not be
#'  utilized if strictly using the main function)
#' 
#' @importFrom tcltk tclvalue tkchooseDirectory
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual theme_bw xlab ylab
#' @import vctrs
#' @import gridExtra
#' @import cowplot
#' @import flowTime
#' 
#' @return either a single line graph for the samples or a PDF file of samples
#' plotted in a grid
#' @export
#'
#' @examples
#' library(BiocFileCache)
#' original_wd <- getwd()
#' cachePath <- tempfile(pattern = "bioc_cache_")
#' dir.create(cachePath, showWarnings = FALSE)
#' bfc <- BiocFileCache(cachePath, ask = FALSE)
#' on.exit({
#'  if(normalizePath(cachePath) == normalizePath(getwd())){
#'    setwd(tempdir())  
#'  }
#'  unlink(cachePath, recursive = TRUE, force = TRUE)
#' }, add = TRUE)
#' gatedFolder <- file.path(cachePath, "fcs_data")
#' dir.create(gatedFolder, showWarnings = FALSE)
#' urls <- c(
#'   "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_1.fcs",
#'   "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_2.fcs",
#'   "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_3.fcs",
#'   "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_4.fcs",
#'   "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_5.fcs",
#'   "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_6.fcs"
#' )
#' rids <- vapply(urls, function(url){
#'   res <- bfcquery(bfc, url)
#'   if(nrow(res) == 0){
#'     bfcadd(bfc, basename(url), url)
#'   }else{
#'     res$rid[1]
#'   }
#' }, character(1))
#' fpaths <- bfcrpath(bfc, urls)
#' newFPaths <- file.path(gatedFolder, basename(urls))
#' file.copy(fpaths, newFPaths, overwrite = TRUE)
#' downloadedFiles <- newFPaths
#' flowLineGraph(
#'  flowControl = "A1_1.fcs",
#'  flowSamples = c("A1_2.fcs", "A1_3.fcs", "A1_4.fcs", "A1_5.fcs", "A1_6.fcs"),
#'  flowColours = NA,
#'  xVariable = "FITC-A",
#'  flowDir = gatedFolder,
#'  grid = TRUE,
#'  annotations = c("RSE", 12.41, 13.16, 13.45, 8.15, 4.21),
#'  vertLine1 = c(150, 155, 150, 145, 140),
#'  vertLine2 = c(335, 330, 340, 330, 305),
#'  fileName = "Gated_Grid"
#' )

flowLineGraph = function(
    flowDir = NA,
    flowControl = NA,
    flowSamples = NA,
    xVariable = "FL1-A",
    flowColours = NA,
    grid = FALSE,
    annotations = NA,
    vertLine1 = NA,
    vertLine2 = NA,
    fileName = "Samples",
    samplePeaks = NA
){
    ##Removing NOTE 'no visible binding for global variable'
    x<-y<-Data<-NULL
    if(is.na(flowDir)){
      getwd()
      flowDir <- tclvalue(tkchooseDirectory())
    }
    ##Checking if the folder they selected is empty
    if(purrr::is_empty(flowDir)){
        stop("Your directory is empty")
    }
    ##Checking if samplePeaks was filled & FLG was called directly
    calledByRSE <- any(sapply(sys.calls(), function(call){deparse(call)[[1]] == "RSEOutlierDetection"}))
    if (!is.null(samplePeaks) && "samplePeaks" %in% names(match.call()) && calledByRSE){
      warning("Parameter 'samplePeaks' should not be filled when calling
              flowLineGraph directly. Setting it to NA.")
      samplePeaks <- NA
    }
    ##samples dataset
    if(TRUE %in% !is.na(flowSamples)){   # sample set exists
        if(length(flowSamples) > 1){     # number of samples is greater than 1
            sampleDs <- c()
            for(k in seq_len(length(flowSamples))){
                if(!flowSamples[k] %in% list.files(flowDir)){
                        errorMsg<-paste0("The flow frame ",flowSamples[k]," is not
                        in the folder, check on the spelling of flowName and/or   
                        make sure you selected the proper folder")
                        stop(errorMsg)
                        rm(errorMsg)
                }
                xpectr::suppress_mw(
                        flowName <- flowCore::read.FCS(
                            paste0(flowDir, "/", flowSamples[k]),
                            transformation=FALSE
                        )
                )
                ##Checking to see if the user input the correct X variable
                if(!xVariable %in% flowName@parameters@data$name){
                        stop("Your X variable is not in the dataset")
                }
                ds <- .smoothData( flowName, xVariable, 5)
                ds$Data <- flowSamples[k]
                sampleDs <- rbind(
                    sampleDs,
                    ds
                )
            }
        }else{    # number of samples is 1
            xpectr::suppress_mw(
                    flowName <- flowCore::read.FCS(
                        paste0(flowDir, "/", flowSamples), transformation=FALSE
                    ) 
            )
            ##Checking to see if the user input the correct X variable
            if(!xVariable %in% flowName@parameters@data$name){
                stop("Your X variable is not in the dataset")
            }
            sampleDs <- .smoothData(flowName, xVariable, 5)
            sampleDs$Data <- flowSamples
        }
        if(!is.na(flowControl)){     # control sample is given
            if(!flowControl %in% list.files(flowDir)){
                errorMsg<-paste0(
                  "The flow frame ",flowControl," is not in the folder, 
                  check on the spelling of flowName and/or make sure you
                  selected the proper folder"
                )
                stop(errorMsg)
                rm(errorMsg)
            }
            ##control dataset
            xpectr::suppress_mw(
              flowNameControl <- flowCore::read.FCS(
                paste0(flowDir, "/", flowControl), transformation=FALSE
              ) 
            )
            ##Checking to see if the user input the correct X variable
            if(!xVariable %in% flowNameControl@parameters@data$name){
                stop("Your X variable is not in the dataset")
            }
            controlDs <- .smoothData(flowNameControl, xVariable, 5)
            controlDs$Data <- flowControl
            if(grid == FALSE){
              if(is.na(flowColours[1])){
                  ##plotting
                  flowPlot <- ggplot() +
                    geom_line(data=controlDs, aes(x=x, y=y, group=2), linewidth = 1,
                              color='gray62') +
                    geom_line(data=sampleDs, aes(x=x, y=y, group=Data,
                                                 color = Data)) + 
                    ylab("Counts") + xlab(xVariable) + theme_bw()
              }else{
                  if(length(flowColours) != length(unique(sampleDs$Data))){
                    stop(
                      "'flowColours' vector length does not match the number of
                      unique samples in 'flowSamples'"
                    )
                  }
                  ##plotting
                  flowPlot <- ggplot() +
                    geom_line(data=controlDs, aes(x=x, y=y, group=2), linewidth = 1,
                              color='gray62') + 
                    geom_line(data=sampleDs, aes(x=x, y=y, group=Data, 
                                                 color = Data)) + 
                    ylab("Counts") + xlab(xVariable) + 
                    scale_color_manual(values=flowColours) + theme_bw()
              }
            }else{    # grid function w/ control
              if(TRUE %in% !is.na(annotations)){    # annotations != NA
                if(length(annotations) != (length(flowSamples) + 1)){
                  stop("'annotations' vector length does not match the number of
                       samples in 'flowSamples' + 1")
                }else{   # length is good, check for first string
                  if(!is.character(annotations[1])){
                    stop("First item in 'annotations' vector is not a string")
                  }
                }
              }
              if(is.na(fileName) || !is.character(fileName)){ # NA or not a string
                stop("Invalid fileName (must be a string)")
              }else{
                fileName <- sub(" ", "_", fileName)  # clean up (no white space)
              }
              if(TRUE %in% !is.na(vertLine1)){
                if(length(vertLine1) != length(unique(sampleDs$Data))){
                  stop("'vertLine1' vector length does not match the number of
                       unique samples in 'flowSamples'")
                }
              }
              if(TRUE %in% !is.na(vertLine2)){
                if(length(vertLine2) != length(unique(sampleDs$Data))){
                  stop("'vertLine2' vector length does not match the number of
                       unique samples in 'flowSamples'")
                }
              }
              # plotting
              pdfFile <- paste0(fileName, "_plotGrid")
              plots <- .plotSamples(sampleDs, controlDs, xVariable, annotations, vertLine1, vertLine2, samplePeaks)
              flowPlot <- .createPDF(flowDir, pdfFile, plots)
            }
        }else{    # no control sample is given
            if(grid == FALSE){
              if(length(flowSamples) > 1){  # number of samples is greater than 1
                  if(is.na(flowColours[1])){
                      ##plotting
                      flowPlot <- ggplot() + 
                        geom_line(data=sampleDs,aes(x=x, y=y, group=Data,
                                                    color = Data)) + 
                        ylab("Counts") + xlab(xVariable)+theme_bw()
                  }else{
                      if(length(flowColours) != length(unique(sampleDs$Data))){
                        stop("'flowColours' vector length does not
                                  match the number of unique samples in
                                  'flowSamples'")
                        }
                        ##plotting
                        flowPlot <- ggplot() +
                          geom_line(data=sampleDs, aes(x=x, y=y, group=Data,
                                                       color = Data)) +
                          ylab("Counts") + xlab(xVariable) +
                          scale_color_manual(values=flowColours) + theme_bw()
                      }
              }else if(length(flowSamples) == 1){
                  if(is.na(flowColours[1])){
                      ##plotting
                      flowPlot <- ggplot() +
                              geom_line(data=sampleDs, aes(x=x, y=y), linewidth = 1,
                                        color='black') + ylab("Counts") +
                        xlab(xVariable) + theme_bw()
                  }else{
                      if(length(flowColours) != length(unique(sampleDs$Data))){
                        stop("'flowColours' vector length does not 
                             match the number of unique samples in 'flowSamples'")
                      }
                      ##plotting
                      flowPlot <- ggplot() +
                          geom_line(data=sampleDs, aes(x=x, y=y, group=Data,
                                                       color = Data)) +
                        ylab("Counts") + xlab(xVariable) +
                        scale_color_manual(values=flowColours) + theme_bw()
                  }
              }
            }else{   # grid function w/out control
              if(TRUE %in% !is.na(annotations)){    # annotations != NA
                if(length(annotations) != (length(flowSamples) + 1)){
                  stop("'annotations' vector length does not match the number of
                       samples in 'flowSamples' + 1")
                }else{   # length is good, check for first string
                  if(!is.character(annotations[1])){
                    stop("First item in 'annotations' vector is not a string")
                  }
                }
              }
              if(is.na(fileName) || !is.character(fileName)){ # NA or not a string
                stop("Invalid fileName (must be a string)")
              }else{
                fileName <- sub(" ", "_", fileName)  # clean up (no white space)
              }
              if(TRUE %in% !is.na(vertLine1)){
                if(length(vertLine1) != length(unique(sampleDs$Data))){
                  stop("'vertLine1' vector length does not match the number of
                       unique samples in 'flowSamples'")
                }
              }
              if(TRUE %in% !is.na(vertLine2)){
                if(length(vertLine2) != length(unique(sampleDs$Data))){
                  stop("'vertLine2' vector length does not match the number of
                       unique samples in 'flowSamples'")
                }
              }
              # plotting
              pdfFile <- paste0(fileName, "_plotGrid")
              plots <- .plotSamples(sampleDs, NA, xVariable, annotations, vertLine1, vertLine2, samplePeaks)
              flowPlot <- .createPDF(flowDir, pdfFile, plots)
            }
        }
    }else{    # sample set does not exist 
        if(!flowControl %in% list.files(flowDir)){
            errorMsg<-paste0("The flow frame ",flowControl," is not in
                                the folder, check on the spelling of
                                flowName and/or make sure you selected
                                the proper folder")
            stop(errorMsg)
            rm(errorMsg)
        }
        ##control dataset
        xpectr::suppress_mw(flowNameControl <- flowCore::read.FCS(
          paste0(flowDir, "/", flowControl), transformation=FALSE))
        ##Checking to see if the user input the correct X variable
        if(!xVariable %in% flowNameControl@parameters@data$name){
              stop("Your X variable is not in the dataset")
        }
        controlDs <- .smoothData(flowNameControl, xVariable, 5)
        controlDs$Data <- flowControl
        ##plotting
        if(TRUE %in% !is.na(annotations)){
          if(length(annotations) != 2){
            stop("'annotations' vector length is not equal to 2")
          }else{   # length is good, check for first string
            if(!is.character(annotations[1])){
              stop("First item in 'annotations' vector is not a string")
            }
          }
        }
        if(TRUE %in% !is.na(vertLine1)){
          if(length(vertLine1) != 1){
            stop("'vertLine1' vector length is not equal to 1")
          }
        }
        if(TRUE %in% !is.na(vertLine2)){
          if(length(vertLine2) != 1){
            stop("'vertLine2' vector length is not equal to 1")
          }
        }
        flowPlot <- .plotSamples(NA, controlDs, xVariable, annotations, vertLine1, vertLine2, samplePeaks)
    }
    return(flowPlot)
}
.smoothData = function(flowDs, xVariable, smoothLevel){
    ##Get counts and breaks from histogram
    histData <- hist(flowDs@exprs[, xVariable], breaks=256, plot=FALSE)
    ##Data that will be smoothed
    data <- histData$counts
    ##Apply smoothing to the counts with 'smoothLevel'
    smoothedDs <- zoo::rollmean(data, k=smoothLevel, fill=0)
    ##Create data frame with smoothed data
    ds <- data.frame(
            x=histData$breaks,
            y=c(0, smoothedDs)
    )
    return(ds)
}

.createPlotGrid = function(plots, start, end, plotWidths = NULL, plotHeights = NULL){
  ##Creates a grid of plots
  grid.arrange(grobs = plots[start:end], ncol = 3, widths = plotWidths, heights = plotHeights)
}

.splitString = function(inputString, delimiter){
  ##Splits string if specified delimiter exists, otherwise return original string
  if(grepl(delimiter, inputString)) {
    result <- unlist(strsplit(inputString, delimiter))[1]
  }else{
    result <- inputString
  }
  return(result)
}

.plotSamples = function(samples, control, xVariable, annotations, vertLine1, vertLine2, samplePeaks){
  P <- list()
  slot <- 1
  if((TRUE %in% is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){   # check if vertLine1 is NA but vertLine2 is not
    vertLine1 <- vertLine2
  }
  if(TRUE %in% !is.na(samples)){      # there are samples
    if(TRUE %in% !is.na(control)){    # there is a control sample
      if(TRUE %in% !is.na(annotations)){   # there are annotations
        if((TRUE %in% !is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){   # two vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            yVal <- (max(temp$y, control$y) * 0.9)
            sampName <- .splitString(samp, "[.]")      # clean up name
            p <- ggplot() +
              geom_line(data = control, aes(x = x, y = y, group = 2),
                        linewidth = 1, colour = 'gray62') +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              geom_vline(xintercept = vertLine2[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw() +
              annotate("text", x = (max(temp$x) * 0.85), y = yVal,
                       label = paste0(paste0(annotations[1], ":"),
                                     annotations[slot + 1]), colour = "#DF536B",
                       size = 3) 
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
          }
        }else if(TRUE %in% !is.na(vertLine1)){   # one vertical line
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            yVal <- (max(temp$y, control$y) * 0.9)
            sampName <- .splitString(samp, "[.]")     
            p <- ggplot() +
              geom_line(data = control, aes(x = x, y = y, group = 2),
                        linewidth = 1, colour = 'gray62') +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw() +
              annotate("text", x = (max(temp$x) * 0.85), y = yVal,
                       label = paste0(paste0(annotations[1], ":"),
                                     annotations[slot + 1]), colour = "#DF536B",
                       size = 3) 
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
          }
        }else{   # no vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            yVal <- (max(temp$y, control$y) * 0.9)
            sampName <- .splitString(samp, "[.]")     
            p <- ggplot() +
              geom_line(data = control, aes(x = x, y = y, group = 2),
                        linewidth = 1, colour = 'gray62') +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw() +
              annotate("text", x = (max(temp$x) * 0.85), y = yVal,
                       label = paste0(paste0(annotations[1], ":"),
                                     annotations[slot + 1]), colour = "#DF536B",
                       size = 3) 
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
          }
        }
      }else{    # no annotations
        if((TRUE %in% !is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){   # two vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")     
            p <- ggplot() +
              geom_line(data = control, aes(x = x, y = y, group = 2),
                        linewidth = 1, colour = 'gray62') +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              geom_vline(xintercept = vertLine2[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw()
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
          }
        }else if(TRUE %in% !is.na(vertLine1)){   # one vertical line
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")
            p <- ggplot() +
              geom_line(data = control, aes(x = x, y = y, group = 2),
                        linewidth = 1, colour = 'gray62') +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw()
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
          }
        }else{   # no vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")       # cleans up file name
            p <- ggplot() +
              geom_line(data = control, aes(x = x, y = y, group = 2),
                        linewidth = 1, colour = 'gray62') +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw()
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
          }
        }
      }
    }else{      # no control sample
      if(TRUE %in% !is.na(annotations)){   # there are annotations
        if((TRUE %in% !is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){   # two vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")   
            p <- ggplot() +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              geom_vline(xintercept = vertLine2[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw() +
              annotate("text", x = (max(temp$x) * 0.85), y = (max(temp$y) * 0.9),
                       label = paste0(paste0(annotations[1], ":"),
                                     annotations[slot + 1]), colour = "#DF536B",
                       size = 3) 
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
            #P <- c(P, list(p))
          }
        }else if(TRUE %in% !is.na(vertLine1)){   # one vertical line
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")   
            p <- ggplot() +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw() +
              annotate("text", x = (max(temp$x) * 0.85), y = (max(temp$y) * 0.9),
                       label = paste0(paste0(annotations[1], ":"),
                                     annotations[slot + 1]), colour = "#DF536B",
                       size = 3) 
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
            #P <- c(P, list(p))
          }
        }else{   # no vertical lines
          if(TRUE %in% !is.na(samplePeaks)){   # RSEOutlierDetection called FLG
            spTrans <- data.frame(t(samplePeaks))
            colours <- vec_rep_each(c("#2297E6", "#61D04F"), length(samplePeaks)/2)
            for(samp in unique(samples$Data)){
              temp <- samples[samples$Data == samp, ]
              sampName <- .splitString(samp, "[.]")
              xCol <- spTrans[, slot]
              p <- ggplot() +
                geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                          colour = 'black') +
                geom_point(data = spTrans, aes(x = xCol, y = 0,
                                               colour = colours), size = 2,na.rm = TRUE) +
                labs(x = xVariable, y = "Counts", title = sampName) +
                scale_colour_identity() + 
                theme_bw() +
                annotate("text", x = (max(temp$x) * 0.85), y = (max(temp$y) * 0.9),
                         label = paste0(paste0(annotations[1], ":"),
                                       annotations[slot + 1]), colour = "#DF536B",
                         size = 3) 
              P[[slot]] <- ggplot_gtable(ggplot_build(p))
              slot <- slot + 1
            }
          }else{      # regularly called by FLG
            for(samp in unique(samples$Data)){
              temp <- samples[samples$Data == samp, ]
              sampName <- .splitString(samp, "[.]")     
              p <- ggplot() +
                geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                          colour = 'black') +
                labs(x = xVariable, y = "Counts", title = sampName) +
                theme_bw() +
                annotate("text", x = (max(temp$x) * 0.85), y = (max(temp$y) * 0.9),
                         label = paste0(paste0(annotations[1], ":"),
                                       annotations[slot + 1]), colour = "#DF536B",
                         size = 3) 
              P[[slot]] <- ggplot_gtable(ggplot_build(p))
              slot <- slot + 1
              #P <- c(P, list(p))
            }
          }
        }
      }else{    # no annotations
        if((TRUE %in% !is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){   # two vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")     
            p <- ggplot() +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              geom_vline(xintercept = vertLine2[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw() 
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
            #P <- c(P, list(p))
          }
        }else if(TRUE %in% !is.na(vertLine1)){   # one vertical line
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")     
            p <- ggplot() +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              geom_vline(xintercept = vertLine1[slot], color = "#2297E6",
                         linetype = "dashed") +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw()
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
            #P <- c(P, list(p))
          }
        }else{   # no vertical lines
          for(samp in unique(samples$Data)){
            temp <- samples[samples$Data == samp, ]
            sampName <- .splitString(samp, "[.]")       # cleans up file name
            p <- ggplot() +
              geom_line(data = temp, aes(x = x, y = y), linewidth = 1,
                        colour = 'black') +
              labs(x = xVariable, y = "Counts", title = sampName) +
              theme_bw()
            P[[slot]] <- ggplot_gtable(ggplot_build(p))
            slot <- slot + 1
            #P <- c(P, list(p))
          }
        }
      }
    }
  }else{   # no samples, only control
    if(TRUE %in% !is.na(annotations)){   # annotations
      if((TRUE %in% !is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){     # two vertical lines
        controlName <- .splitString(unique(control$Data), "[.]")       # cleans up file name
        P <- ggplot() +
          geom_line(data=control, aes(x=x, y=y, group=2), linewidth = 1,
                    color='black') +
          geom_vline(xintercept = vertLine1[1], color = "#2297E6",
                     linetype = "dashed") +
          geom_vline(xintercept = vertLine2[1], color = "#2297E6",
                     linetype = "dashed") +
          labs(x = xVariable, y = "Counts", title = controlName) +
          theme_bw() +
          annotate("text", x = (max(control$x) * 0.85),
                   y = (max(control$y) * 0.9),
                   label = paste0(paste0(annotations[1], ":"), annotations[2]),
                   colour = "#DF536B", size = 3)
      }else if(TRUE %in% !is.na(vertLine1)){                    # one vertical line
        controlName <- .splitString(unique(control$Data), "[.]")       # cleans up file name
        P <- ggplot() +
          geom_line(data=control, aes(x=x, y=y, group=2), linewidth = 1,
                    color='black') +
          geom_vline(xintercept = vertLine1[1], color = "#2297E6",
                     linetype = "dashed") +
          labs(x = xVariable, y = "Counts", title = controlName) +
          theme_bw() +
          annotate("text", x = (max(control$x) * 0.85),
                   y = (max(control$y) * 0.9),
                   label = paste0(paste0(annotations[1], ":"), annotations[2]),
                   colour = "#DF536B", size = 3)
      }else{                                          # no vertical lines
        controlName <- .splitString(unique(control$Data), "[.]")       # cleans up file name
        P <- ggplot() +
          geom_line(data=control, aes(x=x, y=y, group=2), linewidth = 1,
                    color='black') +
          labs(x = xVariable, y = "Counts", title = controlName) +
          theme_bw() +
          annotate("text", x = (max(control$x) * 0.85),
                   y = (max(control$y) * 0.9),
                   label = paste0(paste0(annotations[1], ":"), annotations[2]),
                   colour = "#DF536B", size = 3)
      }
    }else{                     # no annotations
      if((TRUE %in% !is.na(vertLine1)) && (TRUE %in% !is.na(vertLine2))){     # two vertical lines
        controlName <- .splitString(unique(control$Data), "[.]")       # cleans up file name
        P <- ggplot() +
          geom_line(data=control, aes(x=x, y=y, group=2), linewidth = 1,
                    color='black') +
          geom_vline(xintercept = vertLine1[1], color = "#2297E6",
                     linetype = "dashed") +
          geom_vline(xintercept = vertLine2[1], color = "#2297E6",
                     linetype = "dashed") +
          labs(x = xVariable, y = "Counts", title = controlName) +
          theme_bw()
      }else if(TRUE %in% !is.na(vertLine1)){                    # one vertical line
        controlName <- .splitString(unique(control$Data), "[.]")       # cleans up file name
        P <- ggplot() +
          geom_line(data=control, aes(x=x, y=y, group=2), linewidth = 1,
                    color='black') +
          geom_vline(xintercept = vertLine1[1], color = "#2297E6",
                     linetype = "dashed") +
          labs(x = xVariable, y = "Counts", title = controlName) +
          theme_bw()
      }else{                                          # no vertical lines
        controlName <- .splitString(unique(control$Data), "[.]")       # cleans up file name
        P <- ggplot() +
          geom_line(data=control, aes(x=x, y=y, group=2), linewidth = 1,
                    color='black') +
          labs(x = xVariable, y = "Counts", title = controlName) +
          theme_bw()
      }
    }
  }
  return(P)
}

.createPDF = function(flowDir, pdfFileName, allPlots){
  plotsPerPage <- 12    # set the # of plots per page
  totalPlots <- length(allPlots)
  totalPages <- ceiling(totalPlots / plotsPerPage)  # calculate number of pages 
  pdf(paste0(dirname(flowDir), pdfFileName, ".pdf"), width = 8.5, height = 11)  # open a pdf file
  for(page in 1:totalPages){     # loop through pages & save to pdf file
    start <- (page - 1) * plotsPerPage + 1
    end <- min(page * plotsPerPage, totalPlots)
    
    # create grid for current page
    currentGrid <- .createPlotGrid(allPlots, start, end,
                                  plotWidths = rep(unit(2.75, "inches"), 3),
                                  plotHeights = rep(unit(2.7, "inches"), 4))
    print(currentGrid)    # plot current page
    if(page < totalPages){    # add page break (if not last page)
      cat("\\newpage\n")
    }
  }
  # close pdf device
  dev.off()
  return("Your PDF file has been created")
}
