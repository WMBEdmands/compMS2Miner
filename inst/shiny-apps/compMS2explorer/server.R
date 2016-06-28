library(CompMS2miner)
library(shiny)
library(igraph)
library(rhandsontable)

shiny::shinyServer(function(input,  output, session){
  # comment lines below if action button is used to commit changes
  values = reactiveValues()
  setHot = function(x) values[["hot"]] = x
  
  
  #options(shiny.trace=T)
  observe({
    if(input$CloseAppBtn > 0){
      if(!is.null(values[["hot"]])){
        object@Comments <- values[["hot"]]
      }
      shiny::stopApp(object)
    }
  })
  # feature selection shiny::reactive
  Featureselection <- shiny::reactive({
    if(input$goButton == 0){ 
      return() } else if (input$goButton > 0){
        
        if(input$specDBMatch == T){
          FeaturesIndx <- rep(F, length(Features.v))
          FeaturesIndx[indxSpectralDb] <- T
        } else {
        FeaturesIndx <- rep(T, nrow(SubStr_types))
        NoFeaturesIndx <- rep(F, nrow(SubStr_types))
        if(any(input$NotSubStrTypes != "")){
          ###escape special characters
          NotSubStrTypes <- gsub("\\+",  "\\\\+",  input$NotSubStrTypes)
          NotSubStrTypes <- gsub("\\[",  "\\\\[",  NotSubStrTypes)
          NotSubStrTypes <- gsub("\\]",  "\\\\]",  NotSubStrTypes)
          NoFeaturesIndx <- grepl(paste(NotSubStrTypes,  collapse = "|"),  SubStr_types[,  1]) | grepl(paste(NotSubStrTypes,  collapse = "|"),  SubStr_types[,  2])
        }
        
        if (any(input$SubStrTypes != "")){
          ###escape special characters
          SubStrTypes <- gsub("\\+",  "\\\\+",  input$SubStrTypes)
          SubStrTypes <- gsub("\\[",  "\\\\[",  SubStrTypes)
          SubStrTypes <- gsub("\\]",  "\\\\]",  SubStrTypes)
          FeaturesIndx <- grepl(paste(SubStrTypes,  collapse = "|"),  SubStr_types[,  1]) | grepl(paste(SubStrTypes,  collapse = "|"),  SubStr_types[,  2])
        }
        
        FeaturesIndx[NoFeaturesIndx == T] <- F 
        
        if(any(input$subStrAnnoTypes != '')){
          if(input$subStrAnnoThresh != ''){
            likelySubStrIndx  <- sapply(subStrAnno.list, function(x) any(grepl(paste(input$subStrAnnoTypes,  collapse = "|"), x$SubStrType) & as.numeric(x$SumRelInt) >= as.numeric(input$subStrAnnoThresh)))
          } else {
        likelySubStrIndx  <- sapply(subStrAnno.list, function(x) any(grepl(paste(input$subStrAnnoTypes,  collapse = "|"), x$SubStrType)))
          }
          FeaturesIndx[likelySubStrIndx == F] <- F 
        }

        # filter by mz and Rt
        if(input$All_Features == F){
          # if values in boxes then process
        if(input$mass_to_charge != '' & input$mass_accuracy != ''){
        FeaturesIndx <-  mass.v < (as.numeric(input$mass_to_charge) + ((as.numeric(input$mass_to_charge) / 1E06) * as.numeric(input$mass_accuracy)))  & mass.v > (as.numeric(input$mass_to_charge) - ((as.numeric(input$mass_to_charge) / 1E06) * as.numeric(input$mass_accuracy))) & FeaturesIndx
        }
        if(input$retentionTime != '' & input$RTtolerance != ''){
        FeaturesIndx <-  RT.v < (as.numeric(input$retentionTime) + as.numeric(input$RTtolerance))  & RT.v > (as.numeric(input$retentionTime) - as.numeric(input$RTtolerance)) & FeaturesIndx 
        }
        }
      } 
        Feature.v.sub <- Features.v[FeaturesIndx]
       
        if(length(Feature.v.sub) == 0){
          return("No MS2 features found")
        } else {
          return(Feature.v.sub)}
      } 
  })
  # DB name match shiny::reactive
  DBFeatureselection <- shiny::reactive({
    
    if(input$DBbutton  ==  0){ return()
    } else if (input$DBbutton > 0){
      if(input$DB_match_table == 'DB Annotations'){
      FeaturesIndx <- grepl(input$DB_match_name,  DBmatches[,  1],  ignore.case=T)
      } else if (input$DB_match_table == 'Best Annotations'){
        FeaturesIndx <- grepl(input$DB_match_name,   DBBestMatches[,  1],  ignore.case=T)   
      }
      # FeaturesIndx <- ifelse(FeaturesIndx == T | grepl(input$DB_match_name,  DBmatches[,  2],  ignore.case=T), T, F)
      # 
      # DBnamesTMP <- unique(c(DBnamesTMP,  grep(input$DB_match_name,  DBmatches[,  2], 
      #                                         ignore.case=T)))
      
      NoFeaturesIndx <- rep(F, nrow(SubStr_types))
      if(any(input$NotSubStrTypes != "")){
        ###escape special characters
        NotSubStrTypes <- gsub("\\+",  "\\\\+",  input$NotSubStrTypes)
        NotSubStrTypes <- gsub("\\[",  "\\\\[",  NotSubStrTypes)
        NotSubStrTypes <- gsub("\\]",  "\\\\]",  NotSubStrTypes)
        NoFeaturesIndx <- grepl(paste(NotSubStrTypes,  collapse = "|"),  SubStr_types[,  1]) | grepl(paste(NotSubStrTypes,  collapse = "|"),  SubStr_types[,  2])
      }
      
      if (any(input$SubStrTypes != "")){
        ###escape special characters
        SubStrTypes <- gsub("\\+",  "\\\\+",  input$SubStrTypes)
        SubStrTypes <- gsub("\\[",  "\\\\[",  SubStrTypes)
        SubStrTypes <- gsub("\\]",  "\\\\]",  SubStrTypes)
        FeaturesIndx <- grepl(paste(SubStrTypes,  collapse = "|"),  SubStr_types[,  1]) | grepl(paste(SubStrTypes,  collapse = "|"),  SubStr_types[,  2])
      }
      
      FeaturesIndx[NoFeaturesIndx == T] <- F 
      
      if(any(input$subStrAnnoTypes != '')){
        if(input$subStrAnnoThresh != ''){
          likelySubStrIndx  <- sapply(subStrAnno.list, function(x) any(grepl(paste(input$subStrAnnoTypes,  collapse = "|"), x$SubStrType) & as.numeric(x$SumRelInt) >= as.numeric(input$subStrAnnoThresh)))
        } else {
          likelySubStrIndx  <- sapply(subStrAnno.list, function(x) any(grepl(paste(input$subStrAnnoTypes,  collapse = "|"), x$SubStrType)))
        }
        FeaturesIndx[likelySubStrIndx == F] <- F 
      }
      
      # filter by mz and Rt
      if(input$All_Features == F){
        # if values in boxes then process
        if(input$mass_to_charge != '' & input$mass_accuracy != ''){
          FeaturesIndx <-  mass.v < (as.numeric(input$mass_to_charge) + ((as.numeric(input$mass_to_charge) / 1E06) * as.numeric(input$mass_accuracy)))  & mass.v > (as.numeric(input$mass_to_charge) - ((as.numeric(input$mass_to_charge) / 1E06) * as.numeric(input$mass_accuracy))) & FeaturesIndx
        }
        if(input$retentionTime != '' & input$RTtolerance != ''){
          FeaturesIndx <-  RT.v < (as.numeric(input$retentionTime) + as.numeric(input$RTtolerance))  & RT.v > (as.numeric(input$retentionTime) - as.numeric(input$RTtolerance)) & FeaturesIndx 
        }
      }
      
      Feature.v.sub <- Features.v[FeaturesIndx]
      if(length(Feature.v.sub) == 0)
      {
        return("No MS2 features found")
      } else {
        return(Feature.v.sub)}
    }
  })
  # feature selection observer
  observe({ if(input$goButton){
    output$FeatureNames = shiny::renderUI({
      Featurenames <- shiny::isolate({Featureselection()})
      shiny::selectizeInput('FeatureNames',  'Choose a feature to plot :',  
                            choices=Featurenames, options=list(maxOptions=10000))
    })
    
    output$matchSummaryText <- shiny::renderText({
      return("Match summary table")
    })
    
    output$matchSummary <- shiny::renderTable({
      Featurenames <- shiny::isolate({Featureselection()})
      if(Featurenames[1]!="No MS2 features found"){
        matchSummary <- data.frame(nMatches=length(unique(gsub(".+_", "", Featurenames))), 
                                 nCompositeSpectra=length(Featurenames))
        return(matchSummary)
      } else {
        return()
      }
    })
    output$tabbedPanelButton = shiny::renderUI({
      shiny::actionButton("tabbedPanelButton", "ViewData")
    })
  }
  })
  
  # DB matches observer
  observe({ if(input$DBbutton){
    output$FeatureNames = shiny::renderUI({
      Featurenames <- shiny::isolate({DBFeatureselection()})
      shiny::selectizeInput('FeatureNames',  'Choose a feature to plot :',  
                            choices=Featurenames, options=list(maxOptions=10000))
    })
    
    
    
    output$matchSummary <- shiny::renderTable({
      Featurenames <- shiny::isolate({DBFeatureselection()})
      if(Featurenames[1]!="No MS2 features found"){
        matchSummary <- data.frame(nMatches=length(unique(gsub(".+_", "", Featurenames))), nCompositeSpectra=length(Featurenames))
        return(matchSummary)
      } else {
        return()
      }
    })
    
    output$tabbedPanelButton = shiny::renderUI({
      shiny::actionButton("tabbedPanelButton", "ViewData")
    })
  }
  })
  # main body observer tabbed panel
  shiny::observe({if(!is.null(input$tabbedPanelButton)){ 
    if(input$tabbedPanelButton > 0){
      if(input$FeatureNames %in% Features.v){
        
        # if(!is.null(input$CommentButton))
        # {
        #   if(input$CommentButton == 0)
        #   {
        #     vars <- shiny::reactiveValues(actionCounter=0)
        #   } else {
        #     vars <- shiny::reactiveValues(actionCounter=0)
        #     vars$actionCounter <- input$CommentButton 
        #   }
        # }
        # isolate feature index
        feat.indx <- which(Features.v %in% shiny::isolate(input$FeatureNames))
        
        ###########################
        ##### 1. Raw data plot ####
        ###########################
        
        plotDf <- reactive({
          MS2_data <- data.frame(composite_spectra[[feat.indx]],  stringsAsFactors = F)
          MS2_data[] <- lapply(MS2_data,  as.character)
          if("interfrag.diff" %in% colnames(MS2_data)){
            MS2_data[, c(1:5)] <- apply(MS2_data[, c(1:5)], 2, function(x) round(as.numeric(x), digits=4))
            SMILESindx <- grep("SMILES$", colnames(MS2_data))
            IDindx <- sapply(paste0(gsub("\\.SMILES", "", colnames(MS2_data)[SMILESindx]), "$"), grep, colnames(MS2_data))
            
            for(i in 1:ncol(MS2_data[, SMILESindx]))
            { 
              SMILESsubIndx <- which(MS2_data[, IDindx[i]]!="")
              if(length(SMILESsubIndx)>0)
              {
                MS2_data[SMILESsubIndx, IDindx[i]] <- sapply(c(1:length(MS2_data[SMILESsubIndx, SMILESindx[i]])),  function(x){
                  Smiles <- unlist(strsplit(MS2_data[SMILESsubIndx[x], SMILESindx[i]], ";"))
                  SmilesInc <- which(Smiles!="")
                  SMILEShtml <- unlist(strsplit(MS2_data[SMILESsubIndx[x], IDindx[i]], ";"))
                  SMILEShtml[SmilesInc] <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", Smiles[SmilesInc], "/PNG", "' target='_blank'>", SMILEShtml[SmilesInc], "</a>")
                  SMILEShtml <- paste(SMILEShtml, collapse=" ")
                  return(SMILEShtml)
                })
              }
            }
            MS2_data <- MS2_data[, -SMILESindx, drop=F]
          } else {
            MS2_data[, c(1:2)] <- apply(MS2_data[, c(1:2)], 2, function(x) round(as.numeric(x), digits=4))
            
          }
          colnames(MS2_data) <- gsub("\\.", "_", colnames(MS2_data))
          MS2_data$Precursorfrag_diff <- as.numeric(MS2_data$Precursorfrag_diff)
          MS2_data$Fragment_Assigned <- ifelse(apply(MS2_data[, c("Frag_ID", "Neutral_loss", "interfrag_loss")], 1, function(x) any(x!="")), "Fragment_identified", "No_Fragment_identified")
          MS2_data <- MS2_data[, grepl('_SMILES', colnames(MS2_data))  ==  F, drop=F]
          # if spectral DB matches
          if(feat.indx %in% indxSpectralDb){
          if(!is.null(input$spectralDBtable_rows_selected)){
            specDBtableTmp <- specDBmatches[[feat.indx]]
            indxTmp <- duplicated(specDBtableTmp$dbSpectra[, 'compound_msp']) == F
            indivDBentries <- specDBtableTmp$dbSpectra[indxTmp, , drop=F]
            
            selectedDBentry <- indivDBentries[input$spectralDBtable_rows_selected, 'compound_msp', drop=F]
          dbMassIntensities <- specDBtableTmp$dbSpectra[specDBtableTmp$dbSpectra[, 'compound_msp'] %in% selectedDBentry, 1:2, drop=F]
          dbMassIntensities <- apply(dbMassIntensities, 2, as.numeric)
          if(!is.matrix(dbMassIntensities)){
            dbMassIntensities <- matrix(dbMassIntensities, ncol=2, byrow=T)
          }
          colnames(dbMassIntensities) <- c('mass', 'intensity')
          dbMassIntensities <- cbind(dbMassIntensities, Rel_Intensity=-(dbMassIntensities[, 'intensity']/max(dbMassIntensities[, 'intensity']) * 100), dbData=rep(1, nrow(dbMassIntensities)))
          MS2_data$dbData <- 0
          MS2_data <- rbind(MS2_data[, c('mass', 'Rel_Intensity', 'dbData')],
                            dbMassIntensities[, c('mass', 'Rel_Intensity', 'dbData')])
            }
          }
          return(MS2_data)
        })
        
        output$MS2_plot <- shiny::renderPlot({
         
              plotDfTmp <- plotDf()
              if(!is.null(input$compMS2_brush)){
                xlimTmp <- c(input$compMS2_brush$xmin, input$compMS2_brush$xmax)
                ylimTmp <- c(input$compMS2_brush$ymin, input$compMS2_brush$ymax)
              } else if('dbData' %in% colnames(plotDfTmp)){
                xlimTmp <- c(0, max(plotDfTmp$mass) * 1.1)
                ylimTmp <- c(-100, max(plotDfTmp$Rel_Intensity))  
              } else {
                xlimTmp <- c(0, max(plotDfTmp$mass) * 1.1)
                ylimTmp <- c(0, max(plotDfTmp$intensity))
              }
              
              if('dbData' %in% colnames(plotDfTmp)){
                colsTmp <- ifelse(plotDfTmp[, 'dbData'] == 1, "#009E73", '#000000')
                plot(plotDfTmp[, c('mass', 'Rel_Intensity'), drop=F], xlim=xlimTmp, ylim=ylimTmp, yaxt='n', ylab='relative intensity',
                     type='h', col=colsTmp, cex.axis=1.5, cex.lab=1.5)
                axis(2, at=seq(-100, 100, 20), labels=c(abs(seq(-100, -20, 20)), seq(0, 100, 20)), las=2, cex.axis=1.3)
                abline(h=0)
                legend('topright', c("composite spectrum", 'database spectrum'), lty=c(1, 1), lwd=c(4, 4), col=c('#000000', "#009E73"), pt.cex=4, cex=1.8, ncol=1)
              } else {
              plot(plotDfTmp[, c('mass', 'intensity')], xlim=xlimTmp, ylim=ylimTmp,
                   type='h', col=ifelse(plotDfTmp$Fragment_Assigned  ==  'Fragment_identified', 'red', 'black'), cex.axis=1.5, cex.lab=1.5)
              }
          })
        
        # spectral db table
        output$spectralDBtable <- DT::renderDataTable({
          if(feat.indx %in% indxSpectralDb){
          specDBtableTmp <- specDBmatches[[feat.indx]]
          indxTmp <- duplicated(specDBtableTmp$dbSpectra[, 'compound_msp']) == F
          dbSpectraTmp <- do.call(rbind, strsplit(specDBtableTmp$dbSpectra[indxTmp, 'compound_msp'], '__'))
          dbSpectraTmp <- cbind(dbSpectraTmp, round(as.numeric(specDBtableTmp$dbSpectra[indxTmp, 'dotProductScore']), 2))
          colnames(dbSpectraTmp) <- c('compound', 'mspFile', 'dotProdScore')
          return(dbSpectraTmp)
          } else {
          return(data.frame(spectral_database_match='no spectral database match'))  
          }
        }, selection='single', rownames=FALSE)
        
        
        # ui tab 1 text plot hover and main plot
        output$compMS2Hover <- shiny::renderText({
          if(!is.null(input$compMS2_hover)){
            paste0('m/z = ', round(input$compMS2_hover$x, 3), ' intensity = ', round(input$compMS2_hover$y, 3))
          } else {
            return('m/z = NULL intensity = NULL')
          }
        })
        # ui tab 1 nearpoints plot click 
        output$compMS2tableInfo <- DT::renderDataTable({
          plotDfTmp <- plotDf()
          
          brushedPoints(plotDfTmp, input$compMS2_brush, xvar = "mass", yvar = "intensity")
        }, rownames=FALSE,  escape = F)
          
        ###########################
        ##### 9. overview plot ####
        ########################### 
      
       
        output$overview_plot <- shiny::renderPlot({
          
          if(!is.null(input$overview_brush)){
            xlimTmp <- c(input$overview_brush$xmin, input$overview_brush$xmax)
            ylimTmp <- c(input$overview_brush$ymin, input$overview_brush$ymax)
          } else {
            xlimTmp <- c(min(allFeatTable$rt), max(allFeatTable$rt))
            ylimTmp <- c(min(allFeatTable$mass), max(allFeatTable$mass))
          }
          if(input$goButton){
          Featurenames <- shiny::isolate({Featureselection()})
          }
          if(input$DBbutton){
          Featurenames <- shiny::isolate({DBFeatureselection()})
          }
          subFeatTable <- allFeatTable[allFeatTable$specNames %in% Featurenames, , drop=F]
          if(nrow(subFeatTable) ==  0){
          subFeatTable <- allFeatTable  
          }
          colsTmp <- rep("darkblue", nrow(subFeatTable))
          # if any commented then change colour
          if (!is.null(input$hot)){
            metIDcomments <- hot_to_r(input$hot)
          } 
          
          compMSCommented <- apply(metIDcomments[, 3:ncol(metIDcomments)], 1, function(x) any(x != '')) 
          overviewMatchCommented <- match(subFeatTable$specNames, metIDcomments$compSpectrum)
          compMSCommented <- overviewMatchCommented %in% which(compMSCommented)
          colsTmp[compMSCommented] <- "#C77CFF"  
          selectedFeat <- allFeatTable$specNames[feat.indx]
          colsTmp[subFeatTable$specNames %in% selectedFeat] <- 'red'
          with(subFeatTable, symbols(x=rt, y=mass, circles=precursorInt_group, inches=1/8, bg=colsTmp, fg=NULL, ylab='m/z', xlab='retentionTime', cex.axis=1.5, cex.lab=1.5, ylim = ylimTmp, xlim = xlimTmp))
          legend('topleft', c('already commented', 'currently selected'),
                 pch=c(NA, NA), lty=c(1, 1), lwd=c(4, 4),
                 col=c("#C77CFF", 'red'), cex=1.4, ncol=1)
      })
        
        
        # ui tab 9 nearpoints plot click 
        output$overviewtableInfo <- DT::renderDataTable({
          if(input$goButton){
            Featurenames <- shiny::isolate({Featureselection()})
          }
          if(input$DBbutton){
            Featurenames <- shiny::isolate({DBFeatureselection()})
          }
          subFeatTable <- allFeatTable[allFeatTable$specNames %in% Featurenames, , drop=F]
          if(nrow(subFeatTable) ==  0){
            subFeatTable <- allFeatTable  
          }
          
          brushedPoints(subFeatTable, input$overview_brush, yvar = "mass", xvar = "rt")}, rownames=FALSE)
        
        #######################################
        ##### 10. correlation network plot ####
        ####################################### 
        output$corrNodesEdges <- shiny::renderText({
          paste0('nodes: ', length(igraph::V(corrNetTmp)), ' edges: ', length(igraph::E(corrNetTmp)), ' (N.B. large numbers of nodes e.g. >= 300 may not display properly)')
        })
      
      output$corrNetworkTableBrush <- DT::renderDataTable({
      bpDfTmp  <- brushedPoints(corrScaledLayout, input$corr_network_brush, xvar='xvar', yvar='yvar')
      bpDfTmp <- bpDfTmp[, 3:6]
      }, rownames=FALSE, options = list(pageLength = 20))
               
      output$corr_network_plot <- shiny::renderPlot({
          if(length(object@network) > 0){
            if(!is.null(input$corr_network_brush)){
              xlimTmp <- c(input$corr_network_brush$xmin, input$corr_network_brush$xmax)
              ylimTmp <- c(input$corr_network_brush$ymin, input$corr_network_brush$ymax)
            } else {
              xlimTmp <- c(-1, 1)
              ylimTmp <- c(-1, 1)
            } 
             
            vertexSizeSub <- igraph::V(corrNetTmp)$vertexSize
            MS2netColsSub <- igraph::V(corrNetTmp)$MS2netColours
            vertexShapesSub <- igraph::V(corrNetTmp)$vertexShapes
            # if any commented then change colour
            if(!is.null(input$hot)){
              metIDcomments <- hot_to_r(input$hot)
            } 
            
            compMSCommented <- apply(metIDcomments[, 3:ncol(metIDcomments)], 1, function(x) any(x != '')) 
            corrNetMatchCommented <- match(paste0('CC_', igraph::V(corrNetTmp)$name), metIDcomments$compSpectrum)
            compMSCommented <- corrNetMatchCommented %in% which(compMSCommented)
            if(any(compMSCommented)){
            MS2netColsSub[compMSCommented] <- "#C77CFF"  
            }
            # colour selected features
            selFeatIndx <- corrNetMatchIndx %in% feat.indx
            if(any(selFeatIndx)){
            vertexSizeSub[selFeatIndx] <- 6
            MS2netColsSub[selFeatIndx] <- "#7CAE00"
            # id first neighbours
            neighSel <- neighbors(corrNetTmp, which(selFeatIndx), mode='all')
            MS2netColsSub[neighSel] <- "#7CAE00"
            }
            # highlight subset features as triangles
            if(input$goButton){
              Featurenames <- shiny::isolate({Featureselection()})
            }
            if(input$DBbutton){
              Featurenames <- shiny::isolate({DBFeatureselection()})
            }
            subsetFeatures <- which(Features.v %in% Featurenames)
            subsetFeatures <- corrNetMatchIndx %in% subsetFeatures
            if(any(subsetFeatures)){
            vertexShapesSub[subsetFeatures] <- 'csquare'  
            }
           
            # black background igraph
            par(bg = "black")
             
            plot(corrNetTmp, layout=corrLayoutTmp[, 1:2], edge.arrow.size=.1, edge.color="gray33", vertex.color=MS2netColsSub, vertex.label.font=2, vertex.label.color= "gray83", vertex.label=V(corrNetTmp)$name, vertex.shape=vertexShapesSub, vertex.size=vertexSizeSub, vertex.label.cex=1.5, xlim=xlimTmp, ylim=ylimTmp) #layout=layout.circle,
            legend('topleft', c("currently selected spectrum and 1st neighbours (if present)", 'currently subset EIC', "MS2 matched EIC", "unmatched EIC", 'already commented', paste0('edge corrCoeff >= ', round(object@Parameters$corrThresh, 2))), pch=c(NA, 22, 21, 21, 21, NA), lty=c(1, NA, NA, NA, NA, 1), lwd=c(4, 1, 1, 1, 1, 4),
                   col=c("#7CAE00", "black", "black", "black", "black", 'gray33'), pt.bg=c("#7CAE00", "#D55E00", "#D55E00", "#0072B2", "#C77CFF", "gray33"), pt.cex=4, cex=1.8, bg='gray79', ncol=1)#text.col='white', bty="n",
          }
        })
        
      ###############################################
      ##### 12. spectral similarity network plot ####
      ############################################### 
      output$specSimNodesEdges <- shiny::renderText({
        paste0('nodes: ', length(igraph::V(specSimNetTmp)), ' edges: ', length(igraph::E(specSimNetTmp)), ' (N.B. large numbers of nodes e.g. >= 300 may not display properly)')
      })
      
      output$specSimNetworkTableBrush <- DT::renderDataTable({
        bpDfTmp  <- brushedPoints(specSimScaledLayout, input$specSim_network_brush, xvar='xvar', yvar='yvar')
        bpDfTmp <- bpDfTmp[, 3:6, drop=F]
      }, rownames=FALSE, options = list(pageLength = 20))
      
      output$specSim_network_plot <- shiny::renderPlot({
        if(length(object@network) > 0){
          if(!is.null(input$specSim_network_brush)){
            xlimTmp <- c(input$specSim_network_brush$xmin, input$specSim_network_brush$xmax)
            ylimTmp <- c(input$specSim_network_brush$ymin, input$specSim_network_brush$ymax)
          } else {
            xlimTmp <- c(-1, 1)
            ylimTmp <- c(-1, 1)
          } 
          
          vertexSizeSub <- igraph::V(specSimNetTmp)$vertexSize
          MS2netColsSub <- igraph::V(specSimNetTmp)$MS2netColours
          vertexShapesSub <- igraph::V(specSimNetTmp)$vertexShapes
          # if any commented then change colour
          if (!is.null(input$hot)){
            metIDcomments <- hot_to_r(input$hot)
          } 
          
          compMSCommented <- apply(metIDcomments[, 3:ncol(metIDcomments)], 1, function(x) any(x != '')) 
          specSimMatchCommented <- match(igraph::V(specSimNetTmp)$name, metIDcomments$compSpectrum)
          compMSCommented <- specSimMatchCommented %in% which(compMSCommented)
          if(any(compMSCommented)){
            MS2netColsSub[compMSCommented] <- "#C77CFF"  
          }
          # colour selected features
          selFeatIndx <- specSimMatchIndx %in% feat.indx
          if(any(selFeatIndx)){
            vertexSizeSub[selFeatIndx] <- 6
            MS2netColsSub[selFeatIndx] <- "#7CAE00"
            # id first neighbours
            neighSel <- neighbors(specSimNetTmp, which(selFeatIndx), mode='all')
            MS2netColsSub[neighSel] <- "#7CAE00"
          }
          # highlight subset features as triangles
          if(input$goButton){
            Featurenames <- shiny::isolate({Featureselection()})
          }
          if(input$DBbutton){
            Featurenames <- shiny::isolate({DBFeatureselection()})
          }
          subsetFeatures <- which(Features.v %in% Featurenames)
          subsetFeatures <- specSimMatchIndx %in% subsetFeatures
          if(any(subsetFeatures)){
            vertexShapesSub[subsetFeatures] <- 'csquare'  
          }
          
          # black background igraph
          par(bg = "black")
          
          plot(specSimNetTmp, layout=specSimLayoutTmp[, 1:2], edge.arrow.size=.1, edge.color=igraph::E(specSimNetTmp)$color, vertex.color=MS2netColsSub, vertex.label.font=2, vertex.label.color= "gray83", vertex.label=gsub('CC_', '', V(specSimNetTmp)$name), vertex.shape=vertexShapesSub, vertex.size=vertexSizeSub, vertex.label.cex=1.5, xlim=xlimTmp, ylim=ylimTmp) #layout=layout.circle,
          legend('topleft', c("currently selected spectrum and 1st neighbours (if present)", 'currently subset EIC', "MS2 matched EIC", 'already commented', paste0('edge fragment ions (dot product >= ', round(object@Parameters$minDotProdThresh, 2), ')'), paste0('edge neutral losses (dot product >= ', round(object@Parameters$minDotProdThresh, 2), ')')), pch=c(NA, 22, 21, 21, NA, NA), lty=c(1, NA, NA, NA, 1, 1), lwd=c(4, 1, 1, 1, 4, 4),
                 col=c("#7CAE00", "black", "black", "black", "#CC79A7", "#56B4E9"), pt.bg=c("#7CAE00", "#D55E00", "#D55E00", "#C77CFF", "#CC79A7", "#56B4E9"), pt.cex=4, cex=1.8, bg='gray79', ncol=1)#text.col='white', bty="n",
        }
      })
      
        ##################################
        ##### 11. metID comments table ###
        ##################################
        output$hot = renderRHandsontable({
          if (!is.null(input$hot)) {
            metIDcomments = hot_to_r(input$hot)
          } 
          
          if(input$goButton){
            Featurenames <- shiny::isolate({Featureselection()})
          }
          if(input$DBbutton){
            Featurenames <- shiny::isolate({DBFeatureselection()})
          }
          currSubset <- metIDcomments$compSpectrum %in% Featurenames
          currSubset <- ifelse(currSubset, 'Yes', 'No')
          metIDcomments$currently_subset <- currSubset
          metIDcomments <- metIDcomments[, c("compSpectrum", 'currently_subset', "possible_identity", "compound_class", "user_comments")]
          # order by currently subset
          metIDcomments <- metIDcomments[order(currSubset, decreasing = T), , drop=F]
          if(all(currSubset == 'Yes')){
          # resort to Features.v order
          metIDcomments <- metIDcomments[order(as.numeric(row.names(metIDcomments))), ]  
          }
          # currently selected to top of table
          metIDcomments <- metIDcomments[order(metIDcomments$compSpectrum %in% Features.v[feat.indx], decreasing = T), ]
          setHot(metIDcomments)
          rhandsontable(metIDcomments, readOnly = object@Parameters$readOnly) %>%
            # hot_col('compSpectrum', readOnly = T) %>%
            hot_table(highlightCol = TRUE, highlightRow = TRUE)
        })
        ###########################
        ##### 2. metadata table ###
        ###########################
        
        output$metadata <- DT::renderDataTable({
          metadata.tmp.sub <- metaData.tmp[[feat.indx]]
          row.names.tmp <- gsub(".+\\.mzXML_[0-9]*_",  "",  names(metadata.tmp.sub))
          #   gsub("_[:alpha:].+",  "",  names(metadata))
          column.names.tmp <- gsub(gsub("\\.",  "\\\\.",  
                                        paste(unique(row.names.tmp),  collapse = "$|")), 
                                   "",  names(metadata.tmp.sub))
          column.names.tmp <- gsub("_$",  "",  column.names.tmp)
          metadata.tmp.sub <- sapply(metadata.tmp.sub,  function(x) {
            if(nchar(gsub("(.*\\.)",  "",  as.character(x[1]))) > 4){
              tmp <- round(x,  digits = 4) 
              tmp <- paste(tmp,  collapse = "; ")
              return(tmp)
            } else {
              tmp <- paste(x,  collapse = "; ")
              return(tmp)
            }
          })
          metadata.tmp.sub <- sapply(unique(column.names.tmp),  function(x) 
            metadata.tmp.sub[grep(x,  names(metadata.tmp.sub))])
          row.names(metadata.tmp.sub) <- unique(row.names.tmp)
          return(metadata.tmp.sub)    
        },  options = list(pageLength = 21))
        
        ###########################
        ##### 3. MS2 data table ###
        ###########################
        
#         output$MS2_data <- DT::renderDataTable({
#           MS2_data <- data.frame(composite_spectra[[feat.indx]],  stringsAsFactors = F)
#           MS2_data[] <- lapply(MS2_data,  as.character)
#           if("interfrag.diff" %in% colnames(MS2_data)){
#             MS2_data[, c(1:5)] <- apply(MS2_data[, c(1:5)], 2, function(x) round(as.numeric(x), digits=4))
#             SMILESindx <- grep("SMILES$", colnames(MS2_data))
#             IDindx <- sapply(paste0(gsub("\\.SMILES", "", colnames(MS2_data)[SMILESindx]), "$"), grep, colnames(MS2_data))
#             
#             for(i in 1:ncol(MS2_data[, SMILESindx]))
#             { 
#               SMILESsubIndx <- which(MS2_data[, IDindx[i]]!="")
#               if(length(SMILESsubIndx)>0)
#               {
#                 MS2_data[SMILESsubIndx, IDindx[i]] <- sapply(c(1:length(MS2_data[SMILESsubIndx, SMILESindx[i]])),  function(x){
#                   Smiles <- unlist(strsplit(MS2_data[SMILESsubIndx[x], SMILESindx[i]], ";"))
#                   SmilesInc <- which(Smiles!="")
#                   SMILEShtml <- unlist(strsplit(MS2_data[SMILESsubIndx[x], IDindx[i]], ";"))
#                   SMILEShtml[SmilesInc] <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", Smiles[SmilesInc], "/PNG", "' target='_blank'>", SMILEShtml[SmilesInc], "</a>")
#                   SMILEShtml <- paste(SMILEShtml, collapse=" ")
#                   return(SMILEShtml)
#                 })
#               }
#             }
#             MS2_data <- MS2_data[, -SMILESindx, drop=F]
#           } else {
#             MS2_data[, c(1:2)] <- apply(MS2_data[, c(1:2)], 2, function(x) round(as.numeric(x), digits=4))
#             
#           }
#           return(MS2_data)    
#         },  options = list(pageLength = 25),  escape = F)
        
        ###########################
        ##### 4. DB results #####
        ###########################
        
        output$DB.results <- DT::renderDataTable({
          if(length(object@DBanno)  ==  0){
            DB.results <- data.frame("metID.dbAnnotate function has not yet been run",  stringsAsFactors=F)
            colnames(DB.results) <- "Result"
            return(DB.results) 
          } else {
            DB.results <- tmp.DBanno.res[[feat.indx]]
            ###create clickable url links to DB in table
            DB.results$DBid  <- paste0("<a href='http://", DB.results$WebAddress,  DB.results$DBid,  
                                     "' target='_blank'>",  DB.results$DBid,  "</a>")
            # clickable SMILES pugRest
            SMILESindx <- grep("SMILES$", colnames(DB.results))
            DB.results[, SMILESindx] <- sapply(SMILESindx, function(x){
              NoFuncIndx <- DB.results[,  x]  ==  ""
              ifelse(NoFuncIndx, "", 
                     paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", 
                            gsub("/|\\\\", "",  DB.results[, x]), "/PNG", "' target='_blank'>",  
                            paste0(substring(DB.results[, x],  1,  6),  "..."), "</a>"))
            })
            DB.results$WebAddress <- NULL 
            DB.results[,  c("expectedMass",  "candidateMass",  "ppmMatch")] <- lapply(DB.results[,  c("expectedMass",  "candidateMass",  "ppmMatch")],  as.numeric)
            return(DB.results)
          }}, rownames=FALSE,  escape=F,  options = list(pageLength = 20))#, digits=4,  sanitize.text.function = function(x) x)
        
        ##############################
        ##### 5. Best annotations #####
        ##############################
        
        output$BestCandidates <- DT::renderDataTable({
          if(length(object@BestAnno)  ==  0){
            DB.results <- data.frame("metID.dbProb function has not yet been run", stringsAsFactors=F)
            colnames(DB.results) <- "Result"
            return(DB.results) 
          } else {
            DB.results <- tmp.BestAnno[[feat.indx]]
          
          if(nrow(DB.results) == 0){
              DB.results <- data.frame("no best/ most likely annotations based on substructures identified", stringsAsFactors=F)
              colnames(DB.results) <- "Result"
              return(DB.results) 
            } else {
            ###create clickable url links to DB in table
            DB.results$DBid  <- paste0("<a href='http://", DB.results$WebAddress,  DB.results$DBid,  
                                     "' target='_blank'>",  DB.results$DBid,  "</a>")
            # clickable SMILES pugRest
            SMILESindx <- grep("SMILES$", colnames(DB.results))
            DB.results[, SMILESindx] <- sapply(SMILESindx, function(x){
              NoFuncIndx <- DB.results[,  x]  ==  ""
              ifelse(NoFuncIndx, "", 
                     paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", 
                            gsub("/|\\\\", "",  DB.results[, x]), "/PNG", "' target='_blank'>",  
                            paste0(substring(DB.results[, x],  1,  6),  "..."), "</a>"))
            })
            DB.results$WebAddress <- NULL 
            DB.results[,  c("expectedMass",  "candidateMass",  "ppmMatch")] <- lapply(DB.results[,  c("expectedMass",  "candidateMass",  "ppmMatch")],  as.numeric)
            return(DB.results)
          }}}, rownames=FALSE,  escape=F,  options = list(pageLength = 20))#, digits=4,  sanitize.text.function = function(x) x)
        
        
        ############################################
        ##### 6. Best substructure annotations #####
        ############################################
        
        output$BestSubStrAnno <- DT::renderDataTable({
          if(nrow(subStrAnno.df)  ==  0){
            SbStrResults <- data.frame("subStructure.prob function has not yet been run",  stringsAsFactors=F)
            colnames(SbStrResults) <- "Result"
            return(SbStrResults) 
          } else {
            SbStrResults <-  subStrAnno.df[subStrAnno.df$compSpecName %in% Features.v[feat.indx],  ,  drop = F]
            return(SbStrResults)
          }}, rownames=FALSE,  escape=F,  options = list(pageLength = 8))
        # text entry for custom searches
        output$customPubMedSearch <- shiny::renderUI({
          shiny::textInput('customPubMedSearch', 'Enter custom search text:')
        })
        #ui for max number of HMDB abstracts to return
        output$nPMIDAbstracts <- shiny::renderUI({
          shiny::numericInput('nPMIDAbstracts', 'Number of PubMed abstracts used to calculate word cloud (max = 10000) : ',  value=100, min=10, max=10000, step=10)
        }) 
        # #ui for max number of random articles to return
        # output$nRandomArticles <- shiny::renderUI({
        #   shiny::numericInput('nRandomArticles', 'Number of randomly selected abstracts to display (max = 20) : ',  value=5, min=1, max=20, step=1)
        # }) 
        
        ##ui for word cloud drop-down
        output$wordCloudSelect <- renderUI({
          if(length(object@DBanno)  ==  0){
            DBmatches <- "metID.dbAnnotate function has not yet been run/ No matches to DB"
          } else {
            DBmatches <- tmp.DBanno.res[[feat.indx]]$DBname
            shiny::selectizeInput('wordCloudSelect', 'Or DB annotation match name : ',  choices = DBmatches)  
          }}) 
        # actionButton to start pubmed search
        output$pubMedSearchButton <- shiny::renderUI({
          shiny::actionButton('pubMedSearchButton', 'Search PubMed')
        })
        
        ###word cloud text output  
        output$WordCloudText <- shiny::renderText({"If a large number of abstracts (>100) from PubMed are returned the word cloud will take longer to update or may fail to load...Please Wait"})
        
        # mine PubMed reactive
        pubMedMine <- reactive({
          if(length(object@DBanno)  ==  0){
            ClAbs <- data.frame(word="metID.dbAnnotate function not run",  freq=1)
            wordCloudDf <- data.frame("No PMIDs returned")
            colnames(wordCloudDf) <- input$FeatureNames
            PMIDs <- 'none'
            names(PMIDs) <-  ifelse(searchTerm == '', 'noSearchTerm', searchTerm)
            Abs <- 'none'
          } else if (nrow(tmp.DBanno.res[[feat.indx]])   ==  0){
            ClAbs <- data.frame(word="No matches to DB",  freq=1)
            wordCloudDf <- data.frame("No PMIDs returned")
            colnames(wordCloudDf) <- input$FeatureNames
            PMIDs <- 'none'
            names(PMIDs) <-  ifelse(searchTerm == '', 'noSearchTerm', searchTerm)
            Abs <- 'none'
          } else {
            
            if(input$customPubMedSearch == ''){
              #subset -1 to remove count value obtained from eutils XML file
              PMIDs <- PMIDsearch(keys=input$wordCloudSelect,  n=99999)
              Count <- PMIDs[1]
              PMIDs <- PMIDs[-1]
              searchTerm <- input$wordCloudSelect
              } else {
              # else custom text search
              PMIDs <- PMIDsearch(keys=input$customPubMedSearch,  n=99999)
              Count <- PMIDs[1]
              PMIDs <- PMIDs[-1]
              searchTerm <- input$customPubMedSearch
              }
           
          if(length(PMIDs) > 0){
            # random sample of all PMIDs to sample
           randIndx <- sample(1:length(PMIDs), 
                              ifelse(length(PMIDs) >= input$nPMIDAbstracts, 
                                     input$nPMIDAbstracts, length(PMIDs)))
           randPMIDs <- PMIDs[randIndx]
           articleTitles <- getTitles(randPMIDs)
           Abs <- PubMedWordcloud::getAbstracts(randPMIDs)
           wordCloudDf <- data.frame(randPMIDs, articleTitles, stringsAsFactors=F)
           wordCloudDf$url <- paste0('http://www.ncbi.nlm.nih.gov/pubmed/', wordCloudDf$randPMIDs)   
           row.names(wordCloudDf) <- paste0("article ", 
                                            seq(1, nrow(wordCloudDf), 1))
           colnames(wordCloudDf)[1] <- paste0('Search term: ', searchTerm, " total PubMedIds returned: ", Count)
                    
            if(length(Abs) > 0){
                  ClAbs <- PubMedWordcloud::cleanAbstracts(Abs)
                  ###remove compound names from word frequency table
                  subsName <- strsplit(input$wordCloudSelect, " ")[[1]]
                  subsName <- unlist(lapply(subsName, function(x) grep(paste0("\\b", x, "\\b"), ClAbs$word, ignore.case=T)))
                  if(length(subsName) > 0){
                    ClAbs <- ClAbs[-subsName, ]
                  }
                  ###only keep word which are less than 50 characters
                  ClAbs <- ClAbs[which(sapply(as.character(ClAbs$word), nchar) < 50), ]
                  PMIDs <- paste0(PMIDs, collapse = '; ')
                  names(PMIDs) <-  searchTerm
                  Abs <- data.frame(Abs)
                  colnames(Abs) <- searchTerm
                  # row.names(Abs) <- randPMIDs
                } else {
                  ###if no PMIDs returned then plot 
                  ClAbs <- data.frame(word="No Abstract text available", freq=1)
                  wordCloudDf <- data.frame("No PMIDs returned")
                  colnames(wordCloudDf) <- input$FeatureNames
                  PMIDs <- 'none'
                  names(PMIDs) <-  ifelse(searchTerm == '', 'noSearchTerm', searchTerm)
                  Abs <- 'none'
                }
              } else {
                ###if no PMIDs returned then plot 
                ClAbs <- data.frame(word="No PubMedIDs returned", freq=1)
                wordCloudDf <- data.frame("No PMIDs returned")
                colnames(wordCloudDf) <- input$FeatureNames
                PMIDs <- 'none'
                names(PMIDs) <-  ifelse(searchTerm == '', 'noSearchTerm', searchTerm)
                Abs <- 'none'
              }
          }
          return(list(ClAbs=ClAbs, wordCloudDf=wordCloudDf, PMIDs=PMIDs, Abs=Abs))
        })
        ####################################
        ##### 6. PubMed wordcloud plot #####
        ####################################
        
        output$WordCloud <- shiny::renderPlot({
        if(input$pubMedSearchButton == 0){ 
         return() } else if(input$pubMedSearchButton > 0){
        ClAbs <- shiny::isolate(pubMedMine()[[1]])
        if(nrow(ClAbs) == 1){
          PubMedWordcloud::plotWordCloud(ClAbs,  min.freq=1,  max.words=100,  rot.per=0)     } else {
          suppressWarnings(PubMedWordcloud::plotWordCloud(ClAbs, max.words=100, scale=c(4, 0.5)))
          }
          }})
        
        ###########################################
        ##### 7. PubMed random article table  #####
        ###########################################
        
        output$WordCloudTable <- DT::renderDataTable({
          if(input$pubMedSearchButton == 0){ 
            wordCloudDf <- data.frame("No Search of PubMed performed")
            colnames(wordCloudDf) <- input$FeatureNames
            return(wordCloudDf)
             } else if(input$pubMedSearchButton > 0){
              wordCloudDf <- shiny::isolate(pubMedMine()[[2]])
              
              wordCloudDf$htmlUrl <- paste0("<a href='http://www.ncbi.nlm.nih.gov/pubmed/", 
                     wordCloudDf[, 1], "' target='_blank'> PMID : ", 
                     wordCloudDf[, 1], " Title : ", 
                     wordCloudDf$articleTitles, "</a>")
              htmlUrlOnly <- wordCloudDf[, 'htmlUrl', drop=F]
              colnames(htmlUrlOnly) <- colnames(wordCloudDf)[1]
              return(htmlUrlOnly)
             }}, escape=F,  options = list(pageLength = 10))
          # ,  sanitize.text.function = function(x) x)
        
        # download pubmed data
        output$downloadPubMedData <- shiny::downloadHandler(
          filename = function(){
            resFiles <- isolate(pubMedMine())
            fileNameTmp <- paste0("PubMedSearch_", names(resFiles[[3]]), "_",
                                         Sys.Date(), '.zip')},
          content = function(file){
            resFiles <- isolate(pubMedMine())
            tmpdir <- tempdir()
            setwd(tempdir())
            print(tempdir())
            
            fileNames <- paste0("PubMedSearch_", names(resFiles[[3]]),
                                c("_wordFreq_", "_randArticles_", '_allPMIDs_', '_Abstracts_'), Sys.Date(), '.txt')
            
            write.table(resFiles['ClAbs'], fileNames[1], sep='\t', row.names=F)
            write.table(resFiles['wordCloudDf'], fileNames[2], sep='\t', row.names=F)
            write.table(resFiles['PMIDs'], fileNames[3], sep='')
            # write.table(resFiles['Abs'], fileNames[4], sep='\t')
  
            print(fileNames)
            
            zip(zipfile=file, files=fileNames[1:3])
          }, contentType = "application/zip")
        
        #####################################
        ##### 8. output MetFrag results #####
        #####################################
        
        output$MetFragTable <- shiny::renderTable({
          if(length(tmp.metFrag)  ==  0){
            metFrag.df.tmp <- "metID.MetFrag function has not yet been run"
            metFrag.df.tmp <- data.frame(metFrag.df.tmp)
            colnames(metFrag.df.tmp) <- "Result"
            return(metFrag.df.tmp)
          } else {
            metFrag.df.tmp <- tmp.metFrag[[feat.indx]]
            if(is.null(metFrag.df.tmp)){
              metFrag.df.tmp <- "metID.MetFrag function has not yet been run"
              metFrag.df.tmp <- data.frame(metFrag.df.tmp)
              colnames(metFrag.df.tmp) <- "Result" 
              return(metFrag.df.tmp)
            } else if (is.character(metFrag.df.tmp)){
              metFrag.df.tmp <- "MetFrag returned no results"
              metFrag.df.tmp <- data.frame(metFrag.df.tmp)
              colnames(metFrag.df.tmp) <- "Result" 
              return(metFrag.df.tmp)
            } else {  
              # reduce n decimal places to 4
              metFrag.df.tmp[,  c("PeakScore",  "BondEnergyScore",  "Score")] <- sapply(metFrag.df.tmp[,  c("PeakScore",  "BondEnergyScore",  "Score")],  
                                                                                     function(x) {
                                                                                       if(nchar(gsub("(.*\\.)",  "",  as.character(x[1]))) > 4){
                                                                                         x <- round(as.numeric(x),  digits = 4) 
                                                                                       }
                                                                                     })
              ###create clickable url links to DB in table
              metFrag.df.tmp$DBid  <- paste0("<a href='http://",  metFrag.df.tmp$WebAddress,  metFrag.df.tmp$DBid,  
                                           "' target='_blank'>",  metFrag.df.tmp$DBid,  "</a>")
              metFrag.df.tmp$WebAddress <- NULL 
              metFrag.df.tmp$SMILES <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", 
                                              metFrag.df.tmp$SMILES,  "/PNG",  "' target='_blank'>", 
                                              paste0(substring(metFrag.df.tmp$SMILES,  1,  6),  "..."),  "</a>")
              return(metFrag.df.tmp)
            }
          }
        },  sanitize.text.function = function(x) x)
        
        ##load in any previous user comments
        ###add text area for storing user comments
        
        output$CommentSectionHeader <- shiny::renderUI({shiny::tags$b("Add metabolite identification notes including any useful html links here...\n
                                                                      NB. Do not forget to save your comments to the features results directory using the button on the left")})
        output$UserComments <- shiny::renderUI({
          if(length(Comments(compMS2))  ==  0){
            PrevComment <- "No previous comments"
          } else if(is.null(Comments(object)[[feat.indx]])){
            PrevComment <- "No previous comments"
          } else {
            PrevComment <- Comments(object)[[feat.indx]]
            PrevComment <- as.character(PrevComment[nrow(PrevComment),  2])
          }
          
          shiny::tags$textarea(id="UserComments",  rows=4,  cols=60,  PrevComment)
        })
        
        
        #########################
        ##### 13. Comments  #####
        #########################
        ####save  
        # observe({if(!is.null(input$CommentButton))
        # {  
        #   if(input$CommentButton > 0)
        #   {
        #     if(input$CommentButton > vars$actionCounter)
        #     {
        #       vars$actionCounter <- input$CommentButton
        #       UserComments.v[[feat.indx]] <- input$UserComments
        # 
        #     }
        #   }
        # }
        # })
        
        
        }
  }
  }
  })
  
  session$onSessionEnded(function(){
    observe({
    if(!is.null(values[["hot"]])){
    metIDcomments  <- values[["hot"]]
    metIDcomments$currently_subset <- NULL
    metIDcomments <- metIDcomments[order(as.numeric(row.names(metIDcomments))), ]
    object@Comments <- metIDcomments
    }
    shiny::stopApp(object)})
    })
  
}) # END compMS2server