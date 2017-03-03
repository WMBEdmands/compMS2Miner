library(compMS2Miner)
library(DT)
library(shiny)
library(igraph)   
library(rhandsontable)

shiny::shinyServer(function(input,  output, session){
  # comment lines below if action button is used to commit changes
  values <- reactiveValues()
  setHot <- function(x) values[["hot"]] = x
  currEICRv <- reactiveValues(currEIC='')
  output$pdfviewer <- renderText({
    return(paste('<iframe style="height:900px; width:100%" src="', pdfFile, '"></iframe>', sep = ""))
  })
  
  output$currEIC <- shiny::renderText({
    paste0('Current EIC:', gsub('.+_', '', currEICRv$currEIC))
  })
  output$sessionTxt <- shiny::renderUI({shiny::HTML(seshTxtTmp)})
  #options(shiny.trace=T)
  observe({
    if(input$CloseAppBtn > 0){
      if(!is.null(values[["hot"]])){
        metIDcomments  <- values[["hot"]]
        metIDcomments$currently_subset <- NULL
        metIDcomments <- metIDcomments[order(as.numeric(row.names(metIDcomments))), ]       
        object@Comments <- metIDcomments
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
        } else if(input$inSilicoMatch == T){
          FeaturesIndx <- rep(F, length(Features.v))
          FeaturesIndx[inSilicoIndx] <- T
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
    
    currEICRv$currEIC <- ''
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
    currEICRv$currEIC <- ''
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
        currEICRv$currEIC <- input$FeatureNames
        ###########################
        ##### 1. Raw data plot ####
        ###########################
        
        plotDf <- reactive({
          if(!is.null(composite_spectra[[feat.indx]])){
          MS2_data <- data.frame(composite_spectra[[feat.indx]],  stringsAsFactors = F)
          MS2_data[] <- lapply(MS2_data,  as.character)
          if("interfrag.diff" %in% colnames(MS2_data)){
            if(nrow(MS2_data) > 1){
            MS2_data <- as.data.frame(apply(MS2_data, 2, function(x) gsub('noSMILES|noSMILES;|noID|noID;', '', x)), stringsAsFactors=FALSE)
            } else {
            MS2_data[1, ] <- gsub('noSMILES|noSMILES;|noID|noID;', '', MS2_data)  
            }
            MS2_data[, c(1:5)] <- apply(MS2_data[, c(1:5)], 2, function(x) round(as.numeric(x), digits=4))
            SMILESindx <- grep("SMILES$", colnames(MS2_data))
            IDindx <- sapply(paste0(gsub("\\.SMILES", "", colnames(MS2_data)[SMILESindx]), "$"), grep, colnames(MS2_data))
            
            for(i in 1:ncol(MS2_data[, SMILESindx])){ 
              SMILESsubIndx <- which(MS2_data[, IDindx[i]]!="")
              if(length(SMILESsubIndx)>0){
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
          if(!is.null(MS2_data$Precursorfrag_diff)){
          MS2_data$Precursorfrag_diff <- as.numeric(MS2_data$Precursorfrag_diff)
          MS2_data$Fragment_Assigned <- ifelse(apply(MS2_data[, c("Frag_ID", "Neutral_loss", "interfrag_loss")], 1, function(x) any(x!="")), "Fragment_identified", "No_Fragment_identified")
          
          MS2_data <- MS2_data[, grepl('_SMILES', colnames(MS2_data))  ==  F, drop=F]
          }
          # if spectral DB matches
          if(feat.indx %in% indxSpectralDb){
          if(is.numeric(input$spectralDBtable_rows_selected)){
            specDBtableTmp <- specDBmatches[[feat.indx]]
            selectedDBentry <- unique(specDBtableTmp$dbSpectra[, 'compound_msp'])[input$spectralDBtable_rows_selected]
           
          dbMassIntensities <- specDBtableTmp$dbSpectra[specDBtableTmp$dbSpectra[, 'compound_msp'] %in% selectedDBentry, 1:2, drop=F]
          dbMassIntensities <- apply(dbMassIntensities, 2, as.numeric)
          if(!is.matrix(dbMassIntensities)){
            dbMassIntensities <- matrix(dbMassIntensities, ncol=2, byrow=T)
          }
          colnames(dbMassIntensities) <- c('mass', 'intensity')
          dbMassIntensities <- cbind(dbMassIntensities, Rel_Intensity=-(dbMassIntensities[, 'intensity']/max(dbMassIntensities[, 'intensity']) * 100), dbData=rep(1, nrow(dbMassIntensities)))
          MS2_data$dbData <- 0
          if(!'Rel_Intensity' %in% colnames(MS2_data)){
            MS2_data[, 2] <- {MS2_data[, 2]/max(MS2_data[, 2])} * 100
            colnames(MS2_data)[2] <- 'Rel_Intensity'
          }
          MS2_data <- rbind(MS2_data[, c('mass', 'Rel_Intensity', 'dbData')],
                            dbMassIntensities[, c('mass', 'Rel_Intensity', 'dbData')])
            }
          }
        } else {
        MS2_data <- data.frame()
        }  
          return(MS2_data)
        })
        
        output$MS2_plot <- shiny::renderPlot({
              plotDfTmp <- plotDf()
              if(ncol(plotDfTmp) > 0){
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
                if(!is.null(plotDfTmp$Fragment_Assigned)){
                  colsTmp <- ifelse(plotDfTmp$Fragment_Assigned  ==  'Fragment_identified', 'red', 'black')
                } else{
                  colsTmp <- 'black'
                }
              plot(plotDfTmp[, c('mass', 'intensity')], xlim=xlimTmp, ylim=ylimTmp,
                   type='h', col=colsTmp, cex.axis=1.5, cex.lab=1.5)
              }
              } else {
                plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
                text(x = 0.5, y = 0.5, paste('no MS2 data matched to this feature'), cex = 1.6, col = "black")  
              }
          })
        
        # spectral db table
        output$spectralDBtable <- DT::renderDataTable({
          if(feat.indx %in% indxSpectralDb){
          specDBtableTmp <- specDBmatches[[feat.indx]]
          indxTmp <- duplicated(specDBtableTmp$dbSpectra[, 'compound_msp']) == F
          dbSpectraTmp <- do.call(rbind, strsplit(specDBtableTmp$dbSpectra[indxTmp, 'compound_msp'], '__'))
          dbSpectraTmp <- cbind(dbSpectraTmp, round(as.numeric(specDBtableTmp$dbSpectra[indxTmp, 'dotProductScore']), 2), round(as.numeric(specDBtableTmp$dbSpectra[indxTmp, 'propTIC_explained']), 2))
          colnames(dbSpectraTmp) <- c('compound', 'mspFile', 'dotProdScore', 'propTIC_explained')
          return(dbSpectraTmp)
          } else {
          return(data.frame(spectral_database_match='no spectral database match'))  
          }
        }, selection='single', rownames=F)
        
        
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
          if(feat.indx %in% indxSpectralDb){
            if(is.numeric(input$spectralDBtable_rows_selected)){
              specDBtableTmp <- specDBmatches[[feat.indx]]
              selectedDBentry <- unique(specDBtableTmp$dbSpectra[, 'compound_msp'])[input$spectralDBtable_rows_selected]
              indxTmp <- which(names(specDBtableTmp[[2]]) %in% selectedDBentry)
              if(length(indxTmp) > 1){
              dbInfoDfTmp <- do.call(cbind, specDBtableTmp[[2]][indxTmp])  
              dbInfoDfTmp <- cbind(row.names(dbInfoDfTmp), dbInfoDfTmp)
              colnames(dbInfoDfTmp)[1] <- 'EntryNo_EntryName'
              } else {
              dbInfoDfTmp <- as.data.frame(specDBtableTmp[[2]][[indxTmp]], stringsAsFactors = F)
              dbInfoDfTmp <- cbind(row.names(dbInfoDfTmp), dbInfoDfTmp[, 1])
              colnames(dbInfoDfTmp) <- c('EntryNo_EntryName', 'Information')
              }
              return(dbInfoDfTmp)
          } else {
          plotDfTmp <- plotDf()  
          brushedPoints(plotDfTmp, input$compMS2_brush, xvar = "mass", yvar = "intensity")  
          }
          } else {
          plotDfTmp <- plotDf()  
          brushedPoints(plotDfTmp, input$compMS2_brush, xvar = "mass", yvar = "intensity")
          }
        }, rownames=F,  escape = F)
          
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
          
          compMSCommented <- apply(metIDcomments[, 'possible_identity', drop=FALSE], 1, function(x) any(x != '')) 
          overviewMatchCommented <- match(subFeatTable$specNames, metIDcomments$compSpectrum)
          compMSCommented <- overviewMatchCommented %in% which(compMSCommented)
          colsTmp[compMSCommented] <- "#C77CFF"  
          selectedFeat <- allFeatTable$specNames[feat.indx]
          tmpIdx <- subFeatTable$specNames %in% selectedFeat
          colsTmp[tmpIdx] <- 'red'
          coordsTmp <- c(subFeatTable$mass[tmpIdx], subFeatTable$rt[tmpIdx])
          with(subFeatTable, symbols(x=rt, y=mass, circles=precursorInt_group, inches=1/8, bg=colsTmp, fg=NULL, ylab='m/z', xlab='retentionTime', cex.axis=1.5, cex.lab=1.5, ylim = ylimTmp, xlim = xlimTmp))
          legend('topleft', c('uncommented', 'metID Commented', 'currently selected'),
                 pch=c(21, 21, 21), pt.cex=3, #lty=c(1, 1, 1), lwd=c(4, 4, 4),
                 pt.bg=c("darkblue", "#C77CFF", 'red'), cex=1.4, ncol=1)
          abline(h=rep(coordsTmp[1], nrow(subFeatTable)), col="red", 
                 lwd=0.8, lty=2)
          abline(v=rep(coordsTmp[2], nrow(subFeatTable)), col="red", 
                 lwd=0.8, lty=2)
          
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
          
          brushedPoints(subFeatTable, input$overview_brush, yvar = "mass", xvar = "rt")}, rownames=F)
        
        #######################################
        ##### 10. correlation network plot ####
        ####################################### 
        output$corrNodesEdges <- shiny::renderText({
          if(!is.null(object@network$corrNetworkGraph)){
          paste0('nodes: ', length(igraph::V(corrNetTmp)), ' edges: ', length(igraph::E(corrNetTmp)), ' (N.B. large numbers of nodes e.g. >= 300 may not display properly)')
          } else {
          paste0('?metID.corrNetwork function has not yet been run.')  
          }
        })
        
        output$corrNetworkTableBrush  <- DT::renderDataTable({
          if(!is.null(object@network$corrNetworkGraph)){
            if(!is.null(input$reconSub_network_brush)){
              bpDfTmp  <- brushedPoints(corrScaledLayout, input$corr_network_brush, xvar='xvar', yvar='yvar')
            } else {
              bpDfTmp <- corrScaledLayout
            }
            
            # if any commented then change colour
            if(!is.null(input$hot)){
              metIDcomments <- hot_to_r(input$hot)
            } 
            corrNetMatch <- match(bpDfTmp[, 3],  gsub('.+_',  '', metIDcomments$compSpectrum))
            bpDfTmp <- as.data.frame(bpDfTmp[, 3:ncol(bpDfTmp), drop=FALSE]) 
            bpDfTmp$possible_identity <- as.character(metIDcomments$possible_identity[corrNetMatch])
            bpDfTmp$compound_class <- as.character(metIDcomments$compound_class[corrNetMatch])
            selFeatIndx <- corrNetMatchIndx %in% feat.indx
            if(any(selFeatIndx) & input$firstNeigh == TRUE){
            # id first neighbours
            neighSel <- neighbors(corrNetTmp, which(selFeatIndx), mode='all')
            neighSel <- c(Features.v[feat.indx], names(neighSel))
            neighIdx <- bpDfTmp[, 1] %in% gsub('.+_',  '', neighSel)
            bpDfTmp <- bpDfTmp[neighIdx, , drop=FALSE]
            }
            return(bpDfTmp)
          } else {
            data.frame(result='?metID.corrNetwork function has not yet been run.')
          }
        }, rownames=FALSE, options = list(pageLength = 20))
               
      output$corr_network_plot <- shiny::renderPlot({
        
        if(!is.null(object@network$corrNetworkGraph)){
            if(!is.null(input$corr_network_brush)){
              xlimTmp <- c(input$corr_network_brush$xmin, input$corr_network_brush$xmax)
              ylimTmp <- c(input$corr_network_brush$ymin, input$corr_network_brush$ymax)
            } else {
              xlimTmp <- c(-1, 1)
              ylimTmp <- c(-1, 1)
            } 
            coordTmp <- NULL
            vertexSizeSub <- igraph::V(corrNetTmp)$vertexSize
            MS2netColsSub <- igraph::V(corrNetTmp)$MS2netColours
            vertexShapesSub <- igraph::V(corrNetTmp)$vertexShapes
            edgeNames <- attr(igraph::E(corrNetTmp), 'vnames')
            edgeColours <- rep("gray33", length(edgeNames))
            firstNeighSubIdx <- vector('logical', length(vertexShapesSub))
            vertexLabCols <- rep("gray83", length(firstNeighSubIdx))
            vertFrameCols <- rep("black", length(firstNeighSubIdx))
            # if any commented then change colour
            if(!is.null(input$hot)){
              metIDcomments <- hot_to_r(input$hot)
            } 
            
            compMSCommented <- apply(metIDcomments[, 'possible_identity', drop=FALSE], 1, function(x) any(x != '')) 
            corrNetMatchCommented <- match(igraph::V(corrNetTmp)$name, metIDcomments$compSpectrum)
            compMSCommented <- corrNetMatchCommented %in% which(compMSCommented)
            if(any(compMSCommented)){
            MS2netColsSub[compMSCommented] <- "#C77CFF"  
            }
            # colour selected features
            selFeatIndx <- corrNetMatchIndx %in% feat.indx
            if(any(selFeatIndx)){
            vertexSizeSub[selFeatIndx] <- 6
            MS2netColsSub[selFeatIndx] <- "#7CAE00"
            firstNeighSubIdx[selFeatIndx] <- TRUE
            coordTmp <- corrScaledLayout$xvar[selFeatIndx]
            coordTmp <- c(coordTmp, corrScaledLayout$yvar[selFeatIndx])
            # id first neighbours
            neighSel <- neighbors(corrNetTmp, which(selFeatIndx), mode='all')
            MS2netColsSub[neighSel] <- "#7CAE00"
            if(input$firstNeigh == TRUE){
            firstNeighSubIdx[neighSel] <- TRUE
            subGraphTmp <- induced_subgraph(corrNetTmp, c(Features.v[feat.indx], names(neighSel)))
            namesSubGraph <- attr(igraph::E(subGraphTmp), 'vnames')
            edgeColours[{edgeNames %in% namesSubGraph} == FALSE] <- '#00000000' 
            } else {
            firstNeighSubIdx[1:length(firstNeighSubIdx)] <- TRUE
            }
            } else {
            firstNeighSubIdx[1:length(firstNeighSubIdx)] <- TRUE
            }
            
            # change colour of nodes to match background if first neighbor only
            MS2netColsSub[firstNeighSubIdx == FALSE] <- '#00000000' 
            vertexLabCols[firstNeighSubIdx == FALSE] <- '#00000000'
            vertFrameCols[firstNeighSubIdx == FALSE] <- '#00000000'
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
             
            plot(corrNetTmp, layout=as.matrix(corrLayoutTmp[, 1:2]), edge.arrow.size=.1, edge.color=edgeColours, vertex.color=MS2netColsSub, vertex.label.font=2, vertex.label.color= vertexLabCols, vertex.label=gsub('CC_', '', V(corrNetTmp)$name), vertex.shape=vertexShapesSub, vertex.frame.color=vertFrameCols, vertex.size=vertexSizeSub, vertex.label.cex=1.5, xlim=xlimTmp, ylim=ylimTmp)
            # Now set the plot region to grey
            
            #layout=layout.circle,
            legend('topleft', c("currently selected spectrum and 1st neighbours (if present)", 'currently subset EIC', "MS2 matched EIC", "unmatched EIC", 'already commented', paste0('edge corrCoeff >= ', round(object@Parameters$corrThresh, 2))), pch=c(NA, 22, 21, 21, 21, NA), lty=c(1, NA, NA, NA, NA, 1), lwd=c(4, 1, 1, 1, 1, 4),
                   col=c("#7CAE00", "black", "black", "black", "black", 'gray33'), pt.bg=c("#7CAE00", "#D55E00", "#D55E00", "#0072B2", "#C77CFF", "gray33"), pt.cex=4, cex=1.8, bg='gray79', ncol=1)#text.col='white', bty="n",
            # if necc add node location
            if(!is.null(coordTmp)){
              abline(h=rep(coordTmp[2], nrow(corrScaledLayout)), col="#7CAE00", 
                     lwd=1.5, lty=2)
              abline(v=rep(coordTmp[1], nrow(corrScaledLayout)), col="#7CAE00", 
                     lwd=1.5, lty=2)
            }
            
      } else {
        plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
        text(x = 0.5, y = 0.5, paste('?metID.corrNetwork function has not yet been run.'), cex = 1.6, col = "black")
      }
      })
        
      ###############################################
      ##### 12. spectral similarity network plot ####
      ############################################### 
      output$specSimNodesEdges <- shiny::renderText({
        if(!is.null(object@network$specSimGraph)){
        paste0('nodes: ', length(igraph::V(specSimNetTmp)), ' edges: ', length(igraph::E(specSimNetTmp)), ' (N.B. large numbers of nodes e.g. >= 300 may not display properly)')
        } else {
        return('?metID.specSimNetwork function has not yet been run.')  
        }
      })
      
      output$specSimNetworkTableBrush <- DT::renderDataTable({
        if(!is.null(object@network$specSimGraph)){
          if(!is.null(input$reconSub_network_brush)){
            bpDfTmp  <- brushedPoints(specSimScaledLayout, input$specSim_network_brush, xvar='xvar', yvar='yvar')
          } else {
            bpDfTmp <- specSimScaledLayout
          }
        # if any commented then change colour
        if(!is.null(input$hot)){
          metIDcomments <- hot_to_r(input$hot)
        } 
        specSimNetMatch <- match(bpDfTmp[, 3],  metIDcomments$compSpectrum)
        bpDfTmp <- as.data.frame(bpDfTmp[, 3:5, drop=FALSE]) 
        bpDfTmp$possible_identity <- as.character(metIDcomments$possible_identity[specSimNetMatch])
        bpDfTmp$compound_class <- as.character(metIDcomments$compound_class[specSimNetMatch])
        selFeatIndx <- scIdx %in% feat.indx
        if(any(selFeatIndx) & input$firstNeigh == TRUE){
          # id first neighbours
          neighSel <- neighbors(specSimNetTmp, which(selFeatIndx), mode='all')
          neighSel <- c(Features.v[feat.indx], names(neighSel))
          neighIdx <- bpDfTmp[, 1] %in% neighSel
          bpDfTmp <- bpDfTmp[neighIdx, , drop=FALSE]
        }
        return(bpDfTmp)
        } else {
          data.frame(result='?metID.specSimNetwork function has not yet been run.')
        }
      }, rownames=FALSE, options = list(pageLength = 20))
      
      output$specSim_network_plot <- shiny::renderPlot({
        if(!is.null(object@network$specSimGraph)){
         if(!is.null(input$specSim_network_brush)){
            xlimTmp <- c(input$specSim_network_brush$xmin, input$specSim_network_brush$xmax)
            ylimTmp <- c(input$specSim_network_brush$ymin, input$specSim_network_brush$ymax)
          } else {
            xlimTmp <- c(-1, 1)
            ylimTmp <- c(-1, 1)
          } 
          coordTmp <- NULL
          vertexSizeSub <- igraph::V(specSimNetTmp)$vertexSize
          MS2netColsSub <- igraph::V(specSimNetTmp)$MS2netColours
          vertexShapesSub <- igraph::V(specSimNetTmp)$vertexShapes
          edgeNames <- attr(igraph::E(specSimNetTmp), 'vnames')
          edgeColours <- igraph::E(specSimNetTmp)$color
          firstNeighSubIdx <- vector('logical', length(vertexShapesSub))
          vertexLabCols <- rep("gray83", length(firstNeighSubIdx))
          vertFrameCols <- rep("black", length(firstNeighSubIdx))
          # if any commented then change colour
          if (!is.null(input$hot)){
            metIDcomments <- hot_to_r(input$hot)
          } 
          
          compMSCommented <- apply(metIDcomments[, 'possible_identity', drop=FALSE], 1, function(x) any(x != '')) 
          specSimMatchCommented <- match(igraph::V(specSimNetTmp)$name, metIDcomments$compSpectrum)
          compMSCommented <- specSimMatchCommented %in% which(compMSCommented)
          if(any(compMSCommented)){
            MS2netColsSub[compMSCommented] <- "#C77CFF"  
          }
          # colour selected features
          selFeatIndx <- scIdx %in% feat.indx
          if(any(selFeatIndx)){
            vertexSizeSub[selFeatIndx] <- 6
            MS2netColsSub[selFeatIndx] <- "#7CAE00"
            firstNeighSubIdx[selFeatIndx] <- TRUE
            coordTmp <- specSimScaledLayout$xvar[selFeatIndx]
            coordTmp <- c(coordTmp, specSimScaledLayout$yvar[selFeatIndx])
            # id first neighbours
            neighSel <- neighbors(specSimNetTmp, which(selFeatIndx), mode='all')
            MS2netColsSub[neighSel] <- "#7CAE00"
            if(input$firstNeigh == TRUE){
              firstNeighSubIdx[neighSel] <- TRUE
              subGraphTmp <- induced_subgraph(specSimNetTmp, c(Features.v[feat.indx], names(neighSel)))
              namesSubGraph <- attr(igraph::E(subGraphTmp), 'vnames')
              edgeColours[{edgeNames %in% namesSubGraph} == FALSE] <- '#00000000' 
            } else {
              firstNeighSubIdx[1:length(firstNeighSubIdx)] <- TRUE
            }
          } else {
            firstNeighSubIdx[1:length(firstNeighSubIdx)] <- TRUE
          }
          
          # change colour of nodes to match background if first neighbor only
          MS2netColsSub[firstNeighSubIdx == FALSE] <- '#00000000' 
          vertexLabCols[firstNeighSubIdx == FALSE] <- '#00000000'
          vertFrameCols[firstNeighSubIdx == FALSE] <- '#00000000'
          # highlight subset features as triangles
          if(input$goButton){
            Featurenames <- shiny::isolate({Featureselection()})
          }
          if(input$DBbutton){
            Featurenames <- shiny::isolate({DBFeatureselection()})
          }
          subsetFeatures <- which(Features.v %in% Featurenames)
          subsetFeatures <- scIdx %in% subsetFeatures
          if(any(subsetFeatures)){
            vertexShapesSub[subsetFeatures] <- 'csquare'  
          }
          
          plot(specSimNetTmp, layout=specSimLayoutTmp[, 1:2], edge.arrow.size=.1, edge.color=edgeColours, vertex.color=MS2netColsSub, vertex.label.font=2, vertex.label.color= vertexLabCols, vertex.label=gsub('CC_', '', V(specSimNetTmp)$name), vertex.frame.color=vertFrameCols, vertex.shape=vertexShapesSub, vertex.size=vertexSizeSub, vertex.label.cex=1.5, xlim=xlimTmp, ylim=ylimTmp) #layout=layout.circle,
          legend('topleft', c("currently selected spectrum and 1st neighbours (if present)", 'currently subset EIC', "MS2 matched EIC", 'already commented', paste0('edge fragment ions (dot product >= ', round(object@Parameters$minDotProdThresh, 2), ')'), paste0('edge neutral losses (dot product >= ', round(object@Parameters$minDotProdThresh, 2), ')')), pch=c(NA, 22, 21, 21, NA, NA), lty=c(1, NA, NA, NA, 1, 1), lwd=c(4, 1, 1, 1, 4, 4),
                 col=c("#7CAE00", "black", "black", "black", "#CC79A7", "#56B4E9"), pt.bg=c("#7CAE00", "#D55E00", "#D55E00", "#C77CFF", "#CC79A7", "#56B4E9"), pt.cex=4, cex=1.8, bg='gray79', ncol=1)#text.col='white', bty="n",
          # if necc add node location
          if(!is.null(coordTmp)){
            abline(h=rep(coordTmp[2], nrow(specSimScaledLayout)), col="#7CAE00", 
                   lwd=1.5, lty=2)
            abline(v=rep(coordTmp[1], nrow(specSimScaledLayout)), col="#7CAE00", 
                   lwd=1.5, lty=2)
          }
        } else {
          plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
          text(x = 0.5, y = 0.5, paste('?metID.specSimNetwork function has not yet been run.'), cex = 1.6, col = "black")
        }
      })
      
      ######################################################
      ##### 12. reconstructed substructure network plot ####
      ###################################################### 
      output$reconSubNodesEdges <- shiny::renderText({
        if(!is.null(object@network$reconSubGraph)){
          paste0('nodes: ', length(igraph::V(reconSubNetTmp)), ' edges: ', length(igraph::E(reconSubNetTmp)), ' (N.B. large numbers of nodes e.g. >= 300 may not display properly)')
        } else {
          return('?metID.reconSubNetwork function has not yet been run.')  
        }
      })
      
      output$reconSubNetworkTableBrush <- DT::renderDataTable({
        if(!is.null(object@network$reconSubGraph)){
          if(!is.null(input$reconSub_network_brush)){
          bpDfTmp  <- brushedPoints(reconSubScaledLayout, input$reconSub_network_brush, xvar='xvar', yvar='yvar')
          } else {
            bpDfTmp <- reconSubScaledLayout
          }
          # if any commented then change colour
          if(!is.null(input$hot)){
            metIDcomments <- hot_to_r(input$hot)
          } 
          reconSubNetMatch <- match(bpDfTmp[, 3],  metIDcomments$compSpectrum)
          bpDfTmp <- as.data.frame(bpDfTmp[, 3:5, drop=FALSE]) 
          bpDfTmp$possible_identity <- as.character(metIDcomments$possible_identity[reconSubNetMatch])
          bpDfTmp$compound_class <- as.character(metIDcomments$compound_class[reconSubNetMatch])
          selFeatIndx <- rcIdx %in% feat.indx
          if(any(selFeatIndx) & input$firstNeigh == TRUE){
            # id first neighbours
            neighSel <- neighbors(reconSubNetTmp, which(selFeatIndx), mode='all')
            fragsNLTmp <- names(neighSel)
            # all neighbours of the frags/neutral losses
            for(i in 1:length(neighSel)){
            fragsNLTmp <- c(fragsNLTmp, names(neighbors(reconSubNetTmp, neighSel[i], mode = 'all')))
            }
            fragsNLTmp <- unique(fragsNLTmp)
            fragsNLTmp <- c(Features.v[feat.indx], fragsNLTmp)
            neighIdx <- bpDfTmp[, 1] %in% fragsNLTmp
            bpDfTmp <- bpDfTmp[neighIdx, , drop=FALSE]
          }
          return(bpDfTmp)
        } else {
          data.frame(result='?metID.reconSubNetwork function has not yet been run.')
        }
      }, rownames=FALSE, options = list(pageLength = 20))
      
      # reactive plot test
      # reconSubNetPlot <- reactive({
      # })
      # 
      output$reconSub_network_plot <- shiny::renderPlot({
        if(!is.null(object@network$reconSubGraph)){
          if(!is.null(input$reconSub_network_brush)){
            xlimTmp <- c(input$reconSub_network_brush$xmin, input$reconSub_network_brush$xmax)
            ylimTmp <- c(input$reconSub_network_brush$ymin, input$reconSub_network_brush$ymax)
          } else {
            xlimTmp <- c(-1, 1)
            ylimTmp <- c(-1, 1)
          } 
        
          coordTmp <- NULL
          vertexSizeSub <- igraph::V(reconSubNetTmp)$sizeNDetects 
          MS2netColsSub <- igraph::V(reconSubNetTmp)$color
          vertexShapesSub <- igraph::V(reconSubNetTmp)$vertexShapes
          edgeNames <- attr(igraph::E(reconSubNetTmp), 'vnames')
          edgeColours <- rep('gray83')
          firstNeighSubIdx <- vector('logical', length(vertexShapesSub))
          vertexLabCols <- rep("gray83", length(firstNeighSubIdx))
          vertFrameCols <- rep("black", length(firstNeighSubIdx))
          
          # colour selected features
          selFeatIndx <- rcIdx %in% feat.indx
          if(any(selFeatIndx)){
            MS2netColsSub[selFeatIndx] <- "#7CAE00"
            firstNeighSubIdx[selFeatIndx] <- TRUE
            coordTmp <- reconSubScaledLayout$xvar[selFeatIndx]
            coordTmp <- c(coordTmp, reconSubScaledLayout$yvar[selFeatIndx])
            # id first neighbours
            neighSel <- neighbors(reconSubNetTmp, which(selFeatIndx), mode='all')
            firstNeighSubIdx[neighSel] <- TRUE
            fragsNLTmp <- as.numeric()
            # all neighbours of the frags/neutral losses
            for(i in 1:length(neighSel)){
              fragsNLTmp <- c(fragsNLTmp, neighbors(reconSubNetTmp, neighSel[i], mode = 'all'))
            }
            fragsNLTmp <- unique(fragsNLTmp)
            MS2netColsSub[fragsNLTmp] <- "#7CAE00"
            
            if(input$firstNeigh == TRUE){
              firstNeighSubIdx[fragsNLTmp] <- TRUE
              
              subGraphTmp <- induced_subgraph(reconSubNetTmp, c(Features.v[feat.indx], names(fragsNLTmp)))
              namesSubGraph <- attr(igraph::E(subGraphTmp), 'vnames')
              edgeColours[{edgeNames %in% namesSubGraph} == FALSE] <- '#00000000' 
            } else {
              firstNeighSubIdx[1:length(firstNeighSubIdx)] <- TRUE
            }
          } else {
            firstNeighSubIdx[1:length(firstNeighSubIdx)] <- TRUE
          }
          
          # change colour of nodes to match background if first neighbor only
          MS2netColsSub[firstNeighSubIdx == FALSE] <- '#00000000' 
          vertexLabCols[firstNeighSubIdx == FALSE] <- '#00000000'
          vertFrameCols[firstNeighSubIdx == FALSE] <- '#00000000'
          
          par(bg='black')
          plot(reconSubNetTmp, layout=reconSubLayoutTmp[, 1:2], edge.arrow.size=.1, edge.color=edgeColours, vertex.color=MS2netColsSub, vertex.label.font=2, vertex.label.color= vertexLabCols, vertex.label=gsub('CC_', '', V(reconSubNetTmp)$name), vertex.frame.color=vertFrameCols, vertex.shape=vertexShapesSub, vertex.size=vertexSizeSub, vertex.label.cex=1.5, xlim=xlimTmp, ylim=ylimTmp) #layout=layout.circle,
          legend('topleft', c("currently selected spectrum and 1st neighbours (if present)", 'spectrum', "fragment", 'neutral loss'), pch=c(21, 21, 21, 21), lwd=c(1, 1, 1, 1),
                 col=c("black", "black", "black", 'black'), pt.bg=c("#7CAE00", "#56B4E9", "#009E73", "#D55E00"), pt.cex=4, cex=1.8, bg='gray79', ncol=1)
          # if necc add node location
          if(!is.null(coordTmp)){
            abline(h=rep(coordTmp[2], nrow(reconSubScaledLayout)), col="#7CAE00", 
                   lwd=1.5, lty=2)
            abline(v=rep(coordTmp[1], nrow(reconSubScaledLayout)), col="#7CAE00", 
                   lwd=1.5, lty=2)
          }
        } else {
          plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
          text(x = 0.5, y = 0.5, paste('?metID.reconSubNetwork function has not yet been run.'), cex = 1.6, col = "black")
        }
      })
        ##################################
        ##### 11. metID comments table ###
        ##################################
        output$hot = renderRHandsontable({
          # add anno to comments
         
            if(!is.null(input$hot)){
            metIDcomments <- hot_to_r(input$hot)
            }
            # if(input$addAnno){
              if(!is.null(input$BestCandidates_rows_selected)){
                DB.results <- tmp.BestAnno[[feat.indx]]
                # sort by bc score if present
                if(!is.null(DB.results$BC_meanScore)){
                  DB.results <- DB.results[order(DB.results$BC_meanScore, decreasing = TRUE), , drop=FALSE]
                }
                selectedMetab <- as.character(DB.results[input$BestCandidates_rows_selected, 'DBname'])
                selectedESI <- as.character(DB.results[input$BestCandidates_rows_selected, 'ESI_type'])
                commIndx <- metIDcomments$compSpectrum %in% Features.v[feat.indx]
                # empty cell
                commIndx <- commIndx & {metIDcomments$possible_identity == ''}
                metIDcomments$possible_identity[commIndx] <- selectedMetab
                metIDcomments$ESI_type[commIndx] <- selectedESI
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
              metIDcomments <- metIDcomments[, c("compSpectrum", 'currently_subset', "possible_identity", "ESI_type", "compound_class", "user_comments")]
              # order by currently subset
              metIDcomments <- metIDcomments[order(currSubset, decreasing = T), , drop=F]
              if(all(currSubset == 'Yes')){
                # resort to Features.v order
                metIDcomments <- metIDcomments[order(as.numeric(gsub('.+_', '', metIDcomments$compSpectrum))), ]  
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
          splashCodeIndx <- grep('splash', metadata.tmp.sub)
          if(length(splashCodeIndx) > 0){
          splashCodeTmp <- metadata.tmp.sub[[splashCodeIndx]]
          metadata.tmp.sub[[splashCodeIndx]] <- NULL
          }
          
         metadata.tmp.sub <- tapply(metadata.tmp.sub, names(metadata.tmp.sub), function(x){
            x <- unlist(x)
            if(any(grepl(';', x))){
              x <- paste(x, collapse = '; ') 
            } else { 
              if(any(grepl('(.*\\.)', x))){
                if(any(nchar(gsub("(.*\\.)",  "",  as.character(x))) > 4)){
                  x <- as.character(round(as.numeric(x),  digits = 4))
                }
              }
              x <- paste(x,  collapse = "; ")
            }   
          })
         
          row.names.tmp <- gsub(".+\\.mzXML_[0-9]*_",  "",  names(metadata.tmp.sub))
          #   gsub("_[:alpha:].+",  "",  names(metadata))
          column.names.tmp <- gsub(gsub("\\.",  "\\\\.",  
                                        paste(unique(row.names.tmp),  collapse = "$|")), 
                                   "",  names(metadata.tmp.sub))
          column.names.tmp <- gsub("_$",  "",  column.names.tmp)
          
          metadata.tmp.sub <- do.call(cbind, lapply(unique(column.names.tmp),  function(x) 
            metadata.tmp.sub[grep(x,  names(metadata.tmp.sub))]))
          row.names(metadata.tmp.sub) <- unique(row.names.tmp)
          colnames(metadata.tmp.sub) <- unique(column.names.tmp)
          if(all(colnames(metadata.tmp.sub) %in% 'V1')){
            colnames(metadata.tmp.sub) <- 'Parameter'
          }
          # remove unnecessary rows
          remIndx <- grep('MS\\.scantype|TICabovefilter|MS2TICfilt.Indx', row.names(metadata.tmp.sub), 
                          ignore.case = TRUE)
          metadata.tmp.sub <- metadata.tmp.sub[-remIndx, , drop=FALSE]
          # only 1 value needed
          oneOnlyIndx <- grepl('MS1_|ppmDiff', row.names(metadata.tmp.sub))
          singleValues <- t(apply(metadata.tmp.sub[oneOnlyIndx, , drop=FALSE], 1, function(x){
            splitTmp <- strsplit(x, ';')
            singVal <- vector('character', length(splitTmp))
            for(sp in 1:length(splitTmp)){
            singVal[sp] <- splitTmp[[sp]][1]
            }
            return(singVal)
          }))
          metadata.tmp.sub[oneOnlyIndx, ] <- singleValues 
          if(length(splashCodeIndx) > 0){
          splashCodeTmp <-  matrix(splashCodeTmp, ncol=ncol(metadata.tmp.sub), nrow=1)
          colnames(splashCodeTmp) <- colnames(metadata.tmp.sub)
          row.names(splashCodeTmp) <- 'compSpectrum_splashCode'
          metadata.tmp.sub <- rbind(metadata.tmp.sub, splashCodeTmp)
          }
          metadata.tmp.sub <- metadata.tmp.sub[-grep('isoMass|isoInt', row.names(metadata.tmp.sub)), , drop=FALSE]
          return(metadata.tmp.sub)    
        },  options = list(pageLength = 20))
        
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
          if(is.null(DB.results)){
            DB.results <- data.frame("no annotations have yet been made for this feature", stringsAsFactors=F)
            colnames(DB.results) <- "Result"
            return(DB.results)  
          } else if(nrow(DB.results) == 0){
          DB.results <- data.frame("no annotations were made to any database", stringsAsFactors=F)
          colnames(DB.results) <- "Result"
          return(DB.results)
            } else {
              # retention time prediction model remove molecular descriptors if necessary
              DB.results <- DB.results[, grepl('^MD_|trainingSet|chemFP', colnames(DB.results)) == FALSE, drop=FALSE]
            ###create clickable url links to DB in table
            DB.results$DBid  <- paste0("<a href='http://", DB.results$WebAddress,  gsub('_.+', '', DB.results$DBid), ifelse(grepl('chemspider', DB.results$WebAddress), ".html", ""), "' target='_blank'>",  DB.results$DBid,  "</a>")
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
            DB.results$BestAnno <- NULL
            DB.results$chemFP <- NULL
            numIndxTmp <- grep('observedMass', colnames(DB.results)):ncol(DB.results)
            DB.results[, numIndxTmp] <- lapply(DB.results[, numIndxTmp, drop=F],  as.numeric)
            DB.results[, numIndxTmp] <- apply(DB.results[, numIndxTmp, drop=F], 2,  round, 4)
            return(DB.results)
          }}}, rownames=F,  escape=F,  options = list(pageLength = 20))#, digits=4,  sanitize.text.function = function(x) x)
        
        ##############################
        ##### 5. Best annotations #####
        ##############################
        
        output$BestCandidates <- DT::renderDataTable({
          if(length(BestAnno(object))  ==  0){
            DB.results <- data.frame("metID.dbProb function has not yet been run", stringsAsFactors=F)
            colnames(DB.results) <- "Result"
            return(DB.results) 
          } else {
            DB.results <- tmp.BestAnno[[feat.indx]]
            if(is.null(DB.results)){
              DB.results <- data.frame("no annotations have yet been made for this feature", stringsAsFactors=FALSE)
              colnames(DB.results) <- "Result"
              return(DB.results)  
            } else if(nrow(DB.results) == 0){
              DB.results <- data.frame("no best/ most likely annotations based on substructures identified", stringsAsFactors=F)
              colnames(DB.results) <- "Result"
              return(DB.results) 
            } else {
              # retention time prediction model remove molecular descriptors if necessary
              DB.results <- DB.results[, grepl('^MD_|trainingSet|chemFP', colnames(DB.results)) == FALSE, drop=FALSE]
            ###create clickable url links to DB in table
            DB.results$DBid  <- paste0("<a href='http://", DB.results$WebAddress,  gsub('_.+', '', DB.results$DBid), ifelse(grepl('chemspider', DB.results$WebAddress), ".html", ""), 
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
            DB.results$BestAnno <- NULL
            numIndxTmp <- grep('observedMass', colnames(DB.results)):ncol(DB.results)
            DB.results[, numIndxTmp] <- lapply(DB.results[, numIndxTmp, drop=F],  as.numeric)
            DB.results[, numIndxTmp] <- apply(DB.results[, numIndxTmp, drop=F], 2,  round, 4)
            # sort by bc score if present
            if(!is.null(DB.results$BC_meanScore)){
              DB.results <- DB.results[order(DB.results$BC_meanScore, decreasing = TRUE), , drop=FALSE]
            }
            return(DB.results)
          }}}, selection='single', rownames=F,  escape=F,  options = list(pageLength = 20))#, digits=4,  sanitize.text.function = function(x) x)
        
        ################################################
        ##### 12. random forest rt prediction plot #####
        ################################################
        output$rtPredPlot <- shiny::renderPlot({
          # retention time prediction model
          if(any(rtPredIndx)){
            if(!is.null(input$rtPredPlot_brush)){
              xlimTmp <- c(input$rtPredPlot_brush$xmin, input$rtPredPlot_brush$xmax)
              ylimTmp <- c(input$rtPredPlot_brush$ymin, input$rtPredPlot_brush$ymax)
            } else {
              xlimTmp <- c(min(bestAnnoDf$ms1Rt) * 0.99, max(bestAnnoDf$ms1Rt) * 1.01)
              ylimTmp <- c(min(bestAnnoDf$predRts) * 0.99, max(bestAnnoDf$predRts) * 1.01)
            }  
          # colour currently selected spectrum
          indxTmp <- bestAnnoDf$specNames %in% Features.v[feat.indx]
              
          plot(bestAnnoDf$ms1Rt, bestAnnoDf$predRts, xlab='MS1 feature retentionTime', ylab='randomForest predicted retentionTime', xlim = xlimTmp,
               ylim=ylimTmp, col='#0072B2', cex=3, pch=19, cex.axis=1.3, 
               cex.lab=1.3)
          text(bestAnnoDf$ms1Rt, bestAnnoDf$predRts, labels = bestAnnoDf$specNames,
               pos=3, col="#999999", cex=0.8)
          points(bestAnnoDf$ms1Rt[indxTmp], bestAnnoDf$predRts[indxTmp], col='red', cex=3, pch=19)
          points(bestAnnoDf$ms1Rt[as.logical(bestAnnoDf$trainingSet)], bestAnnoDf$predRts[as.logical(bestAnnoDf$trainingSet)], col='#E69F00', cex=3, pch=19)
          
          abline(lmPredRt, col='black')
          legend('topleft', c('trainingSet', 'unknown', 'currently selected'), 
                 pch=21, 
                 pt.bg = c("#E69F00", "#0072B2", "red"), pt.cex = 3, cex=1.3)
          legend("topright", bty="n", legend=paste("r2:", format(summary(lmPredRt)$adj.r.squared, digits=3)), cex=1.3)
          } else {
            plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
            text(x = 0.5, y = 0.5, paste("?metID.rtPred function has not been run.\n"), cex = 1.6, col = "black")
          }  
        })
        
        output$rtPredTable <- DT::renderDataTable({
          if(any(rtPredIndx)){
          brushedPoints(bestAnnoDf, input$rtPredPlot_brush, xvar = "ms1Rt", yvar = "predRts") 
          } else {
          return(data.frame(rtPredTable="?metID.rtPred function has not been run."))           }
        }, rownames=F, escape=F)
        
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
          }}, rownames=F,  escape=F,  options = list(pageLength = 8))
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
        if(!is.null(input$pubMedSearchButton)){ 
        if(input$pubMedSearchButton == 0){ 
         return() } else if(input$pubMedSearchButton > 0){
        ClAbs <- shiny::isolate(pubMedMine()[[1]])
        if(nrow(ClAbs) == 1){
          PubMedWordcloud::plotWordCloud(ClAbs,  min.freq=1,  max.words=100,  rot.per=0)     } else {
          suppressWarnings(PubMedWordcloud::plotWordCloud(ClAbs, max.words=100, scale=c(4, 0.5)))
          }
        }
        }})
        
        ###########################################
        ##### 7. PubMed random article table  #####
        ###########################################
        
        output$WordCloudTable <- DT::renderDataTable({
          if(!is.null(input$pubMedSearchButton)){ 
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
             }} else {
               wordCloudDf <- data.frame("No Search of PubMed performed")
               colnames(wordCloudDf) <- input$FeatureNames
               return(wordCloudDf)   
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
        
        #######################################
        ##### 8. output in silico results #####
        #######################################
        # reactive metfragtable
        metFragTable <- reactive({
          if(length(tmp.metFrag)  ==  0){
            metFrag.df.tmp <- "metID.MetFrag function has not yet been run"
            metFrag.df.tmp <- data.frame(metFrag.df.tmp)
            colnames(metFrag.df.tmp) <- "Result"
            return(list(metFragSelectDf=metFrag.df.tmp))
          } else {
            metFrag.df.tmp <- tmp.metFrag[[feat.indx]]
            if(is.null(metFrag.df.tmp)){
              metFrag.df.tmp <- "metID.MetFrag returned no results"
              metFrag.df.tmp <- data.frame(metFrag.df.tmp)
              colnames(metFrag.df.tmp) <- "Result" 
              return(list(metFragSelectDf=metFrag.df.tmp))
            } else if(ncol(metFrag.df.tmp) == 1){
              metFrag.df.tmp <- "metID.MetFrag returned no results"
              metFrag.df.tmp <- data.frame(metFrag.df.tmp)
              colnames(metFrag.df.tmp) <- "Result" 
              return(list(metFragSelectDf=metFrag.df.tmp))
            } else {  
              # reduce n decimal places to 4
              # metFrag.df.tmp[,  c("PeakScore",  "BondEnergyScore",  "Score")] <- sapply(metFrag.df.tmp[,  c("PeakScore",  "BondEnergyScore",  "Score")], function(x) {if(nchar(gsub("(.*\\.)",  "",  as.character(x[1]))) > 4){x <- round(as.numeric(x),  digits = 4)}})
              ###create clickable url links to DB in table
              metFrag.df.tmp$DBid  <- paste0("<a href='http://",  metFrag.df.tmp$WebAddress,  gsub('_.+', '', metFrag.df.tmp$DBid), ifelse(grepl('chemspider', metFrag.df.tmp$WebAddress), ".html", ""), 
                                             "' target='_blank'>",  metFrag.df.tmp$DBid,  "</a>")
              metFrag.df.tmp$WebAddress <- NULL 
              metFrag.df.tmp$SMILES <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", 
                                              metFrag.df.tmp$SMILES,  "/PNG",  "' target='_blank'>", 
                                              paste0(substring(metFrag.df.tmp$SMILES,  1,  6),  "..."),  "</a>")
              names(metFrag.df.tmp)[names(metFrag.df.tmp) == 'NoExplPeaks'] <- 'nPeaksEx'
              names(metFrag.df.tmp)[names(metFrag.df.tmp) == 'propIntEx'] <- 'metFrag_totPropEx'
              metFrag.df.tmp$Score <- round(metFrag.df.tmp$Score, 2)
              return(list(metFragSelectDf=metFrag.df.tmp[, c("DBid", "DBname", "metFrag_totPropEx", "nPeaksEx", 'SMILES', 'ESI_type', 'Score', 'mass', 'intensity', 'frag_smiles')]))
            }
          }
        })
        
        output$metFragTable <- DT::renderDataTable({
        metFragSelectDf <-  metFragTable()$metFragSelectDf
        metFragSelectDf$mass <- NULL
        metFragSelectDf$intensity <- NULL
        metFragSelectDf$frag_smiles <- NULL
        return(metFragSelectDf) 
        }, selection='single', escape=F, options = list(pageLength = 5), rownames=F)
        
        
        # MetFrag head-to-tail plot
        output$metFragPlot <- shiny::renderPlot({
          metFragTableTmp <- metFragTable()$metFragSelectDf
          if(ncol(metFragTableTmp) > 1){
          if(!is.null(input$metFragTable_rows_selected)){
            plotDfTmp <- data.frame(composite_spectra[[feat.indx]],  stringsAsFactors = F)[, c('mass', 'Rel_Intensity', 'intensity')]
            plotDfTmp$dbData <- 0
            metFragEntryTmp <- metFragTableTmp[input$metFragTable_rows_selected, c('mass', 'intensity'), drop=F]
            if(!is.na(metFragEntryTmp$mass)){
            massTmp <- as.numeric(strsplit(as.character(metFragEntryTmp$mass), '; ')[[1]])
            intTmp <- as.numeric(strsplit(as.character(metFragEntryTmp$intensity), '; ')[[1]])
            metFragEntryTmp <- data.frame(mass=massTmp, Rel_Intensity=intTmp, stringsAsFactors = F)  
            if(nrow(metFragEntryTmp) > 0){
            metFragEntryTmp$Rel_Intensity <- -(metFragEntryTmp$Rel_Intensity / max(plotDfTmp$intensity) * 100)
            metFragEntryTmp$dbData <- 1
            plotDfTmp$intensity <- NULL
            # rbind
            plotDfTmp <- rbind(plotDfTmp, metFragEntryTmp)
            colsTmp <- ifelse(plotDfTmp[, 'dbData'] == 1, "#009E73", '#000000')
            
            if(!is.null(input$metFrag_brush)){
              xlimTmp <- c(input$metFrag_brush$xmin, input$metFrag_brush$xmax)
              ylimTmp <- c(input$metFrag_brush$ymin, input$metFrag_brush$ymax)
            } else {
              xlimTmp <- c(0, max(plotDfTmp$mass) * 1.1)
              ylimTmp <- c(min(plotDfTmp$Rel_Intensity), max(plotDfTmp$Rel_Intensity))
            }
            
            plot(plotDfTmp[, c('mass', 'Rel_Intensity'), drop=F], xlim=xlimTmp, ylim=ylimTmp, yaxt='n', ylab='relative intensity',
                 type='h', col=colsTmp, cex.axis=1.5, cex.lab=1.5)
            axis(2, at=seq(-100, 100, 20), labels=c(abs(seq(-100, -20, 20)), seq(0, 100, 20)), las=2, cex.axis=1.3)
            abline(h=0)
            legend('topright', c("composite spectrum", 'metFrag result'), lty=c(1, 1), lwd=c(3, 3), col=c('#000000', "#009E73"), pt.cex=4, cex=1.2, ncol=1)
          } else {
            plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
            text(x = 0.5, y = 0.5, 'no peaks were explained', cex = 1.6, col = "black")   
          }
          } else {
            plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
            text(x = 0.5, y = 0.5, 'no peaks were explained', cex = 1.6, col = "black")   
          }
          } else {
            plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
            text(x = 0.5, y = 0.5, 'click a row in the table above to view the plot', cex = 1.6, col = "black") 
          }
          } else {
            plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
            text(x = 0.5, y = 0.5, metFragTableTmp[, 1], cex = 1.6, col = "black") 
          }
        })
        # metfragment table
        output$metFragFragments <- DT::renderDataTable({
          metFragTableTmp <-  metFragTable()$metFragSelectDf
          metFragEntryTmp <- data.frame(result='select a row in the table above')
          if(ncol(metFragTableTmp) > 1){
            if(!is.null(input$metFragTable_rows_selected)){
             
              metFragEntryTmp <- metFragTableTmp[input$metFragTable_rows_selected, c('mass', 'intensity', 'frag_smiles'), drop=F]
              if(!is.na(metFragEntryTmp$mass)){
                massTmp <- as.numeric(strsplit(as.character(metFragEntryTmp$mass), '; ')[[1]])
                massTmp <- round(massTmp, 4)#[-length(massTmp)]
                # intTmp <- as.numeric(strsplit(as.character(metFragEntryTmp$intensity), '; ')[[1]])
                fragSmilesTmp <- strsplit(as.character(metFragEntryTmp$frag_smiles), '; ')[[1]]
                if(length(massTmp) > length(fragSmilesTmp)){
                  massTmp <- massTmp[1:length(fragSmilesTmp)]  
                }
                metFragEntryTmp <- data.frame(mass=massTmp, frag_formula=fragSmilesTmp, stringsAsFactors = F) 
                metFragEntryTmp <- metFragEntryTmp[order(metFragEntryTmp$mass, decreasing = T), , drop=F]
                # brushedPoints(metFragEntryTmp, input$metFrag_brush, xvar = "mass") 
               }
            }
          }
          return(metFragEntryTmp)
        }, rownames=F)
        # reactive CFM table
        cfmTable <- reactive({
          if(length(cfmSelectTable)  ==  0){
            cfmSelectDf <- "metID.CFM function has not yet been run"
            cfmSelectDf <- data.frame(cfmSelectDf)
            colnames(cfmSelectDf) <- "Result"
            return(list(cfmSelectDf=cfmSelectDf))
          } else {
            cfmSelectDf <- cfmSelectTable[[feat.indx]]
            cfmMatchesDf <- tmp.CFM[[feat.indx]]
            if(is.null(cfmSelectDf)){
              cfmSelectDf <- "metID.CFM function returned no results"
              cfmSelectDf <- data.frame(cfmSelectDf)
              colnames(cfmSelectDf) <- "Result" 
              return(list(cfmSelectDf=cfmSelectDf))
            } else {
            ###create clickable url links to DB in table
            cfmSelectDf$DBid  <- paste0("<a href='http://",  cfmSelectDf$WebAddress,  gsub('_.+', '', cfmSelectDf$DBid), ifelse(grepl('chemspider', cfmSelectDf$WebAddress), ".html", ""), 
                                             "' target='_blank'>",  cfmSelectDf$DBid,  "</a>")
              cfmSelectDf$WebAddress <- NULL 
              cfmSelectDf$CFM_totPropEx <- round(as.numeric(cfmSelectDf$CFM_totPropEx),  digits = 2)
              cfmSelectDf <- cfmSelectDf[order(cfmSelectDf$CFM_totPropEx, decreasing = T), , drop=F]
              cfmSelectDf$DBname <- ifelse(nchar(cfmSelectDf$DBname) >= 35, 
                                           paste0(substring(cfmSelectDf$DBname,  1,  35),  "..."), 
                                           cfmSelectDf$DBname)
             
              return(list(cfmSelectDf=cfmSelectDf, cfmMatchesDf=cfmMatchesDf))
            }
          }
        })
        
        output$cfmTable <- DT::renderDataTable({
          return(cfmTable()$cfmSelectDf) 
        }, selection='single', escape=F, rownames=F, options = list(pageLength = 5))
        
        
        # MetFrag head-to-tail plot
        output$cfmPlot <- shiny::renderPlot({
          cfmTableTmp <- cfmTable()
          if(ncol(cfmTableTmp$cfmSelectDf) > 1){
            if(!is.null(input$cfmTable_rows_selected)){
              plotDfTmp <- data.frame(composite_spectra[[feat.indx]],  stringsAsFactors = F)[, c('mass', 'intensity', 'Rel_Intensity')]
              plotDfTmp$dbData <- 0
             cfmEntryTmp <- cfmTableTmp$cfmSelectDf$DBid[input$cfmTable_rows_selected]
             cfmMatchesTmp <- cfmTableTmp$cfmMatchesDf[cfmTableTmp$cfmMatchesDf$DBid %in% gsub(".+'>|</a>", '', cfmEntryTmp), c('mass', 'intensity')]
             # remove duplicates
             cfmMatchesTmp <- cfmMatchesTmp[duplicated(cfmMatchesTmp$mass) == F, , drop=F]
             cfmMatchesTmp$Rel_Intensity <- -(cfmMatchesTmp$intensity / max(plotDfTmp$intensity) * 100)
             cfmMatchesTmp$dbData <- 1
                # rbind
                plotDfTmp <- rbind(plotDfTmp, cfmMatchesTmp)
                colsTmp <- ifelse(plotDfTmp[, 'dbData'] == 1, "#009E73", '#000000')
                
                if(!is.null(input$cfm_brush)){
                  xlimTmp <- c(input$cfm_brush$xmin, input$cfm_brush$xmax)
                  ylimTmp <- c(input$cfm_brush$ymin, input$cfm_brush$ymax)
                } else {
                  xlimTmp <- c(0, max(plotDfTmp$mass) * 1.1)
                  ylimTmp <- c(min(plotDfTmp$Rel_Intensity), max(plotDfTmp$Rel_Intensity))
                }
                
                plot(plotDfTmp[, c('mass', 'Rel_Intensity'), drop=F], xlim=xlimTmp, ylim=ylimTmp, yaxt='n', ylab='relative intensity',
                     type='h', col=colsTmp, cex.axis=1.5, cex.lab=1.5)
                axis(2, at=seq(-100, 100, 20), labels=c(abs(seq(-100, -20, 20)), seq(0, 100, 20)), las=2, cex.axis=1.3)
                abline(h=0)
                legend('topright', c("composite spectrum", 'CFM result'), lty=c(1, 1), lwd=c(3, 3), col=c('#000000', "#009E73"), pt.cex=4, cex=1.2, ncol=1)
            } else {
              plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
              text(x = 0.5, y = 0.5, 'click a row in the table above to view the plot', cex = 1.6, col = "black") 
            }
           
          } else {
            plot(c(0, 1), c(0, 1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
            text(x = 0.5, y = 0.5, cfmTableTmp$cfmSelectDf[, 1], cex = 1.6, col = "black") 
          }
        })
        
        # cfm fragments table
        output$cfmFragments <- DT::renderDataTable({
          cfmEntryTmp <- data.frame(result='select a row in the table above')
          cfmTableTmp <- cfmTable()
          if(ncol(cfmTableTmp$cfmSelectDf) > 1){
            if(!is.null(input$cfmTable_rows_selected)){
              
              cfmEntryTmp <- cfmTableTmp$cfmSelectDf$DBname[input$cfmTable_rows_selected]
              cfmEntryTmp <- cfmTableTmp$cfmMatchesDf[cfmTableTmp$cfmMatchesDf$DBname %in% cfmEntryTmp, c('CFM_rank', 'CFM_mass', 'CFM_fragSMILES'), drop=F]
              if(nrow(cfmEntryTmp) > 0){
                cfmEntryTmp$CFM_fragSMILES <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", 
                                                     cfmEntryTmp$CFM_fragSMILES,  "/PNG",  "' target='_blank'>", 
                       paste0(substring(cfmEntryTmp$CFM_fragSMILES,  1, 25),  "..."),  "</a>")
                cfmEntryTmp <- cfmEntryTmp[duplicated(cfmEntryTmp$CFM_rank) == F, , drop=F]
              }
            }
          }
          return(cfmEntryTmp)
        }, escape=F, rownames=F)
        ##load in any previous user comments
        ###add text area for storing user comments
        
        # output$CommentSectionHeader <- shiny::renderUI({shiny::tags$b("Add metabolite identification notes including any useful html links here...\n
        #                                                               NB. Do not forget to save your comments to the features results directory using the button on the left")})
        # output$UserComments <- shiny::renderUI({
        #   if(length(Comments(compMS2))  ==  0){
        #     PrevComment <- "No previous comments"
        #   } else if(is.null(Comments(object)[[feat.indx]])){
        #     PrevComment <- "No previous comments"
        #   } else {
        #     PrevComment <- Comments(object)[[feat.indx]]
        #     PrevComment <- as.character(PrevComment[nrow(PrevComment),  2])
        #   }
        #   
        #   shiny::tags$textarea(id="UserComments",  rows=4,  cols=60,  PrevComment)
        # })
        # 
        
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
        
        
      } else {
        currEICRv$currEIC <- ''    
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
