#' shiny viz
#' @export
setGeneric("compMS2explorer", function(object, ...) standardGeneric("compMS2explorer"))

setMethod("compMS2explorer", signature = "CompMS2", function(object, browserLaunch = TRUE){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(object@file.paths) == 0) {
    stop("The CompMS2 class file is empty")
  } else {
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    # Function to be removed
    compMS2shinyW <- function(object){
      ####GLOBAL.R
      #    library(XML)
      #    library(RCurl)
      ###feature names from compMS2 object
      
      composite_spectra <- compSpectra(object)
      Features.v <- names(composite_spectra)
      ###DB search names
      tmp.DBanno.res <- DBanno(object)
      ### best anno
      tmp.BestAnno <- BestAnno(object)
      metaData.tmp <- metaData(object)
      UserComments.v <- vector("list", length(composite_spectra))
      # best substructure anno
      subStrAnno.df <- subStrAnno(object)
      # metFrag 
      tmp.metFrag <- MetFrag(object)
      ###order by EIC number
      EICorderIndx <- order(as.numeric(gsub(".+_","",Features.v)))
      Features.v <- Features.v[EICorderIndx]
      composite_spectra <- composite_spectra[EICorderIndx]
      metaData.tmp <- metaData.tmp[EICorderIndx]
      
      names(UserComments.v) <- names(composite_spectra)
      if(length(tmp.DBanno.res) > 0){
        tmp.DBanno.res <- tmp.DBanno.res[EICorderIndx]
      }
      if(length(tmp.BestAnno) > 0){
        tmp.BestAnno <- tmp.BestAnno[EICorderIndx]
      }
      if(length(tmp.metFrag) > 0){
        tmp.metFrag<- tmp.metFrag[EICorderIndx]
      }
      # DB match df to string match
      if(length(DBanno(object)) > 0){
        DBmatches <-  t(sapply(tmp.DBanno.res, function(x){
          names.tmp <- x$DBname
          ESI_type.tmp <- x$ESI_type
          SubStr_type.tmp <- x$SubStr_type
          names.tmp <- data.frame(names.tmp, ESI_type.tmp, SubStr_type.tmp, 
                                  stringsAsFactors = F)
        }))
      } else {
        DBmatches <- matrix("The metID.dbAnnotate function has not yet been run", 
                            ncol = 2, nrow = length(composite_spectra))
      }
      # substructures id'd
      if(!is.null(composite_spectra[[1]]$Frag.ID)){
        
        SubStr_types <- t(sapply(composite_spectra, function(x){
          Frag.ID.tmp <- paste0("Frag_", as.character(x$Frag.ID))
          Neutral.loss.tmp <- paste0("NeutLoss_", as.character(x$Neutral.loss))
          SubStrType.tmp <- data.frame(Frag.ID.tmp, Neutral.loss.tmp, 
                                       stringsAsFactors = F)
        }))
      } else {
        SubStr_types <- matrix("The subStructure.Annotate function has not yet been run", 
                               ncol = 2, nrow = length(composite_spectra))
      }
      ###obtain SubStr types for filtration
      SubStrType.inputs <- c(unique(unlist(SubStr_types[, 1])), 
                             unique(unlist(SubStr_types[, 2])))
      
      SubStrType.inputs <- SubStrType.inputs[-grep("Frag_$|NeutLoss_$", SubStrType.inputs)]
      
      if(all(SubStrType.inputs == "The subStructure.Annotate function has not yet been run")){
        SubStrType.inputs <- "The subStructure.Annotate function has not yet been run"
      }
      ###extract mass and RT values
      mass.v <- sapply(metaData(object),function(x) as.numeric(x[grep("MS1_mz", names(x))][[1]][1]))
      RT.v <- sapply(metaData(object),function(x) as.numeric(x[grep("MS1_RT", names(x))][[1]][1]))
      TotalFeatures<-length(unique(gsub(".+_","",Features.v)))
      TotalCompSpectra<-length(Features.v)
     
      
      compMS2server <- shiny::shinyServer(function(input, output,session){ 
        
        observe({
          if(input$CloseAppBtn > 0){
            
            shiny::stopApp(UserComments.v)
          }
        })
        # feature selection shiny::reactive
        Featureselection<-shiny::reactive({
          if(input$goButton==0){ 
            return() } else if (input$goButton > 0 & input$All_Features==TRUE) {
              
              NoFeaturesIndx <- ""
              FeaturesIndx <- ""  
              
              if (any(input$NotSubStrTypes!="")) 
              {
                ###escape special characters
                NotSubStrTypes <- gsub("\\+", "\\\\+", input$NotSubStrTypes)
                NotSubStrTypes <- gsub("\\[", "\\\\[", NotSubStrTypes)
                NotSubStrTypes <- gsub("\\]", "\\\\]", NotSubStrTypes)
                NoFeaturesIndx <- grep(paste(NotSubStrTypes, collapse = "|"), SubStr_types[, 1])
                NoFeaturesIndx <- unique(c(NoFeaturesIndx, grep(paste(NotSubStrTypes, collapse = "|"), SubStr_types[, 2])))
                NoSubStrmatchFreq <- length(NoFeaturesIndx)
                #               NoFeaturesIndx <- Features.v[NoFeaturesIndx]
                if(length(NoFeaturesIndx) == 0){NoFeaturesIndx<-""}
              }
              
              if (any(input$SubStrTypes!=""))
              {
                ###escape special characters
                SubStrTypes <- gsub("\\+", "\\\\+", input$SubStrTypes)
                SubStrTypes <- gsub("\\[", "\\\\[", SubStrTypes)
                SubStrTypes <- gsub("\\]", "\\\\]", SubStrTypes)
                FeaturesIndx <- grep(paste(SubStrTypes, collapse = "|"), SubStr_types[, 1])
                FeaturesIndx <- unique(c(FeaturesIndx, grep(paste(SubStrTypes, collapse = "|"), SubStr_types[, 2])))
                SubStrmatchFreq <- length(NoFeaturesIndx)
                if(length(FeaturesIndx) == 0){FeaturesIndx <- ""}
              }
              
              # if features indices not equal 
              if(setequal(FeaturesIndx, NoFeaturesIndx)==F)
              { 
                if(!is.character(FeaturesIndx))
                {
                  FeaturesIndx <- setdiff(FeaturesIndx, NoFeaturesIndx)
                  
                  Feature.v.sub <- Features.v[FeaturesIndx]
                  
                } else {
                  
                  Feature.v.sub<-Features.v[-NoFeaturesIndx]
                }
                
              } else {
                Feature.v.sub<-Features.v
              }
              
              if(length(Feature.v.sub)==0)
              {
                return("No MS2 features found")
              } else {
                return(Feature.v.sub)}
            } 
        })
        # DB name match shiny::reactive
        DBFeatureselection<-shiny::reactive({
          
          if(input$DBbutton == 0){ return()
          } else if (input$DBbutton > 0){
            DBnamesTMP <- grep(input$DB_match_name, DBmatches[, 1], ignore.case=T)
            DBnamesTMP <- unique(c(DBnamesTMP, grep(input$DB_match_name, DBmatches[, 2],
                                                    ignore.case=T)))
            
            Feature.v.sub<-Features.v[DBnamesTMP]
            if(length(Feature.v.sub)==0)
            {
              return("No MS2 features found")
            } else {
              return(Feature.v.sub)}
          }
        })
        # feature selection observer
        observe({ if(input$goButton)
        {
          output$FeatureNames = shiny::renderUI({
            Featurenames<-shiny::isolate({Featureselection()})
            shiny::selectizeInput('FeatureNames', 'Choose a feature to plot :', 
                                  choices=Featurenames,options=list(maxOptions=10000))
          })
          
          output$matchSummaryText<-shiny::renderText({
            return("Match summary table")
          })
          
          output$matchSummary<-shiny::renderTable({
            Featurenames<-shiny::isolate({Featureselection()})
            if(Featurenames[1]!="No MS2 features found"){
              matchSummary<-data.frame(nMatches=length(unique(gsub(".+_","",Featurenames))),
                                       nCompositeSpectra=length(Featurenames))
              return(matchSummary)
            } else {
              return()
            }
          })
          output$tabbedPanelButton = shiny::renderUI({
            shiny::actionButton("tabbedPanelButton","ViewData")
          })
        }
        })
        
        # DB matches observer
        observe({ if(input$DBbutton)
        {
          output$FeatureNames = shiny::renderUI({
            Featurenames<-shiny::isolate({DBFeatureselection()})
            shiny::selectizeInput('FeatureNames', 'Choose a feature to plot :', 
                                  choices=Featurenames,options=list(maxOptions=10000))
          })
          
          
          
          output$matchSummary<-shiny::renderTable({
            Featurenames<-shiny::isolate({DBFeatureselection()})
            if(Featurenames[1]!="No MS2 features found"){
              matchSummary<-data.frame(nMatches=length(unique(gsub(".+_","",Featurenames))),
                                       nCompositeSpectra=length(Featurenames))
              return(matchSummary)
            } else {
              return()
            }
          })
          
          output$tabbedPanelButton = shiny::renderUI({
            shiny::actionButton("tabbedPanelButton","ViewData")
          })
        }
        })
        # main body observer tabbed panel
        shiny::observe({if(!is.null(input$tabbedPanelButton)){ 
          if(input$tabbedPanelButton > 0){
            if(input$FeatureNames %in% Features.v){
              
              if(!is.null(input$CommentButton))
              {
                if(input$CommentButton==0)
                {
                  vars <- shiny::reactiveValues(actionCounter=0)
                } else {
                  vars <- shiny::reactiveValues(actionCounter=0)
                  vars$actionCounter <- input$CommentButton 
                }
              }
              # isolate feature index
              feat.indx <- which(Features.v %in% shiny::isolate(input$FeatureNames))
              
              ###########################
              ##### 1. Raw data plot ####
              ###########################
              
              # output$Raw_data_plot <- rCharts::renderChart({  
              #   
              #   plot.df<-data.frame(composite_spectra[[feat.indx]], stringsAsFactors = F)
              #   colnames(plot.df)<-gsub("\\.","_",colnames(plot.df))
              #   if(!is.null(plot.df$Frag_ID))
              #   {
              #     
              #     #           plot.df[, c("mass", "intensity", "Rel_Intensity", "Precursorfrag_diff")] <- apply(plot.df[, c("mass", "intensity",  "Rel_Intensity", "Precursorfrag_diff")], 2, as.numeric)
              #     plot.df$Precursorfrag_diff<-as.numeric(plot.df$Precursorfrag_diff)
              #     plot.df$Fragment_Assigned<-ifelse(apply(plot.df[,c("Frag_ID","Neutral_loss","interfrag_loss")],1,function(x) any(x!="")),"Fragment_identified","No_Fragment_identified")
              #     r1<-rCharts::rPlot(x = "bin(mass,2)", y = "intensity", data=plot.df , type = "bar",size = list(const = 5),color="Fragment_Assigned",
              #                        tooltip = "#!function(item){ return ' intensity: ' + item.intensity + ' mass: ' + item.mass + ' Relative Intensity: ' + item.Rel_Intensity + ' Precursor Frag diff: ' + item.Precursorfrag_diff + 
              #                        ' interfrag diff: ' + item.interfrag_diff + ' Frag ID: ' + item.Frag_ID + ' Neutral Loss ID: ' + item.Neutral_loss + ' interfrag loss:  ' + item.interfrag_loss }!#")
              #     #           
              #     r1$guides(
              #       color=list(scale = "#! function(value){
              #                  color_mapping = {No_Fragment_identified: '#6495ED',Fragment_identified:'#FF3030'}
              #                  return color_mapping[value];                  
              # } !#"),
              #       x = list(min=0,max=max(plot.df$mass)+50,title = 'm/z'),
              #       y = list(title = 'intensity')
              #       )
              #     r1$set(width = 1000, height = 700)
              #     r1$set(dom="Raw_data_plot")
              #     
              #     return(r1)
              #     
              # } else {
              #   #           plot.df[, c("mass", "intensity", "Rel_Intensity", "Precursorfrag_diff")] <- apply(plot.df[, c("mass", "intensity",  "Rel_Intensity", "Precursorfrag_diff")], 2, as.numeric)
              #   r1<-rCharts::rPlot(x = "bin(mass,2)", y = "intensity", data=plot.df , type = "bar",size = list(const = 5),#,color="Fragment_Assigned",
              #                      tooltip = "#!function(item){ return ' intensity: ' + item.intensity + ' mass: ' + item.mass }!#")
              #   #           
              #   r1$guides(
              #     x = list(min=0,max=max(plot.df$mass)+50,title = 'm/z'),
              #     y = list(title = 'intensity')
              #   )
              #   r1$set(width = 1000, height = 700)
              #   r1$set(dom="Raw_data_plot")
              #   
              #   return(r1)
              # }
              # #         plot.df<-as.data.frame(composite_spectrum$composite_spectrum,stringsAsFactors=F)
              # #         colnames(plot.df)<-gsub("\\.","_",colnames(plot.df))
              # #         plot.df$MetFrag<-ifelse(plot.df[,input$plotChoice]=="","No_Fragment_identified","MetFrag_assigned")
              # #         r1<-rPlot(x = "bin(mass,2)", y = "intensity", data=plot.df , type = "bar",size = list(const = 5),color="MetFrag",
              # #                   tooltip = "#!function(item){ return ' intensity: ' + item.intensity + ' mass: ' + item.mass + ' Relative Intensity: ' + item.Rel_Intensity + ' Precursor Frag diff: ' + item.Precursorfrag_diff + 
              # #                   ' interfrag diff: ' + item.interfrag_diff + ' Frag ID: ' + item.Frag_ID + ' Neutral Loss ID: ' + item.Neutral_loss + ' interfrag loss:  ' + item.interfrag_loss }!#")
              # #         r1$guides(
              # #           color=list(scale = "#! function(value){
              # #                      color_mapping = {No_Fragment_identified: '#6495ED',MetFrag_assigned:'#BF3EFF'}
              # #                      return color_mapping[value];                  
              # #                      } !#"),
              # #           x = list(min=0,max=max(plot.df$mass)+50,title = 'm/z'),
              # #           y = list(title = 'intensity')
              # #         )
              # #         r1$set(width = 1000, height = 700)
              # #         r1$set(dom="Raw_data_plot")
              # #         
              # #         return(r1)
              # #       }
              # })
              
              ###########################
              ##### 2. metadata table ###
              ###########################
              
              output$metadata <- DT::renderDataTable({
                metadata.tmp.sub <- metaData.tmp[[feat.indx]]
                row.names.tmp <- gsub(".+\\.mzXML_[0-9]*_", "", names(metadata.tmp.sub))
                #   gsub("_[:alpha:].+", "", names(metadata))
                column.names.tmp <- gsub(gsub("\\.", "\\\\.", 
                                              paste(unique(row.names.tmp), collapse = "$|")),
                                         "", names(metadata.tmp.sub))
                column.names.tmp <- gsub("_$", "", column.names.tmp)
                metadata.tmp.sub <- sapply(metadata.tmp.sub, function(x) {
                  if(nchar(gsub("(.*\\.)", "", as.character(x[1]))) > 4){
                    tmp <- round(x, digits = 4) 
                    tmp <- paste(tmp, collapse = "; ")
                    return(tmp)
                  } else {
                    tmp <- paste(x, collapse = "; ")
                    return(tmp)
                  }
                })
                metadata.tmp.sub <- sapply(unique(column.names.tmp), function(x) 
                  metadata.tmp.sub[grep(x, names(metadata.tmp.sub))])
                row.names(metadata.tmp.sub) <- unique(row.names.tmp)
                return(metadata.tmp.sub)    
              }, options = list(pageLength = 21))
              
              ###########################
              ##### 3. MS2 data table ###
              ###########################
              
              output$MS2_data <- DT::renderDataTable({
                MS2_data<-data.frame(composite_spectra[[feat.indx]], stringsAsFactors = F)
                MS2_data[] <- lapply(MS2_data, as.character)
                if("interfrag.diff" %in% colnames(MS2_data)){
                MS2_data[,c(1:5)]<-apply(MS2_data[,c(1:5)],2,function(x) round(as.numeric(x),digits=4))
                SMILESindx<-grep("SMILES$",colnames(MS2_data))
                IDindx<-sapply(paste0(gsub("\\.SMILES","",colnames(MS2_data)[SMILESindx]),"$"),grep,colnames(MS2_data))
                
                for(i in 1:ncol(MS2_data[,SMILESindx]))
                { 
                  SMILESsubIndx<-which(MS2_data[,IDindx[i]]!="")
                  if(length(SMILESsubIndx)>0)
                  {
                    MS2_data[SMILESsubIndx,IDindx[i]]<-sapply(c(1:length(MS2_data[SMILESsubIndx,SMILESindx[i]])), function(x){
                      Smiles<-unlist(strsplit(MS2_data[SMILESsubIndx[x],SMILESindx[i]],";"))
                      SmilesInc<-which(Smiles!="")
                      SMILEShtml<-unlist(strsplit(MS2_data[SMILESsubIndx[x],IDindx[i]],";"))
                      SMILEShtml[SmilesInc]<-paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",Smiles[SmilesInc],"/PNG","' target='_blank'>",SMILEShtml[SmilesInc],"</a>")
                      SMILEShtml<-paste(SMILEShtml,collapse=" ")
                      return(SMILEShtml)
                    })
                  }
                }
                MS2_data<-MS2_data[,-SMILESindx,drop=F]
                } else {
                  MS2_data[,c(1:2)]<-apply(MS2_data[,c(1:2)],2,function(x) round(as.numeric(x),digits=4))
                  
                }
                return(MS2_data)    
              }, options = list(pageLength = 25), escape = F)
              
              ###########################
              ##### 4. DB results #####
              ###########################
              
              output$DB.results <- DT::renderDataTable({
                if(length(DBanno(object)) == 0){
                  DB.results<-data.frame("metID.dbAnnotate function has not yet been run", stringsAsFactors=F)
                  colnames(DB.results)<-"Result"
                  return(DB.results) 
                } else {
                  DB.results<-tmp.DBanno.res[[feat.indx]]
                  ###create clickable url links to DB in table
                  DB.results$DBid <-paste0("<a href='http://",DB.results$WebAddress, DB.results$DBid, 
                                           "' target='_blank'>", DB.results$DBid, "</a>")
                  # clickable SMILES pugRest
                  SMILESindx<-grep("SMILES$",colnames(DB.results))
                  DB.results[,SMILESindx]<-sapply(SMILESindx,function(x){
                    NoFuncIndx<-DB.results[, x] == ""
                    ifelse(NoFuncIndx,"",
                           paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                                  gsub("/|\\\\","", DB.results[,x]),"/PNG","' target='_blank'>", 
                                  paste0(substring(DB.results[,x], 1, 6), "..."),"</a>"))
                  })
                  DB.results$WebAddress <- NULL 
                  DB.results[, c("expectedMass", "candidateMass", "ppmMatch")] <- lapply(DB.results[, c("expectedMass", "candidateMass", "ppmMatch")], as.numeric)
                  return(DB.results)
                }}, escape=F, options = list(pageLength = 20))#,digits=4, sanitize.text.function = function(x) x)
              
              ##############################
              ##### 5. Best annotations #####
              ##############################
              
              output$BestCandidates <- DT::renderDataTable({
                if(length(BestAnno(object)) == 0){
                  DB.results<-data.frame("metID.dbProb function has not yet been run",stringsAsFactors=F)
                  colnames(DB.results)<-"Result"
                  return(DB.results) 
                } else {
                  DB.results<-tmp.BestAnno[[feat.indx]]
                  ###create clickable url links to DB in table
                  DB.results$DBid <-paste0("<a href='http://",DB.results$WebAddress, DB.results$DBid, 
                                           "' target='_blank'>", DB.results$DBid, "</a>")
                  # clickable SMILES pugRest
                  SMILESindx<-grep("SMILES$",colnames(DB.results))
                  DB.results[,SMILESindx]<-sapply(SMILESindx,function(x){
                    NoFuncIndx<-DB.results[, x] == ""
                    ifelse(NoFuncIndx,"",
                           paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                                  gsub("/|\\\\","", DB.results[,x]),"/PNG","' target='_blank'>", 
                                  paste0(substring(DB.results[,x], 1, 6), "..."),"</a>"))
                  })
                  DB.results$WebAddress <- NULL 
                  DB.results[, c("expectedMass", "candidateMass", "ppmMatch")] <- lapply(DB.results[, c("expectedMass", "candidateMass", "ppmMatch")], as.numeric)
                  return(DB.results)
                }}, escape=F, options = list(pageLength = 20))#,digits=4, sanitize.text.function = function(x) x)
             
              
              ############################################
              ##### 6. Best substructure annotations #####
              ############################################
              
              output$BestSubStrAnno <- DT::renderDataTable({
                if(nrow(subStrAnno.df) == 0){
                  SbStrResults <- data.frame("subStructure.prob function has not yet been run", stringsAsFactors=F)
                  colnames(SbStrResults)<-"Result"
                  return(SbStrResults) 
                } else {
             SbStrResults <-  subStrAnno.df[subStrAnno.df$compSpecName %in% Features.v[feat.indx], , drop = F]
             return(SbStrResults)
                }}, escape=F, options = list(pageLength = 8))
              
              
              ##ui for word cloud
              output$wordCloudSelect<-renderUI({
                if(length(DBanno(object)) == 0){
                  DBmatches <- "metID.dbAnnotate function has not yet been run/ No matches to DB"
                } else {
                  DBmatches <- tmp.DBanno.res[[feat.indx]]$DBname
                  shiny::selectizeInput('wordCloudSelect','DB match name : ', choices = DBmatches)
                }}) 
              
              #ui for max number of HMDB abstracts to return
              output$nPMIDAbstracts <- shiny::renderUI({
                shiny::numericInput('nPMIDAbstracts','Number of PubMed abstracts used to calculate word cloud (max = 10000) : ', value=100,min=10,max=10000,step=10)
              }) 
              #ui for max number of random articles to return
              output$nRandomArticles <- shiny::renderUI({
                shiny::numericInput('nRandomArticles','Number of randomly selected abstracts to display (max = 20) : ', value=5,min=1,max=20,step=1)
              }) 
              ###word cloud text output  
              output$WordCloudText <- shiny::renderText({"If a large number of abstracts (>100) from PubMed are returned the word cloud will take longer to update or may fail to load...Please Wait"})
              
              ####################################
              ##### 6. PubMed wordcloud plot #####
              ####################################
              
              output$WordCloud<-shiny::renderPlot({
                if(length(DBanno(object)) == 0){
                  ClAbs <- data.frame(word="metID.dbAnnotate function not run", freq=1)
                  PubMedWordcloud::plotWordCloud(ClAbs, min.freq=1, max.words=100, rot.per=0)
                } else if (nrow(tmp.DBanno.res[[feat.indx]])  == 0) {
                  ClAbs<-data.frame(word="No matches to DB", freq=1)
                  PubMedWordcloud::plotWordCloud(ClAbs, min.freq=1, max.words=100, rot.per=0)
                } else {
                  
                  if(!is.null(input$wordCloudSelect) & !is.null(input$nPMIDAbstracts))
                  {
                    #subset -1 to remove count value obtained from eutils XML file
                    PMIDs <- PMIDsearch(keys=input$wordCloudSelect, n=input$nPMIDAbstracts)[-1]      
                    if(length(PMIDs) > 0)
                    {
                      Abs <- PubMedWordcloud::getAbstracts(PMIDs)
                      if(length(Abs) > 0)
                      {
                        ClAbs <- PubMedWordcloud::cleanAbstracts(Abs)
                        ###remove compound names from word frequency table
                        SubsName<-strsplit(input$wordCloudSelect," ")[[1]]
                        SubsName<-unlist(lapply(SubsName,function(x) grep(paste0("\\b",x,"\\b"),ClAbs$word,ignore.case=T)))
                        if(length(SubsName)>0)
                        {
                          ClAbs<-ClAbs[-SubsName,]
                        }
                        ###only keep word which are less than 50 characters
                        ClAbs<-ClAbs[which(sapply(as.character(ClAbs$word),nchar)<50),]
                        suppressWarnings(PubMedWordcloud::plotWordCloud(ClAbs,max.words=100,scale=c(4,0.5)))
                      } else {
                        ###if no PMIDs returned then plot 
                        ClAbs<-data.frame(word="No Abstract text available",freq=1)
                        PubMedWordcloud::plotWordCloud(ClAbs,min.freq=1,max.words=100,rot.per=0)
                      }
                    } else {
                      ###if no PMIDs returned then plot 
                      ClAbs<-data.frame(word="No PubMedIDs returned",freq=1)
                      PubMedWordcloud::plotWordCloud(ClAbs,min.freq=1,max.words=100,rot.per=0)
                    }
                  }
                }
              })
              
              ###########################################
              ##### 7. PubMed random article table  #####
              ###########################################
              
              output$WordCloudTable <- shiny::renderTable({
                if(length(DBanno(object)) == 0){
                  WordCloud.df <- data.frame("No PMIDs returned")
                  colnames(WordCloud.df) <- input$FeatureNames
                  return(WordCloud.df)
                } else if (nrow(tmp.DBanno.res[[feat.indx]])  == 0) {
                  WordCloud.df <- data.frame("No PMIDs returned")
                  colnames(WordCloud.df) <- input$FeatureNames
                  return(WordCloud.df)
                } else if(!is.null(input$wordCloudSelect) & !is.null(input$nPMIDAbstracts)){
                  PMIDs <- PMIDsearch(keys = input$wordCloudSelect, n = input$nPMIDAbstracts)
                  Count <- PMIDs[1]
                  PMIDs <- PMIDs[-1]
                  RandArticleSample <- sample(PMIDs,
                                              ifelse(length(PMIDs)>=input$nRandomArticles,
                                                     input$nRandomArticles,length(PMIDs)))
                  ##obtain titles for randomly selected articles
                  RandArticleTitles <- getTitles(RandArticleSample)
                  WordCloud.df <- data.frame(c(Count,
                                             paste0("<a href='http://www.ncbi.nlm.nih.gov/pubmed/",
                                                    RandArticleSample,"' target='_blank'> PMID : ",
                                                    RandArticleSample," Title : ",
                                                    RandArticleTitles,"</a>")))
                  row.names(WordCloud.df)<-c("number PubMedIds returned ",
                                             paste0("random article ", 
                                                    seq(1,nrow(WordCloud.df)-1,1)))
                  colnames(WordCloud.df)<-input$FeatureNames
                  return(WordCloud.df)
                } else {
                  WordCloud.df<-data.frame("No PMIDs returned")
                  colnames(WordCloud.df) <- input$FeatureNames
                  return(WordCloud.df)
                }
              }, sanitize.text.function = function(x) x)
              
              #####################################
              ##### 8. output MetFrag results #####
              #####################################
              
              output$MetFragTable <- shiny::renderTable({
                if(length(tmp.metFrag) == 0)
                {
                  metFrag.df.tmp <- "metID.MetFrag function has not yet been run"
                  metFrag.df.tmp <- data.frame(metFrag.df.tmp)
                  colnames(metFrag.df.tmp)<-"Result"
                  return(metFrag.df.tmp)
                } else {
                  metFrag.df.tmp <- tmp.metFrag[[feat.indx]]
                  if(is.null(metFrag.df.tmp)){
                    metFrag.df.tmp <- "metID.MetFrag function has not yet been run"
                    metFrag.df.tmp <- data.frame(metFrag.df.tmp)
                    colnames(metFrag.df.tmp)<-"Result" 
                    return(metFrag.df.tmp)
                  } else if (is.character(metFrag.df.tmp)){
                    metFrag.df.tmp <- "MetFrag returned no results"
                    metFrag.df.tmp <- data.frame(metFrag.df.tmp)
                    colnames(metFrag.df.tmp)<-"Result" 
                    return(metFrag.df.tmp)
                  } else {  
                    # reduce n decimal places to 4
                    metFrag.df.tmp[, c("PeakScore", "BondEnergyScore", "Score")] <- sapply(metFrag.df.tmp[, c("PeakScore", "BondEnergyScore", "Score")], 
                                                                                           function(x) {
                                                                                             if(nchar(gsub("(.*\\.)", "", as.character(x[1]))) > 4){
                                                                                               x <- round(as.numeric(x), digits = 4) 
                                                                                             }
                                                                                           })
                    ###create clickable url links to DB in table
                    metFrag.df.tmp$DBid <-paste0("<a href='http://", metFrag.df.tmp$WebAddress, metFrag.df.tmp$DBid, 
                                                 "' target='_blank'>", metFrag.df.tmp$DBid, "</a>")
                    metFrag.df.tmp$WebAddress <- NULL 
                    metFrag.df.tmp$SMILES <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                                                    metFrag.df.tmp$SMILES, "/PNG", "' target='_blank'>",
                                                    paste0(substring(metFrag.df.tmp$SMILES, 1, 6), "..."), "</a>")
                    return(metFrag.df.tmp)
                  }
                }
              }, sanitize.text.function = function(x) x)
              
              ##load in any previous user comments
              ###add text area for storing user comments
              
              output$CommentSectionHeader <- shiny::renderUI({shiny::tags$b("Add metabolite identification notes including any useful html links here...\n
                                                                            NB. Do not forget to save your comments to the features results directory using the button on the left")})
              output$UserComments <- shiny::renderUI({
                if(length(Comments(compMS2)) == 0){
                  PrevComment <- "No previous comments"
                } else if(is.null(Comments(object)[[feat.indx]])){
                  PrevComment <- "No previous comments"
                } else {
                  PrevComment <- Comments(object)[[feat.indx]]
                  PrevComment <- as.character(PrevComment[nrow(PrevComment), 2])
                }
                
                shiny::tags$textarea(id="UserComments", rows=4, cols=60, PrevComment)
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
      }) # END compMS2server
      
      # Shiny UI    
      compMS2ui <- shiny::shinyUI(shiny::fluidPage(
        shiny::titlePanel("Composite MS2 Explorer "),
        shiny::fluidRow(shiny::column(shiny::h4("Options"),width=2,
                                      shiny::textInput("mass_to_charge",label=shiny::tags$b("mass to charge ratio : "),value=275.0223),
                                      shiny::numericInput("mass_accuracy",label=shiny::tags$b("mass accuracy (ppm) : "),value=10,min=0.1,max=10000,step=0.1),
                                      shiny::textInput("retentionTime",label=shiny::tags$b("retention time (seconds) : "),value=410.8),
                                      shiny::numericInput("RTtolerance",label=shiny::tags$b("retention time tolerance (+/- seconds) : "),value=5,min=1,max=700,step=1),
                                      shiny::selectInput("SubStrTypes",label=shiny::tags$b("Substructure types to include (multiple combinations) :"),choices=c(" ",SubStrType.inputs),multiple=T),
                                      shiny::selectInput("NotSubStrTypes",label=shiny::tags$b("Substructure types to exclude (multiple combinations) :"),choices=c(" ",SubStrType.inputs),multiple=T),
                                      #selectInput("Possible_contaminants",label=tags$b("Possible contaminant/false positive annotations to exclude (name_frequency) :"),choices=c(" ",as.character(AnnoFreq$name)),multiple=T),
                                      shiny::checkboxInput("All_Features","All features",value=TRUE),
                                      shiny::checkboxInput("DB_matches","Database matches",value=FALSE),
                                      shiny::actionButton("goButton","Submit")
        ),
        shiny::column(width=2,
                      shiny::textInput("DB_match_name","Search for database matches by name",value="e.g. p-cresol"),
                      shiny::actionButton("DBbutton","Submit"),
                      shiny::br(), 
                      shiny::br(),
                      shiny::tags$b("Save comments and best candidates"),
                      shiny::actionButton("CommentButton","SAVE"),
                      shiny::br(),
                      shiny::tags$b("Close app and return comments"),
                      shiny::br(),
                      shiny::actionButton("CloseAppBtn","Close"),
                      #                br(),br(),sliderInput("Loess_smooth","Select the Loess smoother span :",min=0.01,max=1,value=0.66),
                      #                numericInput("CorrCoefMin","Minimum correlation coefficient :",min=0.3,max=0.99,value=0.3,step=0.01),
                      shiny::br(),
                      shiny::br(),
                      shiny::uiOutput("FeatureNames"),
                      shiny::br(),
                      shiny::h5(shiny::textOutput("matchSummaryText")),
                      shiny::tableOutput("matchSummary"),
                      shiny::uiOutput("tabbedPanelButton"),
                      shiny::br(),
                      shiny::br(),
                      shiny::uiOutput("Plot_select")
        ),
        shiny::column(width=8,
                      shiny::conditionalPanel(condition="input.tabbedPanelButton >'0' & input.FeatureNames!='No MS2 features found'",
                                              shiny::tabsetPanel(
                                                # shiny::tabPanel("Composite MS2 plot",rCharts::chartOutput("Raw_data_plot","PolyCharts")),#,#plotOutput("MS2_plot",width = "800px", height = "600px")), 
                                                shiny::tabPanel("MS2 spectrum table",DT::dataTableOutput(outputId="MS2_data")),
                                                shiny::tabPanel("MS1 MS2 match summary",DT::dataTableOutput(outputId="metadata")),
                                                shiny::tabPanel("DB Annotations", DT::dataTableOutput(outputId="DB.results")),
                                                shiny::tabPanel("Best Annotations", DT::dataTableOutput("BestCandidates")),
                                                shiny::tabPanel("Substructure Annotations", DT::dataTableOutput("BestSubStrAnno")),
                                                  #                                   tabPanel("Best candidates",dataTableOutput("BestCandidate")),
                                                shiny::tabPanel("PubMed Word Cloud", shiny::uiOutput("wordCloudSelect"), shiny::uiOutput("nPMIDAbstracts"), shiny::verbatimTextOutput("WordCloudText"),
                                                                shiny::uiOutput("nRandomArticles"), shiny::tableOutput(outputId="WordCloudTable"), shiny::plotOutput("WordCloud",width = "800px", height = "600px")),
                                                shiny::tabPanel("MetFrag results", shiny::tableOutput(outputId="MetFragTable")),
                                                shiny::tabPanel("Substructure type table"),#,tableOutput(outputId="SubStr_type")),
                                                shiny::tabPanel("InterFeature Correlation"),#,chartOutput("InterFeatureCorr","morris")),#uiOutput("CorrCoefMin"),
                                                shiny::tabPanel("Chemical Similarity Scores"),#uiOutput("Compound"),uiOutput("IntraInter"),chartOutput("ChemicalSimilarity","morris")),
                                                 shiny::tabPanel("Comments"),#, shiny::uiOutput("CommentSectionHeader"), shiny::uiOutput("UserComments"))#,uiOutput("CommentButtonUI")
                                              shiny::tabPanel("LogD RT plot")#,uiOutput("RTwindow"),chartOutput("LogD_RT_plot","PolyCharts"))
                                              ))
        ))
        
        #         fluidRow(column(h4(paste0("Results summary ","(nMatches=",TotalFeatures,", nCompositeSpectra=",TotalCompSpectra,")")),width=6,
        #                         dataTableOutput("MS2_features_detected")))
      )) # end CompMS2 shiny UI
      #      
      object <-  shiny::runApp(list(ui = compMS2ui, server = compMS2server), launch.browser = browserLaunch)
      return(object)
    }
    object <- compMS2shinyW(object)
    return(object)
    
  }
  
})
