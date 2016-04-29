#' Monoisotopic mass matching unknown to data base entry
monoMassMatch <- function(unknowns = NULL, metMasses.df = NULL, 
                          metMasses = NULL, DB_unique_IDs = NULL,
                          DB_entry_names = NULL, DB_SMILES = NULL,
                          dbEntryWebAddress = NULL, 
                          subStrMasses = NULL, subStrNames = NULL, 
                          nSlaves = NULL, mode="pos", ppm = 10){
  # error handling
  if(is.null(unknowns)){
    stop("unknowns argument is missing with no default")    
  } else if (is.null(metMasses.df)){
    if (is.null(metMasses)){
      stop("metMasses argument is missing with no default")    
    }
    if (is.null(DB_unique_IDs)){
      stop("DB_unique_IDs argument is missing with no default")  
    }
    if (is.null(DB_entry_names)){
      stop("DB_entry_names argument is missing with no default")    
    }  
    
  } else {
    if(!all(c("Unique_DB_ID", "monoisotopic_weight", "name", "SMILES") %in% colnames(metMasses.df)))
    {
      stop('metMasses.df colnames for unique database IDs, monoisotopic weights, names and canonical SMILES codes of data base
           entries must be named "Unique_DB_ID", "monoisotopic_weight", "name" and "SMILES" respectively' )
    }
    metMasses <- metMasses.df$monoisotopic_weight
    DB_unique_IDs <- metMasses.df$Unique_DB_ID
    DB_entry_names <- metMasses.df$name
    DB_SMILES <- metMasses.df$SMILES
    
    }  
  if (!grepl("pos|neg",mode)){
    stop('mode must be one of either "pos" or "neg"') 
  } else if (!is.null(subStrMasses) & is.null(subStrNames)){
    stop("Argument subStrMasses is supplied but corresponding substructure
         names have not (i.e. argument subStrNames)")
  } else if (!is.null(subStrNames) & is.null(subStrMasses)){
    stop("Argument subStrNames is supplied but corresponding substructure
         mass shifts have not (i.e. argument subStrMasses)")
  } else if (is.null(dbEntryWebAddress)) {
    stop("dbEntryWebAddress is missing with no default (e.g. www.hmdb.ca/metabolites/),
         this address is used for generating links to the database for matches")
  } else {  
    
    ###calculate parent ions from metMasses supplied 
    # electrospray adducts taken from 
    # http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/
    if(mode == "pos"){
      
      metMasses <- data.frame(matrix(c(metMasses + 1.007276,  
                                       metMasses / 3 + 1.007276,  
                                       metMasses / 3 + 8.334590,  
                                       metMasses / 3 + 15.7661904,  
                                       metMasses / 3 + 22.989218,  
                                       metMasses / 2 + 1.007276,  
                                       metMasses / 2 + 9.520550,	
                                       metMasses / 2 + 11.998247,	
                                       metMasses / 2 + 19.985217,	
                                       metMasses / 2 + 22.989218,	
                                       metMasses + 18.033823,
                                       metMasses + 22.989218,	
                                       metMasses + 33.033489,
                                       metMasses + 38.963158,	
                                       metMasses + 44.971160,
                                       metMasses + 76.919040,	
                                       2 * metMasses + 1.007276,	
                                       2 * metMasses + 83.060373,	
                                       2 * metMasses + 22.989218,	
                                       2 * metMasses + 28.02312,	
                                       2 * metMasses + 38.963158), 
                                     nrow = length(metMasses)))
      colnames(metMasses)<- c("M+H", "M+3H", "M+2H+Na", "M+H+2Na", "M+3Na", "M+2H",
                              "M+H+NH4;", "M+H+Na", "M+H+K", "M+2Na", "M+NH4;", 
                              "M+Na", "M+CH3OH+H", "M+K", "M+2Na-H", "M+2K-H", 
                              "M2+H", "M2+NH4;", "M2+Na", "M2+3H2O+2H", "M2+K")
    } else if(mode == "neg"){
      
      metMasses <- data.frame(matrix(c(metMasses - 1.007276,
                                       metMasses - 19.01839,
                                       metMasses + 44.998201,  
                                       metMasses + 20.974666,	
                                       metMasses + 36.948606,	
                                       metMasses + 59.013851,
                                       metMasses + 82.00307,
                                       2 * metMasses + 44.998201,	
                                       metMasses /2 - 1.007276,	
                                       metMasses /3 - 1.007276,	
                                       2 * metMasses - 1.007276,	
                                       3 * metMasses - 1.007276,	
                                       2 * metMasses + 59.013851), nrow = length(metMasses)))
      colnames(metMasses)<- c("M-H", "M-H2O-H", "M+FA-H", "M+Na-2H", "M+K-2H", 
                              "M+Hac-H", "M+Hac+Na-H", "M2+FA-H", "M-2H", "M-3H", "M2-H", "M3-H",
                              "M2+Hac-H")
    }
    # add substr mass shifts
    if(!is.null(subStrMasses)){
      colnames.tmp <- paste0(colnames(metMasses), ";", 
                             rep(subStrNames, each = ncol(metMasses)))
      subStr.df.tmp <- do.call(cbind, lapply(subStrMasses, function(ssM) 
        apply(metMasses, 2, function(x) x + ssM)))
      colnames(subStr.df.tmp) <- colnames.tmp
      metMasses <- cbind(metMasses, subStr.df.tmp)
      
    }
    metMasses <- unlist(metMasses)
    # max min range to limit search space
    uk.range <- range(unknowns) 
    metMasses <- metMasses[metMasses > (uk.range[1]-1) & metMasses < (uk.range[2]+1)]
    
    message("matching ", length(unknowns), " unknowns to ", 
            format(length(metMasses), big.mark=",", scientific=FALSE),
            " possible ESI artefacts and substructure mass shifts...")
    flush.console()
    
    names(metMasses) <- paste(metMasses, names(metMasses), sep = "/")
    # split masses into seperate chunks of mass range to speed up computation 
    #        metMasses.indx <- cut(metMasses, 1000)
    #        levels(metMasses.indx) <- unique(ave(metMasses, metMasses.indx))
    #        meanIntervals <- tapply(metMasses, metMasses, mean)
    #        metMasses <- split(metMasses, metMasses.indx)
    #        metMasses.indx <- unlist(sapply(1:length(metMasses), function(x) 
    #                                 data.frame(range(metMasses[x]))))
    #        names(metMasses.indx) <- rep(1:1000, each =2)
    #        ###create function to match with mapply
    #        exactMassMatch <- function(unknown.mz){
    #          unk.mass.diff <- abs(as.numeric(levels(metMasses.indx)) - unknown.mz)
    #          closest.indx <- sort(unk.mass.diff, index.return = TRUE)$ix[1:3]
    #          metMasses.sub.indx <- which(metMasses.indx %in% levels(metMasses.indx)[closest.indx])
    #          ##if a match is made then add the information into the MSfeatures dataframe
    #          metMasses.sub.indx <- metMasses.sub.indx[which((
    #                       ##match MS1 mass by ppm tolerance to precursor m/z
    #                       metMasses[metMasses.sub.indx] < unknown.mz + (unknown.mz/1E06) * ppm &
    #                       metMasses[metMasses.sub.indx] > unknown.mz - (unknown.mz/1E06) * ppm) == T )]
    #          MS1.match.indx <- metMasses[metMasses.sub.indx]
    #          ##if a match is made then add the information into the MSfeatures dataframe
    #    MS1.match.indx <- metMasses[which((
    #    ##match MS1 mass by ppm tolerance to precursor m/z
    #    metMasses < unknown.mz + (unknown.mz/1E06) * ppm &
    #    metMasses > unknown.mz - (unknown.mz/1E06) * ppm) == T )]
    #                   
    #          if(length(MS1.match.indx) > 0){
    #            ppm.tmp <- ((MS1.match.indx - unknown.mz) / MS1.match.indx) * 1E06
    #            MS1.match.indx <- paste(MS1.match.indx, ppm.tmp, sep = ";")
    #            names(MS1.match.indx) <- gsub(".+\\.", "", names(ppm.tmp))
    #          }
    #            
    #          return(MS1.match.indx)   
    #        }
    # compile function using compiler package
    #    exactMassMatch.c <- compiler::cmpfun(exactMassMatch)
    # matching unknowns to possible ESI artefacts and substructure masses shifts
    
    #    if(!is.null(nSlaves)){
    #      if(require(Rcpp) == T){
    #        ##############################################################################
    #        ##############################################################################
    #        ##############################################################################
    #        ##############################################################################
    #        
    #        exactMassMatchC_parallel <- function(chunk, unkns = unknowns, metMss = metMasses, 
    #                                             metMssNames = names(metMasses), ppmTol = ppm){
    #          Rcpp::cppFunction('List exactMassMatchC(List Res, NumericVector unknowns, 
    #                            NumericVector metMasses, CharacterVector metMassesNames, 
    #                            double ppm, double mill) {
    #                            int n = unknowns.size();
    #                            for(int i = 0; i < n; ++i) {
    #                            Res[i] = metMassesNames[metMasses < unknowns[i] + (unknowns[i] / mill) * ppm & metMasses > unknowns[i] - (unknowns[i] / mill) * ppm];
    #                            
    #                            }
    #                            return Res;
    #                            
    #                            }')
    #   #chunkIndx <- as.numeric(unlist(chunkIndx))
    #   matches.tmp <- vector("list", length(chunk))
    #   matches.tmp <- exactMassMatchC(matches.tmp, unkns[chunk], metMss, metMssNames, ppmTol, 1E06)
    #   return(matches.tmp)
    #        }
    #   ##############################################################################
    #   ##############################################################################
    #   ##############################################################################
    #   ##############################################################################
    #   
    #   message(paste0("Starting SNOW cluster with ", nSlaves,
    #                  " local sockets..."))
    #   flush.console()
    #   pmt <- proc.time()
    #   cl <- parallel::makeCluster(nSlaves) 
    #   doSNOW::registerDoSNOW(cl)
    #   
    #   chunkIndx <- split(seq_along(unknowns), cut(seq_along(unknowns), nSlaves))
    #   # foreach and dopar from foreach package
    #   matches.tmp <- foreach(j = 1:length(chunkIndx), .packages = "Rcpp") %dopar% {
    #     exactMassMatchC_parallel(chunk = chunkIndx[[j]])}
    #   # stop SNOW cluster
    #   parallel::stopCluster(cl) 
    #   proc.time() - pmt
    #   matches.tmp <- unlist(matches.tmp, recursive = F)
    #   
    #      } else {
    #      # create a cluster using the doSNOW package
    #      message(paste0("Starting SNOW cluster with ", nSlaves,
    #                     " local sockets..."))
    #      flush.console()
    #      
    #      cl <- parallel::makeCluster(nSlaves) 
    #      doSNOW::registerDoSNOW(cl)
    #      
    #      # foreach and dopar from foreach package
    #      matches.tmp <- foreach(j = 1:length(unknowns)) %dopar% {
    #                              exactMassMatchC(unknowns[j])}
    #      # stop SNOW cluster
    #      parallel::stopCluster(cl) 
    #      } 
    #    } else {
    #      if(require(Rcpp) == T){
    ##############################################################################
    ##############################################################################
    ##############################################################################
    ##############################################################################
    Rcpp::cppFunction('List exactMassMatchC(List Res, NumericVector unknowns, 
                      NumericVector metMasses, CharacterVector metMassesNames, 
                      double ppm, double mill) {
                      int n = unknowns.size();
                      for(int i = 0; i < n; ++i) {
                      Res[i] = metMassesNames[metMasses < unknowns[i] + (unknowns[i] / mill) * ppm & metMasses > unknowns[i] - (unknowns[i] / mill) * ppm];
                      
                      }
                      return Res;
                      
                      }')
       
    ##############################################################################
    ##############################################################################
    ##############################################################################
    ##############################################################################
    
    matches.tmp <- vector("list", length(unknowns))
    pmt <- proc.time()
    matches.tmp <- exactMassMatchC(matches.tmp, unknowns, metMasses, names(metMasses), ppm, 1E06)
    proc.time() - pmt
    #      }
    #    }
    #      } else {
    #      # create list to store results
    #      pmt <- proc.time()
    #      matches.tmp <- vector("list", length(unknowns))
    #      #progress bar
    #      pb <- txtProgressBar(min=0, max=length(matches.tmp), style=3)
    #      
    #      for(i in 1:length(matches.tmp)){
    #        #progress bar
    #        Sys.sleep(0.01)
    #        setTxtProgressBar(pb, i)
    #        flush.console()
    #        
    # #        match.tmp <- exactMassMatch(unknowns[i])
    # match.tmp <- exactMassMatchC(unknowns[i], metMasses, names(metMasses), ppm, 1E06)
    # 
    #        matches.tmp[[i]] <- match.tmp
    #      }
    #      proc.time() -  pmt
    #      }
    #      test <- rep(list(list(ppm = 0, match.indx = 0)), 11)
    #create C++ function to match mS1 features to database more rapidly 
    # 
    # 
    # Rcpp::cppFunction('CharacterVector exactMassMatchC(NumericVec unknownMz, 
    # NumericVector metMasses, CharacterVector metMassesNames, double mill) {
    #                       return metMassesNames[metMasses < (unknownMz + mill) & metMasses > (unknownMz - mill)];
    #                      }')
    # 
    # 
    #    Rcpp::cppFunction('CharacterVector exactMassMatchC(double unknownMz, 
    # NumericVector metMasses, CharacterVector metMassesNames, 
    # double ppm, double mill) {
    #                       return metMassesNames[metMasses < unknownMz + (unknownMz / mill) * ppm & metMasses > unknownMz - (unknownMz/mill) * ppm];
    #                      }')
    # 
    # op <- microbenchmark(test.cppfunc <- lapply(unknowns, function(unknownMz) exactMassMatchC(unknownMz, metMasses, names(metMasses), 10, 1E06)), 
    #                      test.Rfunc <- lapply(unknowns, function(unknown.mz) exactMassMatch(unknown.mz)))
    
    # double mpe(List mod) {
    #   if (!mod.inherits("lm")) stop("Input must be a linear model");
    #   
    #   NumericVector resid = as<NumericVector>(mod["residuals"]);
    #   NumericVector fitted = as<NumericVector>(mod["fitted.values"]);
    #   
    #   int n = resid.size();
    #   double err = 0;
    #   for(int i = 0; i < n; ++i) {
    #     err += resid[i] / (fitted[i] + resid[i]);
    #   }
    #   return err / n;
    # }
    # 
    #    cppFunction('double sumC(NumericVector x) {
    #   int n = x.size();
    #                double total = 0;
    #                for(int i = 0; i < n; ++i) {
    #                total += x[i];
    #                }
    #                return total;
    #    }')
    # #        ###use mapply to match the unknowns
    # #       system.time(  <- mapply(exactMassMatch, unknowns))
    # #             
    
    
    matches.tmp <- lapply(1:length(matches.tmp), function(x) {
      monoMass <- gsub("/.+", "", matches.tmp[[x]])
      ESI_type <- gsub(".+/|[0-9]*$", "", matches.tmp[[x]])
      SubStr_type <- gsub(".+;", "", ESI_type)
      ESI_type <- gsub(";.+", "", ESI_type)
      SubStr_type <- ifelse(SubStr_type == ESI_type, "", SubStr_type)
      indx.tmp <- gsub(".+([[:alpha:]])", "", matches.tmp[[x]])
      indx.tmp <- as.numeric(gsub('.+;', '', indx.tmp))
      tmp.df <- data.frame(WebAddress = dbEntryWebAddress[indx.tmp],
                           DBid = DB_unique_IDs[indx.tmp], 
                           DBname = DB_entry_names[indx.tmp],
                           SMILES = DB_SMILES[indx.tmp],
                           ESI_type = gsub(';', '', ESI_type), 
                           SubStr_type = SubStr_type,
                           expectedMass =  rep(unknowns[x], length(monoMass)),
                           candidateMass = as.numeric(monoMass), stringsAsFactors = F)
      tmp.df$ppmMatch <- ((tmp.df$expectedMass - tmp.df$candidateMass) / tmp.df$expectedMass) *1E06
      return(tmp.df)})
    
    return(matches.tmp)
  }
}