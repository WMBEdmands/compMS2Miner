#' CouchDB login
#' send CompMS2 data set to couchDB
#' @param couchDBname New or existing CouchDB database name (must be all lower case,  can contain underscores) 
#' @return CouchDB records : All Results from the current stage of the compMS2Miner
#'  are sent to the already established/ newly created couchDB database. 
#'  The following documents are sent to couchDB :
#'  
#' @export
setGeneric("couchDBcurate", function(object, ...) standardGeneric("couchDBcurate"))

setMethod("couchDBcurate", signature = "compMS2", function(object,
                                                           couchDBname=NULL, 
                                                           nSlaves=NULL, 
                                                           Username=NULL, 
                                                           Password=NULL, 
                                                           Host = NULL){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else  if(is.null(couchDBname)){
    stop("argument couchDBname is missing with no default")
  }
  
  # if couch DB credentials already available
  if(length(couchDBconn(compMS2)) == 0){

  # if username or password null get login details for couchDB
  if(is.null(Username) | is.null(Password) | is.null(Host)){    
    credentials <- getLoginDetails()
  } else {
    credentials <- c(Host=Host,Username=Username,Password=Password) 
  }
  
  message("Establishing connection with CouchDB...")
  flush.console()
  # login to couchDB and send ping
  pingReq <- couchDBpingReq(credentials = credentials)
  # loop through and break if connection established to give user 3 chances 
  # to get host,  username and password correct
  for(i in 1:2){
    #headers and cache control do not appear in 
    if(any((pingReq == "Error in response from CouchDB") == TRUE)){ 
      # if incorrect credentials then give user a message
      tcltk::tkmessageBox(message = 'Connection could not be made with CouchDB make sure the localhost server is running, and username/password are correct see futon interface "http://localhost:5984/_utils"')
      # run getLoginDetails if wrong details supplied
      credentials <- getLoginDetails()
      # login to couchDB and send ping
      pingReq <- couchDBpingReq(credentials = credentials)
    } else {
      message("Connection with CouchDB made...")
      flush.console()
      break
    }
  }
  
  # final error if not able to login
  if(any((pingReq=="Error in response from CouchDB")==TRUE))
  {
    tcltk::tkmessageBox(message = 'Too many failed attempts to connect with 
                        CouchDB try again')
    stop('Too many failed attempts to connect with CouchDB try again')
  }
  
  # check to see if database name supplied is found if not then create new 
  # database and tell user
  if(is.null(couchDBname))
  {
    stop("a new or existing CouchDB database name must be supplied in order to 
         store the results output")
  } else {
    # create base_url
  auth <- ifelse(credentials["Username"]=="", "", 
                 paste0(credentials["Username"], ":", credentials["Password"],
                        "@"))
  base_url  <-  paste0(proto="http",  "://",  auth,  myConn$couch_http_host, 
                         ":", myConn$couch_http_port)
  # data bases
  path  <-  paste(base_url,  "_all_dbs",  sep = "/")
  DBs <- rjson::fromJSON(file=path)
  # convert to lower case
  couchDBname <- tolower(couchDBname)
  # if database name is not in the database list then create new 
  if(!any(DBs == couchDBname)){
  DBcreate <- couchDB::couch_create_database(myConn, couchDBname)
  }
  }
  # add couchDB connection to compMS2 object
  couchDBconn(object) <- list(base_url = base_url, myConn = myConn)
  Parameters(object) <- data.frame(Username = credentials["Username"], 
                                   couchDBname = couchDBname, Parameters(object)) 
  }
# send CompMS2 class object to couchDB 


return(object)

})
