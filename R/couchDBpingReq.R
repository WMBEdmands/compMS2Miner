#' CouchDBpingReq
#' send ping to couch db using log in credentials from GetLoginDetails or
#' named character vector
#' 
#' @param credentials either output from GetLoginDetails or named character vector "host", "Username", "Password" containing user Login parameters (e.g. c( host = "localhost", Username = "", Password = ""))
#' @return ping request result 
#' @export
couchDBpingReq<-function(credentials){
  # couchDB http connection 
  myConn <<- couchDB::couch_http_connection(host=credentials['Host'], 
                                            user=credentials['Username'], 
                                            password=credentials['Password'])
  # authorization
  auth <<- paste0(credentials['Username'], ":", credentials['Password'], "@")
  # send ping
  pingReq <- couchDB::couch_ping(myConn)
  
  return(pingReq)
  
}