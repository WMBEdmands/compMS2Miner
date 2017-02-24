#' Get Login details for CouchDB
#' 
#' tcltk GUI to get login details for couchDB
#' 
#' @param Host couchDB host name, defaults to localhost. Can be online repository.
#' @param Name couchDB administrator username (used for http commands, not stored).
#' @param Password couchDB administrator password (used for http commands, not stored).
#' @return login details for couchDB
#' @export
getLoginDetails  <-  function(){
  # Based on code by Barry Rowlingson
  # http://r.789695.n4.nabble.com/tkentry-that-exits-after-RETURN-tt854721.html
  #none
  tt  <-  tcltk::tktoplevel()
  tcltk::tkwm.title(tt,  "Get login details")
  Host <- tcltk::tclVar("localhost")
  Name  <-  tcltk::tclVar("")
  Password  <-  tcltk::tclVar("")
  entry.Host  <-  tcltk::tkentry(tt, width="20",  textvariable=Host)
  entry.Name  <-  tcltk::tkentry(tt, width="20",  textvariable=Name)
  entry.Password  <-  tcltk::tkentry(tt,  width="20",  show="*",  
                              textvariable=Password)
  tcltk::tkgrid(tcltk::tklabel(tt,  text="Please enter your login details."))
  tcltk::tkgrid(tcltk::tklabel(tt,  text="couchDB host: "), entry.Host)
  tcltk::tkgrid(tcltk::tklabel(tt,  text="username: "), entry.Name)
  tcltk::tkgrid(tcltk::tklabel(tt,  text="password: "), entry.Password)
  
  OnOK  <-  function()
  { 
    tcltk::tkdestroy(tt) 
  }
  OK.but  <- tcltk::tkbutton(tt, text=" OK ",  command=OnOK)
  tcltk::tkbind(entry.Password,  "<Return>",  OnOK)
  tcltk::tkgrid(OK.but)
  tcltk::tkfocus(tt)
  tcltk::tkwait.window(tt)
  
  invisible(c(Host=tcltk::tclvalue(Host), Username=tcltk::tclvalue(Name),  
              Password=tcltk::tclvalue(Password)))
}
