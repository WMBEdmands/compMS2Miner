# ' deconvolute noise according to a maximum allowed intensity 
#   
#   setGeneric("deconvNoise.maxInt", function(object, ...) standardGeneric("deconvNoise.maxInt"))
#   
#   setMethod("deconvNoise.maxInt", signature = "CompMS2", function(object) {
#     # error handling
#     if(class(object) != "CompMS2"){
#       stop("argument object is not an CompMS2 class object")
#     } else {
#       
#     }
#   })