setClass("Domain",representation(xlim="numeric",ylim="numeric"),validity=function(object)
{
	if (length(object@xlim) != 2) return("length(xlim) != 2")
	if (length(object@ylim) != 2) return("length(xlim) != 2")

	if (any(is.na(object@xlim))) return("xlim NA")
	if (any(is.na(object@ylim))) return("ylim NA")

	if (any(is.infinite(object@xlim))) return("xlim Inf")
	if (any(is.infinite(object@ylim))) return("ylim Inf")

	if (any(object@xlim%%360 != object@xlim)) return("xlim out of [0,360[")
	if (any(abs(object@ylim) > 90)) return("ylim out of [-90,90]")
	if (diff(object@ylim) <= 0) return("ylim not strictly increasing")

	return(TRUE)
}
)

setMethod("initialize","Domain",def=function(.Object,long="numeric",lat="numeric",...)
{
	if (! missing(long)) .Object@xlim = rangeLong(long)
	if (! missing(lat)) .Object@ylim = range(lat)

	callNextMethod(.Object,...)
}
)

setGeneric("shift",valueClass="Domain",def=function(domain,x,y,...)
{
	standardGeneric("shift")
}
)

setMethod("shift","Domain",def=function(domain,x,y,...)
{
	new("Domain",xlim=domain@xlim+x,ylim=domain@ylim+y,...)
}
)

rangeLong = function(long)
{
	dlong1 = range(long%%360)
	dlong2 = range((long+180)%%360-180)

	if (diff(dlong1) > 180 && diff(dlong2) > 180)
		stop("longitudes span more than 180dg (case not supported)")

	if (diff(dlong1) <= diff(dlong2)) {
		return(dlong1)
	} else {
		return(dlong2)
	}
}
