setClass("GaussGrid",representation(nlong="integer",pole="numeric"),
	validity=function(object)
{
	nlat = length(object@nlong)
	if (nlat > 1) return("length(nlong) <= 1")
	if (nlat %% 2 != 0) return("length(nlong) is odd (uneven)")

	if (any(is.na(object@nlong))) return("nlong NA")
	if (any(is.infinite(object@nlong))) return("nlong Inf")
	if (any(object@nlong <= 0)) return("nlong <= 0")
	#if (any(object@nlong %% 2 != 0))) return("nlong has odd (uneven) values")

	if (any(object@nlong[seq(nlat)] != rev(object@nlong[-seq(nlat)]))) {
		return("nlong not symmetric")
	}

	if (diff(object@nlong[seq(nlat)]) >= 0) return("nlong not increasing")

	if (length(object@pole) != 2) return("length(pole) != 2")

	if (any(is.na(object@pole))) return("pole NA")
	if (any(is.infinite(object@pole))) return("pole Inf")

	if (abs(object@pole[1]) > 90) return("pole[1] (ie lat) out of [-90,90]")

	return(TRUE)
}
)
