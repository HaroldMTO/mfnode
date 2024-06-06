Gnum = "-?(\\d+\\.\\d*|\\d*\\.\\d+([eE]?[-+]?\\d+)?\\>)"
#Gnum = "-?\\d*\\.\\d+([eE]?[-+]?\\d+)?\\>"
Gint = "-?\\d+\\>"
gpfre1 = "[UVW] VELOCITY|(SURFACE )?PRESSURE|TEMPERATURE|GRAD[LM]_\\w+|GEOPOTENTIAL"
gpfre2 = "MOIST AIR SPECIF|ISOBARE CAPACITY|SURFACE DIV|d\\(DIV\\)\\*dP"
gpfre3 = "(ATND|ADIAB|CTY|GNH|SISL|SI|SL)_\\w+"
gpfre = paste(gpfre1,gpfre2,gpfre3,sep="|")

getvar = function(var,nd,sep="=")
{
	re = sprintf("^ *\\<%s *%s *(%s|%s).*",var,sep,Gint,Gnum)
	ndvar = grep(re,nd,value=TRUE)
	if (length(ndvar) == 0) {
		re = sprintf(" *\\<%s *%s *(%s|%s).*",var,sep,Gint,Gnum)
		ndvar = grep(re,nd,value=TRUE)
		if (length(ndvar) == 0) return(NULL)

		ndvar = sub(sprintf("^.+(\\<%s\\> *%s)",var,sep),"\\1",ndvar)
	}

	unique(as.numeric(gsub(re,"\\1",ndvar)))
}

line2num = function(nd)
{
	lre = regmatches(nd,gregexpr(sprintf("(%s|\\<NaN\\>)",Gnum),nd))
	lre = lapply(lre,function(x) gsub("(\\d+)([-+]\\d+)","\\1E\\2",x))
	sapply(lre,as.numeric)
}

intlines = function(nd)
{
	as.integer(unlist(regmatches(nd,gregexpr(Gint,nd))))
}

numlines = function(nd)
{
	as.numeric(unlist(regmatches(nd,gregexpr(Gnum,nd))))
}

sulines = function(nd,strict=FALSE)
{
	iend = grep("\\<END OF SETUPS\\>",nd)
	if (length(iend) == 0) {
		if (strict) stop("no line matching 'END OF SETUPS'")
		return(nd)
	}

	nd[seq(iend)]
}

initlines = function(nd,TLAD.ok=TRUE,strict=FALSE)
{
	istart = grep("^ *START +CNT3\\>",nd)
	iend = grep("^ *START +CNT4\\>",nd)

	if (strict) {
		if (length(istart) == 0) stop("no line matching 'START CNT3'")
		if (length(iend) == 0) stop("no line matching 'START CNT4'")
	}

	if (length(istart) == 0) istart = 1
	if (length(iend) == 0) iend = length(nd)
	nd[istart[1]:iend[length(iend)]]
}

fclines = function(nd,TLAD.ok=TRUE,strict=FALSE)
{
	if (TLAD.ok) {
		icnt4 = grep("^ *START +CNT4(TL|AD)?",nd)
		icnt3 = grep("END +CNT3(TL|AD)?",nd)
	} else {
		icnt4 = grep("^ *START +CNT4\\>",nd)
		icnt3 = grep("END +CNT3\\>",nd)
	}

	if (strict) {
		if (length(icnt4) == 0) stop("no line matching 'START CNT4'")
		if (length(icnt3) == 0) stop("no line matching 'END CNT3'")
	}

	if (length(icnt4) == 0) icnt4 = 1
	if (length(icnt3) == 0) icnt3 = length(nd)
	ndfc = nd[icnt4[1]:icnt3[length(icnt3)]]

	# filter out empty lines (important for conceptual search of norm heading lines)
	grep("^ *$",ndfc,invert=TRUE,value=TRUE)
}

checklev = function(lev,has.levels,nflevg)
{
	if (has.levels) {
		stopifnot(length(lev) == 1 && lev %in% 0:nflevg || all(lev %in% 1:nflevg))
	} else {
		stopifnot(identical(lev,0L))
	}
}

normtags = function(nd,type)
{
	tags = character()

	if (type == "SP") {
		ind = grep("NORMS AT (START|NSTEP|END) CNT4",nd,ignore.case=TRUE)
		if (length(ind) > 0) tags = "NORMS AT (START|NSTEP|END) CNT4"
		ndsp = grep("^ *spnorm +",nd,ignore.case=TRUE,value=TRUE)
		if (length(ndsp) > 0) tags = c(tags,unique(sort(ndsp)))
	} else if (type == "GP") {
		ndgp = grep("^ *gpnorm +",nd,ignore.case=TRUE,value=TRUE)
		ndgp = grep("AVERAGE +MIN(IMUM)? +MAX",ndgp,invert=TRUE,value=TRUE)
		if (length(ndgp) > 0) tags = c(tags,unique(sort(ndgp)))
	}

	tags
}

spnorm = function(nd,lev,tag="NORMS AT (START|NSTEP|END) CNT4",abbrev=TRUE)
{
	lev = as.integer(lev)

	ndsu = sulines(nd)
	has.levels = getvar("NSPPR",ndsu) > 0
	nflevg = getvar("NFLEVG",ndsu)
	checklev(lev,has.levels,nflevg)

	# exclude lines before CNT4 (avoid norms for DFI)
	nd = fclines(nd)

	ind = grep("SPECTRAL NORMS",nd)

	# normal (ie predictor) is 1 line before, corrector (if present) is 2 lines before
	indg = grep(tag,nd[ind-1],ignore.case=TRUE)
	lpc = any(regexpr("\\<LPC_FULL *= *T",ndsu) > 0)
	if (tag == "NORMS AT (START|NSTEP|END) CNT4" && lpc) {
		indg2 = grep(tag,nd[ind-2],ignore.case=TRUE)
		indg = unique(sort(c(indg,indg2)))
	}

	if (length(indg) == 0) {
		cat("--> no group found for tag",tag,"\n")
		return(NULL)
	}

	ind = ind[indg]

	spsp = gsub("SPECTRAL NORMS - +LOG\\(PREHYDS\\) +([-0-9.E+]+|NaN)( +OROGRAPHY.+)?$",
		"\\1",nd[ind])
	spsp = as.numeric(gsub("(\\d+)([-+]\\d+)","\\1E\\2",spsp))

	noms = strsplit(nd[ind[1]+1]," {2,}")[[1]][-1]

	indi = rep(ind,each=length(lev))+lev+2

	spn = line2num(nd[indi])

	nval = length(lev)*length(noms)
	nt = length(spn)%/%nval
	lnhdyn = any(regexpr("\\<LNHDYN *= *T",ndsu) > 0)
	cat("--> found",nt,"time-steps of SP norms - NHdyn/PC scheme:",lnhdyn,"/",lpc,"\n")
	if (length(spn) > nt*nval) {
		cat("--> truncated last time-step, delete it\n")
		length(spn) = nt*nval
	}

	dim(spn) = c(length(noms),length(lev),nt)
	spn = aperm(spn,c(3,2,1))

	if (identical(lev,0L)) {
		spn = c(spn,spsp)
		noms = c(noms,"SP")
		dim(spn) = c(nt,length(lev),length(noms))
	}

	if (lnhdyn) {
		nl2 = 2+has.levels*nflevg
		indnh = ind[1]+nl2+1
		nomsnh = strsplit(nd[indnh]," {2,}")[[1]][-1]

		indi = indi+nl2
		spnh = line2num(nd[indi])
		nval = length(lev)*length(nomsnh)
		if (length(spnh) > nt*nval) length(spnh) = nt*nval
		dim(spnh) = c(length(nomsnh),length(lev),nt)
		spnh = aperm(spnh,c(3,2,1))

		spn = c(spn,spnh)
		noms = c(noms,nomsnh)
		dim(spn) = c(nt,length(lev),length(noms))
	}

	prestepo = tag == "NORMS AT (START|NSTEP|END) CNT4"
	spstep = stepindex(nd,ind,prestepo=TRUE)
	if (! prestepo) {
		spstepf = stepindex(nd,ind)
		ist = (is.na(spstep) | duplicated(spstep) | regexpr("^X",spstepf) > 0) &
			! is.na(spstepf)
		if (any(ist)) sptep[ist] = spstepf[ist]
	}

	if (abbrev) {
		noms = shortnames(noms)
		ip = grep("LOG\\(P/P_hyd\\)|d4|d5",noms,invert=TRUE)
		noms[ip] = abbreviate(noms[ip])
	}

	dimnames(spn) = list(spstep,lev,noms)

	# account for STEPX
	if (dim(spn)[1] > 1 && identical(all.equal(spn[1,,],spn[2,,]),TRUE)) {
		cat("--> duplicated 1st step (ie STEPX), remove it\n")
		spn = spn[-1,,,drop=FALSE]
	}

	spn
}

gpnorm = function(nd,lev,tag="NORMS AT (START|NSTEP|END) CNT4",gpin="\\w+.*",
	gpout=character(),abbrev=TRUE)
{
	lev = as.integer(lev)

	# limit lines (for performance)
	ndsu = sulines(nd)
	has.levels = getvar("NSPPR",ndsu) > 0
	nflevg = getvar("NFLEVG",ndsu)
	checklev(lev,has.levels,nflevg)

	# exclude lines before CNT4 (avoid norms for DFI)
	nd = fclines(nd)

	ind = grep(sprintf("GPNORM +(%s) +AVERAGE",gpin),nd)
	gpout = paste(c(gpout,"OUTPUT"),collapse="|")
	indo = grep(sprintf("GPNORM +(%s) +AVERAGE",gpout),nd[ind],invert=TRUE)
	ind = ind[indo]
	if (length(ind) == 0) return(NULL)

	nl2 = 2+has.levels*nflevg

	# normal (ie predictor) is 1 line before, corrector (if present) is 2 lines before
	lpc = any(regexpr("\\<LPC_FULL *= *T",ndsu) > 0)
	lnhdyn = any(regexpr("\\<LNHDYN *= *T",ndsu) > 0)
	if (tag == "NORMS AT (START|NSTEP|END) CNT4") {
		up = nl2
		if (lnhdyn) up = 2*nl2
		indg = grep(tag,nd[ind-up-2],ignore.case=TRUE)
		if (lpc) {
			indg2 = grep(tag,nd[ind-up-3],ignore.case=TRUE)
			indg = sort(unique(indg,indg2))
		}

		if (length(indg) == 0) {
			cat("--> no group found, try 1 more line before\n")
			indi = ind-up-4
			indg = grep(tag,nd[indi[indi > 0]],ignore.case=TRUE)
		}
	} else {
		indg = grep(tag,nd[ind-1],ignore.case=TRUE)
	}

	if (length(indg) == 0) {
		cat("--> no group found for tag",tag,"\n")
		return(NULL)
	}

	indh = indexpand(nd,ind,indg,nl2,average=identical(lev,0L))
	if (length(indh) == 0) return(NULL)

	# in TL/AD, norms for GFL are only printed for some fields (q and some others)
	if (any(! indh %in% ind)) {
		cat("--> non-constant group, reduce fields to constant ones\n")
		indh = indh[indh %in% ind]
		noms = sub(" *GPNORM +(\\w+.+?) +AVERAGE.+","\\1",nd[indh])
		frf = table(noms)
		ix = which(frf == max(frf))
		i = grep(paste(names(frf)[ix],collapse="|"),nd[ind])
		ind = ind[i]
		if (tag == "NORMS AT (START|NSTEP|END) CNT4") {
			indg = grep(tag,nd[ind-up-2],ignore.case=TRUE)
		} else {
			indg = grep(tag,nd[ind-1],ignore.case=TRUE)
		}
		indh = indexpand(nd,ind,indg,nl2,average=identical(lev,0L))
	}

	noms = unique(sub(" *GPNORM +(\\w+.+?) +AVERAGE.+","\\1",nd[indh]))

	indi = rep(indh,each=length(lev))+lev+1

	gpn = line2num(nd[indi])

	nval = 3*length(lev)*length(noms)
	nt = length(gpn)%/%nval
	cat("--> found",nt,"groups of GP norms - NH dyn/PC scheme:",lnhdyn,"/",lpc,"\n")
	if (length(gpn) > nt*nval) {
		cat("--> truncated last time-step, delete it\n")
		length(gpn) = nt*nval
	}

	dim(gpn) = c(3,length(lev),length(noms),nt)
	gpn = aperm(gpn,c(4,2,1,3))

	prestepo = tag == "NORMS AT (START|NSTEP|END) CNT4"
	gpstep = stepindex(nd,ind[indg],prestepo=TRUE)
	if (! prestepo) {
		gpstepf = stepindex(nd,ind[indg])
		ist = (is.na(gpstep) | duplicated(gpstep) | regexpr("^X",gpstepf) > 0) &
			! is.na(gpstepf)
		if (any(ist)) gpstep[ist] = gpstepf[ist]
	}

	if (abbrev) noms = shortnames(noms)
	dimnames(gpn) = list(gpstep,lev,c("ave","min","max"),noms)

	# account for STEPX
	if (dim(gpn)[1] > 1 && identical(all.equal(gpn[1,,,],gpn[2,,,]),TRUE)) {
		cat("--> duplicated 1st step (ie STEPX), remove it\n")
		gpn = gpn[-1,,,,drop=FALSE]
	}

	gpn
}

gpnorm2D = function(nd)
{
	ind = grep("^ *NUMFLDS=",nd)
	nfg = as.integer(sub(" *NUMFLDS= *(\\d+) .+","\\1",nd[ind]))
	ind = ind[nfg > 0]
	nfg = nfg[nfg > 0]
	# in CANARI, several prints of the setup of surface fields
	ind = ind[! duplicated(nd[ind-1])]

	if (length(ind) == 0) {
		cat("--> no GP norms found for surface fields\n")
		return(NULL)
	}

	surf = list()

	indo = grep("^ *GPNORM OUTPUT +AVERAGE",nd)

	for (i in seq(along=ind)) {
		# group name is (always?) printed 1 line before
		group = sub("^.+ (\\w+) +- +.+","\\1",nd[ind[i]-1])
		gnames = sub("^ *\\w+( +\\d+)+ +(\\w+(\\.\\w+)?).+","\\2",nd[ind[i]+seq(nfg[i])])
		gnames = substr(gnames,1,16)

		indh = grep(sprintf("\\<%s\\>",group),nd[indo-1])
		if (length(indh) == 0) {
			# group name is also always printed 1 line after
			group = sub("^ *(\\w+) +.+","\\1",nd[ind[i]+1])
			indh = grep(sprintf("\\<%s\\>",group),nd[indo-1])
			# group may not have GP norms
			if (length(indh) == 0) next
		}

		while (length(indh) > 0) {
			if (regexpr(", +FIELD +\\d+",nd[indo[indh[1]]-1]) > 0) {
				# nfg lines AVE, every 4 lines (group, GPNORM, AVE, 1)
				indi = indo[indh[1]]+(seq(nfg[i])-1)*4+1
				indh = indh[-seq(nfg[i])]
			} else if (regexpr(" \\d+ +FIELDS\\>",nd[indo[indh[1]]-1]) > 0) {
				# nfg lines AVE, every 3 lines (GPNORM, AVE, 1)
				indi = indo[indh[1]]+(seq(nfg[i])-1)*3+1
				indh = indh[-1]
			} else {
				# nfg lines after GPNORM and AVE
				indi = indo[indh[1]]+seq(nfg[i])+1
				indh = indh[-1]
			}

			gpn = t(line2num(nd[indi]))
			dimnames(gpn) = list(gnames,c("ave","min","max"))
			if (length(surf) < i) {
				surf[[i]] = gpn
			} else {
				# add this group, managing 2D and pseudo-2D fields
				n = length(surf[[i]])
				stopifnot(n %% length(gpn) == 0)
				surf[[i]] = array(c(surf[[i]],gpn),c(dim(gpn),n/length(gpn)+1))
				dimnames(surf[[i]])[1:2] = dimnames(gpn)
			}
		}

		names(surf)[i] = group
	}

	surf[sapply(surf,length) > 0]
}

fpgpnorm = function(nd,lev,tag=NULL,quiet=FALSE)
{
	ind = grep("^ *(FULL-POS +GPNORMS|GPNORMS +OF FIELDS)",nd)

	if (length(ind) == 0) {
		cat("--> no FP norms found\n")
		return(NULL)
	}

	if (! is.null(tag)) {
		indg = grep(tag,nd[ind-2],ignore.case=TRUE)
		indg = c(indg,grep(tag,nd[ind-1],ignore.case=TRUE))
		ind = ind[indg]
	}

	fpstep = unique(stepindex(nd,ind))
	if (length(fpstep) == 0) {
		cat("--> no FP norms found\n")
		return(NULL)
	}

	nt = length(fpstep)
	cat("-->",nt,"FP events found\n")

	nflevg = getvar("NFLEVG",nd)
	nfp3s = getvar("NFP3S",nd)
	indf = grep("^ *(\\w+|[. '/])+\\w+ *: +[+-]?\\d*\\.\\d+",nd)
	#indf = grep(sprintf("^ *(\\w+|[. '/])+\\w+ *:( +%s){1,3}$",Gnum),nd)
	for (i in seq(along=ind)) {
		ii = indf[indf >= ind[i]+1]
		ind2 = which(diff(ii) > 2)
		if (length(ind2) > 0) ii = ii[1:ind2[1]]

		i3d = grep("^ *S\\d+\\w+",nd[ii])
		if (length(i3d) == 0) {
			gpl = list()
		} else {
			noms = unique(sub("^ *S\\d+(\\w+.*\\>) *:.+$","\\1",nd[ii[i3d]]))
			# '/000': whole domain, leave it
			noms = gsub(" */000","",noms)
			gpl = vector("list",length(noms))
			names(gpl) = noms

			for (j in seq(along=noms)) {
				inds = grep(sprintf(" *S\\d+%s\\>.*? *:",noms[j]),nd[ii])
				if (length(inds) < 3 && nflevg > 2) {
					if (! quiet) cat("--> field",noms[j],"not 3D (maybe only top/bottom PP)\n")
					next
				}

				stopifnot(length(inds) %in% c(nfp3s,nflevg))
				fp = line2num(nd[ii[inds]])
				re = sprintf(" *S0*(\\d+)%s\\>.* *:.+",noms[j])
				levs = as.integer(sub(re,"\\1",nd[ii[inds]]))
				gpl[[j]] = fp[,order(levs)]
				attr(gpl[[j]],"nlev") = length(levs)
				ii = ii[-inds]
			}

			gpl = gpl[! sapply(gpl,is.null)]
		}

		if (length(ii) > 0) {
			gpl2 = lapply(strsplit(sub(".+: +","",nd[ii])," +"),as.numeric)
			noms = sub("^ *((\\w+|[. '/])+\\w+) *:.+","\\1",nd[ii])
			# '/000': whole domain, leave it
			noms = gsub(" */000","",noms)
			stopifnot(all(sapply(gpl2,length) == 3))
			for (j in seq(along=gpl2)) attr(gpl2[[j]],"nlev") = 1
			names(gpl2) = noms
			gpl = c(gpl,gpl2)
		}

		if (i == 1) {
			fpl = gpl
			next
		}

		indv = match(names(gpl),names(fpl))
		for (j in which(! is.na(indv))) {
			stopifnot(attr(fpl[[indv[j]]],"nlev") == attr(gpl[[j]],"nlev"))
			fpl[[indv[j]]] = c(fpl[[indv[j]]],gpl[[j]])
			attr(fpl[[indv[j]]],"nlev") = attr(gpl[[j]],"nlev")
		}

		iv = names(gpl) %in% names(fpl)
		if (! quiet && any(! iv)) cat("adding fields:",names(gpl)[! iv],"\n")
		if (any(! iv)) fpl = c(fpl,gpl[! iv])
	}

	indi = integer()

	for (i in seq(along=fpl)) {
		fpi = fpl[[i]]
		nlev = attr(fpi,"nlev")
		if (nlev == 1) {
			nti = length(fpi)%/%3
			if (length(fpi) != 3*nti) {
				if (! quiet) {
					cat("--> truncated events for",names(fpl)[i],nti,length(fpi)/3,"\n")
				}

				indi = c(indi,i)
            next
         }

			fp = array(dim=c(3,1,nt))
			#fp = matrix(nrow=3,ncol=nt)
			fp[,,seq(nti)] = fpi
			fpl[[i]] = aperm(fp)
			dimnames(fpl[[i]]) = list(seq(nt),1,c("ave","min","max"))
		} else {
			nti = length(fpi)%/%(3*nlev)
			if (length(fpi) != 3*nlev*nti) {
				if (! quiet) {
					cat("--> truncated events or levels for",names(fpl)[i],nti,length(fpi)/3,"\n")
				}

				indi = c(indi,i)
				next
			}

			fp = array(dim=c(3,nlev,nt))
			fp[,,seq(nti)] = fpi
			dim(fp) = c(3,nlev,nt)
			fpl[[i]] = aperm(fp)
			dimnames(fpl[[i]]) = list(seq(nt),seq(nlev),c("ave","min","max"))
		}

		attr(fpl[[i]],"nlev") = nlev
	}

	if (length(indi) == length(fpl)) return(NULL)

	if (length(indi) > 0) fpl = fpl[-indi]

	it = seq(along=fpstep)
	lev = as.integer(lev)
	if (identical(lev,0L)) {
		#fp = sapply(fpl,vmean,it,simplify="array")
		fp = sapply(fpl,function(x) apply(x,c(1,3),mean),simplify="array")
	} else {
		nlev = sapply(fpl,attr,"nlev")
		nlevx = max(nlev)
		stopifnot(all(lev %in% seq(nlevx)))
		fp3 = fpl[nlev == nlevx]
		if (length(fp3) == 0) return(NULL)

		fp = sapply(fp3,function(x) x[it,lev,,drop=FALSE],simplify="array")
	}

	# try setting names for steps, not sure it succeeds
	nt = dim(fp)[1]
	if (length(fpstep) < nt) {
		warning("fpstep not used")
	} else {
		dimnames(fp)[[1]] = fpstep[1:nt]
	}

	dn = dimnames(fp)
	if (length(dn) == 3) {
		d = dim(fp)
		dim(fp) = c(d[1],1,d[2:3])
		dimnames(fp) = c(dn[1],list(lev),dn[2:3])
	}

	fp
}

vmean = function(x,it)
{
	if (attr(x,"nlev") == 1) {
		x[it,,drop=FALSE]
	} else {
		apply(x[it,,,drop=FALSE],c(1,3),mean)
	}
}

fpspnorm = function(nd,lev)
{
	ind = grep("^ *(FULL-POS +)?SPNORMS",nd)

	if (length(ind) == 0) {
		cat("--> no FP norms found\n")
		return(NULL)
	}

	nflevg = getvar("NFLEVG",nd)
	nfp3s = getvar("NFP3S",nd)
	nt = length(ind)

	cat("-->",nt,"FP events found\n")

	resp = "(S\\d+)?(\\w+(\\.\\w+)*)"
	indf = grep(sprintf("^ *%s *: %s$",resp,Gnum),nd)
	indf = indf[indf > ind[1]]

	ndfp = nd[indf]

	noms = unique(sub(sprintf("^ *%s.+",resp),"\\2",ndfp))
	i3d = grep("^ *S\\d+\\w+",ndfp)
	noms3d = unique(sub("^ *S\\d+(\\w+(\\.\\w+)*).+","\\1",ndfp[i3d]))
	fp3d = array(dim=c(nt,nfp3s,length(noms3d)),dimnames=list(1:nt,1:nfp3s,noms3d))
	noms2d = noms[! noms %in% noms3d]
	fp2d = array(dim=c(nt,length(noms2d)),dimnames=list(1:nt,noms2d))

	for (i in seq(along=ind)) {
		ii = indf[indf >= ind[i]+1]
		ind2 = which(diff(ii) > 2)
		if (length(ind2) > 0) ii = ii[1:ind2[1]]

		i3d = grep("^ *S\\d+\\w+",nd[ii])
		if (length(i3d) == 0) next

		for (j in seq(along=noms3d)) {
			inds = grep(sprintf(" *S\\d+%s\\>.*? *:",noms3d[j]),nd[ii])
			if (length(inds) < 3 && nflevg > 2) {
				if (! quiet) cat("--> field",noms[j],"not 3D (maybe only top/bottom PP)\n")
				next
			}

			stopifnot(length(inds) %in% c(nfp3s,nflevg))
			fp = line2num(nd[ii[inds]])
			re = sprintf(" *S0*(\\d+)%s\\>.* *:.+",noms3d[j])
			levs = as.integer(sub(re,"\\1",nd[ii[inds]]))
			fp3d[i,,j] = fp[order(levs)]
			ii = ii[-inds]
		}

		if (length(ii) > 0) {
			fp2d[i,] = sapply(strsplit(sub(".+: +","",nd[ii])," +"),as.numeric)
		}
	}

	if (length(lev) > 1) {
		stopifnot(all(lev %in% seq(nfp3s)))
		fp = fp3d[,lev,,drop=FALSE]
	} else if (lev == 0) {
		fp3m = apply(fp3d,c(1,3),mean)
		fp = cbind(fp2d,fp3m)
	} else {
		stopifnot(lev %in% seq(nfp3s))
		fp = fp3d[,lev,]
	}

	fp
}

stepindex = function(nd,ind,prestepo=FALSE)
{
	icnt4 = grep("START CNT4",nd)
	if (length(icnt4) == 0) {
		warning("no line matching 'START CNT4'")
		return(ind)
	}

	if (prestepo) {
		indt = grep("^ *NORMS AT (START|NSTEP|END) CNT4(TL|AD)?",nd)
		nstep = sub(" *NORMS AT NSTEP CNT4(TL|AD)?( \\((PREDICTOR|(C)ORRECTOR)\\))? +(\\d+)",
			"\\1\\4\\5",nd[indt])
		ix = grep("^ *NORMS AT (START|END) CNT4",nd[indt])
		nstep[ix] = sub(" *NORMS AT (START|END) CNT4(TL|AD)?.*","\\2\\1",nd[indt[ix]])
		if (length(ix) > 0) {
			inds = grep("^ *START CNT4(TL|AD)",nd)
			ixs = sapply(ix,function(i) rev(which(inds < indt[i]))[1])
			nsim = sub("^ *START CNT4(TL|AD), +NSIM4D *= *(\\d+).*","\\2",nd[inds[ixs]])
			nstep[ix] = paste(nstep[ix],nsim,sep="")
		}
	} else {
		indt = grep("^ *NSTEP = *\\d+ +STEP[OX](TL|AD)? +",nd)
		nstep = sub("^ *NSTEP = *(\\d+) +STEP(O|(X))(TL|AD)? +.*","\\3\\4\\1",nd[indt])
	}

	jstep = rep(NA_character_,length(ind))
	for (i in seq(along=ind)) {
		it = which(indt <= ind[i])
		if (ind[i] < icnt4[1]) {
			jstep[i] = "-1"
		} else if (length(it) > 0) {
			jstep[i] = nstep[it[length(it)]]
		}
	}

	jstep
}

diagstep = function(nd,type,nstop)
{
	stopifnot(type %in% c("SP","GP","FP","HIS","SFXHIS","RES","AN"))

	if (missing(nstop)) nstop = getvar("NSTOP",nd)

	if (type == "FP") {
		tag = "NFRPOS"
		nfr = tag
	} else {
		tag = switch(EXPR=type,"GP"="NFRGDI","SP"="NFRSDI","HIS"=,"RES"=,"AN"="NFRHIS",
			"SFXHIS"="NFRSFXHIS")
		nfr = sprintf(".+ %s",tag)
	}

	# don't know steps for RES and AN...
	if (type == "SP") {
		i1 = grep("^ *SPECTRUM, *ISDITS\\>",nd)
		i2 = grep("^ *IDHFGTS\\>",nd)
	} else if (type == "HIS") {
		i1 = grep("^ *HISTORY WRITE-UP, *IHISTS\\>",nd)
		i2 = grep("^ *HISTORY SURFACE WRITE-UP, *ISFXHISTS\\>",nd)
	} else if (type == "SFXHIS") {
		i1 = grep("^ *HISTORY SURFACE WRITE-UP, *ISFXHISTS\\>",nd)
		i2 = grep("^ *MASS CONSERVATION FIXUP, IMASSCONS\\>",nd)
	} else if (type == "FP") {
		i1 = grep("^ *POST-PROCESSING EVENTS, *IPOSTS\\>",nd)
		i2 = grep("^ *FULLPOS SUB-OBJECT EVENTS, *IFPOSTS\\>",nd)
	} else if (type == "RES") {
		i1 = grep("^ *RESTART WRITE-UP, *IRESTS\\>",nd)
		i2 = grep("^ *SPECTRUM, *ISDITS\\>",nd)
	} else if (type == "AN") {
		i1 = grep("^ *AN WRITE-OUT EVENTS, *IANATS\\>",nd)
		i2 = grep("^ *\\d+",nd,invert=TRUE)
		i2 = i2[which(i2 > i1)[1]]
	}

	freq = getvar(nfr,nd)
	step = seq(0,nstop,by=freq)

	ind = seq(i1+1,i2-1)
	events = intlines(nd[ind])
	step = step[which(events == 1)]
	attr(step,"nstop") = nstop
	attr(step,"tag") = tag
	attr(step,"freq") = freq
	attr(step,"events") = events
	step
}

countfield = function(ind,indg,nl2)
{
	if (length(ind) == 1) return(length(ind))

	# indt: 1st and 2nd occurrences of GPNORM heading lines (ie steps 1 and 2)
	if (length(indg) == 1) {
		indt = ind[indg:length(ind)]
	} else {
		indt = ind[seq(indg[1],indg[2])]
	}

	# find the 1st line among ind separated by nl2, if any
	if (all(diff(indt) <= nl2)) return(length(indt))

	which(diff(indt) > nl2)[1]
}

indexpand = function(nd,ind,indg,nl2,average=TRUE)
{
	# ind: GPNORM heading lines
	# indg: among ind, 1 line befone GPNORM heading lines, lines for a given group name

	# nf: nb of fields in the group
	nf = countfield(ind,indg,nl2)

	# offset: offsets in ind of GPNORM heading lines for nf fields
	indt1 = ind[indg[1]+seq(nf)-1]
	dist = diff(indt1)
	offset = c(0,cumsum(dist))
	if (! average) {
		indl1 = ind[indg[1]+seq(nf)-1]+2
		i3d = which(regexpr(sprintf("^ *1 +%s",Gnum),nd[indl1]) > 0)
		offset = offset[i3d]
	}

	rep(ind[indg],each=length(offset))+offset
}

shortnames = function(x)
{
	x[x == "KINETIC ENERGY"] = "KIN ENERGY"
	x[x == "LOG(PRE/PREHYD)"] = "LOG(P/P_hyd)"
	x[grep("d4 *= *VERT DIV *\\+ *X",x)] = "d4 (= vdiv+X)"
	x[grep("d5 *= *VERT DIV *\\+ *X",x)] = "d5 (= vdiv+X+Xs)"
	x[x == "SURFACE PRESSURE"] = "SURF P"
	x[x == "TEMPERATURE"] = "TEMP"
	x[x == "U VELOCITY"] = "U WIND"
	x[x == "V VELOCITY"] = "V WIND"

	x
}

varqc = function(nd)
{
	indo = grep("^ *VAR *= *\\d+",nd)
	iv = as.integer(sub("VAR *= *(\\d+).+","\\1",nd[indo]))
	qc = vector(length(unique(iv)),mode="list")
	for (i in iv) {
		ind = which(iv == i)
		notvar = intlines(sub(".+ NOTVAR *=","",nd[indo[ind[1]]]))
		v = numlines(gsub("\\*+","-9999.",nd[indo[ind[-1]]]))
		v[v==-9999] = NA
		qc[[i]] = matrix(c(v,notvar),nc=length(ind))
	}

	qc
}

jotable = function(nd)
{
	iobst = grep("Obstype +\\d+ +=+",nd,ignore.case=TRUE)
	ijog = grep("Jo Global",nd)
	ndo = nd[iobst[1]:ijog[1]]
	ijoh = grep("Codetype +\\d+ +=+",ndo)
	ijot = grep("^ +\\w+ +\\d+( +\\d+\\.\\d+){3}",ndo)
	lj = strsplit(ndo[ijot],split=" +")
	jot = t(sapply(lj,function(x) as.numeric(x[3:6])))
	jot = as.data.frame(jot)
	names(jot) = c("DataCount","Jo_Costfunction","Jo/n","ObsErr")
	vars = sapply(lj,"[",2)

	indi = findInterval(ijot,ijoh)
	code = sub(" +Codetype +(\\d+) +.+","\\1",ndo[ijoh[indi]])

	jot = cbind(Codetype=as.integer(code),Variable=vars,jot)

	jot
}

runtime = function(nd)
{
   indw = grep("STEP +\\d+ +H=.+\\+CPU=",nd)
   walls = as.difftime(gsub("^ *([[:digit:]:]+) .+","\\1",nd[indw]),units="secs")
   cpus = as.difftime(as.numeric(gsub(".+\\+CPU= *","",nd[indw])),units="secs")

	dwalls = diff(walls)

	# in case of change of date (time goes to 00:00)
	ind = which(dwalls < 0)
	for (i in ind) dwalls[-(1:i)] = dwalls[-(1:i)]+86400

	# small escalating over steps within 1s
	i1 = 1
	ind = which(c(dwalls,dwalls[length(dwalls)]+1) > 0)
	for (i in ind) {
		if (i > i1) {
			n = i-i1+1
			dt = seq(0,1,length.out=n+1)[-(n+1)]
			walls[i1:i] = walls[i1]+dt
		}

		i1 = i+1
	}

	rt = data.frame(wall=walls,dwall=c(0,dwalls),cpu=cpus)

	i1 = grep("TIME OF START *=",nd)
	attr(rt,"start") = as.difftime(gsub("^ *TIME OF START *= *","",nd[i1]),units="secs")

	rt
}
