
getFrame = function(file)
{
	Gficbin = tempfile(fileext=".bin")
	system(sprintf("epy_dump.py %s -f frame -o %s",file,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	dims = readBin(con,"integer",5,size=8)

	# global (Gauss grid) or LAM (Cartesian grid)
	if (dims[5] == 0) {
		nlong = readBin(con,"integer",dims[1],size=8)
		nwave = readBin(con,"integer",dims[2],size=8)
		ngp = dims[4]
		nl = dims[3]
	} else {
		ngp = prod(dims[1:2])
		nl = dims[5]
	}

	base = readBin(con,"numeric",1)
	base = as.POSIXct(base,origin="1970-01-01")
	step = readBin(con,"integer",1,size=8)
	if (regexpr(".+\\+0*([0-9]+)",file) > 0) {
		ech = as.integer(gsub(".+\\+0*([0-9]+)","\\1",file))
		if (step != ech*3600) step = ech*3600
	}

	lats = readBin(con,"numeric",ngp)
	longs = readBin(con,"numeric",ngp)

	Gvp0 = 101325
	Ah = readBin(con,"numeric",nl+1)
	Bh = readBin(con,"numeric",nl+1)
	eta = ((Ah[-1]+Ah[-(nl+1)])/Gvp0+(Bh[-1]+Bh[-(nl+1)]))/2

	close(con)

	if (dims[5] == 0) {
		nlat = length(nlong)
		theta = 90-180*(seq(nlat)-.5)/nlat

		frame = list(nlat=nlat,nwave=dims[2],nlevel=nl,npdg=ngp,nlong=nlong,theta=theta,
			lat=lats,long=longs,eta=eta,base=base,step=step,gauss=TRUE,lam=FALSE)
	} else {
		nlong = rep(dims[2],dims[1])
		frame = list(nlat=dims[1],nlong=nlong,nwavex=dims[3],nwavey=dims[4],nlevel=nl,
			npdg=ngp,lat=lats,long=longs,eta=eta,base=base,step=step,gauss=FALSE,lam=TRUE)
	}

	frame
}

getField = function(fic,param,frame,symbol,frlow=frame,interp=FALSE)
{
	fsave = sprintf("%s/e%d/%s.RData",dirname(fic),frame$step/3600,symbol)
	if (file.exists(fsave) && ! is.null(frlow$ilev)) {
		ilev = 0
		load(fsave)
		if (dim(data)[1] == frlow$npdg && identical(ilev,frlow$ilev)) return(data)

		cat("--> different grid or levels, read data again\n")
	}

	data = getGPFields(fic,param,frame)

	if (frame$npdg != frlow$npdg) {
		if (interp) {
			data = interpGauss(data,frame,frlow)
		} else {
			stopifnot("ind" %in% names(frlow))
			data = data[frlow$ind,,drop=FALSE]
		}
	}

	ilev = frlow$ilev

	if (dim(data)[2] > 1) {
		if (frame$nlevel != frlow$nlevel) data = interpAB(data,frame,frlow)

		data = data[,ilev,drop=FALSE]
	}

	if (! file.exists(dirname(fsave))) dir.create(dirname(fsave))
	save("data","ilev",file=fsave)

	data
}

getGPFields = function(file,field,frame)
{
	Gficbin = tempfile(fileext=".bin")
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	nlevel = readBin(con,"integer",1,size=8)
	stopifnot(length(nlevel) == 1)

	data = matrix(nrow=frame$npdg,ncol=nlevel)
	for (j in seq(nlevel)) data[,j] = readBin(con,"numeric",frame$npdg)
	stopifnot(nlevel == readBin(con,"integer",1,size=8))
	close(con)

	data
}

getSPFields = function(file,field,frame)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	nsp2 = readBin(con,"integer",1,size=8)
	data = readBin(con,"numeric",nsp2)
	close(con)

	stopifnot(nsp2 == frame$nwave*(frame$nwave-1))

	matrix(data,nrow=frame$nwave)
}
