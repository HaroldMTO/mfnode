csLat = function(frame)
{
	mucen = sin(pi/180*latPole(frame))
	locen = pi/180*lonPole(frame)
	C = 2.4
	c2 = C^2

	sqm2 = sqrt(1-mucen^2)
	sinlocen = sin(locen)
	coslocen = cos(locen)

	lat = numeric(frame$nlat)

	ind = frame$lat[cumsum(frame$nlong)-frame$nlong/2]
	cslon = pi

	for (ilat in seq(frame$nlat)) {
		mu = sin(pi/180*frame$theta[ilat])
		sqmu2 = sqrt(1-mu^2)
		a = 1/(c2+1+(c2-1)*mu)
		b = c2-1+(c2+1)*mu
		ca1 = a*(2*C*mucen*sqmu2-b*sqm2*cos(cslon))
		lat[ilat] = a*(2*C*sqm2*sqmu2*cos(cslon)+b*mucen)
	}

	lat
}

csLong = function(frame)
{
	longs = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (i in seq(frame$nlat)) {
		nl = frame$nlong[i]
		longs[seq(nl),i] = 360*(seq(nl)-1)/nl
	}

	longs
}

latPole = function(frame)
{
	n1 = frame$nlong[1]

	if (all(duplicated(frame$lat[1:n1])[-1])) {
		latcen = 90
	} else if (diff(range(frame$long[1:n1])) < 180) {
		# lat of Pole is so that Pole is out of 1st lat circle
		latcen = mean(frame$lat[c(1,n1/2+1)])
	} else {
		# lat of Pole is close to geo North Pole and is inside 1st lat circle
		# latcen = (latx+2*(90-latx)+latn)/2 = 90+(latn-latx)/2 = 90-(latx-latn)/2
		latcen = 90-abs(diff(frame$lat[c(1,n1/2+1)]))/2
	}

	latcen
}

lonPole = function(frame)
{
	n1 = frame$nlong[1]
	frame$long[n1/2+1]
}

compass = function(frame)
{
	mucen = sin(pi/180*latPole(frame))
	locen = pi/180*lonPole(frame)
	C = 2.4
	c2 = C^2

	sqm2 = sqrt(1-mucen^2)
	sinlocen = sin(locen)
	coslocen = cos(locen)

	l = m = numeric(frame$npdg)

	off = 0
	for (ilat in seq(frame$nlat)) {
		nlon = frame$nlong[ilat]
		cslon = 2*pi*seq(0,nlon-1)/nlon
		mu = sin(pi/180*frame$theta[ilat])
		sqmu2 = sqrt(1-mu^2)
		a = 1/(c2+1+(c2-1)*mu)
		b = c2-1+(c2+1)*mu
		ca1 = a*(2*C*mucen*sqmu2-b*sqm2*cos(cslon))
		gemu = a*(2*C*sqm2*sqmu2*cos(cslon)+b*mucen)
		rcoslat = 1/sqrt(1-gemu^2)
		l[off+1:nlon] = -sqm2*sin(cslon)*rcoslat
		m[off+1:nlon] = ca1*rcoslat
		off = off+nlon
	}

	data.frame(l=l,m=m)
}

interpGauss = function(datao,framo,frame)
{
	stopifnot(frame$gauss)

	# only for (Gauss) grids nested inside
	stopifnot(length(framo$nlong) >= length(frame$nlong))
	stopifnot(max(framo$nlong) >= max(frame$nlong))

	longs = csLong(frame)

	data = matrix(nrow=frame$npdg,ncol=dim(datao)[2])

	clats = c(0,cumsum(frame$nlong))
	clato = c(0,cumsum(framo$nlong))
	dlat = 180/framo$nlat
	dlon = 360/framo$nlong

	for (ip in seq(frame$npdg)) {
		ilat = which.max(ip <= clats[-1])
		ilon = ip-clats[ilat]
		e0 = (90-dlat/2-frame$theta[ilat])/dlat
		ilato = floor(e0)+1
		e0 = e0-(ilato-1)
		stopifnot(0 < ilato && ilato < length(framo$nlong))
		stopifnot(0 <= e0 && e0 < 1)

		e1 = longs[ilon,ilat]/dlon[ilato]
		e2 = longs[ilon,ilat]/dlon[ilato+1]
		ilono1 = floor(e1)+1
		ilono2 = floor(e2)+1
		stopifnot(0 < ilono1 && ilono1 <= framo$nlong[ilato])
		stopifnot(0 < ilono2 && ilono2 <= framo$nlong[ilato+1])
		e1 = e1-(ilono1-1)
		e2 = e2-(ilono2-1)
		stopifnot(0 <= e1 && e1 < 1)
		stopifnot(0 <= e2 && e2 < 1)
		data[ip,] = interp(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
	}

	data
}

interp = function(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
{
	ip1 = clato[ilato]+ilono1
	ip2 = ip1%%clato[ilato+1]+1
	d1 = datao[ip1,]+e1*(datao[ip2,]-datao[ip1,])
	ip1 = clato[ilato+1]+ilono2
	ip2 = ip1%%clato[ilato+2]+1
	d2 = datao[ip1,]+e2*(datao[ip2,]-datao[ip1,])
	d1+e0*(d2-d1)
}

interpAB = function(datao,framo,frame)
{
	eta = frame$A/Gvp0+frame$B
	etao = framo$A/Gvp0+framo$B

	ind = findInterval(eta,etao)
	e = (eta-etao[ind])/(etao[ind+1]-etao[ind])
	stopifnot(all(0 < ind & ind < length(etao)))
	stopifnot(all(0 <= e & e < 1))

	data = matrix(nrow=dim(datao)[1],ncol=frame$nlevel)

	for (i in seq(frame$nlevel)) {
		d1 = datao[,ind[i]]
		d2 = datao[,ind[i]+1]
		data[,i] = d1+e[i]*(d2-d1)
	}

	data
}

makeDomain = function(frame,ind,xlim,ylim)
{
	stopifnot(length(ind) == 1)
	stopifnot(length(ylim) == 2 && length(xlim) == 2)
	stopifnot(diff(ylim) > 0 && diff(xlim) > 0)

	xlim = xlim+frame$long[ind]
	ylim = ylim+frame$lat[ind]
}

getDomain = function(frame)
{
	xlim = longRange(frame$long)
	ylim = range(frame$lat)

	list(xlim=xlim,ylim=ylim)
}

longRange = function(long)
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

inDomain = function(points,domain)
{
	ilat = domain$ylim[1] <= points$lat & points$lat <= domain$ylim[2]

	# case xlim=[0,0] (whole globe) is in alternative
	xlim = (domain$xlim+180)%%360-180
	if (diff(xlim) > 0) {
		ilon = xlim[1] <= points$long & points$long <= xlim[2]
	} else {
		ilon = points$long <= xlim[2] | points$long >= xlim[1]
	}

	ilat & ilon
}

area = function(dom1,dom2)
{
	diff(dom1$xlim)*diff(dom1$ylim)/(diff(dom2$xlim)*diff(dom2$ylim))
}

select = function(frame,ind)
{
	list(lat=frame$lat[ind],long=frame$long[ind])
}

zoom = function(data,frame,domain)
{
	mask = inDomain(frame,domain)
	ind = which(mask)

	lam = frame$lam || diff(domain$xlim) < 360
	list(lat=frame$lat[ind],long=frame$long[ind],lam=lam,data=data[ind,,drop=FALSE])
}

sectiongeo = function(data,frame,long=c(0,360),lat=c(-90,90))
{
	stopifnot(all(-90 <= lat & lat <= 90))

	# frame longs belong to [-180,180[
	if (length(long) == 1) {
		stopifnot(diff(lat) > 0)
		xf = frame$lat
		yf = frame$long
		x = lat
		y = (long+180)%%360-180
	} else if (length(lat) == 1) {
		stopifnot(diff(long) > 0)
		xf = frame$long
		yf = frame$lat
		x = (long+180)%%360-180
		y = lat
	} else {
		stop("lat or long must be of length 1")
	}

	clats = c(0,cumsum(frame$nlong))

	data1 = array(dim=c(2,frame$nlat,dim(data)[2]))
	x1 = array(dim=c(2,frame$nlat))

	for (ilat in seq(frame$nlat)) {
		off = clats[ilat]
		nlon = frame$nlong[ilat]
		yfn = yf[off+1:nlon]
		yi = y
		if (frame$gauss) {
			# 1 more point in longitude (1st one) since frame is a Gaussian global grid
			yfn = c(yfn,yfn[1])
			il = which(abs(diff(yfn)) > 180)
			stopifnot(length(il) <= 2)
			for (i in rev(il)) yfn[-(1:i)] = yfn[-(1:i)]-360*sign(diff(yfn)[i])
			if (yi < min(yfn)) yi = yi+360

			stopifnot(all(diff(yfn) > -180))
		}

		ind = yfn[-nlon] <= yi & yi < yfn[-1] | yfn[-1] <= yi & yi < yfn[-nlon]
		if (all(! ind)) next

		xfn = xf[off+1:nlon]
		i1 = which(ind)
		stopifnot(length(i1) <= 2)
		i2 = i1+1
		i2[i1 == nlon] = 1
		e = (y-yfn[i1])/(yfn[i2]-yfn[i1])
		stopifnot(all(0 <= e & e < 1))

		for (i in seq(along=i1)) {
			data1[i,ilat,] = (1-e[i])*data[off+i1[i],]+e[i]*data[off+i2[i],]
			x1[i,ilat] = (1-e[i])*xfn[i1[i]]+e[i]*xfn[i2[i]]
		}
	}

	x1 = as.vector(x1)
	if (all(is.na(x1))) stop("no lat/long crossing given parameter")

	ll = sort(x1,index.return=TRUE,na.last=TRUE)
	ind = ll$ix[1:length(na.omit(x1))]
	if (diff(x) > 0) {
		ii = which(x[1] <= x1[ind] & x1[ind] <= x[2])
	} else {
		ii = which(x[1] <= x1[ind] | x1[ind] <= x[2])
	}

	ind = ind[ii]

	dim(data1) = c(2*frame$nlat,dim(data)[2])

	if (length(long) == 1) {
		list(lats=x1[ind],long=long,data=data1[ind,])
	} else {
		list(longs=x1[ind],lat=lat,data=data1[ind,])
	}
}

DOPfill.legend = function(levels,col,...)
{
	u = par("usr")
	p = par("plt")

	width = (1-p[2])/6*diff(u[1:2])/diff(p[1:2])
	x = u[2]+width/3

	nl = length(levels)
	height = diff(u[3:4])
	dy = height/(nl-1)
	y = u[3]
	ybas = y + dy*(seq(nl-1)-1)
	yhaut = ybas + dy

	rect(x,ybas,x+width,yhaut,col=col,border=NA,xpd=TRUE)

	op = par(las=2,yaxt="s")
	axis(4,c(ybas[1],yhaut),levels,tick=FALSE,pos=x,mgp=c(1,.5,0),...)
	par(op)
}

mapxy = function(xy,axes=TRUE,frame.plot=TRUE,...)
{
	l = map(xlim=xy$xlim,ylim=xy$ylim,...)

	if (axes) {
		axis(1)
		axis(2)
	}

	if (frame.plot) box()

	l
}

plotxy = function(x,y,data,breaks="Sturges",palette=terrain.colors,pch=15,ppi=50,
	cex.leg=.8,...)
{
	u = par("usr")
	f = par("fin")

	if (ppi > 144) stop("ppi > 144\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(data) > npmax) {
		cat("--> reducing xy plot from",length(data),"to",npmax,"points\n")
		ind = seq(1,length(data),length.out=npmax)
		dn = diff(ind[1:2])%/%3
		if (dn > 2) {
			ind[seq(2,length(ind)-1,by=3)] = ind[seq(2,length(ind)-1,by=3)]+dn%/%2
			ind[seq(3,length(ind),by=3)] = ind[seq(3,length(ind),by=3)]-dn%/%2
		}

		x = x[ind]
		y = y[ind]
		data = data[ind]
	}

	h = hist(data,breaks,plot=FALSE)
	br = h$breaks
	if (is.character(breaks) && length(br) > 8) {
		h = hist(data,breaks=8,plot=FALSE)
		br = h$breaks
	}

	ind = findInterval(data,br)
	cols = palette(length(br))

	points(x,y,col=cols[ind],pch=pch,...)

	levels = sprintf("% .3g",br)
	DOPfill.legend(levels,col=cols,cex.axis=cex.leg)
}

plotxy2 = function(x,y,zx,zy,breaks="Sturges",colvec=TRUE,length=.05,angle=15,ppi=8,
	cex.leg=.8,...)
{
	f = par("fin")

	if (ppi > 96) stop("ppi > 96\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(zx) > npmax) {
		cat("--> reducing xy2 plot from",length(zx),"to",npmax,"points\n")
		ind = as.integer(seq(1,length(x),length.out=npmax))
		dn = diff(ind[1:2])%/%3
		if (dn > 2) {
			ind[seq(2,length(ind)-1,by=3)] = ind[seq(2,length(ind)-1,by=3)]+dn%/%2
			ind[seq(3,length(ind),by=3)] = ind[seq(3,length(ind),by=3)]-dn%/%2
		}

		x = x[ind]
		y = y[ind]
		zx = zx[ind]
		zy = zy[ind]
	}

	ffz = sqrt(zx^2+zy^2)

	u = par("usr")
	ux = diff(u[1:2])
	uy = diff(u[3:4])

	# scale of wind vectors: a fraction of the smallest side
	asp = min(ux,uy)/max(ux,uy)
	ffu = asp*sqrt((ux^2+uy^2)/length(zx))

	if (colvec) {
		fx = zx/ffz*ffu
		fy = zy/ffz*ffu

		h = hist(ffz,breaks,plot=FALSE)
		br = h$breaks
		if (is.character(breaks) && length(br) > 8) {
			h = hist(ffz,breaks=8,plot=FALSE)
			br = h$breaks
		}

		ic = findInterval(ffz,br)
		palette = rainbow(length(br),start=2/3,end=0)
		cols = palette[ic]
	} else {
		ff1 = ffu*scale(ffz,center=FALSE)
		fx = zx/ffz*ff1
		fy = zy/ffz*ff1
		cols = rep(1,length(ff1))
	}

	x2 = x+fx
	y2 = y+fy

	ffi = sqrt((fx/ux*f[1])^2+(fy/uy*f[2])^2)

	ind = ffi > 1.e-3
	arrows(x[ind],y[ind],x2[ind],y2[ind],length,angle,col=cols[ind],...)
	if (any(! ind)) points(x[! ind],y[! ind],pch=1,col=1)

	if (colvec) {
		levels = sprintf("% .3g",br)
		DOPfill.legend(levels,col=palette,cex.axis=cex.leg)
	}
}

mapdom = function(dom,points,data,main=NULL,mar=par("mar"),...)
{
	l = mapxy(dom,mar=mar)
	plotxy(points$long,points$lat,data,...)
	lines(l)
	title(main)
}


mapdom2 = function(dom,points,datax,datay,main=NULL,mar=par("mar"),...)
{
	l = mapxy(dom,mar=mar)
	plotxy2(points$long,points$lat,datax,datay,...)
	lines(l)
	title(main)
}
