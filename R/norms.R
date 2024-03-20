Gndigits = round(log10(.Machine$double.eps))

digitsdiff = function(x,y)
{
	x0 = pmax(abs(x),abs(y))
	ndiff = pmax(Gndigits,1+log10(abs(y-x))-log10(x0))-Gndigits
	ndiff[which(x0 == 0)] = 0
	#if (any(is.nan(ndiff))) warning("NaN in ndiff")
	ndiff
}

plotmean = function(x,y,main="Norm",legend,xunit,col,lty,scale=1,...)
{
	if (is.vector(y)) y = as.matrix(y)

	if (missing(col)) col = seq(dim(y)[2])
	if (missing(lty)) {
		lty = seq(dim(y)[2])
		lty[-1] = 2
	}

	y = y*scale
	matplot(x,y,type="l",lty=lty,col=col,main=main,...)

	if (! missing(legend)) legend("topleft",legend,col=col,lty=lty,bg="transparent")

	if (all(y == 0,na.rm=TRUE)) {
		text(sum(range(x))/2,.5,"all values = 0")
		return()
	}

	if (! missing(xunit)) {
		reg = line(x,y[,1])
		abline(reg,col="darkgrey",lty=2)
		tend = coef(reg)[2]*86400/xunit
		tt = sprintf("Tend: %+.3e [unit]/day",tend)
		if (xunit == 86400) tt = sprintf("%s, variation: %.2g",tt,tend*diff(range(x)))
		if (scale != 1) tt = sprintf("%s, scaling: *%.0e",tt,scale)

		mtext(tt,cex=par("cex"))
	}
}

plotmnx = function(x,y,main="GP Norm",imnx=1:3,col,lty=c(1,3,3),legend,ylim,...)
{
	tt = c("ave","min","max")
	tt = paste(tt[imnx],collapse="/")
	if (length(main) == 1) main = c(main,tt)

	if (length(dim(y)) < 3) dim(y)[3] = 1
	if (missing(col)) col = seq(dim(y)[3])
	if (missing(ylim)) ylim = range(y[,imnx,],na.rm=TRUE)

	matplot(x,y[,imnx,1],type="l",col=col[1],lty=lty[imnx],lwd=1.05,main=main,ylim=ylim,
		...)

	if (! missing(legend)) legend("topleft",legend,col=col,lty=1,bg="transparent")

	for (k in seq(dim(y)[3])[-1]) {
		if (all(is.na(y[,imnx,k]))) next

		matlines(x,y[,imnx,k],lty=lty[imnx],col=col[k])
	}
}

plotvmean = function(x,y,main="Norm",col=1,lty,ylab="Level",ylim=rev(range(y)),legend,...)
{
	if (is.vector(x)) x = as.matrix(x)
	if (missing(lty)) lty = seq(dim(x)[2])

	matplot(x,y,type="l",col=col,lty=lty,main=main,ylab=ylab,ylim=ylim,...)
	abline(v=0,col="darkgrey",lty=2)

	if (! missing(legend)) legend("topleft",legend,col=col,lty=lty,bg="transparent")
}

plotvmnx = function(x,y,main="Norm",col=1,lty,ylab="Level",legend,ylim=rev(range(y)),...)
{
	tt = c("ave","min","max")
	if (length(main) > 1) tt = sprintf("%s, %s",main[2],tt)

	if (missing(lty)) lty = seq(dim(x)[2])

	for (i in c(2,1,3)) {
		matplot(x[,,i],y,type="l",col=col,lty=lty,main=c(main[1],tt[i]),ylab=ylab,ylim=ylim,
			...)
		abline(v=0,col="darkgrey",lty=2)

		if (! missing(legend)) legend("topleft",legend,col=col,lty=lty,bg="transparent")
	}
}
