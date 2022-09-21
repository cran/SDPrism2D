#' Standard Deviation Visualization
#' @export
sdprism2d <- function(data, hlim=NULL, xyscale=NULL) {
  # Find pd
  if(is.null(hlim)){
    hlim <- 4
  }
  if(is.null(xyscale)){
    xyscale <- 4
  }
  r <- xyscale
  pd <- c()
  cnt = 1
  for (i in data[1:length(data)-1]) {
    for (j in data[0:-cnt]){
      pd <- append(pd, abs(j-i))
    }
    cnt <- cnt + 1
  }

  # Function to draw a prism
  rprism <- function(hyp, xstart, y, t, ratio) {
    x1<-hyp/2
    x2<-x1/ratio
    segments(xstart,y,hyp,y,lty=1,lwd = 1.75)
    segments(hyp,y,hyp,0,lty = 2, lwd = 1, col = "gray")
    segments(hyp,y,x1,y+x2,lty=1,lwd = 2)
    segments(x1,y+x2,0,y,lty=1,lwd = 2)
    segments(0,y,0,y+t,lty=1,lwd = 1.75)
    segments(0,y+t,x1,y+x2+t,lty=1,lwd = 1.75)
    for(i in c(seq(0,x1,.01))) {
      segments(i,y+i/ratio+.015,i,y+i/ratio+t-.015,lty=1,col = "#778899")
    }
    for(i in c(seq(x1,hyp,.01))){
      segments(i,y+(hyp-i)/ratio+.015,i,y+(hyp-i)/ratio+t-.015,lty = 1, col = "#C0C0C0")
    }
    segments(x1,y+x2+t,hyp,y+t,lty=1,lwd = 1.75)
    segments(hyp,y+t,hyp,y,lty=1,lwd = 1.75)
    segments(x1,y+x2,x1,y+x2+t,lty=1,lwd = 1.75)
  }

  # Function to draw the sd prism
  rmspdprism <- function(hyp,ratio) {
    x1<-hyp/2
    x2<-x1/ratio
    segments(0,0,hyp,0,lty=1,lwd = 2, col = "red")
    segments(x1,x2,(hyp+x1)/2-.05*ratio,x2/2+.05,lty=1,lwd = 2, col = "red")
    segments((hyp+x1)/2+.05*ratio,x2/2-.05,hyp,0,lty=1,lwd = 2, col = "red")
    segments(x1,x2,x1/2+.05*ratio,x2/2+.05,lty=1,lwd = 2, col = "red")
    segments(x1/2-.05*ratio,x2/2-.05,0,0,lty=1,lwd = 2, col = "red")
    segments(0,0,0,1,lty=1,lwd = 2, col = "red")
    segments(0,1,x1,x2+1,lty=1,lwd = 2, col = "red")
    segments(x1,x2+1,hyp,1,lty=2,lwd = 2, col = "red")
    segments(hyp,1,hyp,0,lty=1,lwd = 2, col = "red")
    segments(x1,x2,x1,x2+1,lty=1,lwd = 2, col = "red")
    text(x1/2,x2/2, expression("s"), xpd=TRUE, col="red", cex = 1.3, srt=45)
    text((hyp+x1)/2,x2/2, expression("s"), xpd=TRUE, col="red", cex = 1.3, srt=-45)
    text(hyp/(sqrt(2)),-.27, expression("s"), xpd=TRUE, col="red", cex = 1)
    segments(hyp/(sqrt(2)),0,hyp/(sqrt(2)),-.2,lty=1,lwd = 2, col = "red")
    text(hyp,-.27,expression(sqrt(2)*"s"),xpd = TRUE, col = "red",cex = 1)
  }

  # Calculate the cumulative sum to draw step 3 and 4
  pd<-sort(pd)
  n<-length(pd)
  f<-rep(1/n,n)
  F<-cumsum(f)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  dev.new(width=6,height=8)
  par(mfrow=c(4,1),
      mai = c(.15,2,0,2))

  # Step 4
  plot(pd,F,type="s", xlab="", ylab="", las=1, frame.plot = FALSE,
       xaxs="r",yaxs="r",xlim=c(0,max(pd)),ylim=c(0,hlim), xaxt = "n", yaxt = "n")
  segments(0,1,max(pd)+2,1,lty=2,xpd=T)
  segments(max(pd),1, max(pd)+1,1,lty=1,xpd=T)
  axis(1, pos=0, at=c(seq(0,max(pd),1)))
  arrows(0, 0, max(pd)+2, 0, code=2, length=.07, xpd=TRUE)
  text(max(pd)+3,0, expression("PD"), xpd=TRUE, cex = 1.2)
  axis(2, pos=0, at=c(round(seq(0,1,1/5),digits = 1)), las=1)
  arrows(0,-1/5,0,1.4, code=2,length=.07,xpd=TRUE)
  text(0,1.6,expression("CDF\nof PD"),xpd=TRUE, cex=1.1)
  yinit <- 0
  prev <- 0
  for (i in pd) {

    rprism(i,prev,yinit,1/n,r)
    yinit <- yinit + 1/n
    prev = i
  }
  spd <- 0
  for(i in pd) {
    spd <- spd + i * i
  }
  mspd <- spd/n
  rmspd<-sqrt(mspd)
  rmspdprism(rmspd,r)

  # Step 3
  par(mai=c(0,2,0,2))
  plot(pd,F,type="s", xlab="", ylab="",frame.plot = FALSE, las=1,
       xaxs="r",yaxs="r",xlim=c(0,max(pd)),ylim=c(0,hlim),xaxt="n", yaxt="n")
  axis(1, pos=0, at=c(seq(0,max(pd),1)))
  arrows(0, 0, max(pd)+2, 0, code=2, length=.07, xpd=TRUE)
  text(max(pd)+3,0, expression("PD"), xpd=TRUE, cex = 1.2)
  axis(2, pos=0, at=c(round(seq(0,1,1/5),digits = 1)), las=1)
  arrows(0,-1/5,0,1.4, code=2,length=.07,xpd=TRUE)
  text(0,1.6,expression("CDF\nof PD"),xpd=TRUE, cex=1.1)
  #
  segments(0,1,max(pd)+2,1,lty=2,xpd=T)
  segments(max(pd),1, max(pd)+1,1,lty=1,xpd=T)
  segments(min(pd),F[1],min(pd),0,lty = 1)
  #
  m<-mean(pd)

  for(i in c(1:length(pd))){
    if(pd[i] <= m){
      for (j in c(seq(pd[i],m,.02))) {
        segments(j,F[i],j,0,lty = 1,col = "red")
      }
      segments(pd[i],F[i],pd[i],hlim,lty = 2, col = "gray")
    }
    else {
      for (j in c(seq(m,pd[i],.02))) {
        segments(j,F[i-1],j,1,lty = 1,col = "blue")
      }
      segments(pd[i],F[i],pd[i],0,lty = 2, col = "gray")
      segments(pd[i],1,pd[i],hlim,lty = 2, col = "gray")
    }
  }
  for(i in c(1:(length(pd)-1))){
    if(pd[i] < pd[i+1]){
      if(pd[i] <= m){
        segments(pd[i],F[i],0,F[i],lty = 2, col = "gray")
      }
      else {
        segments(m,F[i],0,F[i],lty = 2, col = "gray")
      }
    }
  }
  segments(m,0,m,1,lty = 1, lwd = 2,xpd=T)
  text(m,1.3,expression("MPD"),xpd=TRUE, cex=1.1)
  # Step 2
  par(mai=c(0,2,1,2))
  stripchart(pd,method="stack",offset=.2, pch=19,frame.plot = FALSE,xaxt="n",cex = 1.4,xlim = c(0,max(pd)))
  axis(1, pos=.87, at=c(seq(0,max(pd),1)))
  arrows(-1,.87, max(pd)+2,.87, code=2, length=.07, xpd=TRUE)
  text(max(pd)+3,.87, expression("PD"), xpd=TRUE, cex = 1.2)

  # Step 1
  par(mai=c(1,2,0,2))
  stripchart(data,method="stack",offset=.2, pch=19,frame.plot = FALSE,xaxt="n", cex = 1.4)
  axis(1, pos=.87, at=c(seq(min(data),max(data),1)))
  arrows(min(data)-1, .87, max(data)+2,.87, code=2, length=.07, xpd=TRUE)
  text(max(data)+3,.87, expression("x"), xpd=TRUE, cex = 1.2)
}
