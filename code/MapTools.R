library(maps)
library(plotrix)

cpslines <- function(cps){
    if(cps!=F){
        for(cp in 2:(length(cps)-1)){
            lines(c(cps[[cp]],cps[[cp]]),c(-20,20))
        }
    }
}

plotcloud <- function(qs){    
    years=1940:2014
    polygon(c(years,rev(years)),c(qs[[1]],rev(qs[[5]])), col='lightgrey', border = NA)
    polygon(c(years,rev(years)),c(qs[[2]],rev(qs[[4]])), col='darkgrey', border = NA)
    
    lines(c(1920,2030),c(0,0),lty=2)
    lines(years,qs[[3]], lwd = 2)
}

plotpeak <- function(region,folder,inv,cps,varname='Peak',qtiles = c(.025,.16,.5,.84,.975)){
    filenames = inv$ID[inv$REGION==region]
    if(region=='All') filenames = inv$ID
    years=1940:2014
    vals = vector('list', length(years))
    for(i in 1:length(filenames)){
        data = read.table(paste(folder, filenames[[i]], sep = '/'), 
                           row.names=NULL, fill=TRUE)
        if(length(data)==4){
            data[[1]] = as.integer(data[[1]])+(data[[2]]>9)
            #vyears = unique(data[[1]])
            dm = mean(data[[4]])   
            dstd = sd(data[[4]])
            for(y in 1:length(years)){
                v=F
                v = (data[[4]][data[[1]] == years[y]]-dm)/dstd
                if(length(v)==1){
                    vals[[y]] = append(vals[[y]],v)
                } else{
                    vals[[y]] = append(vals[[y]],NA)
                }
            }
        }
    }
    qs = vector('list', 5)
    for(q in 1:5){
        qs[[q]] = rep(0, length(years))
        for(y in 1:length(years)){
            qs[[q]][[y]] = quantile(vals[[y]],qtiles[[q]], na.rm=T)
        }
    }
    plot(years,qs[[3]],ylim=c(-2,2), type='n', main = paste(region,varname), xlab=F, ylab=expression(sigma))
    monthcolors = c('gray22', 'deepskyblue', 'springgreen','firebrick', 'orange')
    plotcloud(qs)
    
    #polygon(c(years,rev(years)),c(qs[[1]][[1]],rev(qs[[1]][[5]])), col='lightgrey', border = NA)
    #polygon(c(years,rev(years)),c(qs[[1]][[2]],rev(qs[[1]][[4]])), col='darkgrey', border = NA)
    
    #lines(c(1920,2030),c(0,0),lty=2)
    #lines(years,qs[[1]][[3]], col=monthcolors[[1]], lwd = 2)
    if(cps!=F){
        for(i in 1:(length(cps)-1)){
            curmed = mean(qs[[3]][(cps[[i]]-1939):(cps[[i+1]]-1939)])
            medstring = as.character(round(curmed,2))
            text(mean(cps[[i]]:cps[[i+1]]), 2, bquote(.(medstring)~sigma))
            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
         }
    }
    cpslines(cps)
    return(cps)
}

plotregionvariable <- function(region, folder, inv, cps, varname, qtiles = c(.025,.16,.5,.84,.975)){
    filenames = inv$ID[inv$REGION==region]
    if(region=='All') filenames = inv$ID
    years=1940:2014
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    vals = vector('list', 5)
    for(s in 1:5){
        vals[[s]] = vector('list', length(years))
    }
    for(i in 1:length(filenames)){
        data = read.table(paste(folder, filenames[[i]], sep = '/'), 
                           row.names=NULL, fill=TRUE)
        #vyears = unique(data[[1]])
        dm = c( mean(data[[3]][is.element(data[[2]],okmonths[[1]])]),
                mean(data[[3]][is.element(data[[2]],okmonths[[2]])]),
                mean(data[[3]][is.element(data[[2]],okmonths[[3]])]),
                mean(data[[3]][is.element(data[[2]],okmonths[[4]])]),
                mean(data[[3]][is.element(data[[2]],okmonths[[5]])]))
        dstd = c(sd(data[[3]][is.element(data[[2]],okmonths[[1]])]),
                 sd(data[[3]][is.element(data[[2]],okmonths[[2]])]),
                 sd(data[[3]][is.element(data[[2]],okmonths[[3]])]),
                 sd(data[[3]][is.element(data[[2]],okmonths[[4]])]),
                 sd(data[[3]][is.element(data[[2]],okmonths[[5]])]))
        for(s in 1:5){
            for(y in 1:length(years)){
                v = (mean(data[[3]][data[[1]] == years[y] &
                                    is.element(data[[2]], okmonths[[s]])])-dm[[s]])/dstd[[s]]
                cnt = length(data[[3]][data[[1]] == years[y] &
                                    is.element(data[[2]], okmonths[[s]])])
                if(s==1){
                    okyr = cnt >= 8
                } else{
                    okyr = cnt==3
                }
                if(okyr){
                    vals[[s]][[y]] = append(vals[[s]][[y]], v)
                }
            }
        }
    }
    qs = vector('list', 5)
    for(s in 1:5){
        qs[[s]] = vector('list', 5)
        for(q in 1:5){
            qs[[s]][[q]] = rep(0, length(years))
            for(y in 1:length(years)){
                qs[[s]][[q]][[y]] = quantile(vals[[s]][[y]],qtiles[[q]])
            }
        }
    }
    plot(years,qs[[1]][[3]],ylim=c(-2,2), type='n', main = paste(region,varname), xlab=F, ylab=expression(sigma))
    monthcolors = c('gray22', 'deepskyblue', 'springgreen','firebrick', 'orange')
    plotcloud(qs[[1]])
    
    #polygon(c(years,rev(years)),c(qs[[1]][[1]],rev(qs[[1]][[5]])), col='lightgrey', border = NA)
    #polygon(c(years,rev(years)),c(qs[[1]][[2]],rev(qs[[1]][[4]])), col='darkgrey', border = NA)
    
    #lines(c(1920,2030),c(0,0),lty=2)
    #lines(years,qs[[1]][[3]], col=monthcolors[[1]], lwd = 2)
    if(cps==F){
        cpout = Lanzante(qs[[1]][[3]], years)
        cps=c()
        for(i in 1:cpout[[1]]){
            cps=append(cps,cpout[[i*3-1]])
        }
        cps=c(1940,cps,2014)
    }
    if(cps!=F){
        for(i in 1:(length(cps)-1)){
            curmed = mean(qs[[1]][[3]][(cps[[i]]-1939):(cps[[i+1]]-1939)])
            medstring = as.character(round(curmed,2))
            text(mean(cps[[i]]:cps[[i+1]]), 2, bquote(.(medstring)~sigma))
            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
         }
    }
    season=c('All','DJF','MAM','JJA','SON')
    cpslines(cps)
    
    for(k in 2:5){
        plot(years,qs[[1]][[3]],ylim=c(-2,2), type='n', main = paste(region,season[[k]],varname), xlab=F, ylab=expression(sigma))
        plotcloud(qs[[k]])
    #    lines(c(1920,2030),c(0,0),lty=2)
    #    lines(years,qs[[k]][[3]])
        if(cps!=F){
            for(i in 1:(length(cps)-1)){
                curmed = mean(qs[[k]][[3]][(cps[[i]]-1939):(cps[[i+1]]-1939)])
                medstring = as.character(round(curmed,2))
                text(mean(cps[[i]]:cps[[i+1]]), 2, bquote(.(medstring)~sigma))
                lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
             }
        }
        cpslines(cps)
    }
    return(cps)
}

teleseries <- function(Tel,years,tname,cps){
        plot(years,rep(0,length(years)),ylim=c(-3,3), type='n',main=tname, xlab=F, ylab=expression(sigma))
    Tel$V[Tel$V==-99.990]=-.001   
    plotunder(Tel$Y+(Tel$M-1)/11, Tel$V/sd(Tel$V),
          y0 =0)
    if(cps!=F){
        for(i in 1:(length(cps)-1)){
            curmed = mean(Tel$V[is.element(Tel$Y,cps[[i]]:cps[[i+1]])])
            medstring = as.character(round(curmed,2))
            text(mean(cps[[i]]:cps[[i+1]]), 3, bquote(.(medstring)~sigma))
            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
        }
    }
    cpslines(cps)
}

plotregionvals <- function(region, folder, inv, precipfolder=F, precipinv=F, peakfolder=F, peakinv=F, CP=F, signum=F, qtiles = c(.025,.16,.5,.84,.975), cps=FALSE, pdfname='generic.pdf',cpflow=F){
    pdf(pdfname, width = 5, height = 75)
    layout(matrix(1:(24+length(cps)),(24+length(cps))))
    par(bty='l',mai=c(.2,.5,.3,.1))
    years=1940:2014
    for(cp in cps){
        FLOW = QuestionCP('monthlyflow', inv, cpflow, cp)
        newmap(flowinv,FLOW[[1]][1,],FLOW[[2]][1,], 'Mean Annual Flow', .05)
    }
    cps=c(1940,cps,2014)
    cps = plotregionvariable(region, folder, inv, cps, varname='Flow', qtiles = qtiles)
    
    cps = plotpeak(region,peakfolder,peakinv,cps,varname='Peakflow')
    
    cps = plotregionvariable(region, precipfolder, precipinv, cps, varname='Precipitation', qtiles = qtiles)
#    
#    filenames = inv$ID[inv$REGION==region]
#    if(region=='All') filenames = inv$ID
#    years=1940:2015
#    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
#    vals = vector('list', 5)
#    for(s in 1:5){
#        vals[[s]] = vector('list', length(years))
#    }
#    for(i in 1:length(filenames)){
#        data = read.table(paste(folder, filenames[[i]], sep = '/'), 
#                           row.names=NULL, fill=TRUE)
#        #vyears = unique(data[[1]])
#        dm = c( mean(data[[3]][is.element(data[[2]],okmonths[[1]])]),
#                mean(data[[3]][is.element(data[[2]],okmonths[[2]])]),
#                mean(data[[3]][is.element(data[[2]],okmonths[[3]])]),
#                mean(data[[3]][is.element(data[[2]],okmonths[[4]])]),
#                mean(data[[3]][is.element(data[[2]],okmonths[[5]])]))
#        dstd = c(sd(data[[3]][is.element(data[[2]],okmonths[[1]])]),
#                 sd(data[[3]][is.element(data[[2]],okmonths[[2]])]),
#                 sd(data[[3]][is.element(data[[2]],okmonths[[3]])]),
#                 sd(data[[3]][is.element(data[[2]],okmonths[[4]])]),
#                 sd(data[[3]][is.element(data[[2]],okmonths[[5]])]))
#        for(s in 1:5){
#            for(y in 1:length(years)){
#                v = (mean(data[[3]][data[[1]] == years[y] &
#                                    is.element(data[[2]], okmonths[[s]])])-dm[[s]])/dstd[[s]]
#                cnt = length(data[[3]][data[[1]] == years[y] &
#                                    is.element(data[[2]], okmonths[[s]])])
#                if(s==1){
#                    okyr = cnt >= 8
#                } else{
#                    okyr = cnt==3
#                }
#                if(okyr){
#                    vals[[s]][[y]] = append(vals[[s]][[y]], v)
#                }
#            }
#        }
#    }
#    qs = vector('list', 5)
#    for(s in 1:5){
#        qs[[s]] = vector('list', 5)
#        for(q in 1:5){
#            qs[[s]][[q]] = rep(0, length(years))
#            for(y in 1:length(years)){
#                qs[[s]][[q]][[y]] = quantile(vals[[s]][[y]],qtiles[[q]])
#            }
#        }
#    }
#    plot(years,qs[[1]][[3]],ylim=c(-2,2), type='n', main = region, xlab=F, ylab=expression(sigma))
#    monthcolors = c('gray22', 'deepskyblue', 'springgreen','firebrick', 'orange')
#    plotcloud(qs[[1]])
#    
#    #polygon(c(years,rev(years)),c(qs[[1]][[1]],rev(qs[[1]][[5]])), col='lightgrey', border = NA)
#    #polygon(c(years,rev(years)),c(qs[[1]][[2]],rev(qs[[1]][[4]])), col='darkgrey', border = NA)
#    
#    #lines(c(1920,2030),c(0,0),lty=2)
#    #lines(years,qs[[1]][[3]], col=monthcolors[[1]], lwd = 2)
#    if(cps!=F){
#        for(i in 1:(length(cps)-1)){
#            curmed = mean(qs[[1]][[3]][(cps[[i]]-1939):(cps[[i+1]]-1939)])
#            medstring = as.character(round(curmed,2))
#            text(mean(cps[[i]]:cps[[i+1]]), 2, bquote(.(medstring)~sigma))
#            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#         }
#    }
#    season=c('All','DJF','MAM','JJA','SON')
#    cpslines(cps)
#    
#    for(k in 2:5){
#        plot(years,qs[[1]][[3]],ylim=c(-2,2), type='n', main = paste(region,season[[k]]), xlab=F, ylab=expression(sigma))
#        plotcloud(qs[[k]])
#    #    lines(c(1920,2030),c(0,0),lty=2)
#    #    lines(years,qs[[k]][[3]])
#        if(cps!=F){
#            for(i in 1:(length(cps)-1)){
#                curmed = mean(qs[[k]][[3]][(cps[[i]]-1939):(cps[[i+1]]-1939)])
#                medstring = as.character(round(curmed,2))
#                text(mean(cps[[i]]:cps[[i+1]]), 2, bquote(.(medstring)~sigma))
#                lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#             }
#        }
#        cpslines(cps)
#    }
    teleseries(AMO,years,'AMO',cps)
    teleseries(ENSO,years,'ENSO',cps)
    teleseries(NAO,years,'NAO',cps)
    teleseries(PNA,years,'PNA',cps)
    teleseries(AO,years,'AO',cps)
    teleseries(EA,years,'EA',cps)
    teleseries(WP,years,'WP',cps)
    teleseries(EPNP,years,'EPNP',cps)
    teleseries(EAWR,years,'EAWR',cps)
    teleseries(SCA,years,'SCA',cps)
    teleseries(TNH,years,'TNH',cps)
    teleseries(POL,years,'POL',cps)
    teleseries(PT,years,'PT',cps)

#    #AMO
#    plot(years,rep(0,length(years)),ylim=c(-3,3), type='n',main='AMO', xlab=F, ylab=expression(sigma))
#    AMO$V[AMO$V==-99.990]=-.001   
#    plotunder(AMO$Y+(AMO$M-1)/11, AMO$V/sd(AMO$V),
#          y0 =0)
#    if(cps!=F){
#        for(i in 1:(length(cps)-1)){
#            curmed = mean(AMO$V[is.element(AMO$Y,cps[[i]]:cps[[i+1]])])
#            medstring = as.character(round(curmed,2))
#            text(mean(cps[[i]]:cps[[i+1]]), 3, bquote(.(medstring)~sigma))
#            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#        }
#    }
#    cpslines(cps)
#     
#    #ENSO
#    plot(years,rep(0,length(years)),ylim=c(-3,3), type='n',main='ENSO', xlab=F, ylab=expression(sigma))
#    ENSO$V[ENSO$V==-99.990]=-.001   
#    plotunder(ENSO$Y+(ENSO$M-1)/11, ENSO$V/sd(ENSO$V),
#          y0 =0)
#    if(cps!=F){
#        for(i in 1:(length(cps)-1)){
#            curmed = mean(ENSO$V[is.element(ENSO$Y,cps[[i]]:cps[[i+1]])])
#            medstring = as.character(round(curmed,2))
#            text(mean(cps[[i]]:cps[[i+1]]), 3, bquote(.(medstring)~sigma))
#            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#        }
#    }
#    cpslines(cps)
#     
#     
#     
#    #NAO
#    plot(years,rep(0,length(years)),ylim=c(-3,3), type='n',main='NAO', xlab=F, ylab=expression(sigma))
#    NAO$V[NAO$V==-99.990]=-.001   
#    plotunder(NAO$Y+(NAO$M-1)/11, NAO$V/sd(NAO$V),
#          y0 =0)
#    if(cps!=F){
#        for(i in 1:(length(cps)-1)){
#            curmed = mean(NAO$V[is.element(NAO$Y,cps[[i]]:cps[[i+1]])])
#            medstring = as.character(round(curmed,2))
#            text(mean(cps[[i]]:cps[[i+1]]), 3, bquote(.(medstring)~sigma))
#            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#        }
#    }
#    cpslines(cps)
#    
#    #PNA
#    plot(years,rep(0,length(years)),ylim=c(-3,3), type='n',main='PNA', xlab=F, ylab=expression(sigma))
#    PNA$V[PNA$V==-99.990]=-.001   
#    plotunder(PNA$Y+(PNA$M-1)/11, PNA$V/sd(PNA$V),
#          y0 =0)
#    if(cps!=F){
#        for(i in 1:(length(cps)-1)){
#            curmed = mean(PNA$V[is.element(PNA$Y,cps[[i]]:cps[[i+1]])])
#            medstring = as.character(round(curmed,2))
#            text(mean(cps[[i]]:cps[[i+1]]), 3, bquote(.(medstring)~sigma))
#            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#        }
#    }
#    cpslines(cps)
#    
#    #AO
#    plot(years,rep(0,length(years)),ylim=c(-3,3), type='n',main='AO', xlab=F, ylab=expression(sigma))
#    AO$V[AO$V==-99.990]=-.001   
#    plotunder(AO$Y+(AO$M-1)/11, AO$V/sd(AO$V),
#          y0 =0)
#    if(cps!=F){
#        for(i in 1:(length(cps)-1)){
#            curmed = mean(AO$V[is.element(AO$Y,cps[[i]]:cps[[i+1]])])
#            medstring = as.character(round(curmed,2))
#            text(mean(cps[[i]]:cps[[i+1]]), 3, bquote(.(medstring)~sigma))
#            lines(c(cps[[i]],cps[[i+1]]),c(curmed,curmed), lwd=2)
#        }
#    }
#    cpslines(cps)
     
    dev.off()
}



plotunder <- function(x,y,upcol='red',downcol='blue',y0 = 0){
    N = length(x)
    y1 = y
    y1[y1 < y0] = y0 
    y2 = y
    y2[y > y0] = y0
    polygon(c(x[[1]], x, x[[N]]), c(y0, y1, y0),  col=upcol)
    polygon(c(x[[1]], x, x[[N]]), c(y0, y2, y0),  col=downcol)
}

    
mapCPdecadal <- function(CP, title, pdfname, inv, frames = list(1940:1949, 
                         1950:1959, 1960:1969, 1970:1979, 1980:1989,
                         1990:1999, 2000:2010),ny=FALSE){
    pdf(pdfname)
    for(f in frames){
        if(ny){
            map('county','new york')
        }
        else{
            map('state')
        }
        title(main=paste(as.character(f[[1]]), '-', 
                         as.character(f[[length(f)]])))
        points(inv$LON, inv$LAT, pch =20)
        SUBCP = CP[is.element(CP$CP,f) & CP$P<=.05,]
        pch = rep(20,length(SUBCP$ID))
        pch[SUBCP$DM>0] <- 24
        pch[SUBCP$DM<0] <- 25
        color=rep('gray', length(SUBCP[,1]))
        color[SUBCP$DM>0] <- 'green' 
        color[SUBCP$DM<0] <- 'red' 
        points(SUBCP$LON, SUBCP$LAT, bg = color, pch = pch, cex = 1.5)
    }
    dev.off()
}

mapCPannual <- function(CP, title, pdfname, inv, frames = 1940:2009, ny=FALSE){
    pdf(pdfname)
    for(f in frames){
        if(ny){
            map('county','new york')
        }
        else{
            map('state')
        }
        title(main=f)
        points(inv$LON, inv$LAT, pch =20)
        SUBCP = CP[CP$CP==f & CP$P<=.05,]
        pch = rep(20,length(SUBCP$ID))
        pch[SUBCP$DM>0] <- 24
        pch[SUBCP$DM<0] <- 25
        color=rep('gray', length(SUBCP[,1]))
        color[SUBCP$DM>0] <- 'green' 
        color[SUBCP$DM<0] <- 'red' 
        points(SUBCP$LON, SUBCP$LAT, bg = color, pch = pch, cex = 1.5)
    }
    dev.off()
}

CPmap <- function(CP, name, inv, years, states){
    map('state',states)
    title(main=name)
    points(inv$LON,inv$LAT, pch = 20, cex = .5, col='gray')
    up = CP$DM>0 & is.element(CP$CP,years)
    dn = CP$DM<0 & is.element(CP$CP,years)
    points(CP$LON[up],CP$LAT[up], pch=24, bg = 'yellowgreen', col='gray27')
    points(CP$LON[dn],CP$LAT[dn], pch=25, bg = 'brown1', col='gray27')
}



LabelSample <- function(lat,lon,sample){
    for(i in 1:length(sample)){
        points(lon[sample[[i]]],lat[sample[[i]]],pch=as.character(i), font = 2)
    }
}

CPtimeseries <- function(index, CP, sourcefolder, ENSO, NAO, PNA, AO, AMO, years){
    width = (years[[1]]-10):(years[[length(years)]]+10)
    #number of telleconnections
    tc=5
    N = length(index)
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    monthcolors = c('gray22', 'deepskyblue', 'springgreen','firebrick', 'orange')
    plot(c(width[[1]],width[[length(width)]]), c(0,N+tc), type='n', axes=F, xlab ='')
    for(i in 0:(tc+N)) lines(c(1930,2016),c(i,i),col='gray')
    axis(1,width[c(T,F)],pos=0,col='gray')
    for(i in 1:N){ 
        data = read.table(paste(sourcefolder, CP$ID[index[i]], sep = '/'))
        siteyears = unique(data[[1]])
        vals = rep(0,length(siteyears))
#        for(s in 1:5){
        for(s in 1:1){
            for(y in 1:length(siteyears)){
                vals[y] = mean(data[[3]][data[[1]] == siteyears[y] &
                                    is.element(data[[2]], okmonths[[s]])])
            }
            vals = vals-median(vals)
            vals = vals/max(abs(vals))/2
            med1 = median(vals[siteyears<CP$CP[[index[[i]]]]])
            med2 = median(vals[siteyears>CP$CP[[index[[i]]]]])
            #lines(siteyears+.5,vals + N + 4.5 - i,col=monthcolors[[s]])
            lines(siteyears+.5,vals + N + tc+.5 - i,col=monthcolors[[s]])
            lines(c(1930,CP$CP[[index[[i]]]]),c(med1,med1)+N+tc+.5-i)
            lines(c(CP$CP[[index[[i]]]],2016),c(med2,med2)+N+tc+.5-i)
            lines(c(CP$CP[[index[[i]]]],CP$CP[[index[[i]]]]),c(N + tc+1 - i,N + tc - i))
        }
        boxed.labels(years[[1]]-10.67,N+tc+.78-i,labels=as.character(i), bg='white')
    }
    plotunder(AMO$Y+(AMO$M-1)/11,
          AMO$V/max(abs(AMO$V[is.element(AMO$Y,width)]))/2+.5+4,
          y0 = .5+4)
    plotunder(ENSO$Y+(ENSO$M-1)/11,
          ENSO$V/max(abs(ENSO$V[is.element(ENSO$Y,width)]))/2+.5+3,
          y0 = .5+3)
    plotunder(NAO$Y+(NAO$M-1)/11,
          NAO$V/max(abs(NAO$V[is.element(NAO$Y,width)]))/2+.5+2,
          y0 = .5+2)
    plotunder(PNA$Y+(PNA$M-1)/11,
          PNA$V/max(abs(PNA$V[is.element(PNA$Y,width)]))/2+.5+1,
          y0 = .5+1)
    plotunder(AO$Y+(AO$M-1)/11,
          AO$V/max(abs(AO$V[is.element(AO$Y,width)]))/2+.5,
          y0 = .5)
    boxed.labels(rep(years[[1]]-10.67,5)+c(.34,.59,.62,.88,.62),c(.78,1.78,2.78,3.78, 4.78),
                 labels = c('AO','PNA','NAO','ENSO','AMO'),
                 bg='white')
}

CPtimeseries <- function(index, CP, sourcefolder, ENSO, NAO, PNA, AO, AMO, years){
    width = (years[[1]]-20):(years[[length(years)]]+20)
    #number of telleconnections
    tc=5
    N = length(index)
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    monthcolors = c('gray22', 'deepskyblue', 'springgreen','firebrick', 'orange')
    plot(c(width[[1]],width[[length(width)]]), c(0,N+tc), type='n', axes=F, xlab ='')
    for(i in 0:(tc+N)) lines(c(1930,2016),c(i,i),col='gray')
    axis(1,width[c(T,F)],pos=0,col='gray')
    for(i in 1:N){ 
        data = read.table(paste(sourcefolder, CP$ID[index[i]], sep = '/'))
        siteyears = unique(data[[1]])
        vals = rep(0,length(siteyears))
#        for(s in 1:5){
        for(s in 1:1){
            for(y in 1:length(siteyears)){
                vals[y] = mean(data[[3]][data[[1]] == siteyears[y] &
                                    is.element(data[[2]], okmonths[[s]])])
            }
            vals = vals-median(vals)
            vals = vals/max(abs(vals))/2
            med1 = median(vals[siteyears<CP$CP[[index[[i]]]]])
            med2 = median(vals[siteyears>CP$CP[[index[[i]]]]])
            #lines(siteyears+.5,vals + N + 4.5 - i,col=monthcolors[[s]])
            lines(siteyears+.5,vals + N + tc+.5 - i,col=monthcolors[[s]])
            lines(c(1930,CP$CP[[index[[i]]]]),c(med1,med1)+N+tc+.5-i)
            lines(c(CP$CP[[index[[i]]]],2016),c(med2,med2)+N+tc+.5-i)
            lines(c(CP$CP[[index[[i]]]],CP$CP[[index[[i]]]]),c(N + tc+1 - i,N + tc - i))
        }
        boxed.labels(years[[1]]-10.67,N+tc+.78-i,labels=as.character(i), bg='white')
    }
    plotunder(AMO$Y+(AMO$M-1)/11,
          AMO$V/max(abs(AMO$V[is.element(AMO$Y,width)]))/2+.5+4,
          y0 = .5+4)
    plotunder(ENSO$Y+(ENSO$M-1)/11,
          ENSO$V/max(abs(ENSO$V[is.element(ENSO$Y,width)]))/2+.5+3,
          y0 = .5+3)
    plotunder(NAO$Y+(NAO$M-1)/11,
          NAO$V/max(abs(NAO$V[is.element(NAO$Y,width)]))/2+.5+2,
          y0 = .5+2)
    plotunder(PNA$Y+(PNA$M-1)/11,
          PNA$V/max(abs(PNA$V[is.element(PNA$Y,width)]))/2+.5+1,
          y0 = .5+1)
    plotunder(AO$Y+(AO$M-1)/11,
          AO$V/max(abs(AO$V[is.element(AO$Y,width)]))/2+.5,
          y0 = .5)
    boxed.labels(rep(years[[1]]-10.67,5)+c(.34,.59,.62,.88,.62),c(.78,1.78,2.78,3.78, 4.78),
                 labels = c('AO','PNA','NAO','ENSO','AMO'),
                 bg='white')
}


timeseries <- function(index, sourcefolder, ENSO, NAO, PNA, AO, AMO, year){
    width = (year-10):(year+10)
    #number of telleconnections
    tc=5
    N = length(index)
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    monthcolors = c('gray22', 'deepskyblue', 'springgreen','firebrick', 'orange')
    plot(c(width[[1]],width[[length(width)]]), c(0,N+tc), type='n', axes=F, xlab ='')
    for(i in 0:(tc+N)) lines(c(1930,2016),c(i,i),col='gray')
    axis(1,width[c(T,F)],pos=0,col='gray')
    for(i in 1:N){ 
        data = read.table(paste(sourcefolder, CP$ID[index[i]], sep = '/'))
        siteyears = unique(data[[1]])
        vals = rep(0,length(siteyears))
#        for(s in 1:5){
        for(s in 1:1){
            for(y in 1:length(siteyears)){
                vals[y] = mean(data[[3]][data[[1]] == siteyears[y] &
                                    is.element(data[[2]], okmonths[[s]])])
            }
            vals = vals-median(vals)
            vals = vals/max(abs(vals))/2
            med1 = median(vals[siteyears<CP$CP[[index[[i]]]]])
            med2 = median(vals[siteyears>CP$CP[[index[[i]]]]])
            #lines(siteyears+.5,vals + N + 4.5 - i,col=monthcolors[[s]])
            lines(siteyears+.5,vals + N + tc+.5 - i,col=monthcolors[[s]])
            lines(c(1930,CP$CP[[index[[i]]]]),c(med1,med1)+N+tc+.5-i)
            lines(c(CP$CP[[index[[i]]]],2016),c(med2,med2)+N+tc+.5-i)
            lines(c(CP$CP[[index[[i]]]],CP$CP[[index[[i]]]]),c(N + tc+1 - i,N + tc - i))
        }
        boxed.labels(years[[1]]-10.67,N+tc+.78-i,labels=as.character(i), bg='white')
    }
    plotunder(AMO$Y+(AMO$M-1)/11,
          AMO$V/max(abs(AMO$V[is.element(AMO$Y,width)]))/2+.5+4,
          y0 = .5+4)
    plotunder(ENSO$Y+(ENSO$M-1)/11,
          ENSO$V/max(abs(ENSO$V[is.element(ENSO$Y,width)]))/2+.5+3,
          y0 = .5+3)
    plotunder(NAO$Y+(NAO$M-1)/11,
          NAO$V/max(abs(NAO$V[is.element(NAO$Y,width)]))/2+.5+2,
          y0 = .5+2)
    plotunder(PNA$Y+(PNA$M-1)/11,
          PNA$V/max(abs(PNA$V[is.element(PNA$Y,width)]))/2+.5+1,
          y0 = .5+1)
    plotunder(AO$Y+(AO$M-1)/11,
          AO$V/max(abs(AO$V[is.element(AO$Y,width)]))/2+.5,
          y0 = .5)
    boxed.labels(rep(years[[1]]-10.67,5)+c(.34,.59,.62,.88,.62),c(.78,1.78,2.78,3.78, 4.78),
                 labels = c('AO','PNA','NAO','ENSO','AMO'),
                 bg='white')
}

mapregion <- function(name, year, inv, cpdir){
    shape = loadregion(name)
    map('state')
    polygon(shape$x,shape$y,density=0)
    pch=rep('.',length(inv[[1]]))
    points(inv$LON,inv$LAT, pch = 20, cex = .5, col='gray')
    up = cpdir==1
    dn = cpdir==-1
    points(inv$LON[up],inv$LAT[up], pch=24, bg = 'yellowgreen', col='gray27')
    points(inv$LON[dn],inv$LAT[dn], pch=25, bg = 'brown1', col='gray27')
}

newmap <- function(inv, p, dir, name, cut){
    map('state')
    title(main=name)
    points(inv$LON,inv$LAT, pch = 20, cex = .5, col='gray')
    dir[p > cut] = 0
    up = sign(dir)==1
    dn = sign(dir)==-1
    points(inv$LON[up],inv$LAT[up], pch=24, bg = 'yellowgreen', col='gray27')
    points(inv$LON[dn],inv$LAT[dn], pch=25, bg = 'brown1', col='gray27')
}

mapsuite <- function(FLOW, POT, PRECIP, flowinv, precipinv, cut, index){
    newmap(flowinv,FLOW[[1]][1,],FLOW[[2]][1,], 'Mean Annual Flow', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,FLOW[[1]][2,],FLOW[[2]][2,], 'Mean MAM Flow', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,FLOW[[1]][3,],FLOW[[2]][3,], 'Mean JJA Flow', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,FLOW[[1]][4,],FLOW[[2]][4,], 'Mean SON Flow', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,FLOW[[1]][4,],FLOW[[2]][4,], 'Mean DJF Flow', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    
    newmap(flowinv,POT[[1]][1,],POT[[2]][1,], 'Mean Annual POT', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,POT[[1]][2,],POT[[2]][2,], 'Mean MAM POT', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,POT[[1]][3,],POT[[2]][3,], 'Mean JJA POT', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,POT[[1]][4,],POT[[2]][4,], 'Mean SON POT', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(flowinv,POT[[1]][4,],POT[[2]][4,], 'Mean DJF POT', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    
    newmap(precipinv,PRECIP[[1]][1,],PRECIP[[2]][1,], 'Mean Annual Precip', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(precipinv,PRECIP[[1]][2,],PRECIP[[2]][2,], 'Mean MAM Precip', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(precipinv,PRECIP[[1]][3,],PRECIP[[2]][3,], 'Mean JJA Precip', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(precipinv,PRECIP[[1]][4,],PRECIP[[2]][4,], 'Mean SON Precip', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
    newmap(precipinv,PRECIP[[1]][4,],PRECIP[[2]][4,], 'Mean DJF Precip', cut)
    LabelSample(flowinv$LAT,flowinv$LON,index)
}

            
mapCPsuite <- function(CPflow, CPflowpot, CPprecip, pdfname, flowinv, prcpinv,
                       years=1970:1972, states=c()){
    prename= paste(as.character(years[[1]]),
                   as.character(years[[length(years)]]), sep = '-')
    pdf(pdfname)
    index = sample(which(is.element(MLflowCP$ALL$CP,years) & MLflowCP$ALL$LON > -128),5)
    CPmap(CPflow$ALL, paste(prename, 'Annual Flow'), flowinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPflow$MAM, paste(prename, 'Annual Flow MAM'), flowinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPflow$JJA, paste(prename, 'Annual Flow JJA'), flowinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPflow$SON, paste(prename, 'Annual Flow SON'), flowinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPflow$DJF, paste(prename, 'Annual Flow DJF'), flowinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPflowpot$ALL, paste(prename, 'Annual POT'), flowinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPprecip$ALL, paste(prename, 'Annual Precip'), precipinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPprecip$MAM, paste(prename, 'Annual Precip MAM'), precipinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPprecip$JJA, paste(prename, 'Annual Precip JJA'), precipinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPprecip$SON, paste(prename, 'Annual Precip SON'), precipinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPmap(CPprecip$DJF, paste(prename, 'Annual Precip DJF'), precipinv, years, states)
    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
    CPtimeseries(index, CPflow$ALL, 'monthlyflow', ENSO, NAO, PNA, AO, AMO, years)
    dev.off()
}


TeleMap <- function(teleob,inv,nameprefix=''){
    N=length(teleob[[1]][1,])
    suffixes = c('Annual', 'MAM', 'JJA', 'SON', 'DJF')
    for(i in 1:5){
        map('state')
        pchar = rep(20, N)
        psize = rep(.01, N)
        pchar[teleob[[2]][i,]<.05] = 1
        psize[teleob[[2]][i,]<.05] = abs(teleob[[1]][1,][teleob[[2]][i,]<.05]*10)
        points(inv$LON,inv$LAT,cex = psize,pch=pchar)
        title(main = paste(nameprefix,suffixes[[i]]))
    }
}
    
CPhist <- function(CP,region,inv,allseason=TRUE,signum=FALSE){
    cpsmoothup = vector('list',3)
    cpsmoothdn = vector('list',3)
    cpmax = vector('list',3)
    N=length(1945:2010)
    if(signum) signum = signum
    if(!signum) signum = sum(is.element(inv$REGION,region))*.025
    if(allseason) irng=1:5
    if(!allseason) irng=c(1)
    for(i in irng){
        cps = CP[[i]]$CP[is.element(CP[[i]]$REGION, region) & CP[[i]]$DM > 0]
        cnt = rep(0,N)
        for(y in 1945:2010){
            cnt[[y-1944]] = sum(cps==y)
        }        
        cpsmoothup[[i]] = (cnt[c(1,1:(N-1))]+cnt[c(1:N)]+cnt[c(2:N,N)])/3
        #cpsmoothup[[i]] = (cnt[c(1,1,1:(N-2))] + cnt[c(1,1:(N-1))] + cnt[c(1:(N))] + cnt[c(2:(N),N)] + cnt[c(3:N,N,N)])/5
        cps = CP[[i]]$CP[is.element(CP[[i]]$REGION,region) & CP[[i]]$DM < 0]
        cnt = rep(0,N)
        for(y in 1945:2010){
            cnt[[y-1944]] = sum(cps==y)
        }    
        cpsmoothdn[[i]] = (cnt[c(1,1:(N-1))]+cnt[c(1:N)]+cnt[c(2:N,N)])/3
        #cpsmoothdn[[i]] =  (cnt[c(1,1,1:(N-2))] + cnt[c(1,1:(N-1))] + cnt[c(1:(N))] + cnt[c(2:(N),N)] + cnt[c(3:N,N,N)])/5  
        cpmax[[i]] = cpsmoothup[[i]]
        dngtup = cpsmoothdn[[i]]>cpsmoothup[[i]]
        cpmax[[i]][dngtup] = -cpsmoothdn[[i]][dngtup]
    }
#    plot(1945:2010, cpsmoothup[[1]],type='l', ylim = c(-70,70), lwd= 2)
#    lines(1945:2010, -cpsmoothdn[[1]], lwd= 2)
    plot(1945:2010, cpmax[[1]],type='l', ylim = c(-70,70), lwd= 2)
    if(allseason){
#        lines(1945:2010, cpsmoothup[[2]],col='green')
#        lines(1945:2010, cpsmoothup[[3]],col='red')
#        lines(1945:2010, cpsmoothup[[4]],col='orange')
#        lines(1945:2010, cpsmoothup[[5]],col='blue')
#        lines(1945:2010, -cpsmoothdn[[2]],col='green')
#        lines(1945:2010, -cpsmoothdn[[3]],col='red')
#        lines(1945:2010, -cpsmoothdn[[4]],col='orange')
#        lines(1945:2010, -cpsmoothdn[[5]],col='blue')
        lines(1945:2010, cpmax[[2]],col='green')
        lines(1945:2010, cpmax[[3]],col='red')
        lines(1945:2010, cpmax[[4]],col='orange')
        lines(1945:2010, cpmax[[5]],col='blue')
    }
    lines(c(1940,2015), c(signum,signum),lty=2)
    lines(c(1940,2015), c(-signum,-signum),lty=2)
    lines(c(1940,2015), c(0,0), lwd= 2)
    #for(y in 1940:2015){
    #    lines(c(y,y),c(-80,80))
    #}
}
#mapCPsuite <- function(CPflow, CPflowpot, CPprecip, pdfname, flowinv, prcpinv,
#                       years=1970:1972, states=c()){
#    prename= paste(as.character(years[[1]]),
#                   as.character(years[[length(years)]]), sep = '-')
#    pdf(pdfname)
#    index = sample(which(is.element(MLflowCP$ALL$CP,years) & MLflowCP$ALL$LON > -128),5)
#    CPmap(CPflow$ALL, paste(prename, 'Annual Flow'), flowinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPflow$MAM, paste(prename, 'Annual Flow MAM'), flowinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPflow$JJA, paste(prename, 'Annual Flow JJA'), flowinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPflow$SON, paste(prename, 'Annual Flow SON'), flowinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPflow$DJF, paste(prename, 'Annual Flow DJF'), flowinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPflowpot$ALL, paste(prename, 'Annual POT'), flowinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPprecip$ALL, paste(prename, 'Annual Precip'), precipinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPprecip$MAM, paste(prename, 'Annual Precip MAM'), precipinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPprecip$JJA, paste(prename, 'Annual Precip JJA'), precipinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPprecip$SON, paste(prename, 'Annual Precip SON'), precipinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPmap(CPprecip$DJF, paste(prename, 'Annual Precip DJF'), precipinv, years, states)
#    LabelSample(CPflow$ALL$LAT,CPflow$ALL$LON,index)
#    CPtimeseries(index, CPflow$ALL, 'monthlyflow', ENSO, NAO, PNA, AO, AMO, years)
#    dev.off()
#}
