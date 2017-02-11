library(maps)
library(geosphere)

GeneratePointsUS <- function(N=10000){
    #generates N points inside the continental US
    pts = matrix(rep(0,2*N),ncol=2)
    for(i in 1:N){
        repeat{
            pts[i,] = c(runif(1,-124.73,-66.95), runif(1,25.11,49.38))
            if(!is.na(map.where('state',pts[i,1],pts[i,2]))) break
        }
    }
    return(pts)
}

GeneratePointsNY <- function(N=10000){
    #generates N points inside the continental US
    pts = matrix(rep(0,2*N),ncol=2)
    for(i in 1:N){
        repeat{
            pts[i,] = c(runif(1,-81,-71), runif(1,40.02,45.3))
            statename = strsplit(map.where('state',pts[i,1],pts[i,2]),':')[[1]][[1]]
            if(!is.na(statename)){
                if(strsplit(map.where('state',pts[i,1],pts[i,2]),':')[[1]][[1]] == 'new york') break
            }
        }
    }
    return(pts)
}

DistanceGrid <- function(rpts,cpts){
    #compute Distance between rpts and cpts and return matrix rXc
    N = dim(rpts)[[1]]
    M = dim(cpts)[[1]]
    dists = matrix(rep(0,N*M),ncol=M)
    for(i in 1:N){ 
        dists[i,] = distGeo(cpts,rpts[i,])
    }
    return(dists)
}

MapAll <- function(pts1,pts2){
    map('state')
    points(pts1,col='gray')
    points(pts2,pch=20,cex=.5)
}

PlotGridInfluence <- function(grid, stepsize=F, lowerlimit=F, upperlimit=F, mincount=50){
    #plot the number of sites associated with each point
    if(!upperlimit) upperlimit = max(grid)
    if(!lowerlimit) lowerlimit = 0
    if(!stepsize) stepsize = (upperlimit-lowerlimit)/100
    limits = seq(lowerlimit, upperlimit, stepsize)
    counts = c()
    mins = c()
    maxs = c()
    meds = c()
    for(i in limits){
        counts = append(counts,sum(grid<i)/1000)
        cvec = apply((grid<i),1,sum)
        mins = append(mins,min(cvec))
        maxs = append(maxs,max(cvec))
        meds = append(meds,median(cvec))
    }
    lowcut = min(which(mins>mincount))
    print(c(counts[[lowcut-1]], mins[[lowcut-1]], maxs[[lowcut-1]], meds[[lowcut-1]],limits[[lowcut-1]]))
    print(c(counts[[lowcut]], mins[[lowcut]], maxs[[lowcut]], meds[[lowcut]],limits[[lowcut]]))
    par(mfrow=c(2,2))
    plot(limits,counts, main = 'mean')
    #lines(c(lowcut,lowcut),c(0,2000))
    plot(limits,mins, main = 'min')
    #lines(c(lowcut,lowcut),c(0,2000))
    plot(limits,maxs, main = 'max')
    #lines(c(lowcut,lowcut),c(0,2000))
    plot(limits,meds, main = 'median')
    #lines(c(lowcut,lowcut),c(0,2000))
    return(limits[[lowcut]])
}

GridInfluence <- function(grid, stepsize=F, lowerlimit=F, upperlimit=F, mincount=50){
    #calculate the number of sites in each point 
    if(!upperlimit) upperlimit = max(grid)
    if(!lowerlimit) lowerlimit = 0
    if(!stepsize) stepsize = (upperlimit-lowerlimit)/100
    limits = seq(lowerlimit, upperlimit, stepsize)
    mins = c()
    for(i in limits){
        cvec = apply((grid<i),1,sum)
        mins = append(mins,min(cvec))
    }
    lowcut = min(which(mins>mincount))
    return(limits[[lowcut]])
}

findannualvals <- function(filename, folder, years, season=F){
    #aggrogate matrix of sitesXyear by seasons or annual if season=F
    if(!season) season=1:12
    data = read.table(paste(folder, filename, sep = '/'), 
                           row.names=NULL, fill=TRUE)
    vout = rep(0,length(years))
    for(y in 1:length(years)){
        vout[[y]] = mean(data[[3]][data[[1]] == years[y] & 
                                   is.element(data[[2]],season)])
    }
    vout = (vout-mean(vout))/sd(vout)
}

GenerateFlowTable <- function(folder,inv,season=F){
    #read all files in inv and store in a sitesXyear matrix
    filenames = inv$ID
    years = 1940:2014
    vals = sapply(filenames, findannualvals, folder=folder, years = years,
                   season=season)
}

CleanInv <- function(folder, inv){
    #remove any sites in inv which have NA anywhere in the record or are in HI 
    #or AK
    flowtab = GenerateFlowTable(folder, inv)
    flowsum = apply(flowtab,2,sum)
    return(inv[!is.na(flowsum) & !is.element(inv$REGION,c('HI', 'AK')),])
}

GenerateMedianTable <- function(flowtab,mask){
    #multiply matrices to get 
    #flowtab[years,sites]
    #mask[pts,sites]
    #output[pts,years]
    M = dim(mask)[[2]]
    N = dim(mask)[[1]]
    meds = matrix(rep(0,N*75),ncol=75)
    for(i in 1:N){
        K = sum(mask[i,])
        if(K==0){
            meds[i,]=0
            next
        }
        for(j in 1:75){
            if(K>1) meds[i,j] = median(flowtab[,mask[i,]][j,],na.rm=T)
            if(K==1) meds[i,j] = flowtab[,mask[i,]][[j]]
        }
    }
    return(meds)
}

FindCPsThorough <- function(meds,pval=.05){
    years = 1940:2014
    CPS = apply(meds,1,Lanzante,years=1940:2014, pcut=pval)
    cpst = matrix(,ncol=4,nrow=0)
    N=length(CPS)
    for(k in 1:N){
        cps = c()
        CP=CPS[[k]]
        
        if(CP[[1]]>0){
            for(i in 1:CP[[1]]){
                cps = append(cps,CP[i*3-1])
            }
            for(y in 1:length(years)){
                ml = ModLan(cps, meds[k,], years, years[[y]])
                if(!is.na(ml[[1]])){
                    if(ml[[1]]<pval){
                        cpst = rbind(cpst,c(years[[y]],ml,k))
                    }
                }
            }
        }
    }
    return(cpst)
}
    
FindCPs <- function(meds,pval=.05){
    CPS = apply(meds,1,Lanzante,years=1940:2014, pcut=pval)
    cps = matrix(,ncol=4,nrow=0)
    N=length(CPS)
    for(k in 1:N){
        CP=CPS[[k]]
        if(CP[[1]]>0){
            for(i in 1:CP[[1]]){
                cps = rbind(cps,c(CP[i*3+(-1:1)],k))
            }
        }
    }
    return(cps)
}

CountByYear <- function(cps,cutoff,window=0,N=10000, pval=.05){
    countup=rep(0,65)
    countdown=rep(0,65)
    for(y in 1945:2009){
        countup[[y-1944]] =sum(is.element(cps[,1],(y-window):(y+window)) &
                             cps[,3]>0)
        countdown[[y-1944]] =sum(is.element(cps[,1],(y-window):(y+window)) &
                             cps[,3]<0)
    }
    plot(1945:2009,countup,type='l',main=as.character(cutoff))
    lines(1945:2009,countdown,lty=2)
    lines(c(1940,2014),c(N*pval/2,N*pval/2))
    for(y in 1945:2009){
        lines(c(y,y),c(0,1000))
    }
}

vectordestination <- function(latlonpoint, dist, angle){
    Rearth <- 6372795
    Dd <- dist / Rearth
    Cc <- angle
    
    lata = latlonpoint[[1]]*pi/180
    lona = latlonpoint[[2]]*pi/180
    latb = asin(cos(Cc)*cos(lata)*sin(Dd)+sin(lata)*cos(Dd))
    dlon = atan2(cos(Dd)-sin(lata)*sin(latb),sin(Cc)*sin(Dd)*cos(lata))
    lonb = lona-dlon+pi/2
    lonb[lonb > pi]=lonb[lonb > pi]-2*pi
    lonb[lonb < -pi]=lonb[lonb < -pi]+2*pi
    latb <- latb*180/pi
    lonb <- lonb*180/pi
    cbind(lonb,latb)
}

plotcircle <- function(lat, lon, radius, dens, shadeangle){
    if(length(lat)>0){
        circlepoints = vectordestination(c(lat,lon),radius,seq(0,2*pi,pi/50))
        polygon(circlepoints[,1], circlepoints[,2], density = dens, 
            angle = shadeangle, border=F)
        #lines(circlepoints[,1], circlepoints[,2])
    }
}



plotperimeter <- function(cluster, pts, cutoff, year=F, angle=45, 
     colorpal = c('gray70','gray60','gray50','gray40','gray30','gray20','gray10')
           , wgt = 1,dottype=1, startyr=1930, filename = 'blank'){
    if(year) cluster = list(cluster[[1]][cluster[[1]]==year-1939],
                            cluster[[2]][cluster[[1]]==year-1939])
    latup = pts[,2][unique(cluster[[2]])]
    lonup = pts[,1][unique(cluster[[2]])]
    yr = cluster[[1]][[1]]+9
    ci = ceiling((yr+1929-startyr)*8/30)
    Nup = length(latup)
    circlepts = cbind(lonup, latup)
    perimpts = c()
    if(Nup){
        for(i in 1:Nup){ 
            temppts = vectordestination(circlepts[i,c(2,1)],
                        cutoff,seq(0,2*pi,pi/50))
            dgrid = DistanceGrid(circlepts,temppts)
            perimpts = rbind(perimpts,temppts[apply(dgrid,2,min)-cutoff > -1000,])
        }
        points(perimpts, pch=20, col='red',cex=wgt)
        write.table(perimpts,filename, row.names=F, col.names=F)
    }
}



plotperimeterfromsource <- function(filename, color){
    points(read.table(filename),col=color,pch=20)
}







plotcps <- function(cps, year, pts, cutoff, upangle=45, downangle=315){
    cpsshort = cps[cps[,1]==year,]
    cpsup = cpsshort[cpsshort[,3]>0,]
    cpsdn = cpsshort[cpsshort[,3]<0,]
    latup = pts[,2][cpsup[,4]]
    lonup = pts[,1][cpsup[,4]]
    latdn = pts[,2][cpsdn[,4]]
    londn = pts[,1][cpsdn[,4]]
    Nup = length(latup)
    Ndn = length(latdn)
    map('state')
    title(main = as.character(year))
    
    if(Nup){
        for(i in 1:Nup){ 
            plotcircle(latup[[i]],lonup[[i]],cutoff,20,upangle)
        }
    }
    if(Ndn){
        for(i in 1:Ndn){ 
            plotcircle(latdn[[i]],londn[[i]],cutoff,20,downangle)
        }
    }
}


plotcluster <- function(cluster, pts, cutoff, year=F, angle=45){
    if(year) cluster = list(cluster[[1]][cluster[[1]]==year-1939],
                            cluster[[2]][cluster[[1]]==year-1939])
    latup = pts[,2][unique(cluster[[2]])]
    lonup = pts[,1][unique(cluster[[2]])]
    print(length(unique(cluster[[2]])))
    Nup = length(latup)
    map('state')
    title(main = as.character(min(cluster[[1]]+1939)))
    
    if(Nup){
        for(i in 1:Nup){ 
            plotcircle(latup[[i]],lonup[[i]],cutoff,10,angle)
        }
    }
}

plotclusterny <- function(cluster, pts, cutoff, year=F, angle=45){
    if(year) cluster = list(cluster[[1]][cluster[[1]]==year-1939],
                            cluster[[2]][cluster[[1]]==year-1939])
    latup = pts[,2][unique(cluster[[2]])]
    lonup = pts[,1][unique(cluster[[2]])]
    print(length(unique(cluster[[2]])))
    Nup = length(latup)
    map('state','new york')
    title(main = as.character(min(cluster[[1]]+1939)))
    
    if(Nup){
        for(i in 1:Nup){ 
            plotcircle(latup[[i]],lonup[[i]],cutoff,10,angle)
        }
    }
}


cluster <- function(ks,mask,usemax=F,crit=250){
    Nk = length(ks)
    if(Nk==0) return(list(0))
    if(Nk==1) return(list(1,ks[[1]]))
    smask = mask[ks,ks]
    kluster=-Nk:-1
    nextid=1
    for(i in 1:Nk){
        k=i
        kcurrs = sort(unique(kluster[smask[k,]]))
        if(sum(kcurrs>0)==0){
            kluster[smask[k,]] = nextid
            nextid = nextid + 1
        } else{
            newval = kcurrs[length(kcurrs)]
            kluster[smask[k,]]=newval
            kluster[is.element(kluster,kcurrs[-length(kcurrs)])]=newval
        }
    }
    ktab = table(kluster)
    if(usemax) kmax = as.integer(names(ktab)[which.max(ktab)])
    if(!usemax) kmax = as.integer(names(ktab)[which(ktab>=crit)])
    outthing = list(rep(0,length(kmax)*2))
    if(length(kmax)>0){
        for(i in 1:length(kmax)){
            outthing[[i*2-1]] = unname(ktab[names(ktab)==kmax[[i]]])
            outthing[[i*2]] = ks[kluster==kmax[[i]]]
        }
    }
    return(outthing)
}


plotmatgrids <- function(upmat, dnmat,cutoffs){
    cutoffskm = cutoffs/1000
    dcut = (cutoffskm[[2]]-cutoffskm[[1]])
    plot(c(),xlim=c(1940,2015),ylim=c(min(cutoffskm),max(cutoffskm)))
    for(yi in 1:75){
        for(ci in 1:10){
            if(upmat[yi,ci]) polygon(yi+c(1939.5, 1939.5, 1940.5, 1940.5),
                             cutoffskm[[ci]]+c(-dcut,dcut,dcut,-dcut),
                             density=10,angle=45, border=F)    
            if(dnmat[yi,ci]) polygon(yi+c(1939.5, 1939.5, 1940.5, 1940.5),
                             cutoffskm[[ci]]+c(-dcut,dcut,dcut,-dcut),
                             density=10,angle=315, border=F)    
        }
    }
}

GenerateClustersByYear <- function(cps,pval,mask,N=10000){
    years = 1940:2014
    upcluster=list()
    dncluster=list()
    isup = rep(0,75)
    isdn = rep(0,75)
    for(yi in 1:75){
        upk = cps[,4][cps[,1]==years[[yi]] & cps[,3]>0]
        dnk = cps[,4][cps[,1]==years[[yi]] & cps[,3]<0]
        upcluster[[yi]] = cluster(upk,mask,crit=N*pval/2)
        dncluster[[yi]] = cluster(dnk,mask,crit=N*pval/2)
    }
    return(list(upcluster,dncluster))
}

MergeClusters <- function(incluster){
    clusterlist = list()
    upcurr=c()
    currclusts=c()
    nextclust=1
    for(yi in 1:75){
        N = as.integer(length(incluster[[yi]])/2)
        if(N>0){
            for(i in 1:as.integer(N)){
                oc=0
                thisclust = incluster[[yi]][[i*2]]
                for(ci in currclusts){
                    if(sum(is.element(thisclust,clusterlist[[ci]][[2]][clusterlist[[ci]][[1]]==(yi-1)]))){
                        if(oc>0){
                            clusterlist[[oc]][[1]] = append(clusterlist[[oc]][[1]], clusterlist[[ci]][[1]])
                            clusterlist[[oc]][[2]] = append(clusterlist[[oc]][[2]], clusterlist[[ci]][[2]])
                            clusterlist[[ci]] = list(c(),c())
                        }
                        if(oc==0){
                            clusterlist[[ci]][[1]] = append(clusterlist[[ci]][[1]], rep(yi,length(thisclust)))
                            clusterlist[[ci]][[2]] = append(clusterlist[[ci]][[2]], thisclust)
                            oc=ci
                            upcurr = append(upcurr,ci)
                        }
                    }
                }
                if(oc==0){
                    clusterlist[[nextclust]] = list(rep(yi,length(thisclust)),thisclust)
                    upcurr = append(upcurr,nextclust)
                    nextclust = nextclust+1
                }
            }
        }
        currclusts = unique(upcurr)
        upcurr=c()
    }
    return(clusterlist)
}

trimclusters <- function(clusterlist){
    newlist=list()
    for(i in 1:length(clusterlist)){
        if(length(clusterlist[[i]][[1]])>1){
            maxyr = as.integer(names(sort(-table(clusterlist[[i]][[1]]))))[[1]]
            keepk = clusterlist[[i]][[1]]==maxyr
            if(is.element(maxyr,6:70)){
                newlist[[length(newlist)+1]] = list( 
                                clusterlist[[i]][[1]][keepk],
                                clusterlist[[i]][[2]][keepk])
            }
        }
    }
    return(newlist)
}

buildclustermask <- function(clusters,mask){
    N = length(clusters)
    M = dim(mask)[[2]]
    newmask = matrix(rep(0,N*M),ncol=M)
    for(i in 1:N){
        newmask[i,] = colSums(mask[clusters[[i]][[2]],])
    }
    return(newmask&TRUE)
}

findoverlapyrs <- function(clusters1,clusters2){
    cfull = append(clusters1,clusters2)
    N1 = length(clusters1)
    N2 = length(clusters2)
    N = N1+N2
    yrs = vector('list',N)
    for(i in 1:N){
        yrs[[i]]=c(cfull[[i]][[1]][[1]])
        for(j in (1:N)[-i]){
            if(sum(is.element(cfull[[i]][[2]],cfull[[j]][[2]]))>0){
                yrs[[i]] = c(yrs[[i]],cfull[[j]][[1]][[1]])
            }
        }
        yrs[[i]]=yrs[[i]][-1]+1939
    }
    return(list(yrs[1:N1],yrs[(N1+1):N]))
}

pullyears <- function(clusters){
    yrs=c()
    for(ci in 1:length(clusters)){
        yrs[[ci]] = clusters[[ci]][[1]][[1]]    
    }
    return(yrs+1939)
}

plotmeds <- function(meds, yrs, overlapyrs,i){
    #dimz = dim(meds)
    #for(i in 1:dimz[[1]]){
        breaks = sort(c(overlapyrs[[i]],yrs[[i]],1940,2014))
        K=yrs[[i]]-1939
        plot(1940:2014,meds[i,],ylim=c(-4,4),type='l')
        for(l in 2:length(breaks)){
            j = breaks[[l-1]]-1939
            k = breaks[[l]]-1939
            m = median(meds[i,j:k])
            lines(c(j,k)+1939,c(m,m))
            medstring = as.character(round(m,2))
            text((j+k)/2+1939, 2.5, bquote(.(medstring)~sigma))
        }
        lines(c(K,K)+1939,c(-20,20))
    
}

plotvariable <- function(ann,DJF,MAM,JJA,SON,yrs,overlapyrs,i,varname){
    breaks = sort(c(overlapyrs[[i]],yrs[[i]],1940,2014))
    K=yrs[[i]]-1939
    seasoncolors=c('blue', 'forestgreen','firebrick', 'orange')
    plot(1940:2014,ann[i,],ylim=c(-2,4),type='l',lwd=2,ylab=expression(sigma),xlab='')
    lines(1940:2014,DJF[i,], col=seasoncolors[[1]])
    lines(1940:2014,MAM[i,], col=seasoncolors[[2]])
    lines(1940:2014,JJA[i,], col=seasoncolors[[3]])
    lines(1940:2014,SON[i,], col=seasoncolors[[4]])
    for(l in 2:length(breaks)){
        j = breaks[[l-1]]-1939
        k = breaks[[l]]-1939
        m = median(ann[i,j:k])
        m1 = median(DJF[i,j:k])
        m2 = median(MAM[i,j:k])
        m3 = median(JJA[i,j:k])
        m4 = median(SON[i,j:k])
        #lines(c(j,k)+1939,c(m1,m1))
        #lines(c(j,k)+1939,c(m2,m2))
        #lines(c(j,k)+1939,c(m3,m3))
        #lines(c(j,k)+1939,c(m4,m4))
        lines(c(j,k)+1939,c(m,m))
        annstring = as.character(format(round(m,2),nsmall=2))
        DJFstring = as.character(format(round(m1,2),nsmall=2))
        MAMstring = as.character(format(round(m2,2),nsmall=2))
        JJAstring = as.character(format(round(m3,2),nsmall=2))
        SONstring = as.character(format(round(m4,2),nsmall=2))
        text((j+k)/2+1941, 3.5, bquote(.(annstring)~sigma),pos=2)
        text((j+k)/2+1941, 3.25, bquote(.(DJFstring)~sigma),pos=2, col=seasoncolors[[1]])
        text((j+k)/2+1941, 3, bquote(.(MAMstring)~sigma),pos=2, col=seasoncolors[[2]])
        text((j+k)/2+1941, 2.75, bquote(.(JJAstring)~sigma),pos=2, col=seasoncolors[[3]])
        text((j+k)/2+1941, 2.5, bquote(.(SONstring)~sigma),pos=2, col=seasoncolors[[4]])
        
    }
    lines(c(K,K)+1939,c(-20,20))
}


plotvariableNY <- function(ann,DJF,MAM,JJA,SON,yrs,overlapyrs,i,varname){
    breaks = sort(c(overlapyrs[[i]],yrs[[i]],1940,2014))
    K=yrs[[i]]-1939
    seasoncolors=c('deepskyblue', 'springgreen','firebrick', 'orange')
    plot(1940:2014,ann[i,],ylim=c(-2,4),type='l',lwd=2,ylab=expression(sigma),xlab='', main = varname)
    lines(1940:2014,DJF[i,], col=seasoncolors[[1]])
    lines(1940:2014,MAM[i,], col=seasoncolors[[2]])
    lines(1940:2014,JJA[i,], col=seasoncolors[[3]])
    lines(1940:2014,SON[i,], col=seasoncolors[[4]])
    for(l in 2:length(breaks)){
        j = breaks[[l-1]]-1939
        k = breaks[[l]]-1939
        m = median(ann[i,j:k])
        m1 = median(DJF[i,j:k])
        m2 = median(MAM[i,j:k])
        m3 = median(JJA[i,j:k])
        m4 = median(SON[i,j:k])
        #lines(c(j,k)+1939,c(m1,m1))
        #lines(c(j,k)+1939,c(m2,m2))
        #lines(c(j,k)+1939,c(m3,m3))
        #lines(c(j,k)+1939,c(m4,m4))
        lines(c(j,k)+1939,c(m,m))
        annstring = as.character(format(round(m,2),nsmall=2))
        DJFstring = as.character(format(round(m1,2),nsmall=2))
        MAMstring = as.character(format(round(m2,2),nsmall=2))
        JJAstring = as.character(format(round(m3,2),nsmall=2))
        SONstring = as.character(format(round(m4,2),nsmall=2))
        text((j+k)/2+1941, 3.5, bquote(.(annstring)~sigma),pos=2)
        text((j+k)/2+1941, 3.25, bquote(.(DJFstring)~sigma),pos=2, col=seasoncolors[[1]])
        text((j+k)/2+1941, 3, bquote(.(MAMstring)~sigma),pos=2, col=seasoncolors[[2]])
        text((j+k)/2+1941, 2.75, bquote(.(JJAstring)~sigma),pos=2, col=seasoncolors[[3]])
        text((j+k)/2+1941, 2.5, bquote(.(SONstring)~sigma),pos=2, col=seasoncolors[[4]])
        
    }
}


plotannual <- function(ann,years,overlapyrs,varname){
    startyear = years[[1]]-1
    breaks = sort(c(overlapyrs,years[[1]],years[[length(years)]]))
    seasoncolors=c('deepskyblue', 'springgreen','firebrick', 'orange')
    plot(years,ann,ylim=c(-2,4),type='l',lwd=2,ylab=expression(sigma),xlab='', main = varname)
    for(l in 2:length(breaks)){
        j = breaks[[l-1]]-startyear
        k = breaks[[l]]-startyear
        m = median(ann[j:k])
        #lines(c(j,k)+1939,c(m1,m1))
        #lines(c(j,k)+1939,c(m2,m2))
        #lines(c(j,k)+1939,c(m3,m3))
        #lines(c(j,k)+1939,c(m4,m4))
        lines(c(j,k)+startyear,c(m,m))
        annstring = as.character(format(round(m,2),nsmall=2))
        text((j+k)/2+startyear+20, 3.5, bquote(.(annstring)~sigma),pos=2,cex=.5)
        
    }
    for(k in overlapyrs){ 
        K=k-startyear
        lines(c(K,K)+startyear,c(-20,20))
    }
}












