library(trend)

set.seed(8675309)
boots = vector('list',1000)
for(i in 1:1000){
    boots[[i]]=sample(1:75)
}

filterinv <- function(inv, newinvfile, folder){
    filenames= inv[[1]]
    N = length(filenames)
    ok=rep(T,N)
    for(i in 1:N){
        data = read.table(paste(folder, filenames[[i]], sep = '/'))
        for(m in 1:12){
            if(sum(data[[2]]==1) < 60){
                ok[[i]]=F
            }
        }
    }
    nin = data.frame(inv[ok,1:4])
    write.table(nin,newinvfile, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}


pettitt_it <- function(v,years){
    ptr = pettitt.test(v)
    dv = mean(v[ptr$estimate:length(v)]) - mean(v[1:ptr$estimate])
    return(c(1, years[[ptr$estimate]], ptr$p.value, dv))
}


Siegel_Castellan <- function(v,years){
    N = length(v)
    SR = rep(NA, N)
    r = rank(v)
    for(i in 1:N){
        SR[[i]] = sum(r[1:i])
    }
    SA = abs(2*SR-(1:N)*N+1)
    N1 = which.max(SA)
    W = SR[[N1]]
    N2 = N-N1
    Wcrit = N1*(N+1)/2
    sw = sqrt(N1*N2*(N+1)/12)
    delta = .5 * sign(Wcrit-W)
    z = (W-Wcrit+delta)/sw
    p = pnorm(-abs(z))*2
    dv = median(v[N1:length(v)]) - median(v[1:N1])
    return(c(1, years[[N1]], p, dv))
}

normalizemedian <- function(v, cps){
    cps=sort(cps)
    segments=list()
    v1=c()
    for(i in 1:length(cps)){
        if(length(segments)==0){
            segments[[i]]=1:cps[[i]]
        } else{ 
            segments[[i]]=cps[[i-1]]:cps[[i]]
        }
    }
    segments[[i+1]] = cps[[i]]:length(v)
    for(seg in segments){
        v1[seg] = v[seg]-median(v[seg])
    }
    return(v1)
}
        

Lanzante <- function(v, years, cp = c(0), pcut = .05, buffer = 5){
    #Iterative version of Siegel_Castellan repeats untils cp is not sig
    #at p_cut, is in the first or last 10 pts.
    #If within 5 of previous CP, a secondary maximum is used
    cptemp = Siegel_Castellan(v,years)
    cpcheck=FALSE
    if(length(cp)>1){
        cpcheck=sum(abs(cptemp[[2]] - cp[(1:cp[[1]])*3-1]) < buffer)
    }
    if(cptemp[[3]] > pcut |
       which(years == cptemp[[2]]) < buffer |
       which(years == cptemp[[2]]) > (length(v)-buffer) |
       cpcheck){
       return(cp)
    }
    else{
        cp[(1:3)+ 1 + cp[[1]]*3] = cptemp[2:4]
        cp[[1]] = cp[[1]]+1
        cpindex = which(is.element(years, cp[3*(1:cp[[1]])-1]))
        v1=normalizemedian(v, cpindex)
        return(Lanzante(v1, years, cp, pcut, buffer))
    }
}


##CP = c(N, cp1, pval1, dv1, cp2 ....)

whiteitup <- function(data){
    yrs = unique(data[[1]])
    d1 = data[[1]]
    d3 = data[[3]]
    mu = mean(d3)
    vals = c()
    for(y in yrs){
        vals = append(vals,mean(d3[d1 == y]))
    }
    autoc = cor(vals,c(0,vals[1:(length(vals)-1)]))
    for(y in 2:length(yrs)){
        data[[3]][d1 == yrs[[y]]] = data[[3]][d1 == yrs[[y]]] - vals[[y-1]]*autoc
    }
    return(data)
}

generateallCP <- function(folder, inv, ally = 1940:2014, test_func = Lanzante,boots=list(c(1:75)),seasons=1:5, outfile=FALSE, iinit=1, ifin=FALSE, binit=1, normout= FALSE, prewhiten=FALSE){
    filenames= inv[[1]]
    lat = inv[[2]]
    lon = inv[[3]]
    regions = inv$REGION
    N = length(filenames)
    if(ifin==FALSE) ifin=N
    CP = list(ALL = data.frame(ID=character(), LAT = numeric(), LON = numeric(), 
                         CP = integer(), P = numeric(), DM = numeric(), 
                         REGION = character(), BOOT = numeric(), stringsAsFactors = FALSE),
              MAM = data.frame(ID=character(), LAT = numeric(), LON = numeric(), 
                         CP = integer(), P = numeric(), DM = numeric(),  
                         REGION = character(), BOOT = numeric(), stringsAsFactors = FALSE),
              JJA = data.frame(ID=character(), LAT = numeric(), LON = numeric(), 
                         CP = integer(), P = numeric(), DM = numeric(),  
                         REGION = character(), BOOT = numeric(), stringsAsFactors = FALSE),
              SON = data.frame(ID=character(), LAT = numeric(), LON = numeric(), 
                         CP = integer(), P = numeric(), DM = numeric(),  
                         REGION = character(), BOOT = numeric(), stringsAsFactors = FALSE),
              DJF = data.frame(ID=character(), LAT = numeric(), LON = numeric(), 
                         CP = integer(), P = numeric(), DM = numeric(),  
                         REGION = character(), BOOT = numeric(), stringsAsFactors = FALSE))
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    for(i in iinit:ifin){
        print((i-iinit)/(ifin-iinit))
        data = read.table(paste(folder, filenames[[i]], sep = '/'))
        if(prewhiten){
            data = whiteitup(data)
        }
        years = unique(data[[1]])
        
        vals = rep(0,length(years))
        count = rep(0,length(years))
        for(b in binit:length(boots)){     
            kyear = ally[boots[[b]]]
            kyear = kyear[is.element(kyear,ally)]
            for(s in seasons){
                for(y in 1:length(kyear)){
                    vals[y] = mean(data[[3]][data[[1]] == kyear[y] &
                                        is.element(data[[2]], okmonths[[s]])])
                    count[y] = length(data[[3]][data[[1]] == kyear[y] &
                                        is.element(data[[2]], okmonths[[s]])])
                }
                if(s==1){
                    okyr = (count >= 8 & !is.na(vals))
                } else{
                    okyr = (count==3 & !is.na(vals))
                }
                y_all = is.element(kyear,ally)
                if(sum(y_all & okyr)<(length(ally)-10)){
                    #print(filenames[[i]])
                    blankvar=1
                }else{
                    test_result = test_func(vals[okyr & y_all],kyear[okyr & y_all])
                    if(test_result[[1]]>=1){
                        lims=c(1)
                        for(j in 1:test_result[[1]]){
                            lims[j+1] = as.integer(test_result[(-1)+3*j])-1939
                            CP[[s]] = rbind(CP[[s]], 
                                           data.frame(ID = filenames[[i]],
                                             LAT = lat[[i]], 
                                             LON = lon[[i]], 
                                             CP = as.integer(test_result[(-1)+3*j]), 
                                             P = as.numeric(test_result[3*j]), 
                                             DM = as.numeric(test_result[1+3*j]), 
                                             REGION = regions[[i]], BOOT = b,
                                             stringsAsFactors=FALSE))
                            if(outfile!=FALSE){
                                write.table(matrix(c(filenames[[i]],
                                           as.integer(test_result[(-1)+3*j]),
                                           as.numeric(test_result[3*j]),
                                           as.numeric(test_result[1+3*j]),
                                           regions[[i]],b), nrow=1),
                                             outfile,sep='\t',quote=F,
                                             row.names=F, col.names=F,
                                             append=T)
                            }
                        }
                        lims = append(lims,2015-1939)
                        lims=lims[order(lims)]
                        #print(lims)
                        v2=c()
                        if(normout & s==1){
                            for(k in 1:(length(lims)-1)){
                                vtemp=vals[(lims[[k]]):(lims[[k+1]]-1)]  
                                v2 = append(v2,vtemp/median(vtemp))
                            }
                            write.table(cbind(ally,v2), paste('misc/normout',filenames[[i]],sep='/'), quote = F, row.names=F, col.names=F)                    
                        }
                    }
                    
                }

            }
        }
        binit=1
    }
    return(CP)
}

generateannualCP <- function(folder, inv, ally = 1940:2014, test_func = Lanzante){
    rmvlist = c()
    filenames= inv[[1]]
    lat = inv[[2]]
    lon = inv[[3]]
    regions = inv$REGION
    N = length(filenames)
    CP = list(ALL = data.frame(ID = character(), LAT = numeric(), 
                         LON = numeric(), 
                         CP = integer(), P = numeric(), DM = numeric(), 
                         REGION = character(), 
                         stringsAsFactors = FALSE))
    for(i in 1:length(filenames)){
        data = read.table(paste(folder, filenames[[i]], sep = '/'), 
                           row.names=NULL, fill=TRUE)
        data[[1]] = as.integer(data[[1]]) + as.integer(data[[2]] > 9)                
        years = unique(data[[1]])
        vals = rep(0,length(years))
        count = rep(0,length(years))
        for(y in 1:length(years)){
            vals[y] = mean(data[[3]][data[[1]] == years[y]])
            count[y] = length(data[[3]][data[[1]] == years[y]])
        }
        
        y_all = is.element(years,ally)
        if(sum(y_all)<60){
            rmvlist = append(rmvlist,filenames[[i]])
            print(filenames[[i]])
        }else{
            test_result = test_func(vals[y_all],years[y_all])
            if(test_result[[1]]>=1){
                for(j in 1:test_result[[1]]){
                    if(is.na(test_result[1+3*j])){break}
                    print(test_result[1+3*j])
                    CP[[1]] = rbind(CP[[1]], 
                                   data.frame(ID = filenames[[i]],
                                    LAT = lat[[i]], 
                                    LON = lon[[i]], 
                                    CP = as.integer(test_result[(-1)+3*j]), 
                                    P = as.numeric(test_result[3*j]), 
                                    DM = as.numeric(test_result[1+3*j]),  
                                    REGION = regions[[i]],
                                    stringsAsFactors=FALSE))
                    }
                }
            }
        }
    if(length(rmvlist)>0){
        write.table(rmvlist,'rmvlist',row.names=F,col.names=F,quote=F,append=TRUE)
    }
    return(CP)
}

ModLan <- function(cps, v, years, year){
    N=length(cps)
    if(N == 0){
        return(c(.5,0))
    }
    cps = cps[!(abs(cps-year) <= 5)]
    cpindex = which(is.element(years,cps))
    if(length(cpindex>0)){
        v = normalizemedian(v, cpindex)
    }
    N = length(v)
    SR = rep(NA, N)
    r = rank(v)
    for(i in 1:N){
        SR[[i]] = sum(r[1:i])
    }
    SA = abs(2*SR-(1:N)*N+1)
    N1 = which(years==year)
    if(length(N1)==0){
        print(years)
        print(year)
        return(c(.5,0))
    }
    W = SR[[N1]]
    N2 = N-N1
    Wcrit = N1*(N+1)/2
    sw = sqrt(N1*N2*(N+1)/12)
    delta = .5 * sign(Wcrit-W)
    z = (W-Wcrit+delta)/sw
    p = pnorm(-abs(z))
    dv = median(v[N1:length(v)]) - median(v[1:N1])
    return(c(p,dv))
}

QuestionCP <- function(folder, inv, CP, year, test_func = ModLan, ally = 1940:2014){
    filenames= inv[[1]]
    lat = inv[[2]]
    lon = inv[[3]]
    N = length(filenames)
    yrp = matrix(data=.5, nrow=5, ncol=length(filenames))
    yrdm = matrix(data=.5, nrow=5, ncol=length(filenames))
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    for(i in 1:length(filenames)){
        data = read.table(paste(folder, filenames[[i]], sep = '/'))
        years = unique(data[[1]])
        vals = rep(0,length(years))
        count = rep(0,length(years))
        for(s in 1:5){
            for(y in 1:length(years)){
                vals[y] = mean(data[[3]][data[[1]] == years[y] &
                                    is.element(data[[2]], okmonths[[s]])])
                count[y] = length(data[[3]][data[[1]] == years[y] &
                                    is.element(data[[2]], okmonths[[s]])])
            }
            if(s==1){
                okyr = count >= 8
            } else{
                okyr = count==3
            }
            y_all = is.element(years,ally)
            if(sum(y_all & okyr)<(length(ally)-10)){
                print(filenames[[i]])
            }else{
                cps = CP[[s]]$CP[which(CP[[s]]$ID==filenames[[i]])]
                test_result = test_func(cps, vals[okyr & y_all], years[okyr & y_all], year)
                yrp[s,i]=test_result[[1]]
                yrdm[s,i]=test_result[[2]]
            }    
        }
    }
    return(list(yrp,yrdm))
}

runavg <-function(v,winrad=1){
    v2=v
    for(i in 1:length(v)){
        if(i<=winrad){
            v2[[i]]=mean(v[1:(i+winrad)])
        } else{
            if(i>=(length(v)-winrad+1)){
                v2[[i]] = mean(v[(i-winrad):length(v)])
            } else{
                v2[[i]] = mean(v[(i-winrad):(i+winrad)])
            }
        }
    }
    return(v2)
}



inputregion <- function(CPflow,flowinv,years=1970:1972, states=c()){
    CPmap(MLflowCP$ALL, 'Annual Flow', flowinv, years, states)
    shape = locator(n=20,type='o')
    polygon(shape$x,shape$y,density=0)
    return(shape)
}

saveregion <- function(shp,name,years){
    write.table(matrix(c(name, years[[1]], years[[length(years)]]), 
                       ncol=3),
                'regions/key', row.names = F, col.names = F, append = T,
                quote = F)
    write.table(cbind(shp$x, shp$y), paste('regions',name, sep = '/'),
                quote = F, row.names = F, col.names=F)
}

loadregion <- function(name){
    data = read.table(paste('regions',name,sep='/'))
    names(data) <- c('x','y')
    return(data)
}

regionyears <- function(name){
    data = read.table('regions/key')
    index = which(data[[1]] == name)
    return(data[[2]][[index]]:data[[3]][[index]])
}

inshp <- function(shp, lat, lon){
    Wpts = shp$x <= lon
    Spts = shp$y <= lat
    N = length(Wpts)
    SNlines = xor(Spts,c(Spts[2:N],Spts[[1]]))
    EWlines = xor(Wpts,c(Wpts[2:N],Wpts[[1]])) & SNlines
    Wlines = Wpts & c(Wpts[2:N],Wpts[[1]]) & SNlines
    if(sum(EWlines)){
        shp$x1 = c(shp$x[2:N],shp$x[[1]])
        shp$y1 = c(shp$y[2:N],shp$y[[1]])
        x0 = (shp$x1-shp$x)/(shp$y1-shp$y)*(lat-shp$y)+shp$x
        Wlines[SNlines & EWlines & x0 <= lon] = TRUE
    }
    return(as.logical(sum(Wlines)%%2))
}

countin <- function(shp, obj){
    N = dim(obj)[[1]]
    which_in = rep(F,N)    
    for(i in 1:N){
        which_in[[i]] = inshp(shp,obj$LAT[[i]],obj$LON[[i]])
    }
    return(sum(which_in))
}

regionstats <- function(name, CPflow, CPflowpot, CPprecip, flowinv, precipinv, year, cut = .025, plotsuite=TRUE, setindex = c(F,F,F,F,F)){
    # returns number of significant cps in flow/by season, pot, and prec/byseason
    shp = loadregion(name)
    TOTflow = countin(shp,flowinv)
    TOTprecip = countin(shp,precipinv)
    
    FLOW = QuestionCP('monthlyflow', flowinv, CPflow, year)
    POT = QuestionCP('monthlypot', flowinv, CPflowpot, year)
    PRECIP = QuestionCP('monthlyprecip', precipinv, CPprecip, year)
    
    flowdir = FLOW[[2]]
    flowdir[FLOW[[1]]>cut]=0
    mapregion(name,year,flowinv,sign(flowdir[1,]))
    
    CNT_ALL_flow = countin(shp, flowinv[FLOW[[1]][1,] < cut, 1:3])
    CNT_ALL_flowpot = countin(shp, flowinv[POT[[1]][1,] < cut, 1:3])
    CNT_ALL_precip = countin(shp, flowinv[PRECIP[[1]][1,] < cut, 1:3])

    CNT_MAM_flow = countin(shp, flowinv[FLOW[[1]][3,] < cut, 1:3])
    CNT_MAM_flowpot = countin(shp, flowinv[POT[[1]][3,] < cut, 1:3])
    CNT_MAM_precip = countin(shp, flowinv[PRECIP[[1]][3,] < cut, 1:3])
    
    CNT_JJA_flow = countin(shp, flowinv[FLOW[[1]][4,] < cut, 1:3])
    CNT_JJA_flowpot = countin(shp, flowinv[POT[[1]][4,] < cut, 1:3])
    CNT_JJA_precip = countin(shp, flowinv[PRECIP[[1]][4,] < cut, 1:3])
    
    CNT_SON_flow = countin(shp, flowinv[FLOW[[1]][5,] < cut, 1:3])
    CNT_SON_flowpot = countin(shp, flowinv[POT[[1]][5,] < cut, 1:3])
    CNT_SON_precip = countin(shp, flowinv[PRECIP[[1]][5,] < cut, 1:3])
    
    CNT_DJF_flow = countin(shp, flowinv[FLOW[[1]][2,] < cut, 1:3])
    CNT_DJF_flowpot = countin(shp, flowinv[POT[[1]][2,] < cut, 1:3])
    CNT_DJF_precip = countin(shp, flowinv[PRECIP[[1]][2,] < cut, 1:3])
    if(plotsuite){
        mapdev = dev.cur()
        pdf(paste(as.character(year),'suite.pdf',sep=''))
        index = sample(which(FLOW[[1]][1,] < .025 & flowinv$LON > -128),5)
        index[setindex&TRUE]=setindex[setindex&TRUE]
        print(index)
        mapsuite(FLOW, POT, PRECIP, flowinv, precipinv, cut, index)
        CPtimeseries(index, CPflow$ALL, 'monthlyflow', ENSO, NAO, PNA, AO, AMO, year)
        dev.off()
        dev.set(mapdev)
    }
    outobj = c(CNT_ALL_flow, CNT_ALL_flowpot, CNT_ALL_precip, 
               CNT_MAM_flow, CNT_MAM_flowpot, CNT_MAM_precip,
               CNT_JJA_flow, CNT_JJA_flowpot, CNT_JJA_precip,
               CNT_SON_flow, CNT_SON_flowpot, CNT_SON_precip,
               CNT_DJF_flow, CNT_DJF_flowpot, CNT_DJF_precip,
               TOTflow, TOTprecip, year)
    names(outobj) = c('CNT_ALL_flow', 'CNT_ALL_flowpot', 'CNT_ALL_precip', 
                      'CNT_MAM_flow', 'CNT_MAM_flowpot', 'CNT_MAM_precip',
                      'CNT_JJA_flow', 'CNT_JJA_flowpot', 'CNT_JJA_precip',
                      'CNT_SON_flow', 'CNT_SON_flowpot', 'CNT_SON_precip',
                      'CNT_DJF_flow', 'CNT_DJF_flowpot', 'CNT_DJF_precip',
                      'TOTflow', 'TOTprecip', 'year')
    return(outobj)
}

TeleCor <- function(TelIndex, inv, folder){
    filenames= inv[[1]]
    lat = inv[[2]]
    lon = inv[[3]]
    N = length(filenames)
    telcor = matrix(data=.0, nrow=5, ncol=length(filenames))
    telp = matrix(data=.5, nrow=5, ncol=length(filenames))
    okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
    for(i in 1:N){
        data = read.table(paste(folder, filenames[[i]], sep = '/'))
        names(data) = c("Y", "M", "Vdata")
        for(j in 1:5){
            comb = merge(data[is.element(data$M,okmonths[[j]]),],
                         TelIndex[is.element(TelIndex$M,okmonths[[j]]),])
            if(sum(comb[[3]])!=0){                
                tres = cor.test(comb$V,comb$Vdata,alternative='two.sided')
                telcor[j,i] = tres$estimate
                telp[j,i] = tres$p.value
            }
        }
    }
    return(list(telcor,telp))
}
        


AllTel <- function(ENSO,AMO,NAO,PNA,AO, EA, WP, EPNP, EAWR, SCA, TNH, POL, PT,precipinv,flowinv,pfolder,ffolder,potfolder){
    tels = list(ENSO,AMO,NAO,PNA,AO,EA,WP,EPNP,EAWR,SCA,TNH,POL,PT)
    tnames = c('ENSO','AMO','NAO','PNA','AO','EA','WP','EPNP','EAWR','SCA','TNH','POL','PT')
    invs = list(precipinv,flowinv,flowinv)
    folders = list(pfolder,ffolder,potfolder)
    varnames= c('Precipitation', 'Flow', 'POT')
    
    Ntels = length(tels)
    Nfolders = length(folders)
    for(t in 1:Ntels){
        for(f in 1:Nfolders){
            print(c(t,f))
            teleob = TeleCor(tels[[t]],invs[[f]],folders[[f]])
            pdf(paste(tnames[[t]],varnames[[f]],'.pdf', sep = ''))
            TeleMap(teleob,invs[[f]],paste(tnames[[t]],varnames[[f]]))
            dev.off()
        }
    }
    return
}











