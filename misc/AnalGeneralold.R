library(trend)
library(maps)

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
    years[[N1]]
    dv = median(v[N1:length(v)]) - median(v[1:N1])
    return(c(1, years[[N1]], p, dv))
}

normalizemedian <- function(v, cps){
    cps=sort(cps)
    segments=list()
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
        cp[[1]] = cp[[1]] + 1
        cpindex = which(is.element(years, cp[3*(1:cp[[1]])-1]))
        v1=normalizemedian(v, cpindex)
        return(Lanzante(v1, years, cp, pcut, buffer))
    }
}
        





generate1CP <- function(folder, filenames, lat, lon, early = 1940:1990, late = 1970:2014, 
                       ally = 1940:2014, test_func = pettitt_it){
    N = length(filenames)
    CP = data.frame(filenames, lat, lon, 
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N))
    names(CP) <- c('ID', 'LAT', 'LON',
                   'ALL', 'ALLp', 'ALLd',
                   'ALLe', 'ALLpe', 'ALLde',
                   'ALLl', 'ALLpl', 'ALLdl',
                   'DJF', 'DJFp', 'DJFd',
                   'DJFe', 'DJFpe', 'DJFde',
                   'DJFl', 'DJFpl', 'DJFdl',
                   'MAM', 'MAMp', 'MAMd',
                   'MAMe', 'MAMpe', 'MAMde',
                   'MAMl', 'MAMpl', 'MAMdl',
                   'JJA', 'JJAp', 'JJAd',
                   'JJAe', 'JJApe', 'JJAde',
                   'JJAl', 'JJApl', 'JJAdl',
                   'SON', 'SONp', 'SONd',
                   'SONe', 'SONpe', 'SONde',
                   'SONl', 'SONpl', 'SONdl')
    for(i in 1:length(filenames)){
        data = read.table(paste(folder, filenames[[i]], sep = '/'))
        years = unique(data[[1]])
        okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
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
                okyr = count==12
            } else{
                okyr = count==3
            }
            y_all = is.element(years,ally)
            test_result = test_func(vals[okyr & y_all],years[okyr & y_all])
            CP[[4+(s-1)*9]][[i]] = test_result[[2]]
            CP[[5+(s-1)*9]][[i]] = test_result[[3]]
            CP[[6+(s-1)*9]][[i]] = test_result[[4]]
            y_early = is.element(years,early)
            test_result = test_func(vals[okyr & y_early],years[okyr & y_early])
            CP[[7+(s-1)*9]][[i]] = test_result[[2]]
            CP[[8+(s-1)*9]][[i]] = test_result[[3]]
            CP[[9+(s-1)*9]][[i]] = test_result[[4]]
            y_late = is.element(years,late)
            test_result = test_func(vals[okyr & y_late],years[okyr & y_late])
            CP[[10+(s-1)*9]][[i]] = test_result[[2]]
            CP[[11+(s-1)*9]][[i]] = test_result[[3]]
            CP[[12+(s-1)*9]][[i]] = test_result[[4]]
        }
    }
    return(CP)
}


generateallCP <- function(folder, filenames, lat, lon, early = 1940:1990, late = 1970:2014, 
                       ally = 1940:2014, test_func = pettitt_it){
    N = length(filenames)
    CP = data.frame(filenames, lat, lon, 
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N),
                    rep(0,N), rep(0,N), rep(0,N))
    names(CP) <- c('ID', 'LAT', 'LON',
                   'ALL', 'ALLp', 'ALLd',
                   'ALLe', 'ALLpe', 'ALLde',
                   'ALLl', 'ALLpl', 'ALLdl',
                   'DJF', 'DJFp', 'DJFd',
                   'DJFe', 'DJFpe', 'DJFde',
                   'DJFl', 'DJFpl', 'DJFdl',
                   'MAM', 'MAMp', 'MAMd',
                   'MAMe', 'MAMpe', 'MAMde',
                   'MAMl', 'MAMpl', 'MAMdl',
                   'JJA', 'JJAp', 'JJAd',
                   'JJAe', 'JJApe', 'JJAde',
                   'JJAl', 'JJApl', 'JJAdl',
                   'SON', 'SONp', 'SONd',
                   'SONe', 'SONpe', 'SONde',
                   'SONl', 'SONpl', 'SONdl')
    for(i in 1:length(filenames)){
        data = read.table(paste(folder, filenames[[i]], sep = '/'))
        years = unique(data[[1]])
        okmonths = list(1:12, c(1,2,12), 3:5, 6:8, 9:11)
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
                okyr = count==12
            } else{
                okyr = count==3
            }
            y_all = is.element(years,ally)
            test_result = test_func(vals[okyr & y_all],years[okyr & y_all])
            CP[[4+(s-1)*9]][[i]] = test_result[[2]]
            CP[[5+(s-1)*9]][[i]] = test_result[[3]]
            CP[[6+(s-1)*9]][[i]] = test_result[[4]]
            y_early = is.element(years,early)
            test_result = test_func(vals[okyr & y_early],years[okyr & y_early])
            CP[[7+(s-1)*9]][[i]] = test_result[[2]]
            CP[[8+(s-1)*9]][[i]] = test_result[[3]]
            CP[[9+(s-1)*9]][[i]] = test_result[[4]]
            y_late = is.element(years,late)
            test_result = test_func(vals[okyr & y_late],years[okyr & y_late])
            CP[[10+(s-1)*9]][[i]] = test_result[[2]]
            CP[[11+(s-1)*9]][[i]] = test_result[[3]]
            CP[[12+(s-1)*9]][[i]] = test_result[[4]]
        }
    }
    return(CP)
}

mapCP <- function(CP, title, pdfname){
    pdf(pdfname)
    S = c('All months', 'DJF', 'MAM', 'JJA', 'SON')
    t = c('1940-2015', '1940-1990', '1970-2015')
    tmin = c(193, 193, 196) 
    pall = c('red', 'orange', 'gold4', 'green', 'turquoise4', 'blue', 'violetred3', 'purple' )
    for(i in 1:3){
        for(j in 1:5){
            map('state')
            pt_color = rep('grey', length(CP$LAT))
            pt_pch  = rep(20, length(CP$LAT))
            sig = CP[[5+(i-1)*3+(j-1)*9]] < .5
            pt_pch[sig & CP[[6+(i-1)*3+(j-1)*9]] > 0]    <- 24
            pt_pch[sig & CP[[6+(i-1)*3+(j-1)*9]] < 0]    <- 25
            pt_color[sig] = pall[floor(CP[[4+(i-1)*3+(j-1)*9]]/10)[sig] - tmin[[i]]]
            points(CP$LON, CP$LAT, pch= pt_pch, bg = pt_color)
            title(main=paste(S[[j]],t[[i]]))
            if(i == 1){
                legend(-80.584,31.983, c("40's", "50's", "60's", "70's", "80's", "90's", "00's"),
                       col=pall, pch=20)
            } else if(i==2){
                legend(-80.584,31.983, c("40's", "50's", "60's", "70's", "80's", "90's"),
                       col=pall, pch=20)
            } else{
                legend(-80.584,31.983, c("70's", "80's", "90's", "00's"),
                       col=pall, pch=20)
            }
        }
    }
    dev.off()
}

maoCPdecadal <- function(CP, title, pdfname){
    pdf(pdfname)
    for i in CP
    dev.off()
    


precipinv = read.table('PrecipInv')
names(precipinv) = c('ID','LAT','LON','STATE') 
precipCP = generateCP('monthlyprecip', precipinv$ID, precipinv$LAT, precipinv$LON)
mapCP(precipCP,'ghdj','flowCP.pdf')

flowinv = read.table('FlowInv', sep = '\t', colClasses=c('character', 'numeric',
               'numeric','character','character'))
names(flowinv) = c('ID','LAT','LON','STATE','REGION') 
flowCP = generateCP('monthlyflow', flowinv$ID, flowinv$LAT, flowinv$LON)


