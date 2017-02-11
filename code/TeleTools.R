##Read Teleconnection data
ensoraw = read.table('ENSOindex')
amoraw = read.table('AMOindex')

ENSO = data.frame(Y=numeric(),M=numeric(),V=numeric())
AMO = data.frame(Y=numeric(),M=numeric(),V=numeric())
for(y in 1:length(ensoraw[[1]])){
    for(m in 1:12){
        ENSO = rbind(ENSO,c(ensoraw[y,1],m,ensoraw[y,m+1]))
    }
}

for(y in 1:length(amoraw[[1]])){
    for(m in 1:12){
        AMO = rbind(AMO,c(amoraw[y,1],m,amoraw[y,m+1]))
    }
}

names(ENSO) = c('Y','M','V')
names(AMO) = c('Y','M','V')
NAO = read.table('NAOindex',col.names=c('Y','M','V'))
PNA = read.table('PNAindex',col.names=c('Y','M','V'))
AO = read.table('AOindex',col.names=c('Y','M','V'))

a = read.table('Tele.txt',header=TRUE)

b = data.frame(a$Y,a$M,a$EA)
names(b) = c('Y','M','V')
EA= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$WP)
names(b) = c('Y','M','V')
WP= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$EPNP)
names(b) = c('Y','M','V')
EPNP= b[b$V != -99.90,]

#b = data.frame(a$Y,a$M,a$PNA)
#names(b) = c('Y','M','V')
#PNA= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$EAWR)
names(b) = c('Y','M','V')
EAWR= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$SCA)
names(b) = c('Y','M','V')
SCA= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$TNH)
names(b) = c('Y','M','V')
TNH= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$POL)
names(b) = c('Y','M','V')
POL= b[b$V != -99.90,]

b = data.frame(a$Y,a$M,a$PT)
names(b) = c('Y','M','V')
PT= b[b$V != -99.90,]

NAO$V = runavg(NAO$V)
PNA$V = runavg(PNA$V)
AO$V = runavg(AO$V)
