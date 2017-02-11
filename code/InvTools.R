
USGSstatecodes <- c('AL','AK','XT','AZ','AR','CA','XT','CO','CT','DE',
                    'DC','FL','GA','XT','HI','ID','IL','IN','IA','KS',
                    'KY','LA','ME','MD','MA','MI','MN','MS','MO','MT',
                    'NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK',
                    'OR','PA','XT','RI','SC','SD','TN','TX','UT','VT',
                    'VA','XT','WA','WV','WI','WY','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT')
                    
StateNames <-     c('alabama', 'alaska', 'xt1', 'arizona', 'arkansas', 
                    'california', 'xt2', 'colorado', 'connecticut', 'delaware',

                    'district of columbia', 'florida', 'georgia', 'xt3', 'hawaii', 
                    'idaho', 'illinois', 'indiana', 'iowa' ,'kansas', 
                    
                    'kentucky' , 'louisiana', 'maine', 'maryland', 'massachusetts',
                    'michigan', 'minnesota', 'mississippi', 'missouri', 'montana', 
                    
                    'nebraska', 'nevada', 'new hampshire', 'new jersey', 'new mexico', 
                    'new york', 'north carolina', 'north dakota', 'ohio', 'oklahoma',
                    
                    'oregon', 'pennsylvania', 'xt4', 'rhode island', 'south carolina', 
                    'south dakota', 'tennessee', 'texas', 'utah', 'vermont', 
                    
                    'virginia', 'xt5', 'washington', 'west virginia', 'wisconsin', 
                    'wyoming', 'xt6', 'xt7', 'xt8', 'xt9',
                    
                    'xt10','xt11','xt12','xt13','xt14', 
                    'xt15','xt16','xt17','xt18','xt19', 
                    
                    'xt20','xt21','xt22','xt23','xt24', 
                    'xt25','xt26','xt27','xt28','xt29', 
                    
                    'xt30','xt31','xt32','xt33','xt34', 
                    'xt35','xt36','xt37','xt38','xt39', 
                    
                    'xt40','xt41','xt42','xt43','xt44', 
                    'xt45','xt46','xt47','xt48','xt49')
 
USGSRegionList <- c('SE','AK','XT','SW','SE','SW','XT','SW','NE','NE',
                    'NE','SE','SE','XT','HI','NW','MW','MW','MW','SC',
                    'SE','SE','NE','NE','NE','MW','MW','SE','MW','NC',
                    'NC','SW','NE','NE','SW','NE','SE','NC','MW','SC',
                    'NW','NE','XT','NE','SE','NC','SE','SC','SW','NE',
                    'SE','XT','NW','NE','MW','NC','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT',
                    'XT','XT','XT','XT','XT','XT','XT','XT','XT','XT')


precipinv = read.table('misc/pinv1', stringsAsFactors=FALSE)
names(precipinv) = c('ID','LAT','LON','STATE') 
pindex = match(precipinv$STATE,USGSstatecodes)
precipinv$REGION = USGSRegionList[pindex]

flowinv = read.table('misc/finv1', sep = '\t', colClasses=c('character', 'numeric',
               'numeric','character'), stringsAsFactors=FALSE)
names(flowinv) = c('ID','LAT','LON','STATE') 
findex = match(flowinv$STATE,StateNames)
flowinv$REGION = USGSRegionList[findex]


peakinv = read.table('misc/fpinv1', sep = '\t', colClasses=c('character', 'numeric',
               'numeric','character'), stringsAsFactors=FALSE)
names(peakinv) = c('ID','LAT','LON','STATE') 
peakinv$REGION = USGSRegionList[as.integer(peakinv$STATE)]

rmvlist = read.table('misc/rmvlist', colClasses=c('character'),stringsAsFactors=FALSE)[[1]]
peakinv = peakinv[!is.element(peakinv$ID,rmvlist),]

