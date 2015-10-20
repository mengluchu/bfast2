# 1) STFDF: covert array to stsdf: 2) STSDF: array store the change point deter spatial polygon list; prodes
# spatial polygon list 3) way 1: manually build spatial and temporal index 4) way 2: convert from array to
# STFDF, then to STIDF, then to STSDF (from STIDF to STSDF the endtime is by default) 5) polygontopoint:
# rasterize polygons to spatialpoint

# arraytoSTFDF2 : 1) deterpolytopointstoSTSDF 2) 3) changearraytoSTSDF: 2) 3) waystfdf2stsdf: 4)
 

arraytoSTFDF2 <- function(array, attr.name = "value", alltime = bt, x = c(58930:59079), y = c(48210:48359), 
    crs, months = 0.3) {
    
    t <- c(1:dim(array)[3])
    tt2 <- unique(as.POSIXct(alltime[t]))  # array index to time
    y1 <- rep(y, each = length(x))
    x1 <- rep(x, length(y))
    xyd <- as.data.frame(cbind(x1, y1))
    xyc <- getxyMatrix(xyd, 231.6564)
    xyc <- as.data.frame(xyc)
    
    names(xyc) <- c("x", "y")
    coordinates(xyc) <- ~x + y
    proj4string(xyc) <- crs
    
    data2 <- as.vector(array)
    
    data2 <- data.frame(data2)
    if (months != 0) {
        ltime <- tt2 + months * 3600 * 24 * 30
        etime <- tt2 - months * 3600 * 24 * 30
        
        stfdf1 <- STFDF(xyc, etime, data2, ltime)
    } else {
        stfdf1 <- STFDF(xyc, tt2, data2)
    }
    return(stfdf1)
}

polygontopoint <- function(x = c(58930:59079), y = c(48210:48359), sppolygon, crs = MODISCRS) {
    # sppolygon<-prodes67808
    
    x1 <- rep(x, each = length(y))
    y1 <- rep(y, length(x))
    
    xyd <- as.data.frame(cbind(x1, y1))
    xyc <- getxyMatrix(xyd, 231.6564)
    xyc <- as.data.frame(xyc)
    coordinates(xyc) <- ~x + y
    SpatialPoints(xyc)
    # as.character(round(unique(xyc@coords[,2])))
    
    proj4string(xyc) <- crs
    
    proj4string(sppolygon) <- CRS("+proj=utm +zone=21 +south  +ellps=WGS84 ")
    sppolygon = spTransform(sppolygon, CRS(crs))
    points <- xyc[sppolygon, ]
    
    return(points)
}
# p1<-polygontopoint(c(58930:59079),c(48210:48359),spatialPolygons)

deterpolytopointstoSTSDF <- function(tdall2 = tdall, pp01 = pp0, rpp2 = rpp, crs = MODISCRS, attr.name = "value", 
    months = 0.3) {
    
    detertime <- tdall2
    dtime <- as.POSIXct(detertime)
    
    uniquetime <- unique(dtime)
    rp1 <- table(dtime)
    tt2 <- uniquetime[order(uniquetime)]
    
    lt <- length(uniquetime)
    
    tn <- unlist(lapply(lt:1, function(i) rep(i, rp1[i])))
    
    tn2 <- unlist(lapply(1:length(tn), function(i) rep(tn[i], rpp2[[i]])))  #rpp2: the number of points
    
    sn2 <- c(1:length(pp01))
    index1 <- as.matrix(cbind(sn2, tn2))
    
    data2 <- data.frame(rep(1, length(pp01)))
    names(data2) <- attr.name
    proj4string(pp01) <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
    
    ltime <- tt2 + months * 3600 * 24 * 30
    etime <- tt2 - months * 3600 * 24 * 30
    stsdf2 <- STSDF(pp01, etime, index = index1, data2, ltime)
    
    # stpolygon<-spTransform(stsdf2, CRS(crs))
    
    return(stsdf2)
}

changearraytoSTSDF <- function(array, crs, attr.name = "value", alltime = bt, x = c(58930:59079), y = c(48210:48359), 
    months = 6) {
    # woud be easier to get the unordered index of time than get the right index of space
    change7 <- which(!is.na(array), arr.ind = TRUE)  #0.05
    change8 <- which(!is.na(array))
    x11 <- change7[, 1] + x[1] - 1  # for the second 150 by 150 array   
    y11 <- change7[, 2] + y[1] - 1
    t <- change7[, 3]
    
    xyd <- data.frame(cbind(x11, y11))
    xyc <- getxyMatrix(xyd, 231.6564)
    xyc <- data.frame(xyc)
    names(xyc) <- c("x", "y")
    xyc$ID <- c(1:length(x11))
    coordinates(xyc) <- ~x + y
    
    tt2 <- unique(as.POSIXct(alltime[t]))  # array index to time
    
    zdxyc <- zerodist(xyc)  # non-unique spatial points   
    
    sel = lapply(2:length(zdxyc[, 2]), function(i) !identical(zdxyc[i - 1, 2], zdxyc[i, 2]))
    sel <- unlist(sel)
    sel2 <- c(TRUE, sel)
    sel2 <- sel2[-length(sel2)]
    zuni = zdxyc[sel2, ]
    length(zdxyc[, 2])
    
    spl <- length(xyc)
    
    uniquexyc <- remove.duplicates(xyc)
    dunixyc <- data.frame(uniquexyc)
    dunixyc$IDnew <- c(1:length(uniquexyc))
    
    noID <- data.frame(cbind(dunixyc$ID, dunixyc$IDnew))
    
    remo <- match(zuni[, 1], dunixyc$ID)  # return the index of ID
    
    oldindex <- c(dunixyc$ID, zuni[, 2])
    newindex <- c(dunixyc$IDnew, noID[remo, 2])
    
    oldandnewindex <- data.frame(cbind(newindex, oldindex))
    names(oldandnewindex) <- c("new", "old")
    oldnew <- oldandnewindex[order(oldindex), ]
    lt <- length(tt2)
    # spl<-length(newone2) runi<-rank(unique())
    tn <- lapply(1:lt, function(i) rep(i, table(t)[i]))
    # sn<-lapply(1:spl, function(i) rep(runi[i],newone2[i]))
    index1 <- as.matrix(cbind(oldnew$new, unlist(tn)))
    
    data2 <- na.omit(as.vector(array))
    names(data2) <- attr.name
    data2 <- as.data.frame(data2)
    
    
    ltime <- tt2 + months * 3600 * 24 * 30
    etime <- tt2 - months * 3600 * 24 * 30
    spxyc <- as(uniquexyc, "SpatialPoints")
    proj4string(spxyc) <- crs
    stsdf1 <- STSDF(spxyc, etime, index = index1, data2, ltime)
    return(stsdf1)
}

# not so computationally efficient, but easiest way, pay attention to the x and y sequence
waystfdf2stsdf <- function(data = t3darrbfamul2, alltime = bt, x = c(58930:59079), y = c(48210:48359), crs = MODISCRS, 
    months = 0.3) {
    changestfdf <- arraytoSTFDF2(array = data, alltime = alltime, x = x, y = y, crs = crs, months = months)  #change array, points
    # which(!is.na(edivisive1[,,1:160]),arr.ind=TRUE)
    stsdfchangepoints1 <- as(changestfdf, "STIDF")
    stsdfchangepoints2 <- as(stsdfchangepoints1, "STSDF")
    return(stsdfchangepoints2)
}
 