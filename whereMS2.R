library(MetabFUN)
library(xcms)
library(CAMERA)

if (dir.exists('C:/Users/dbraas/Dropbox/R_projects')==T) dropbox <- 'C:/Users/dbraas/Dropbox/R_projects/'
if (dir.exists('E:/Dropbox_Metabolomics/Dropbox/projects')) dropbox <- 'E:/Dropbox_Metabolomics/Dropbox/projects/'
if (dir.exists('C:/Dropbox/R_projects')==T) dropbox <- 'C:/Dropbox/R_projects/'
if (dir.exists('D:/Dropbox/projects')==T) dropbox <- 'D:/Dropbox/projects/'
setwd(paste0(dropbox, 'XCMS_fun/MS2'))

inclusion <- read.csv("working inclusion list.CSV", header=T)
inclusion <- inclusion[,c(1,2,6,12)]
names(inclusion) <- c('Ion','Formula','Polarity','Name')
inclusion <- inclusion %>%
  mutate(minIon = Ion - Ion * 5 / 1e6,
         maxIon = Ion + Ion * 5 / 1e6) %>%
  select(Name, Ion, minIon, maxIon, Formula, Polarity) %>%
  arrange(Name)

files <- list.files(recursive=T, pattern = 'mzXML')
targeted <- xcmsRaw(files[3], includeMSn = T)
untargeted <- xcmsRaw(files[4], includeMSn = T)
untargeted_PM <- xcmsRaw(files[2], includeMSn = T)

#Function to look for MS2 for particular compound and spits out table with matching precursors

where.MS2 <- function(object, name){
  inc <- filter(inclusion, Name==name)
  EIC <- data.frame('RT'=object@scantime, 'Intensity'=rawEIC(object, mzrange=inc[3:4])$intensity)
  EIC$RT <- EIC$RT / 60
  result <- data.frame('RT' = object@msnRt / 60,
                       'Precursor'= object@msnPrecursorMz,
                       'Intensity' = object@msnPrecursorIntensity)
  result$Index <- as.numeric(rownames(result))
  result <- filter(result, Precursor > inc$minIon & Precursor < inc$maxIon)
  graph <- ggplot(EIC, aes(RT, Intensity)) +
    geom_line() +
    geom_hline(yintercept=0) +
    geom_point(data=result, aes(RT, Intensity), color='blue') +
    theme_bw(base_size=14) +
    labs(list(x='Retention time (min)', y='Intensity (A.U.)',
              title=paste0('Extracted Ion Chromatogram of ',inc$Name,'\n Scan range: ',round(inc$minIon,4),' to ',round(inc$maxIon, 4))))
  print(graph)
  return(result)
}

targeted_MS2 <- data.frame('Index' = targeted@msnScanindex,
                           'RT' = targeted@msnRt / 60,
                           'Precursor'= targeted@msnPrecursorMz,
                           'Intensity' = log(targeted@msnPrecursorIntensity,10))

getMS2 <- function(object, scan){
  MS2 <- top_n(data.frame(getMsnScan(object, scan)),30)
  MS2$intensity <- round(MS2$intensity / MS2$intensity[which.max(MS2$intensity)] * 100 ,2)
  MS2 <- filter(MS2, intensity > 5)
  plot(MS2, type='h', main=paste('Normalized MS2 spectrum of scan:', scan,
                                 '\nPrecursor m/z:', round(object@msnPrecursorMz[scan], 4)),xlab='m/z')
  return(MS2)
}

where.MS2b <- function(object, ion, ppm=10){
  
  maxIon <- ion + ion * ppm / 1e6
  minIon <- ion - ion * ppm / 1e6
  EIC <- data.frame('RT'=object@scantime, 'Intensity'=rawEIC(object, mzrange=as.matrix(c(minIon,maxIon)))$intensity)
  EIC$RT <- EIC$RT / 60
  result <- data.frame('RT' = object@msnRt / 60,
                       'Precursor'= object@msnPrecursorMz,
                       'Intensity' = object@msnPrecursorIntensity)
  result$Index <- as.numeric(rownames(result))
  result <- filter(result, Precursor >= minIon & Precursor <= maxIon)
  graph <- ggplot(EIC, aes(RT, Intensity)) +
    geom_line() +
    geom_hline(yintercept=0, color='firebrick') +
    geom_point(data=result, aes(RT, Intensity), color='blue') +
    theme_bw(base_size=14) +
    labs(list(x='Retention time (min)', y='Intensity (A.U.)',
              title=paste0('Extracted Ion Chromatogram of: ', ion,'\n Scan range: ',round(minIon,4),' to ',round(maxIon, 4))))
  print(graph)
  return(result)
}

what.MS2 <- function(object, srange = object@mzrange, RTrange = range(object@msnRt)/60, ppm=10){

  prec <- data.frame('Index' = object@msnScanindex,
                     'RT' = object@msnRt / 60,
                     'Precursor'= object@msnPrecursorMz,
                     'Intensity' = log(object@msnPrecursorIntensity,10))
  
  prec$Ion <- as.numeric(NA)
  
  Ion <- binscans(prec$Precursor)
  
  for (i in 1:length(Ion)){
    prec$Ion[Ion[i] >= prec$Precursor - prec$Precursor * ppm/1e6 & Ion[i] <= prec$Precursor + prec$Precursor * ppm/1e6] <- Ion[i]
  }
  
  prec <- filter(prec, 
                 Precursor >= srange[1],
                 Precursor <= srange[2],
                 RT >= RTrange[1],
                 RT <= RTrange[2])
  
  a <- ggplot(prec, aes(RT, Precursor, color=Intensity))+
    geom_point()+
    scale_color_distiller(type='div', palette = 5)+
    theme_bw()+
    ylim(srange)+
    labs(x='Rt (in min)', 
         y='Precursor m/z', 
         title= paste0('File: ', gsub('(.)*/','',object@filepath), '  Mode: ', unique(object@polarity)))
  
  print(a)
  return(prec)
}

#making bins for individual m/z 

binscans <- function(ions, ppm=10){
  
  l <- as.numeric()
  ions <- sort(ions)
  x=as.numeric()
  
  for (i in 2:length(ions)){
    if (length(l) == 0) {
      k=1 
      if (ions[i] > ions[k]*ppm/1e6+ions[k]) l <- append(l, i-1)
    }
    else if (ions[i] > ions[l[length(l)]+1]*ppm/1e6+ions[l[length(l)]+1]) l <- append(l, i-1)
  }
  
  for (k in 1:(length(l))){
    if (l[k]==1) x <- append(x, ions[1])
    else if (k==1 & l[k] > 1) x <- append(x, median(ions[1:l[1]]))
    else {
      x <- append(x, median(ions[(l[k-1]+1):(l[k])]))
    }
  }
  x <- append(x, median(ions[l[length(l)]:length(ions)]))
  return(x)
}
  
  



