library(DBapps)
library(xcms)
library(CAMERA)

setwd("C:/Users/danie/Dropbox/projects/XCMS/MS2")

inclusion <- read.csv("working inclusion list.csv", header=T)
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

targeted

plot(targeted@scantime/60, targeted@tic, type='l')
plotEIC(targeted, mzrange=c(179.055,179.057), scanrange=c(200,400)) #test has to be xcmsRaw

target_peaks=as.data.frame(findPeaks(targeted,method='centWave',
                       ppm=3, snthresh=20,
                       prefilter=c(15,5000), peakwidth=c(10,120)))
plotEIC(targeted, target_peaks[1,2:3], as.numeric(target_peaks[1,5:6]))
plotEIC(targeted, target_peaks[10,2:3])

plotPeaks(targeted, as.matrix(target_peaks[1:16,]), c(4,4))
levelplot(targeted, col.regions=colorRampPalette(c('blue','red'))(256))
image(targeted)

targeted_MS2 <- data.frame('Index' = targeted@msnScanindex,
                           'Scan' = rownames(targeted_MS2),
                           'RT' = targeted@msnRt / 60,
                         'Precursor'= targeted@msnPrecursorMz,
                         'Intensity' = targeted@msnPrecursorIntensity)

#Function to look for MS2 for particular compound and spits out table with matching precursors

where.MS2 <- function(object, name){
  inc <- filter(inclusion, Name==name)
  EIC <- data.frame(getEIC(object, as.matrix(inc[3:4], ncol=2))@eic)
  names(EIC) <- c('RT','Intensity')
  EIC$RT <- EIC$RT / 60
  result <- filter(targeted_MS2, Precursor > inc$minIon & Precursor < inc$maxIon)
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

Glc <- top_n(data.frame(getMsnScan(targeted,337)),30)

standard <- xcmsSet(files=files[3:4], method='centWave', ppm=3, snthresh=20,
                    prefilter=c(10, 100000), peakwidth=c(5,120))

grouped <- group(standard)
annotated <- xsAnnotate(grouped)
annotated <- groupFWHM(annotated)
annotated <- groupCorr(annotated)
annotatedI <- findIsotopes(annotated)
annotatedIA <- findAdducts(annotatedI, polarity='negative')
peaklist <- getPeaklist(annotatedIA)
View(peaklist)

plotEICs(annotatedIA, pspec=1, maxlabel=5)
