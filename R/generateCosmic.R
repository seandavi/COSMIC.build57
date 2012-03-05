generateCosmic <- function(dat) {
  # dat is a data.frame from cosmic
  # meant to work with this download:
  # ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicMutantExport_v57_180112.tsv
  require(GenomicRanges)
  makelocs = function(dat,prefix='chr') {
    loc=dat[,colnames(dat) %in% c('Mutation.NCBI36.genome.position','Mutation.GRCh37.genome.position')]
    tmp = do.call(rbind,strsplit(as.character(loc),':'))
    tmp2 = do.call(rbind,strsplit(as.character(tmp[,2]),'-'))
    locs = data.frame(seqnames=paste(prefix,tmp[,1],sep=""),
      start=as.numeric(tmp2[,1]),
      end=as.numeric(tmp2[,2]))
    return(locs)
  }
  dat37 = dat[dat$Mutation.GRCh37.genome.position!='',]
  toRemove=which(colnames(dat37)=='Mutation.NCBI36.genome.position')
  dat37 = dat37[,-c(toRemove)]
  dat36 = dat[dat$Mutation.NCBI36.genome.position!='',]
  toRemove=which(colnames(dat36)=='Mutation.GRCh37.genome.position')
  dat36 = dat36[,-c(toRemove)]
  gr37locs = makelocs(dat37)
  gr36locs = makelocs(dat36)
  gr37 = with(gr37locs,GRanges(seqnames=seqnames,
    ranges=IRanges(start=start,end=end),elementMetadata=dat37))
  gr36 = with(gr36locs,GRanges(seqnames=seqnames,
    ranges=IRanges(start=start,end=end),elementMetadata=dat36))
  colnames(elementMetadata(gr37))=sub('elementMetadata\\.','',colnames(elementMetadata(gr37)))
  colnames(elementMetadata(gr36))=sub('elementMetadata\\.','',colnames(elementMetadata(gr36)))
  return(list(b37=unique(gr37),b36=unique(gr36)))
}
  
