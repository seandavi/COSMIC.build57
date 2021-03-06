\name{cosmicGRCh37}
\alias{cosmicGRCh37}
\docType{data}
\title{
COSMIC data mapped to GRCh37 (hg19)
}
\description{
COSMIC data mapped to GRCh37 (hg19)
}
\usage{data(cosmicGRCh37)}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 24 levels "chr1","chr10",..: 2 4 15 18 4 18 12 8 21 9 ...
  .. .. ..@ lengths        : int [1:11937] 2 13 3 1 6 3 2 2 2 1 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:65870] 52575931 52576024 9243930 9259164 9251272 9242564 9260237 9258837 9225356 9229975 ...
  .. .. ..@ width          : int [1:65870] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3
  .. .. ..@ lengths        : int 65870
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 65870
  .. .. ..@ listData       :List of 24
  .. .. .. ..$ Gene.name                      : Factor w/ 15198 levels "A1CF","A2M","A2ML1",..: 1 1 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ Accession.Number               : Factor w/ 15158 levels "AB023190","AB051467.1",..: 9114 9114 10964 10964 10964 10964 10964 10964 10964 10964 ...
  .. .. .. ..$ HGNC.ID                        : int [1:65870] 24086 24086 7 7 7 7 7 7 7 7 ...
  .. .. .. ..$ Sample.name                    : Factor w/ 153739 levels " HOG06-122"," HOG06-1221",..: 118608 152091 111369 118525 128147 128152 128178 148952 148957 148996 ...
  .. .. .. ..$ ID_sample                      : int [1:65870] 1312984 1604872 1551882 1312971 1560554 1560559 1560588 1287214 1287216 1287244 ...
  .. .. .. ..$ ID_tumour                      : int [1:65870] 1223893 1524191 1473988 1223880 1482004 1482009 1482038 1198492 1198494 1198522 ...
  .. .. .. ..$ Primary.site                   : Factor w/ 43 levels "Fallopian tube",..: 8 33 33 8 40 40 40 8 8 8 ...
  .. .. .. ..$ Site.subtype                   : Factor w/ 172 levels "NS","abdomen",..: 16 44 1 16 70 70 96 16 16 16 ...
  .. .. .. ..$ Primary.histology              : Factor w/ 84 levels "Ewings_sarcoma-peripheral_primitive_neuroectodermal_tumour",..: 40 23 57 40 23 23 23 40 40 40 ...
  .. .. .. ..$ Histology.subtype              : Factor w/ 399 levels "ACTH","ACTH-FSH",..: 82 373 33 82 373 373 373 82 82 82 ...
  .. .. .. ..$ Genome.wide.screen             : Factor w/ 3 levels "","n","y": 1 3 3 1 3 3 3 1 1 1 ...
  .. .. .. ..$ Mutation.ID                    : int [1:65870] 40879 138756 110440 39653 127216 127215 127217 35965 35963 35964 ...
  .. .. .. ..$ Mutation.CDS                   : Factor w/ 38881 levels "","Intronic",..: 38446 36685 14397 37309 9399 16660 34285 1522 22775 20223 ...
  .. .. .. ..$ Mutation.AA                    : Factor w/ 54163 levels "","p.","p.(1409)ins6",..: 9502 31055 16907 47270 51096 17136 53295 15450 11774 46083 ...
  .. .. .. ..$ Mutation.Description           : Factor w/ 15 levels "Complex - compound substitution",..: 11 11 11 11 13 13 13 11 11 13 ...
  .. .. .. ..$ Mutation.zygosity              : Factor w/ 3 levels "het","hom","unk": 1 3 2 1 3 3 3 1 1 1 ...
  .. .. .. ..$ Mutation.NCBI36.strand         : Factor w/ 3 levels "","+","-": 3 1 1 3 1 1 1 3 3 3 ...
  .. .. .. ..$ Mutation.GRCh37.genome.position: Factor w/ 65871 levels "","10:100011329-100011329",..: 1305 1306 9761 9764 9762 9760 9765 9763 9757 9758 ...
  .. .. .. ..$ Mutation.GRCh37.strand         : Factor w/ 3 levels "","+","-": 3 3 3 3 3 3 3 3 3 3 ...
  .. .. .. ..$ Mutation.somatic.status        : Factor w/ 6 levels "Confirmed germline variant",..: 2 2 2 2 6 6 6 2 2 2 ...
  .. .. .. ..$ Pubmed_PMID                    : int [1:65870] 18772396 21984974 21499247 18772396 21798893 21798893 21798893 18772890 18772890 18772890 ...
  .. .. .. ..$ Sample.source                  : Factor w/ 27 levels "NS","autopsy - NOS",..: 23 23 24 23 24 24 24 24 24 24 ...
  .. .. .. ..$ Tumour.origin                  : Factor w/ 8 levels "NS","adenoma adjacent to primary tumour",..: 7 5 4 5 5 5 5 6 7 7 ...
  .. .. .. ..$ Comments                       : Factor w/ 2723 levels "","Age at tumour recurrence:1 year after diagnosis of TMD",..: 1002 1 1205 1002 1 1 1 1 1 1 ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots
  .. .. ..@ seqnames   : chr [1:24] "chr1" "chr10" "chr11" "chr12" ...
  .. .. ..@ seqlengths : int [1:24] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:24] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:24] NA NA NA NA ...
  ..@ metadata       : list()
}
\details{
The data are in a GRanges object with elementMetadata from the original
download.  Note that some of the ranges are large, so one may need to
filter when looking for overlapping variants.
}
\source{
\url{ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicMutantExport_v57_180112.tsv}
}
\references{
\url{http://www.sanger.ac.uk/genetics/CGP/cosmic/}
}
\examples{
data(cosmicGRCh37)
}
\keyword{datasets}
