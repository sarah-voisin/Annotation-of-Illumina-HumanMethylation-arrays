#Obtain complete annotation for the 3 arrays
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data")
library(tidyverse)
library(GEOquery)

arrays = c("EPIC",
           "HM450",
           "HM27")

for (a in arrays)
{
    #Download annotation of EPIC array from Zhou et al. 2016, Nucleic Acid Research
    #We can download the file in a plain tsv (tab delimited format similar to txt) format instead
    if (!file.exists(paste0(a,".hg38.manifest.gencode.v22.tsv")))
        {
        download.file(url=paste0("http://zwdzwd.github.io/InfiniumAnnotation#current/",a,"/",a,".hg38.manifest.gencode.v22.tsv.gz"),
                      destfile=paste0(a,".hg38.manifest.gencode.v22.tsv")) #you need to unzip the file
        
        gunzip(filename=paste0(a,".hg38.manifest.gencode.v22.tsv.gz"),
               destname=paste0(a,".hg38.manifest.gencode.v22.tsv"),remove=T)
        }
    #Read the file
    annotation <- read.delim(paste0(a,".hg38.manifest.gencode.v22.tsv")) %>%
        drop_na(CpG_chrm) %>%
        select(probeID,
               CpG_chrm,
               CpG_beg,
               CpG_end,
               genesUniq,
               CGIposition)

    #Roadmap Epigenomics Project on chromatin states
    if (!file.exists(paste0(a,".hg19.REMC.chromHMM.tsv")))
        {
        download.file(url=paste0("http://zwdzwd.github.io/InfiniumAnnotation#current/",a,"/",a,".hg19.REMC.chromHMM.tsv.gz"),
                      destfile=paste0(a,".hg19.REMC.chromHMM.tsv.gz"))
        
        gunzip(filename=paste0(a,".hg19.REMC.chromHMM.tsv.gz"),
               destname=paste0(a,".hg19.REMC.chromHMM.tsv"),remove=T)
        }
    #Let's read the file
    chrom_states <- read.delim(paste0(a,".hg19.REMC.chromHMM.tsv")) %>%
        select(probeID,E107,E108)

    #Merge
    annotation <- left_join(annotation,
                            chrom_states,
                            by = "probeID")

    #ENCODE TF binding
    if (!file.exists(paste0(a,".hg19.ENCODE.TFBS.tsv")))
        {
        download.file(url=paste0("http://zwdzwd.github.io/InfiniumAnnotation#current/",a,"/",a,".hg19.ENCODE.TFBS.tsv.gz"),
                      destfile=paste0(a,".hg19.ENCODE.TFBS.tsv.gz"))
        
        gunzip(filename=paste0(a,".hg19.ENCODE.TFBS.tsv.gz"),
               destname=paste0(a,".hg19.ENCODE.TFBS.tsv"),remove=T)
        }
    #Let's read the file
    ENCODE <- read.delim(paste0(a,".hg19.ENCODE.TFBS.tsv")) %>%
        filter(sample=="HSMMtube") %>% #Obtain data from the right cell type (HSMMtube = differentiated myoblasts)
        select(probeID,
               TF) %>%
        distinct()#Remove duplicated rows
    #Add column
    ENCODE <- ENCODE %>%
        mutate(presence = rep("yes",nrow(ENCODE)))
    #Pivot data by TF
    ENCODE <- ENCODE %>%
        pivot_wider(names_from = TF,
                    names_prefix = "HSMMtube_",
                    values_from = presence)
    
    #Merge
    annotation <- left_join(annotation,
                            ENCODE,
                            by = "probeID")
    
    #Write the whole annotation into one single table
    write.table(annotation,
                file=paste0("Annotation_",a,".txt"),
                quote=FALSE,
                row.names=FALSE,
                col.names=TRUE,
                sep='\t')
    }

#The only annotation I would add is GeneHancer,
#a database of human regulatory elements (enhancers and promoters)
#and their inferred target genes.
#The GeneHancer database was created by integrating >1 million
#regulatory elements from multiple genome-wide databases.
#Associations between the regulatory elements and target genes
#were based on multiple sources of linking molecular data,
#along with distance. For more information, visit:
#https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=806048677_vUXhmharpmwDl90vgwwGQaI4PjbV&c=chr1&g=geneHancer

#GeneHancer is available for visualisation on the Genome browser
# (https://genome.ucsc.edu/index.html)
#Usually, one can download tracks from the Genome browser as they are
#open-access, but GeneHancer belongs to the Weizmann Institute, so you
#can only have limited access to it.

#Read all annotation files
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data")
library(tidyverse)
annotation_EPIC <- read_tsv("Annotation_EPIC.txt")
annotation_HM450 <- read_tsv("Annotation_HM450.txt")
annotation_HM27 <- read_tsv("Annotation_HM27.txt")

annotation <- full_join(x = annotation_EPIC,
                        y = annotation_HM450)

annotation <- full_join(x = annotation,
                        y = annotation_HM27)

#Replace orf gene names with most updated GeneCard database
library(GeneBook)
genesUniq_with_updated_names <- annotation$genesUniq
browsing <- genecard_id$subname
for (i in 1:nrow(annotation))
{
    gene <- genesUniq_with_updated_names[i]
    if(!is.na(gene))
    {
    if (str_detect(gene,
                   "C[:digit:]+orf[:digit:]+"))
    {
        pattern <- str_extract(gene,
                              "C[:digit:]+orf[:digit:]+")
        index_in_genebook <- str_which(browsing,
                                       paste0('"',pattern,'"'))
        if (length(index_in_genebook)>0)
        {
            print(i)
            print("Replacement")
            genesUniq_with_updated_names[i] = str_replace(string = gene,
                                                          pattern = pattern,
                                                          replacement = paste0(genecard_id[index_in_genebook,"gene"],collapse=";"))
        }
    }
    }
}

annotation2 <- annotation
annotation2$genesUniq <- genesUniq_with_updated_names

#Replace KIAA gene names with most updated GeneCard database
library(GeneBook)
genesUniq_with_updated_names <- annotation2$genesUniq
browsing <- genecard_id$subname
for (i in 1:nrow(annotation))
{
    gene <- genesUniq_with_updated_names[i]
    if(!is.na(gene))
    {
        if (str_detect(gene,
                       "KIAA[:digit:]+"))
        {
            pattern <- str_extract(gene,
                                   "KIAA[:digit:]+")
            index_in_genebook <- str_which(browsing,
                                           paste0('"',pattern,'"'))
            if (length(index_in_genebook)>0)
            {
                print(i)
                print("Replacement")
                genesUniq_with_updated_names[i] = str_replace(string = gene,
                                                              pattern = pattern,
                                                              replacement = paste0(genecard_id[index_in_genebook,"gene"],collapse=";"))
            }
        }
    }
}

annotation2$genesUniq <- genesUniq_with_updated_names

library(GenomicRanges)
annotation_GR <- makeGRangesFromDataFrame(annotation2,
                                         keep.extra.columns=TRUE,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field="CpG_chrm",
                                         start.field="CpG_beg",
                                         end.field="CpG_end",
                                         strand.field="probe_strand",
                                         starts.in.df.are.0based=FALSE)

#Add GeneHancer data
GH_files <- list.files()[grep("GeneHancer",list.files())]
GH <- read_tsv(GH_files[1])
for (g in GH_files[-1])
{
    gh <- read_tsv(g)
    GH <- full_join(GH,
                    gh)
}

library(GenomicRanges)
GH_GR <- makeGRangesFromDataFrame(GH,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("geneHancerChrom"),
                         start.field="geneHancerStart",
                         end.field=c("geneHancerEnd"),
                         strand.field="geneStrand",
                         starts.in.df.are.0based=FALSE)

genesidx <- as.data.frame(findOverlaps(annotation_GR, GH_GR))
genesover <- tapply(genesidx$subjectHits, genesidx$queryHits, 
                    function(x) GH_GR$geneName[x])
op.A <- sapply(genesover, function(l) paste(l, collapse = ";"))
name.A <- names(genesover)
m.A <- as.numeric(name.A)
M <- length(annotation_GR)
overlapping.genes <- rep("", M)
overlapping.genes[m.A] <- op.A
annotation2$GeneHancer_interaction <- overlapping.genes

#Replace NA with white space
annotation2 <- annotation2 %>%
    mutate(genesUniq = replace_na(genesUniq,""))

#Annotate those CpGs in enhancers and intersected with a GeneHancer enhancer
gene_enhancer_annotate <- function(v)
{
    output <- v["genesUniq"]
    if (grepl("Enh",v["E107"])&grepl("Enh",v["E108"]))
        output <- v["GeneHancer_interaction"]
    return(output)
}
genesUniq_with_enh <- apply(annotation2,1,gene_enhancer_annotate)
annotation2$genesUniq_with_enh <- genesUniq_with_enh

setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Annotation of DNA methylation data")
write.table(annotation2,
            file="Annotation.txt",
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE,
            sep='\t')
