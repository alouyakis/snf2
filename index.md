SnF2 analysis
================
Artemis S. Louyakis
2022-03-01

### Project links

Project description can be found in the publication:
<https://www.mdpi.com/2076-2607/10/9/1691/htm> Project raw reads can be
found on NCBI: <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA855131>

### Metatranscriptome Analysis

``` r
Sys.time()
```

    ## [1] "2022-10-03 17:57:05 EDT"

``` r
library(tidyr)
library(tidyverse)
library(edgeR)
library(trinotateR)
library(NOISeq)
library(manhattanly)
library(cowplot)
library(wesanderson)
library(vegan)
library(corrplot)
library(ggtern)
library(pheatmap)
# library(EnhancedVolcano) ## load separately - conflict in list
# library(beepr) ## alarm added to code chunks that may take longer than others
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 17:57:07 EDT"

``` r
dir.create("data")
dir.create("tables")
dir.create("noiseq")
dir.create("plots")
```

##### Import data and prepare analysis tables

Raw reads were filter for quality with trimmomatic, host reads and rRNA
were removed, assembled with trinity and annotated with trinotate.

``` r
Sys.time()
```

    ## [1] "2022-10-03 17:57:07 EDT"

``` r
## Create map to uniprot Entry ID on https://www.uniprot.org/uniprotkb?query=* with annotations 
uniprot_map <- read.delim("data/uniprot-filtered-reviewed_yes_2022.tab", header = TRUE)
## KEGG heirarchal database
khier <- read.delim("data/khier_2020_fixed.txt", header = FALSE)
# kegg_annos <- read.delim("data/clean_ko_id.txt", header = FALSE, quote = "") delete?

## Data assembled and aligned with Trinity
tmm_counts <- read.delim("data/raw_rmhost_rsem.TMM.EXPR.matrix", header = TRUE)
counts <- read.delim("data/raw_rmhost_rsem.counts.matrix", header = TRUE)

## Data annotated with trinotate pipeline
annots <- read.delim("data/trinotate_annotation_report_1e3_rmhost_max1.xls", header = TRUE)

## Prepare kegg annotation table for mapping 
kegg_annos <- read.delim("data/clean_ko_id.txt", header = FALSE, quote = "")

## gene lengths were pulled from one of the trinity alignment results files
## awk -F"\t" '{OFS="\t"}{print $1,$3}' 10_raw_rmhost/RSEM.genes.results > data/gene_length.txt
len <- read.delim("data/gene_length.txt")

# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 17:58:28 EDT"

``` r
counts_lc <- counts %>% column_to_rownames("X")
len <- len[match(rownames(counts_lc), len$gene_id, nomatch=0),]

counts_lc <- merge(counts_lc, len, by.x = 0, by.y = "gene_id") %>% column_to_rownames("Row.names")
counts_lc <- counts_lc %>% mutate_each(funs(./(length/1000)), starts_with("X")) %>% select(., -length)

# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 17:58:57 EDT"

``` r
## remove zeroes and merge with annotations
counts_noz <- counts[rowSums(counts[, -1]) > 0,]
anno_counts_noz <- merge(counts_noz, annots, by.x = "X", by.y = "X.gene_id")
## filter out genes with fewer than 10 instances
counts_lt10 <- counts[rowSums(counts[, -1]) > 10,]

## prep normalized raw counts table and remove genes with no alignments
tmm_counts_rn <- tmm_counts %>% 
  column_to_rownames("X")
tmm_counts_rn <- tmm_counts_rn[rowSums(tmm_counts_rn) > 0,]

## isolate blastx and blastp uniprot id
## remove periods
## merge the columns
annots %>% separate(sprot_Top_BLASTX_hit, into = c("uniprotx"), sep = "\\^", extra = "drop", remove = FALSE, fill = "right") -> annos
annos %>% separate(sprot_Top_BLASTP_hit, into = c("uniprotp"), sep = "\\^", extra = "drop", remove = FALSE, fill = "right") -> annos
annos$uniprotx <- gsub("\\.", "", annos$uniprotx)
annos$uniprotp <- gsub("\\.", "", annos$uniprotp)
annos$uniprot <- annos$uniprotx
annos <- annos %>% mutate(uniprot = ifelse(uniprot %in% "", uniprotp, uniprot))

## prepare table for noiseq with genes with uniprot annotation - removes unannotated
uniprot <- annos[,c(1,20)]
uniprot[uniprot==""]<-NA
uniprot <- unique(na.omit(uniprot))
uniprot_counts_noz <- merge(uniprot, counts_lt10, by.x = "X.gene_id", by.y = "X")
uniprot_counts_noz$X.gene_id <- NULL
  colSums(uniprot_counts_noz[,c(2:13)])
```

    ##  X1_raw_rmhost  X4_raw_rmhost  X8_raw_rmhost X10_raw_rmhost X11_raw_rmhost 
    ##        2184909        5612415        5657922       12498130       20076352 
    ## X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost X21_raw_rmhost X23_raw_rmhost 
    ##       16238882       14162381       19591376       15104536       15920571 
    ## X24_raw_rmhost X25_raw_rmhost 
    ##       14147106       10622071

``` r
  uniprot_counts_noz <- aggregate(uniprot_counts_noz[,c(2:13)], by=list(uniprot_counts_noz=uniprot_counts_noz$uniprot), FUN=sum)
  colSums(uniprot_counts_noz[,c(2:13)])
```

    ##  X1_raw_rmhost  X4_raw_rmhost  X8_raw_rmhost X10_raw_rmhost X11_raw_rmhost 
    ##        2184909        5612415        5657922       12498130       20076352 
    ## X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost X21_raw_rmhost X23_raw_rmhost 
    ##       16238882       14162381       19591376       15104536       15920571 
    ## X24_raw_rmhost X25_raw_rmhost 
    ##       14147106       10622071

``` r
  rownames(uniprot_counts_noz) <- uniprot_counts_noz$uniprot
  uniprot_counts_noz$uniprot_counts_noz <- NULL

## prepare table for noiseq with with uniprot entry ids
### split blast annotation columns
annots %>% separate(sprot_Top_BLASTX_hit, c("uniprot_x", NA, NA, "percent_id_x", "evalue_x", "recname_x", "taxonomy_x"), sep = "\\^", remove = FALSE) %>% 
  separate(sprot_Top_BLASTP_hit, c("uniprot_p", NA, NA, "percent_id_p", "evalue_p", "recname_p", "taxonomy_p"), sep = "\\^", remove = FALSE) %>% 
  separate(evalue_x, c(NA, "evalue_x"), sep = "\\:", remove = TRUE) %>%
  separate(evalue_p, c(NA, "evalue_p"), sep = "\\:", remove = TRUE) %>%
  separate(percent_id_x, c("percent_id_x", NA), sep = "\\%", remove = TRUE) %>%
  separate(percent_id_p, c("percent_id_p", NA), sep = "\\%", remove = TRUE) %>% 
  separate(Kegg, c(NA, "kegg1"), sep = ":", remove = FALSE, extra = "merge") %>%
  separate(kegg1, c("kegg1", "kegg2"), sep = ";", remove = TRUE) -> annos_clean
### choose top hit
annos_clean %>% arrange(transcript_id, evalue_x) %>% group_by(transcript_id) %>% slice_head(n = 1) %>% ungroup() -> annos_clean
### combine with updated uniprot annotations for downstream use
updated_annots <- annos_clean %>% na_if(".") %>% mutate(., uniprot_id = coalesce(uniprot_x, uniprot_p)) %>%
  mutate(., taxa = coalesce(taxonomy_x, taxonomy_p)) %>%
  merge(., uniprot_map, by.x = "uniprot_id", by.y = "Entry.name", all.x = TRUE)

## prepare table for noiseq with genes with kegg annotation - removes unannotated
uniprot_all <- merge(uniprot, uniprot_map, by.x = "uniprot", by.y = "Entry.name", all.x = TRUE)
  kegg_only <- uniprot_all[,c(1:2,9)]
  kegg_only <- kegg_only[Reduce(`&`, lapply(kegg_only, function(x) !(is.na(x)|x==""))),]
  kegg_only <- kegg_only %>% mutate(Cross.reference..KEGG. = strsplit(as.character(Cross.reference..KEGG.), ";")) %>% unnest(Cross.reference..KEGG.)
  kegg_only <- merge(kegg_only, kegg_annos, by.x = "Cross.reference..KEGG.", by.y = "V3", all.x = TRUE)
  kegg_only <- merge(kegg_only, counts_lt10, by.x = "X.gene_id", by.y = "X")
  ko_noz_rn <- kegg_only[,c(4,6:17)]
  kegg_noz_rn <- kegg_only[,c(2,6:17)]
  ko_noz_rn <- na.omit(ko_noz_rn)
  kegg_noz_rn <- na.omit(kegg_noz_rn)
  colSums(ko_noz_rn[,c(2:13)])
```

    ##  X1_raw_rmhost  X4_raw_rmhost  X8_raw_rmhost X10_raw_rmhost X11_raw_rmhost 
    ##        1558764        3782582        3658410        8323349       13594964 
    ## X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost X21_raw_rmhost X23_raw_rmhost 
    ##       11403593        9761437       13499101       10062131       11130681 
    ## X24_raw_rmhost X25_raw_rmhost 
    ##        9487881        7047771

``` r
  colSums(kegg_noz_rn[,c(2:13)])
```

    ##  X1_raw_rmhost  X4_raw_rmhost  X8_raw_rmhost X10_raw_rmhost X11_raw_rmhost 
    ##        2063014        5200890        5265549       11446203       18150914 
    ## X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost X21_raw_rmhost X23_raw_rmhost 
    ##       14983957       13274717       18239576       14118287       14665228 
    ## X24_raw_rmhost X25_raw_rmhost 
    ##       12799360        9672401

``` r
  ko_noz_rn <- aggregate(ko_noz_rn[,c(2:13)], by=list(kegg_id=ko_noz_rn[,1]), FUN=sum)
  kegg_noz_rn <- aggregate(kegg_noz_rn[,c(2:13)], by=list(kegg_id=kegg_noz_rn[,1]), FUN=sum)
  rownames(ko_noz_rn) <- ko_noz_rn$kegg_id
  rownames(kegg_noz_rn) <- kegg_noz_rn$kegg_id
  ko_noz_rn$kegg_id <- NULL
  kegg_noz_rn$kegg_id <- NULL
## normalize by kegg values - compare to merging annotations with data normalized by transcript id
ko_noz_rn_cpm <- cpm(ko_noz_rn)
ko_noz_rn_tmm <- as.data.frame(NOISeq::tmm(ko_noz_rn_cpm, long = 1000, lc = 0, k = 0))
save(ko_noz_rn_tmm, file = "data/ko_noz_rn_tmm.RData")
kegg_noz_rn_cpm <- cpm(kegg_noz_rn)
kegg_noz_rn_tmm <- as.data.frame(NOISeq::tmm(kegg_noz_rn_cpm, long = 1000, lc = 0, k = 0))
save(kegg_noz_rn_tmm, file = "data/kegg_noz_rn_tmm.RData")
ko_noz_rn_tmm$ko <- rownames(ko_noz_rn_tmm)
kegg_noz_rn_tmm$ko <- rownames(kegg_noz_rn_tmm)

write.table(ko_noz_rn_tmm, file = "data/ko_noz_rn_tmm.txt", quote = FALSE, sep = "\t")

## Prep files for network and/or line plots
ko_annos <- as.data.frame(unique(kegg_annos[,1:2]))
  ko_annos <- dplyr::rename(ko_annos, ko = V1, kegg_gene = V2)
khier <- dplyr::rename(khier, kegg1 = V1, kegg2 = V2, kegg3 = V3, ko = V4, kegg_description = V5)
ko_hier <- merge(ko_annos[!duplicated(ko_annos$ko), ], khier[!duplicated(khier$ko), ], by = "ko")
  rownames(ko_hier) <- ko_hier$ko
  ko_hier <- as.matrix(ko_hier)

# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:10:20 EDT"

``` r
x <- read_trinotate("data/trinotate_annotation_report_1e3_rmhost_max1.xls")
summary_trinotate(x)
```

    ##                       unique   total
    ## gene_id              2237680 2834845
    ## transcript_id        2685947 2834845
    ## sprot_Top_BLASTX_hit 1325855 1523427
    ## prot_id              1511906 1511906
    ## prot_coords           235069 1511906
    ## gene_ontology_BLASTX   24703 1358924
    ## Kegg                   92207 1273459
    ## eggnog                  8594 1143394
    ## sprot_Top_BLASTP_hit  851585  927594
    ## gene_ontology_BLASTP   18631  818658
    ## TmHMM                 170068  217619
    ## RNAMMER                    0       0
    ## Pfam                       0       0
    ## SignalP                    0       0
    ## gene_ontology_Pfam         0       0
    ## transcript                 0       0
    ## peptide                    0       0

``` r
x1 <- split_GO(x, hit = "gene_ontology_BLASTX")
x2 <- summary_GO(x1)
head(x2)
```

    ##            go           ontology                           name  genes
    ## 1: GO:0005737 cellular_component                      cytoplasm 302561
    ## 2: GO:0005524 molecular_function                    ATP binding 282043
    ## 3: GO:0005886 cellular_component                plasma membrane 171243
    ## 4: GO:0016021 cellular_component integral component of membrane 126051
    ## 5: GO:0046872 molecular_function              metal ion binding 120712
    ## 6: GO:0003677 molecular_function                    DNA binding 102688
    ##    transcripts proteins  total
    ## 1:      374786   243296 402638
    ## 2:      350849   235136 377932
    ## 3:      215113   158746 235337
    ## 4:      158957   120083 173796
    ## 5:      150501   101833 162091
    ## 6:      126067    87323 136343

``` r
attr(x2, "count")
```

    ##              unique annotations
    ## GO            17171     6497436
    ## genes        992576     4920668
    ## transcripts 1250783     6025498
    ## proteins     865343     4080548

``` r
colnames(x)
```

    ##  [1] "gene_id"              "transcript_id"        "sprot_Top_BLASTX_hit"
    ##  [4] "RNAMMER"              "prot_id"              "prot_coords"         
    ##  [7] "sprot_Top_BLASTP_hit" "Pfam"                 "SignalP"             
    ## [10] "TmHMM"                "eggnog"               "Kegg"                
    ## [13] "gene_ontology_BLASTX" "gene_ontology_BLASTP" "gene_ontology_Pfam"  
    ## [16] "transcript"           "peptide"

``` r
# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:11:31 EDT"

``` r
## merge annotations with counts tables - object input order is important for downstream (x = annos and y = counts)
anno_counts <- merge(updated_annots, counts_lc, by.x = "X.gene_id", by.y = 0)
### matrix contains 2 id columns (gene id and transcript id), 25 annotation columns (including the columns made above), and sample columns as anno_counts[,28:ncol(anno_counts)]
## remove rows with no reads aligned
no_aligns <- anno_counts[rowSums(anno_counts[,c(62:ncol(anno_counts))]) == 0,]
anno_counts <- anno_counts[rowSums(anno_counts[,c(62:ncol(anno_counts))]) > 0,]
colnames(anno_counts)
```

    ##  [1] "X.gene_id"                          "uniprot_id"                        
    ##  [3] "transcript_id"                      "sprot_Top_BLASTX_hit"              
    ##  [5] "uniprot_x"                          "percent_id_x"                      
    ##  [7] "evalue_x"                           "recname_x"                         
    ##  [9] "taxonomy_x"                         "RNAMMER"                           
    ## [11] "prot_id"                            "prot_coords"                       
    ## [13] "sprot_Top_BLASTP_hit"               "uniprot_p"                         
    ## [15] "percent_id_p"                       "evalue_p"                          
    ## [17] "recname_p"                          "taxonomy_p"                        
    ## [19] "Pfam"                               "SignalP"                           
    ## [21] "TmHMM"                              "eggnog"                            
    ## [23] "Kegg"                               "kegg1"                             
    ## [25] "kegg2"                              "gene_ontology_BLASTX"              
    ## [27] "gene_ontology_BLASTP"               "gene_ontology_Pfam"                
    ## [29] "transcript"                         "peptide"                           
    ## [31] "taxa"                               "Entry"                             
    ## [33] "Status"                             "Protein.names"                     
    ## [35] "Gene.names"                         "Organism"                          
    ## [37] "Length"                             "Cross.reference..KEGG."            
    ## [39] "Ensembl.transcript"                 "Cross.reference..eggNOG."          
    ## [41] "Cross.reference..OrthoDB."          "Taxonomic.lineage..KINGDOM."       
    ## [43] "Taxonomic.lineage..PHYLUM."         "Taxonomic.lineage..CLASS."         
    ## [45] "Taxonomic.lineage..ORDER."          "Taxonomic.lineage..FAMILY."        
    ## [47] "Taxonomic.lineage..GENUS."          "Taxonomic.lineage..SPECIES."       
    ## [49] "Gene.ontology..biological.process." "Gene.ontology..cellular.component."
    ## [51] "Gene.ontology..GO."                 "Gene.ontology..molecular.function."
    ## [53] "Gene.ontology.IDs"                  "Gene.names...primary.."            
    ## [55] "EnsemblBacteria.transcript"         "EnsemblFungi.transcript"           
    ## [57] "EnsemblMetazoa.transcript"          "EnsemblPlants.transcript"          
    ## [59] "EnsemblProtists.transcript"         "Cross.reference..GeneID."          
    ## [61] "Cross.reference..PATRIC."           "X1_raw_rmhost"                     
    ## [63] "X4_raw_rmhost"                      "X8_raw_rmhost"                     
    ## [65] "X10_raw_rmhost"                     "X11_raw_rmhost"                    
    ## [67] "X12_raw_rmhost"                     "X14_raw_rmhost"                    
    ## [69] "X17_raw_rmhost"                     "X21_raw_rmhost"                    
    ## [71] "X23_raw_rmhost"                     "X24_raw_rmhost"                    
    ## [73] "X25_raw_rmhost"

``` r
# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:12:25 EDT"

``` r
## uniprot - taxa and function
uniprot_counts <- anno_counts[,c(32,62:ncol(anno_counts))] %>% na_if(".") %>%
  drop_na(Entry) %>% group_by(Entry) %>% summarise_if(is.numeric, sum) %>%
  column_to_rownames("Entry")

# beepr::beep(sound = "mario")
```

##### Prepare metadata factors

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:12:39 EDT"

``` r
myfactors <- data.frame("samp" = colnames(counts[2:13]),
                        "sample" = c(1,4,8,10,11,12,14,17,21,23,24,25),
                        "treatment" = c(rep("control",6),rep("stannous",6)),
                        "toothpaste" = c(rep("mfp",6),rep("snf",6)),
                        "panelist" = c("p1","p2","p3","p4","p5","p6","p1","p2","p3","p4","p5","p6"),
                        "long" = c(rep("mono fluorophosphate",6),rep("stannous fluoride",6)))
rownames(myfactors) <- colnames(counts[2:13])
save(myfactors, file = "noiseq/myfactors.RData")

samp <- c(colnames(counts))
treatment <- c("control", "stannous")
toothpaste <- c("mfp", "snf")
panelist <- c("p1","p2","p3","p4","p5","p6")
sample <- c("1","4","8","10","11","12",
           "14","17","21","23","24","25")
mfp <- c("1","4","8","10","11","12")
snf <- c("14","17","21","23","24","25")

write.table(myfactors, "data/metadata.txt", quote = FALSE, sep = "\t")

# beepr::beep(sound = "mario")
```

##### Differential Expression analysis

DE was completed using noiseqbio or noiseq for transcripts, as well as
on aggregated counts by uniprot entry id, uniprot id, and kegg orthology
for before and after treatment, before and after for each individual,
and between individuals.

``` r
# Sys.time()
# 
# ## de for all transcripts
# counts_noz_rn <- counts_noz
#   rownames(counts_noz_rn) <- counts_noz_rn$X
#   counts_noz_rn$X <- NULL
# dds <- NOISeq::readData(data = counts_noz_rn, factors = myfactors)
# save(dds, file ="noiseq/dds_transcripts.RData")
# 
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseqbio(dds, norm = "tmm", k = 0.5, factor = "treatment", conditions = c(s1, s2),
#                              r = 50, filter = 1, nclust = 15, a0per= 0.9, random.seed = 1234, plot = TRUE, lc = 1)
#         save(tmm, file = paste0("noiseq/nsb_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results, row.names = NULL)
#         write.table(tmm@results, paste0("noiseq/nsb_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$trinity_id <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/nsb_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# ## resulted in no de genes
# 
# ## de for uniprot ids
# dds <- NOISeq::readData(data = uniprot_counts_noz, factors = myfactors)
# save(dds, file ="noiseq/dds_uniprot.RData")
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseqbio(dds, norm = "tmm", k = 0.5, factor = "treatment", conditions = c(s1, s2),
#                              r = 50, filter = 1, nclust = 15, a0per= 0.9, random.seed = 54321, plot = TRUE, lc = 1)
#         save(tmm, file = paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results)
#         write.table(tmm@results, paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$uniprot <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# 
# ## de for uniprot ids - Entry 2022.05.09 for enrichment analysis
# dds <- NOISeq::readData(data = uniprot_counts, factors = myfactors)
# save(dds, file = "noiseq/dds_uniprot_entry.RData")
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseqbio(dds, norm = "tmm", k = 0.5, factor = "treatment", conditions = c(s1, s2),
#                              r = 50, filter = 1, nclust = 15, a0per= 0.9, random.seed = 54321, plot = TRUE, lc = 0)
#         save(tmm, file = paste0("noiseq/nsb_uniprot_entry_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results)
#         write.table(tmm@results, paste0("noiseq/nsb_uniprot_entry_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$uniprot <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/nsb_uniprot_entry_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# 
# ## de for panelists using uniprot ids
# dds <- NOISeq::readData(data = uniprot_counts_noz, factors = myfactors)
# save(dds, file ="noiseq/dds_uniprot.RData")
# for(s1 in panelist){
#   for(s2 in panelist){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseqbio(dds, norm = "tmm", k = 0.5, factor = "panelist", conditions = c(s1, s2),
#                              r = 50, filter = 1, nclust = 15, a0per= 0.9, random.seed = 54321, plot = TRUE, lc = 1)
#         save(tmm, file = paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results)
#         write.table(tmm@results, paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$uniprot <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# 
# ## de on before and after on individuals without replication
# for(p in panelist){
#       print(c(p))
#       md_test <- filter(myfactors, panelist == p)
#       up_test <- select(uniprot_counts_noz, as.vector(md_test$samp))
#       dds <- readData(data = up_test, factors = md_test)
#   for(s1 in toothpaste){
#     for(s2 in toothpaste){
#       if(s1!=s2){
#         if(s1<s2){
#           print(c(s1,s2))
#           tmm <- noiseq(dds, norm = "tmm", k = 0.5, factor = "toothpaste", conditions = c(s1, s2), lc = 1, replicates = "no")
#           save(tmm, file = paste0("noiseq/nsb_uniprot_",p,"_",s1,"vs",s2,"_tmm.RData"))
#           tmm_r <- as.data.frame(tmm@results)
#           write.table(tmm@results, paste0("noiseq/nsb_uniprot_",p,"_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                       row.names = TRUE, col.names = NA, quote = FALSE)
#           tmm_r <- na.omit(tmm_r)
#           tmm_r$P <- with(tmm_r, 1-prob)
#           tmm_r$id <- row.names(tmm_r)
#           write.table(tmm_r, paste0("noiseq/nsb_uniprot_",p,"_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                       row.names = TRUE, col.names = NA, quote = FALSE)
#           tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#           tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#           tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#     }}}}}
# 
# ## de for ko ids
# dds <- NOISeq::readData(data = ko_noz_rn, factors = myfactors)
# save(dds, file ="noiseq/dds_ko.RData")
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseqbio(dds, norm = "tmm", k = 0.5, factor = "treatment", conditions = c(s1, s2),
#                              r = 50, filter = 1, nclust = 15, a0per= 0.9, random.seed = 2345, plot = TRUE, lc = 1)
#         save(tmm, file = paste0("noiseq/nsb_ko_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results)
#         write.table(tmm@results, paste0("noiseq/nsb_ko_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$id <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/nsb_ko_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# 
# ## de for ko ids - noiseq
# dds <- NOISeq::readData(data = ko_noz_rn, factors = myfactors)
# save(dds, file ="noiseq/dds_ko.RData")
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseq(dds, norm = "tmm", k = 0.5, factor = "treatment", conditions = c(s1, s2), replicates = "technical")
#         save(tmm, file = paste0("noiseq/ns_ko_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results)
#         write.table(tmm@results, paste0("noiseq/ns_ko_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$id <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/ns_ko_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# 
# ## de for kegg
# dds <- NOISeq::readData(data = kegg_noz_rn, factors = myfactors)
# save(dds, file ="noiseq/dds_kegg.RData")
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm <- noiseqbio(dds, norm = "tmm", k = 0.5, factor = "treatment", conditions = c(s1, s2),
#                              r = 50, filter = 1, nclust = 15, a0per= 0.9, random.seed = 2345, plot = TRUE, lc = 1)
#         save(tmm, file = paste0("noiseq/nsb_kegg_",s1,"vs",s2,"_tmm.RData"))
#         tmm_r <- as.data.frame(tmm@results)
#         write.table(tmm@results, paste0("noiseq/nsb_kegg_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm_r <- na.omit(tmm_r)
#         tmm_r$P <- with(tmm_r, 1-prob)
#         tmm_r$id <- row.names(tmm_r)
#         write.table(tmm_r, paste0("noiseq/nsb_kegg_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                     row.names = TRUE, col.names = NA, quote = FALSE)
#         tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#         tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#         tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#       }}}}
# 
# ## ko before and after for each panelist
# for(p in panelist){
#       print(c(p))
#       md_test <- filter(myfactors, panelist == p)
#       ko_test <- select(ko_noz_rn, as.vector(md_test$samp))
#       dds <- readData(data = ko_test, factors = md_test)
#   for(s1 in toothpaste){
#     for(s2 in toothpaste){
#       if(s1!=s2){
#         if(s1<s2){
#           print(c(s1,s2))
#           tmm <- noiseq(dds, norm = "tmm", k = 0.5, factor = "toothpaste", conditions = c(s1, s2), lc = 1, replicates = "no")
#           save(tmm, file = paste0("noiseq/nsb_ko_",p,"_",s1,"vs",s2,"_tmm.RData"))
#           tmm_r <- as.data.frame(tmm@results)
#           write.table(tmm@results, paste0("noiseq/nsb_ko_",p,"_",s1,"vs",s2,"_tmm.txt"), sep = "\t",
#                       row.names = TRUE, col.names = NA, quote = FALSE)
#           tmm_r <- na.omit(tmm_r)
#           tmm_r$P <- with(tmm_r, 1-prob)
#           tmm_r$id <- row.names(tmm_r)
#           write.table(tmm_r, paste0("noiseq/nsb_ko_",p,"_",s1,"vs",s2,"_tmm_filtered.txt"), sep = "\t",
#                       row.names = TRUE, col.names = NA, quote = FALSE)
#           tmm.deg <- degenes(tmm, q = 0.8, M = NULL)
#           tmm.deg <- degenes(tmm, q = 0.8, M = "up")
#           tmm.deg <- degenes(tmm, q = 0.8, M = "down")
#     }}}}}
# 
# beepr::beep(sound = "mario")
```

##### Outputting plots from DE

Summarize the DE genes from analyses

``` r
## working chunk to look at specific up or down expressed genes
## summaries outputted for each DE analysis above
load("noiseq/nsb_ko_controlvsstannous_tmm.RData")
tmm.deg <- degenes(tmm, q = 0.95, M = NULL)
```

    ## [1] "708 differentially expressed features"

``` r
tmm.deg <- degenes(tmm, q = 0.95, M = "up")
```

    ## [1] "708 differentially expressed features (up in first condition)"

``` r
tmm.deg <- degenes(tmm, q = 0.95, M = "down")
```

    ## [1] "0 differentially expressed features (down in first condition)"

Output M-A plots

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:12:39 EDT"

``` r
## MA plots for DE analyses

load("noiseq/nsb_controlvsstannous_tmm.RData")
DE.plot(tmm, graphic = "MD")
```

![](index_files/figure-gfm/checking%20M-A-1.png)<!-- -->

``` r
tmm_r <- as.data.frame(tmm@results)
head(tmm_r)
```

    ##                        control_mean stannous_mean       theta      prob
    ## TRINITY_DN173563_c0_g1    0.8063871      1.396671 -0.23485339 0.9855454
    ## TRINITY_DN615946_c0_g1    0.8063871      2.300782 -0.37518752 0.9763578
    ## TRINITY_DN338562_c0_g2   15.0207918      5.693387  0.70794701 0.9859535
    ## TRINITY_DN139387_c0_g1    0.8063871      1.002702 -0.08285488 0.9913012
    ## TRINITY_DN800837_c0_g1           NA            NA          NA        NA
    ## TRINITY_DN383060_c2_g1    2.2977800      3.063901 -0.12852053 0.9905702
    ##                            log2FC
    ## TRINITY_DN173563_c0_g1 -0.7924477
    ## TRINITY_DN615946_c0_g1 -1.5125800
    ## TRINITY_DN338562_c0_g2  1.3996017
    ## TRINITY_DN139387_c0_g1 -0.3143480
    ## TRINITY_DN800837_c0_g1         NA
    ## TRINITY_DN383060_c2_g1 -0.4151292

``` r
tmp <- data.frame(log2FC=tmm_r[,5], log_c=log(tmm_r[,1]), log_s=log(tmm_r[,2]))
tmp$A <- rowMeans(tmp[,2:3])
head(tmp)
```

    ##       log2FC      log_c       log_s           A
    ## 1 -0.7924477 -0.2151914 0.334091478  0.05945003
    ## 2 -1.5125800 -0.2151914 0.833249168  0.30902887
    ## 3  1.3996017  2.7094354 1.739305410  2.22437039
    ## 4 -0.3143480 -0.2151914 0.002698018 -0.10624670
    ## 5         NA         NA          NA          NA
    ## 6 -0.4151292  0.8319434 1.119689069  0.97581625

``` r
ggplot(tmp, aes(x = log2FC, y = A)) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-2.png)<!-- -->

``` r
ggplot(tmm_r, aes(x = log2FC, y = log(control_mean))) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-3.png)<!-- -->

``` r
ggplot(tmm_r, aes(x = log2FC, y = log(stannous_mean))) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-4.png)<!-- -->

``` r
ggplot(tmp, aes(x = A, y = log2FC)) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-5.png)<!-- -->

``` r
load("noiseq/nsb_ko_controlvsstannous_tmm.RData")
DE.plot(tmm, graphic = "MD")
```

![](index_files/figure-gfm/checking%20M-A-6.png)<!-- -->

``` r
tmm_r <- as.data.frame(tmm@results)
head(tmm_r)
```

    ##        control_mean stannous_mean       theta      prob      log2FC
    ## K00002     9.630666     0.9906229  0.58929988 0.9538570  3.28122767
    ## K00003  2569.168448  2408.3418782  0.08645252 0.9236563  0.09326128
    ## K00004  1561.633747  1959.4398580 -0.23504052 0.8980809 -0.32738516
    ## K00005  1399.187582  2055.1976373 -0.33586824 0.8641539 -0.55468775
    ## K00007    12.361450    21.7820004 -0.16365636 0.9108817 -0.81728844
    ## K00008   138.725490   108.9052298  0.09605435 0.9234482  0.34915967

``` r
tmp <- data.frame(log2FC=tmm_r[,5], log_c=log(tmm_r[,1]), log_s=log(tmm_r[,2]))
tmp$A <- rowMeans(tmp[,2:3])
head(tmp)
```

    ##        log2FC    log_c        log_s        A
    ## 1  3.28122767 2.264952 -0.009421315 1.127766
    ## 2  0.09326128 7.851338  7.786693772 7.819016
    ## 3 -0.32738516 7.353488  7.580413925 7.466951
    ## 4 -0.55468775 7.243647  7.628127296 7.435887
    ## 5 -0.81728844 2.514583  3.081083957 2.797833
    ## 6  0.34915967 4.932497  4.690478052 4.811488

``` r
ggplot(tmp, aes(x = log2FC, y = A)) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-7.png)<!-- -->

``` r
ggplot(tmm_r, aes(x = log2FC, y = log(control_mean))) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-8.png)<!-- -->

``` r
ggplot(tmm_r, aes(x = log2FC, y = log(stannous_mean))) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-9.png)<!-- -->

``` r
ggplot(tmp, aes(x = A, y = log2FC)) +
  geom_point()
```

![](index_files/figure-gfm/checking%20M-A-10.png)<!-- -->

``` r
# beepr::beep(sound = "mario")
```

Output volcano plots - both interactive and publication quality plots
outputted

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:14:32 EDT"

``` r
## comparing panelists to one another
# for(s1 in panelist){
#   for(s2 in panelist){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm_r <- read.delim(paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
#         tmm_rm <- merge(tmm_r, uniprot_map, by.x = "uniprot", by.y = "Entry.name", all.x=TRUE)
#         tmm_rm$EFFECTSIZE <- tmm_rm$log2FC
#         write.table(tmm_rm, paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
#         volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "uniprot", gene = "Protein.names",
#                              annotation1 = "Cross.reference..KEGG.", annotation2 = "Organism")
#         volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
#                   genomewideline = -log10(1e-2), genomewideline_color = "green", title = paste(s1,"vs",s2))
#         # print(volc_plot)
#         htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_uniprot_",s1,"vs",s2,"_tmm_filtered.html"))
#       }}}}

## treatment by uniprot id
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm_r <- read.delim(paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
#         tmm_rm <- merge(tmm_r, uniprot_map, by.x = "uniprot", by.y = "Entry.name", all.x=TRUE)
#         tmm_rm$EFFECTSIZE <- tmm_rm$log2FC
#         write.table(tmm_rm, paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
#         volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "uniprot", gene = "Protein.names",
#                              annotation1 = "Cross.reference..KEGG.", annotation2 = "Organism")
#         volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
#                   genomewideline = -log10(1e-2), genomewideline_color = "green", title = paste(s1,"vs",s2))
#         # print(volc_plot)
#         htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_uniprot_",s1,"vs",s2,"_tmm_filtered.html"))
#       }}}}

## treatment by uniprot entry id
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm_r <- read.delim(paste0("noiseq/nsb_uniprot_entry_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
#         tmm_rm <- merge(tmm_r, uniprot_map, by.x = "uniprot", by.y = "Entry", all.x=TRUE)
#         tmm_rm$EFFECTSIZE <- tmm_rm$log2FC
#         write.table(tmm_rm, paste0("noiseq/nsb_uniprot_entry_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
#         volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "uniprot", gene = "Protein.names",
#                              annotation1 = "Cross.reference..KEGG.", annotation2 = "Organism")
#         volc_plot <- volcanoly(volc_obj, effect_size_line = c(-2,2), effect_size_line_color = "orange",
#                   genomewideline = -log10(5e-2), genomewideline_color = "green", title = paste(s1,"vs",s2))
#         # print(volc_plot)
#         htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_uniprot_entry_",s1,"vs",s2,"_tmm_filtered.html"))
#       }}}}

## treatment by ko
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm_r <- read.delim(paste0("noiseq/nsb_ko_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
#         tmm_rm <- merge(tmm_r, kegg_annos, by.x = "id", by.y = "V1", all.x=TRUE)
#         tmm_rm$EFFECTSIZE <- tmm_rm$log2FC
#         write.table(tmm_rm, paste0("noiseq/nsb_ko_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
#         volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "id", gene = "V2", annotation1 = "V3")
#         volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
#                   genomewideline = -log10(1e-2), genomewideline_color = "green", title = paste(s1,"vs",s2))
#         # print(volc_plot)
#         htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_ko_",s1,"vs",s2,"_tmm_filtered.html"))
#       }}}}

## treatment by ko - noiseq
for(s1 in treatment){
  for(s2 in treatment){
    if(s1!=s2){
      if(s1<s2){
        print(c(s1,s2))
        tmm_r <- read.delim(paste0("noiseq/ns_ko_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
        tmm_rm <- merge(tmm_r, kegg_annos, by.x = "id", by.y = "V1", all.x=TRUE)
        tmm_rm$EFFECTSIZE <- tmm_rm$M
        write.table(tmm_rm, paste0("noiseq/ns_ko_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
        volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "id", gene = "V2", annotation1 = "V3")
        volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
                  genomewideline = -log10(1e-2), genomewideline_color = "green", title = paste(s1,"vs",s2))
        print(volc_plot)
        htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_ko_ns_",s1,"vs",s2,"_tmm_filtered.html"))
      }}}}
```

    ## [1] "control"  "stannous"

``` r
## treatment by kegg
# for(s1 in treatment){
#   for(s2 in treatment){
#     if(s1!=s2){
#       if(s1<s2){
#         print(c(s1,s2))
#         tmm_r <- read.delim(paste0("noiseq/nsb_kegg_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
#         tmm_rm <- merge(tmm_r, kegg_annos, by.x = "id", by.y = "V3", all.x=TRUE)
#         tmm_rm$EFFECTSIZE <- tmm_rm$log2FC
#         write.table(tmm_rm, paste0("noiseq/nsb_kegg_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
#         volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "id", gene = "V1", annotation1 = "V2")
#         volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
#                   genomewideline = -log10(1e-2), genomewideline_color = "green", title = paste(s1,"vs",s2))
#         print(volc_plot)
#         htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_kegg_",s1,"vs",s2,"_tmm_filtered.html"))
#       }}}}

## interactive for before and after for each panelist tested - uniprot
# for(p in panelist){
#   print(c(p))
#   for(s1 in toothpaste){
#     for(s2 in toothpaste){
#       if(s1!=s2){
#         if(s1<s2){
#           print(c(s1,s2))
#           tmm_r <- read.delim(paste0("noiseq/nsb_uniprot_",p,"_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
#           tmm_rm <- merge(tmm_r, uniprot_map, by.x = "id", by.y = "Entry.name", all.x=TRUE)
#           tmm_rm$EFFECTSIZE <- tmm_rm$M
#           write.table(tmm_rm, paste0("noiseq/nsb_uniprot_",p,"_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
#           volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "id", gene = "Protein.names",
#                                annotation1 = "Cross.reference..KEGG.", annotation2 = "Organism")
#           volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
#                     genomewideline = -log10(0.2), genomewideline_color = "green", title = paste(s1,"vs",s2))
#           print(volc_plot)
#           htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_uniprot_",p,"_",s1,"vs",s2,"_tmm_filtered.html"))
#         }}}}}

## interactive for before and after for each panelist tested - ko
for(p in panelist){
  print(c(p))
  for(s1 in toothpaste){
    for(s2 in toothpaste){
      if(s1!=s2){
        if(s1<s2){
          print(c(s1,s2))
          tmm_r <- read.delim(paste0("noiseq/nsb_ko_",p,"_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
          tmm_rm <- merge(tmm_r, ko_hier, by.x = "id", by.y = "ko", all.x=TRUE)
          tmm_rm$EFFECTSIZE <- tmm_rm$M
          write.table(tmm_rm, paste0("noiseq/nsb_ko_",p,"_",s1,"vs",s2,"_tmm_filtered_annot.txt"), sep = "\t", row.names = F, quote = FALSE)
          volc_obj <- volcanor(tmm_rm, p = "P", effect_size = "EFFECTSIZE", snp = "id", gene = "kegg_gene",
                               annotation1 = "kegg2", annotation2 = "kegg3")
          volc_plot <- volcanoly(volc_obj, effect_size_line = c(-5,5), effect_size_line_color = "orange",
                    genomewideline = -log10(0.2), genomewideline_color = "green", title = paste(s1,"vs",s2))
          print(volc_plot)
          htmlwidgets::saveWidget(volc_plot, paste0("plots/volcano_ko_",p,"_",s1,"vs",s2,"_tmm_filtered.html"))
        }}}}}
```

    ## [1] "p1"
    ## [1] "mfp" "snf"
    ## [1] "p2"
    ## [1] "mfp" "snf"
    ## [1] "p3"
    ## [1] "mfp" "snf"
    ## [1] "p4"
    ## [1] "mfp" "snf"
    ## [1] "p5"
    ## [1] "mfp" "snf"
    ## [1] "p6"
    ## [1] "mfp" "snf"

``` r
# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:15:50 EDT"

``` r
# Volcano plots for publishing
# for(p in panelist){
#   print(c(p))
#   for(s1 in toothpaste){
#     for(s2 in toothpaste){
#       if(s1!=s2){
#         if(s1<s2){
#           print(c(s1,s2))
#           tmm_r <- read.delim(paste0("noiseq/nsb_ko_",p,"_",s1,"vs",s2,"_tmm_filtered.txt"), header=TRUE, row.names = NULL, com='', check.names=F)
#           plot(tmm_r$M, -log10(1-tmm_r$prob), pch=20, xlim=c(-9,9), ylim=c(0,6), main = paste(p, s1, " vs ", s2),
#             xlab = "log2FC", ylab = "-log10(p-value)",
#               col = ifelse(abs(tmm_r$M)>5 & tmm_r$prob>0.999,"red",
#                 ifelse(abs(tmm_r$M)>5,"orange",
#                   ifelse(tmm_r$prob>0.99,"green","black"))))
# }}}}}

## package is conflicting with another loaded package
## if error occurs, suggest restarting rstudio
library(EnhancedVolcano)

load("noiseq/ns_ko_controlvsstannous_tmm.RData")
tmm_r <- as.data.frame(tmm@results)
tmm_r$P <- 1-tmm_r$prob
tmm_r$fc_flip <- 0-tmm_r$M
colnames(tmm_r)
```

    ## [1] "control_mean"  "stannous_mean" "M"             "D"            
    ## [5] "prob"          "ranking"       "P"             "fc_flip"

``` r
EnhancedVolcano(tmm_r, lab = rownames(tmm_r), x = 'fc_flip', y = 'P',axisLabSize = 10,
                pCutoff = 0.2, FCcutoff = 3, pointSize = 1.0, colAlpha = 0.3, labSize = 2.0,
                col=c('grey37', '#6699CC', 'orange', 'red3'), ylim = c(0, 2.5),
                legendPosition = "none", ylab = bquote(~-Log[10] ~ '(' ~ italic(1-Prob) ~ ')' ),
                # legendPosition = 'bottom', legendLabSize = 8, legendIconSize = 2,
                # legendLabels = c("NS", expression(Log[2] ~ FC > 3), "p-value < 0.05", "both"),
                drawConnectors = TRUE, widthConnectors = 0.25, gridlines.minor = FALSE, gridlines.major = FALSE,
                vline = 0, vlineType = "solid", vlineCol = "grey27", vlineWidth = 0.2,
                # subtitle = expression(MFP ~ vs ~ SnF[2]), subtitleLabSize = 10, captionLabSize = 10,
                subtitle = paste0("← higher in MFP | higher in SnF2 →"), subtitleLabSize = 12,
                caption = NULL, title = NULL) + theme(plot.subtitle = element_text(hjust = 0.5))
```

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-1.png)<!-- -->

``` r
ggsave(paste0("plots/figure08_volcano_enhanced_ko_mfp-vs-snf_tmm.png"), height = 4.5, width = 4.5, dpi = 600)
ggsave(paste0("plots/figure08_volcano_enhanced_ko_mfp-vs-snf_tmm.pdf"), height = 4.5, width = 4.5, dpi = 600)
ggsave(paste0("plots/figure08_volcano_enhanced_ko_mfp-vs-snf_tmm.tiff"), height = 4.5, width = 4.5, dpi = 600)

# ### for pairs
plots <- list()
for(p in panelist){
  print(c(p))
  tmm_r <- read.delim(paste0("noiseq/nsb_ko_",p,"_mfpvssnf_tmm_filtered.txt"), header=TRUE, row.names = 1, com='', check.names=F)
  tmm_r$fc_flip <- 0-tmm_r$M
  v <- EnhancedVolcano(tmm_r, lab = rownames(tmm_r), x = 'fc_flip', y = 'P', axisLabSize = 10,
                            pCutoff = 0.1, FCcutoff = 3, pointSize = 1.0, colAlpha = 0.3, labSize = 2.0,
                            col=c('grey37', '#6699CC', 'orange', 'red3'), ylim = c(0, 6.2),
                            legendPosition = 'none', # for cowplot
                            # legendPosition = 'bottom', legendLabSize = 8, legendIconSize = 2, # for individual plots
                            legendLabels = c("NS", expression(Log[2] ~ FC > 3), "p-value < 0.05", "both"),
                            drawConnectors = TRUE, widthConnectors = 0.25, gridlines.minor = FALSE, gridlines.major = FALSE,
                            vline = 0, vlineType = "solid", vlineCol = "grey27", vlineWidth = 0.2,
                            subtitle = paste0("← higher in MFP | ",p," | higher in SnF2 →"), subtitleLabSize = 8, caption = NULL, # for cowplot
                            # subtitle = paste0("Panelist ", nm), subtitleLabSize = 10, captionLabSize = 10, # for individual plots
                            title = NULL, ylab = bquote(~-Log[10] ~ italic(1-Prob))) + theme(plot.subtitle = element_text(hjust = 0.5))
  print(v)
  ggsave(paste0("plots/volcano_enhanced_ko_",p,"_mfpvssnf_tmm_filtered_flip.png"), height = 6, width = 4.5, dpi = 600)
  plots[[p]] <- v
}
```

    ## [1] "p1"

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-2.png)<!-- -->

    ## [1] "p2"

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-3.png)<!-- -->

    ## [1] "p3"

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-4.png)<!-- -->

    ## [1] "p4"

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-5.png)<!-- -->

    ## [1] "p5"

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-6.png)<!-- -->

    ## [1] "p6"

![](index_files/figure-gfm/volcano%20plots%20for%20presentation-7.png)<!-- -->

``` r
##
# plots[["p1"]]
# plots[["p2"]]
# plots[["p3"]]
# plots[["p4"]]
# plots[["p5"]]
# plots[["p6"]]

v1 <- cowplot::plot_grid(plots[["p1"]], plots[["p2"]], plots[["p3"]], plots[["p4"]], plots[["p5"]], plots[["p6"]],
                   ncol = 2, labels = "AUTO")
cowplot::save_plot("plots/figureS2_volcano_pairs_flip.png", v1, base_height = 9, base_width = 7)
cowplot::save_plot("plots/figureS2_volcano_pairs_flip.pdf", v1, base_height = 9, base_width = 7)
cowplot::save_plot("plots/figureS2_volcano_pairs_flip.tiff", v1, base_height = 9, base_width = 7)

# beepr::beep(sound = "mario")
```

Output noiseq plots and data on addition comparisons

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:16:45 EDT"

``` r
for(s1 in treatment){
  for(s2 in treatment){
    if(s1!=s2){
      if(s1<s2){
        print(c(s1,s2))
        load(file = paste0("noiseq/nsb_uniprot_",s1,"vs",s2,"_tmm.RData"))
        DE.plot(tmm, q = 0.99, graphic = "expr", log.scale = TRUE)
        DE.plot(tmm, q = 0.95, graphic = "MD")
}}}}
```

    ## [1] "control"  "stannous"

![](index_files/figure-gfm/noiseq%20plots-1.png)<!-- -->

    ## [1] "0 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-2.png)<!-- -->

    ## [1] "41017 differentially expressed features"

``` r
for(p in panelist){
  print(p)
  load(file = paste0("noiseq/nsb_uniprot_",p,"_mfpvssnf_tmm.RData"))
  DE.plot(tmm, q = 0.9, graphic = "expr", log.scale = TRUE)
  DE.plot(tmm, q = 0.8, graphic = "MD")
}
```

    ## [1] "p1"

![](index_files/figure-gfm/noiseq%20plots-3.png)<!-- -->

    ## [1] "8794 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-4.png)<!-- -->

    ## [1] "15590 differentially expressed features"
    ## [1] "p2"

![](index_files/figure-gfm/noiseq%20plots-5.png)<!-- -->

    ## [1] "4733 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-6.png)<!-- -->

    ## [1] "11546 differentially expressed features"
    ## [1] "p3"

![](index_files/figure-gfm/noiseq%20plots-7.png)<!-- -->

    ## [1] "13607 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-8.png)<!-- -->

    ## [1] "23163 differentially expressed features"
    ## [1] "p4"

![](index_files/figure-gfm/noiseq%20plots-9.png)<!-- -->

    ## [1] "9335 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-10.png)<!-- -->

    ## [1] "17374 differentially expressed features"
    ## [1] "p5"

![](index_files/figure-gfm/noiseq%20plots-11.png)<!-- -->

    ## [1] "9788 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-12.png)<!-- -->

    ## [1] "18324 differentially expressed features"
    ## [1] "p6"

![](index_files/figure-gfm/noiseq%20plots-13.png)<!-- -->

    ## [1] "16901 differentially expressed features"

![](index_files/figure-gfm/noiseq%20plots-14.png)<!-- -->

    ## [1] "26059 differentially expressed features"

``` r
counts_noz_rn <- counts_noz
# counts_noz_rn <- counts_lt10
rownames(counts_noz_rn) <- counts_noz_rn$X
counts_noz_rn$X <- NULL
tpmMatrix <- cpm(counts_noz_rn)
rnaseqMatrix <- as.matrix(tpmMatrix)
rnaseqMatrix <- round(rnaseqMatrix)
exp_study <- DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study <- calcNormFactors(exp_study, method = "TMM")
exp_study$samples$eff.lib.size <- exp_study$samples$lib.size * exp_study$samples$norm.factors
tmm_table <- cpm(exp_study)
mydata <- readData(tmm_table, myfactors)

myPCA <- dat(mydata, type = "PCA")
png("plots/pca_toothpaste_noiseq.png", width = 750, height = 750)
explo.plot(myPCA, factor = "toothpaste")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("plots/pca_panelist_noiseq.png", width = 750, height = 750)
explo.plot(myPCA, factor = "panelist")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
explo.plot(myPCA, factor = "toothpaste")
```

![](index_files/figure-gfm/noiseq%20plots-15.png)<!-- -->

``` r
explo.plot(myPCA, factor = "panelist")
```

![](index_files/figure-gfm/noiseq%20plots-16.png)<!-- -->

``` r
explo.plot(myPCA, factor = "treatment")
```

![](index_files/figure-gfm/noiseq%20plots-17.png)<!-- -->

``` r
explo.plot(myPCA, factor = "sample")
```

![](index_files/figure-gfm/noiseq%20plots-18.png)<!-- -->

``` r
#### on ko
ko_noz_rn_tmm$ko <- NULL
mydata <- readData(ko_noz_rn_tmm, myfactors)
myPCA <- dat(mydata, type = "PCA")
png("plots/pca_toothpaste_noiseq_ko.png", width = 750, height = 750)
explo.plot(myPCA, factor = "toothpaste")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("plots/pca_panelist_noiseq_ko.png", width = 750, height = 750)
explo.plot(myPCA, factor = "panelist")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# beepr::beep(sound = "mario")
```

Conduct ordination and correlations on normalized transcripts

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:17:19 EDT"

``` r
ord1 <- metaMDS(t(tmm_counts_rn), iweigh = 1, distance = "bray")
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.08114975 
    ## Run 1 stress 0.2280596 
    ## Run 2 stress 0.08114944 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001424118  max resid 0.002912781 
    ## ... Similar to previous best
    ## Run 3 stress 0.1234227 
    ## Run 4 stress 0.08114949 
    ## ... Procrustes: rmse 0.001226099  max resid 0.002507747 
    ## ... Similar to previous best
    ## Run 5 stress 0.08114955 
    ## ... Procrustes: rmse 0.001284299  max resid 0.002626271 
    ## ... Similar to previous best
    ## Run 6 stress 0.08114978 
    ## ... Procrustes: rmse 0.00145252  max resid 0.002970944 
    ## ... Similar to previous best
    ## Run 7 stress 0.0811496 
    ## ... Procrustes: rmse 0.001323305  max resid 0.002706456 
    ## ... Similar to previous best
    ## Run 8 stress 0.08114963 
    ## ... Procrustes: rmse 0.0001762102  max resid 0.0003598068 
    ## ... Similar to previous best
    ## Run 9 stress 0.08114985 
    ## ... Procrustes: rmse 0.0002932862  max resid 0.0005983464 
    ## ... Similar to previous best
    ## Run 10 stress 0.08114986 
    ## ... Procrustes: rmse 0.00150894  max resid 0.003086276 
    ## ... Similar to previous best
    ## Run 11 stress 0.09222023 
    ## Run 12 stress 0.09221887 
    ## Run 13 stress 0.09221965 
    ## Run 14 stress 0.08114969 
    ## ... Procrustes: rmse 0.000208131  max resid 0.0004254713 
    ## ... Similar to previous best
    ## Run 15 stress 0.08114936 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001086785  max resid 0.002222479 
    ## ... Similar to previous best
    ## Run 16 stress 0.09221901 
    ## Run 17 stress 0.09221948 
    ## Run 18 stress 0.08114965 
    ## ... Procrustes: rmse 0.000273851  max resid 0.0005600595 
    ## ... Similar to previous best
    ## Run 19 stress 0.08114945 
    ## ... Procrustes: rmse 0.0001029609  max resid 0.0002106724 
    ## ... Similar to previous best
    ## Run 20 stress 0.08114929 
    ## ... New best solution
    ## ... Procrustes: rmse 9.002528e-05  max resid 0.0001839844 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
ord1_df <- as.data.frame(ord1$points) %>% cbind(myfactors, .)

ggplot(ord1_df, aes(x = MDS1, y = MDS2, color = toothpaste)) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  geom_point(size = 3) +
  theme_bw() + 
  labs(x = "NMDS1", y = "NMDS2") + 
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c(expression("MFP"), expression("SnF"["2"]))) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
      legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  geom_path(aes(group = panelist), color = "grey27", arrow = arrow(length = unit(0.15, "cm"), type = "closed"))
```

![](index_files/figure-gfm/ordination%20on%20transcripts-1.png)<!-- -->

``` r
ggsave("plots/figure07_metamds_bray.png", dpi = 600, height = 4.5, width = 5, units = "in")
ggsave("plots/figure07_metamds_bray.pdf", dpi = 600, height = 4.5, width = 5, units = "in")
ggsave("plots/figure07_metamds_bray.tiff", dpi = 600, height = 4.5, width = 5, units = "in")

## stats on first dimension
summary(lm(MDS1 ~ treatment, data = ord1_df))
```

    ## 
    ## Call:
    ## lm(formula = MDS1 ~ treatment, data = ord1_df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.37698 -0.25935  0.08308  0.20642  1.05168 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)       -0.09151    0.27951  -0.327    0.750
    ## treatmentstannous  0.18302    0.39529   0.463    0.653
    ## 
    ## Residual standard error: 0.6847 on 10 degrees of freedom
    ## Multiple R-squared:  0.02099,    Adjusted R-squared:  -0.07691 
    ## F-statistic: 0.2144 on 1 and 10 DF,  p-value: 0.6533

``` r
## stats on second dimension
summary(lm(MDS2 ~ treatment, data = ord1_df))
```

    ## 
    ## Call:
    ## lm(formula = MDS2 ~ treatment, data = ord1_df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.3477 -0.1349 -0.1004  0.2285  0.5133 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)       -0.04611    0.11819  -0.390    0.705
    ## treatmentstannous  0.09222    0.16714   0.552    0.593
    ## 
    ## Residual standard error: 0.2895 on 10 degrees of freedom
    ## Multiple R-squared:  0.02955,    Adjusted R-squared:  -0.0675 
    ## F-statistic: 0.3045 on 1 and 10 DF,  p-value: 0.5932

``` r
## pearson and spearman correlations
pear_corr <- cor(tmm_counts_rn, use="all.obs", method="pearson")
corrplot(pear_corr) %>%
  corrRect(index = c(1, 7, 12), col = "black")
```

![](index_files/figure-gfm/ordination%20on%20transcripts-2.png)<!-- -->

``` r
ggsave("plots/pearson_correlation.jpeg", dpi = 600)
spear_corr <- cor(tmm_counts_rn, use="all.obs", method="spearman")
corrplot(spear_corr) %>%
  corrRect(index = c(1, 7, 12), col = "black")
```

![](index_files/figure-gfm/ordination%20on%20transcripts-3.png)<!-- -->

``` r
ggsave("plots/spearman_correlation.jpeg", dpi = 600)
cor.mtest(tmm_counts_rn)
```

    ## $p
    ##                X1_raw_rmhost X4_raw_rmhost X8_raw_rmhost X10_raw_rmhost
    ## X1_raw_rmhost              0             0             0              0
    ## X4_raw_rmhost              0             0             0              0
    ## X8_raw_rmhost              0             0             0              0
    ## X10_raw_rmhost             0             0             0              0
    ## X11_raw_rmhost             0             0             0              0
    ## X12_raw_rmhost             0             0             0              0
    ## X14_raw_rmhost             0             0             0              0
    ## X17_raw_rmhost             0             0             0              0
    ## X21_raw_rmhost             0             0             0              0
    ## X23_raw_rmhost             0             0             0              0
    ## X24_raw_rmhost             0             0             0              0
    ## X25_raw_rmhost             0             0             0              0
    ##                X11_raw_rmhost X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost
    ## X1_raw_rmhost               0              0              0              0
    ## X4_raw_rmhost               0              0              0              0
    ## X8_raw_rmhost               0              0              0              0
    ## X10_raw_rmhost              0              0              0              0
    ## X11_raw_rmhost              0              0              0              0
    ## X12_raw_rmhost              0              0              0              0
    ## X14_raw_rmhost              0              0              0              0
    ## X17_raw_rmhost              0              0              0              0
    ## X21_raw_rmhost              0              0              0              0
    ## X23_raw_rmhost              0              0              0              0
    ## X24_raw_rmhost              0              0              0              0
    ## X25_raw_rmhost              0              0              0              0
    ##                X21_raw_rmhost X23_raw_rmhost X24_raw_rmhost X25_raw_rmhost
    ## X1_raw_rmhost               0              0              0              0
    ## X4_raw_rmhost               0              0              0              0
    ## X8_raw_rmhost               0              0              0              0
    ## X10_raw_rmhost              0              0              0              0
    ## X11_raw_rmhost              0              0              0              0
    ## X12_raw_rmhost              0              0              0              0
    ## X14_raw_rmhost              0              0              0              0
    ## X17_raw_rmhost              0              0              0              0
    ## X21_raw_rmhost              0              0              0              0
    ## X23_raw_rmhost              0              0              0              0
    ## X24_raw_rmhost              0              0              0              0
    ## X25_raw_rmhost              0              0              0              0
    ## 
    ## $lowCI
    ##                X1_raw_rmhost X4_raw_rmhost X8_raw_rmhost X10_raw_rmhost
    ## X1_raw_rmhost      1.0000000     0.4833208     0.4150790      0.5876817
    ## X4_raw_rmhost      0.4833208     1.0000000     0.7176840      0.7224499
    ## X8_raw_rmhost      0.4150790     0.7176840     1.0000000      0.5195249
    ## X10_raw_rmhost     0.5876817     0.7224499     0.5195249      1.0000000
    ## X11_raw_rmhost     0.5181003     0.6219971     0.4062432      0.7009428
    ## X12_raw_rmhost     0.7095752     0.5837356     0.4826925      0.6163320
    ## X14_raw_rmhost     0.7695052     0.5468514     0.4933480      0.6811750
    ## X17_raw_rmhost     0.5134083     0.7202707     0.5243191      0.5612390
    ## X21_raw_rmhost     0.1511281     0.4953446     0.7653304      0.3061133
    ## X23_raw_rmhost     0.4724954     0.6924364     0.4143240      0.7622395
    ## X24_raw_rmhost     0.5058635     0.4510717     0.3120555      0.5211096
    ## X25_raw_rmhost     0.5778434     0.6210954     0.4683456      0.6669526
    ##                X11_raw_rmhost X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost
    ## X1_raw_rmhost       0.5181003      0.7095752      0.7695052      0.5134083
    ## X4_raw_rmhost       0.6219971      0.5837356      0.5468514      0.7202707
    ## X8_raw_rmhost       0.4062432      0.4826925      0.4933480      0.5243191
    ## X10_raw_rmhost      0.7009428      0.6163320      0.6811750      0.5612390
    ## X11_raw_rmhost      1.0000000      0.7607138      0.6438082      0.6709127
    ## X12_raw_rmhost      0.7607138      1.0000000      0.7534785      0.6565772
    ## X14_raw_rmhost      0.6438082      0.7534785      1.0000000      0.6315040
    ## X17_raw_rmhost      0.6709127      0.6565772      0.6315040      1.0000000
    ## X21_raw_rmhost      0.2495167      0.2510800      0.1989576      0.3450359
    ## X23_raw_rmhost      0.8353080      0.7483401      0.5458344      0.5657405
    ## X24_raw_rmhost      0.8991529      0.8000951      0.6233349      0.6921595
    ## X25_raw_rmhost      0.8408957      0.7934247      0.7339306      0.7279383
    ##                X21_raw_rmhost X23_raw_rmhost X24_raw_rmhost X25_raw_rmhost
    ## X1_raw_rmhost       0.1511281      0.4724954      0.5058635      0.5778434
    ## X4_raw_rmhost       0.4953446      0.6924364      0.4510717      0.6210954
    ## X8_raw_rmhost       0.7653304      0.4143240      0.3120555      0.4683456
    ## X10_raw_rmhost      0.3061133      0.7622395      0.5211096      0.6669526
    ## X11_raw_rmhost      0.2495167      0.8353080      0.8991529      0.8408957
    ## X12_raw_rmhost      0.2510800      0.7483401      0.8000951      0.7934247
    ## X14_raw_rmhost      0.1989576      0.5458344      0.6233349      0.7339306
    ## X17_raw_rmhost      0.3450359      0.5657405      0.6921595      0.7279383
    ## X21_raw_rmhost      1.0000000      0.2551312      0.1858287      0.2870133
    ## X23_raw_rmhost      0.2551312      1.0000000      0.7168661      0.6971543
    ## X24_raw_rmhost      0.1858287      0.7168661      1.0000000      0.8436734
    ## X25_raw_rmhost      0.2870133      0.6971543      0.8436734      1.0000000
    ## 
    ## $uppCI
    ##                X1_raw_rmhost X4_raw_rmhost X8_raw_rmhost X10_raw_rmhost
    ## X1_raw_rmhost      1.0000000     0.4854252     0.4173521      0.5894787
    ## X4_raw_rmhost      0.4854252     1.0000000     0.7190147      0.7237617
    ## X8_raw_rmhost      0.4173521     0.7190147     1.0000000      0.5215294
    ## X10_raw_rmhost     0.5894787     0.7237617     0.5215294      1.0000000
    ## X11_raw_rmhost     0.5201089     0.6236799     0.4085364      0.7023387
    ## X12_raw_rmhost     0.7109376     0.5855453     0.4847986      0.6180341
    ## X14_raw_rmhost     0.7706243     0.5487757     0.4954255      0.6826460
    ## X17_raw_rmhost     0.5154301     0.7215912     0.5263098      0.5631195
    ## X21_raw_rmhost     0.1538136     0.4974166     0.7664671      0.3086030
    ## X23_raw_rmhost     0.4746282     0.6938648     0.4165988      0.7633891
    ## X24_raw_rmhost     0.5079065     0.4532590     0.3145351      0.5231095
    ## X25_raw_rmhost     0.5796719     0.6227814     0.4704892      0.6684762
    ##                X11_raw_rmhost X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost
    ## X1_raw_rmhost       0.5201089      0.7109376      0.7706243      0.5154301
    ## X4_raw_rmhost       0.6236799      0.5855453      0.5487757      0.7215912
    ## X8_raw_rmhost       0.4085364      0.4847986      0.4954255      0.5263098
    ## X10_raw_rmhost      0.7023387      0.6180341      0.6826460      0.5631195
    ## X11_raw_rmhost      1.0000000      0.7618698      0.6454151      0.6724218
    ## X12_raw_rmhost      0.7618698      1.0000000      0.7546646      0.6581386
    ## X14_raw_rmhost      0.6454151      0.7546646      1.0000000      0.6331541
    ## X17_raw_rmhost      0.6724218      0.6581386      0.6331541      1.0000000
    ## X21_raw_rmhost      0.2520932      0.2536544      0.2015968      0.3474557
    ## X23_raw_rmhost      0.8361372      0.7495474      0.5477618      0.5676070
    ## X24_raw_rmhost      0.8996782      0.8010823      0.6250131      0.6935890
    ## X25_raw_rmhost      0.8416992      0.7944410      0.7351965      0.7292282
    ##                X21_raw_rmhost X23_raw_rmhost X24_raw_rmhost X25_raw_rmhost
    ## X1_raw_rmhost       0.1538136      0.4746282      0.5079065      0.5796719
    ## X4_raw_rmhost       0.4974166      0.6938648      0.4532590      0.6227814
    ## X8_raw_rmhost       0.7664671      0.4165988      0.3145351      0.4704892
    ## X10_raw_rmhost      0.3086030      0.7633891      0.5231095      0.6684762
    ## X11_raw_rmhost      0.2520932      0.8361372      0.8996782      0.8416992
    ## X12_raw_rmhost      0.2536544      0.7495474      0.8010823      0.7944410
    ## X14_raw_rmhost      0.2015968      0.5477618      0.6250131      0.7351965
    ## X17_raw_rmhost      0.3474557      0.5676070      0.6935890      0.7292282
    ## X21_raw_rmhost      1.0000000      0.2576998      0.1884819      0.2895343
    ## X23_raw_rmhost      0.2576998      1.0000000      0.7182000      0.6985647
    ## X24_raw_rmhost      0.1884819      0.7182000      1.0000000      0.8444640
    ## X25_raw_rmhost      0.2895343      0.6985647      0.8444640      1.0000000

``` r
# beepr::beep(sound = "mario")
```

Conduct ordination and correlations on normalized kegg orthology
aggregated counts

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:19:02 EDT"

``` r
ord1 <- metaMDS(t(ko_noz_rn_tmm), iweigh = 1, distance = "bray")
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.07181471 
    ## Run 1 stress 0.3359236 
    ## Run 2 stress 0.07216544 
    ## ... Procrustes: rmse 0.02153445  max resid 0.05978852 
    ## Run 3 stress 0.08998867 
    ## Run 4 stress 0.07216543 
    ## ... Procrustes: rmse 0.02154708  max resid 0.05983643 
    ## Run 5 stress 0.1019058 
    ## Run 6 stress 0.0899887 
    ## Run 7 stress 0.07181468 
    ## ... New best solution
    ## ... Procrustes: rmse 4.691774e-05  max resid 8.609884e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.07181467 
    ## ... New best solution
    ## ... Procrustes: rmse 1.259845e-05  max resid 1.849436e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.1019057 
    ## Run 10 stress 0.07181471 
    ## ... Procrustes: rmse 3.836855e-05  max resid 6.347807e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.1191731 
    ## Run 12 stress 0.07216544 
    ## ... Procrustes: rmse 0.02151919  max resid 0.05973914 
    ## Run 13 stress 0.07181467 
    ## ... New best solution
    ## ... Procrustes: rmse 7.653031e-06  max resid 1.406975e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.07181465 
    ## ... New best solution
    ## ... Procrustes: rmse 7.36858e-05  max resid 0.0001379981 
    ## ... Similar to previous best
    ## Run 15 stress 0.08998871 
    ## Run 16 stress 0.07216548 
    ## ... Procrustes: rmse 0.02144061  max resid 0.05945237 
    ## Run 17 stress 0.1019057 
    ## Run 18 stress 0.1019057 
    ## Run 19 stress 0.1019057 
    ## Run 20 stress 0.07216547 
    ## ... Procrustes: rmse 0.02144624  max resid 0.0594741 
    ## *** Solution reached

``` r
ord1_df <- as.data.frame(ord1$points) %>% cbind(myfactors, .)

ggplot(ord1_df, aes(x = MDS1, y = MDS2, color = toothpaste)) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  geom_point(size = 3) +
  theme_bw() + 
  labs(x = "NMDS1", y = "NMDS2") + #ggtitle("bray") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
      legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  geom_path(aes(group = panelist), color = "grey27", arrow = arrow(length = unit(0.15, "cm"), type = "closed"))
```

![](index_files/figure-gfm/ordination%20on%20ko-1.png)<!-- -->

``` r
ggsave("plots/figureS1_metamds_bray_ko.png", dpi = 600, height = 4.5, width = 5, units = "in")
ggsave("plots/figureS1_metamds_bray_ko.pdf", dpi = 600, height = 4.5, width = 5, units = "in")
ggsave("plots/figureS1_metamds_bray_ko.tiff", dpi = 600, height = 4.5, width = 5, units = "in")

## stats on first dimension
summary(lm(MDS1 ~ treatment, data = ord1_df))
```

    ## 
    ## Call:
    ## lm(formula = MDS1 ~ treatment, data = ord1_df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.53964 -0.05833  0.02420  0.12680  0.27645 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)       -0.03026    0.09294  -0.326    0.751
    ## treatmentstannous  0.06052    0.13143   0.460    0.655
    ## 
    ## Residual standard error: 0.2276 on 10 degrees of freedom
    ## Multiple R-squared:  0.02076,    Adjusted R-squared:  -0.07716 
    ## F-statistic: 0.212 on 1 and 10 DF,  p-value: 0.655

``` r
## stats on second dimension
summary(lm(MDS2 ~ treatment, data = ord1_df))
```

    ## 
    ## Call:
    ## lm(formula = MDS2 ~ treatment, data = ord1_df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.26027 -0.03609 -0.01135  0.09053  0.16894 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)       -0.01097    0.05167  -0.212    0.836
    ## treatmentstannous  0.02195    0.07307   0.300    0.770
    ## 
    ## Residual standard error: 0.1266 on 10 degrees of freedom
    ## Multiple R-squared:  0.00894,    Adjusted R-squared:  -0.09017 
    ## F-statistic: 0.09021 on 1 and 10 DF,  p-value: 0.7701

``` r
## pearson and spearman correlations
pear_corr <- cor(ko_noz_rn_tmm, use="all.obs", method="pearson")
corrplot(pear_corr) %>%
  corrRect(index = c(1, 7, 12), col = "black")
```

![](index_files/figure-gfm/ordination%20on%20ko-2.png)<!-- -->

``` r
spear_corr <- cor(ko_noz_rn_tmm, use="all.obs", method="spearman")
corrplot(spear_corr) %>%
  corrRect(index = c(1, 7, 12), col = "black")
```

![](index_files/figure-gfm/ordination%20on%20ko-3.png)<!-- -->

``` r
cor.mtest(ko_noz_rn_tmm)
```

    ## $p
    ##                X1_raw_rmhost X4_raw_rmhost X8_raw_rmhost X10_raw_rmhost
    ## X1_raw_rmhost              0             0             0              0
    ## X4_raw_rmhost              0             0             0              0
    ## X8_raw_rmhost              0             0             0              0
    ## X10_raw_rmhost             0             0             0              0
    ## X11_raw_rmhost             0             0             0              0
    ## X12_raw_rmhost             0             0             0              0
    ## X14_raw_rmhost             0             0             0              0
    ## X17_raw_rmhost             0             0             0              0
    ## X21_raw_rmhost             0             0             0              0
    ## X23_raw_rmhost             0             0             0              0
    ## X24_raw_rmhost             0             0             0              0
    ## X25_raw_rmhost             0             0             0              0
    ##                X11_raw_rmhost X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost
    ## X1_raw_rmhost               0              0              0              0
    ## X4_raw_rmhost               0              0              0              0
    ## X8_raw_rmhost               0              0              0              0
    ## X10_raw_rmhost              0              0              0              0
    ## X11_raw_rmhost              0              0              0              0
    ## X12_raw_rmhost              0              0              0              0
    ## X14_raw_rmhost              0              0              0              0
    ## X17_raw_rmhost              0              0              0              0
    ## X21_raw_rmhost              0              0              0              0
    ## X23_raw_rmhost              0              0              0              0
    ## X24_raw_rmhost              0              0              0              0
    ## X25_raw_rmhost              0              0              0              0
    ##                X21_raw_rmhost X23_raw_rmhost X24_raw_rmhost X25_raw_rmhost
    ## X1_raw_rmhost               0              0              0              0
    ## X4_raw_rmhost               0              0              0              0
    ## X8_raw_rmhost               0              0              0              0
    ## X10_raw_rmhost              0              0              0              0
    ## X11_raw_rmhost              0              0              0              0
    ## X12_raw_rmhost              0              0              0              0
    ## X14_raw_rmhost              0              0              0              0
    ## X17_raw_rmhost              0              0              0              0
    ## X21_raw_rmhost              0              0              0              0
    ## X23_raw_rmhost              0              0              0              0
    ## X24_raw_rmhost              0              0              0              0
    ## X25_raw_rmhost              0              0              0              0
    ## 
    ## $lowCI
    ##                X1_raw_rmhost X4_raw_rmhost X8_raw_rmhost X10_raw_rmhost
    ## X1_raw_rmhost      1.0000000     0.8597024     0.8112826      0.7800641
    ## X4_raw_rmhost      0.8597024     1.0000000     0.8879451      0.8203135
    ## X8_raw_rmhost      0.8112826     0.8879451     1.0000000      0.8897190
    ## X10_raw_rmhost     0.7800641     0.8203135     0.8897190      1.0000000
    ## X11_raw_rmhost     0.8399338     0.8892925     0.8325159      0.7574900
    ## X12_raw_rmhost     0.7238406     0.7305001     0.6953608      0.6506036
    ## X14_raw_rmhost     0.8072096     0.8107291     0.7961608      0.7907443
    ## X17_raw_rmhost     0.8843677     0.9711459     0.8649976      0.7906014
    ## X21_raw_rmhost     0.7643464     0.8449357     0.8943342      0.7429785
    ## X23_raw_rmhost     0.9268449     0.9222494     0.8370468      0.8282017
    ## X24_raw_rmhost     0.8050363     0.8905972     0.7988938      0.7419381
    ## X25_raw_rmhost     0.8063505     0.8949933     0.8935737      0.9154746
    ##                X11_raw_rmhost X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost
    ## X1_raw_rmhost       0.8399338      0.7238406      0.8072096      0.8843677
    ## X4_raw_rmhost       0.8892925      0.7305001      0.8107291      0.9711459
    ## X8_raw_rmhost       0.8325159      0.6953608      0.7961608      0.8649976
    ## X10_raw_rmhost      0.7574900      0.6506036      0.7907443      0.7906014
    ## X11_raw_rmhost      1.0000000      0.7466053      0.7828113      0.8872600
    ## X12_raw_rmhost      0.7466053      1.0000000      0.6859117      0.7310199
    ## X14_raw_rmhost      0.7828113      0.6859117      1.0000000      0.8001672
    ## X17_raw_rmhost      0.8872600      0.7310199      0.8001672      1.0000000
    ## X21_raw_rmhost      0.8308887      0.6735431      0.7203332      0.8229369
    ## X23_raw_rmhost      0.9022314      0.7399524      0.8167087      0.9238046
    ## X24_raw_rmhost      0.9514269      0.7465503      0.7656939      0.8858212
    ## X25_raw_rmhost      0.8412224      0.7581453      0.8442965      0.8738350
    ##                X21_raw_rmhost X23_raw_rmhost X24_raw_rmhost X25_raw_rmhost
    ## X1_raw_rmhost       0.7643464      0.9268449      0.8050363      0.8063505
    ## X4_raw_rmhost       0.8449357      0.9222494      0.8905972      0.8949933
    ## X8_raw_rmhost       0.8943342      0.8370468      0.7988938      0.8935737
    ## X10_raw_rmhost      0.7429785      0.8282017      0.7419381      0.9154746
    ## X11_raw_rmhost      0.8308887      0.9022314      0.9514269      0.8412224
    ## X12_raw_rmhost      0.6735431      0.7399524      0.7465503      0.7581453
    ## X14_raw_rmhost      0.7203332      0.8167087      0.7656939      0.8442965
    ## X17_raw_rmhost      0.8229369      0.9238046      0.8858212      0.8738350
    ## X21_raw_rmhost      1.0000000      0.8041085      0.7764145      0.7870871
    ## X23_raw_rmhost      0.8041085      1.0000000      0.8763082      0.8635352
    ## X24_raw_rmhost      0.7764145      0.8763082      1.0000000      0.8515513
    ## X25_raw_rmhost      0.7870871      0.8635352      0.8515513      1.0000000
    ## 
    ## $uppCI
    ##                X1_raw_rmhost X4_raw_rmhost X8_raw_rmhost X10_raw_rmhost
    ## X1_raw_rmhost      1.0000000     0.8709202     0.8260098      0.7969543
    ## X4_raw_rmhost      0.8709202     1.0000000     0.8970298      0.8344003
    ## X8_raw_rmhost      0.8260098     0.8970298     1.0000000      0.8986675
    ## X10_raw_rmhost     0.7969543     0.8344003     0.8986675      1.0000000
    ## X11_raw_rmhost     0.8526070     0.8982738     0.8457272      0.7758957
    ## X12_raw_rmhost     0.7444287     0.7506635     0.7177244      0.6756243
    ## X14_raw_rmhost     0.8222234     0.8254953     0.8119455      0.8069033
    ## X17_raw_rmhost     0.8937260     0.9735796     0.8758203      0.8067703
    ## X21_raw_rmhost     0.7822961     0.8572436     0.9029275      0.7623365
    ## X23_raw_rmhost     0.9328879     0.9286580     0.8499299      0.8417239
    ## X24_raw_rmhost     0.8202025     0.8994783     0.8144887      0.7613637
    ## X25_raw_rmhost     0.8214246     0.9035357     0.9022256      0.9224190
    ##                X11_raw_rmhost X12_raw_rmhost X14_raw_rmhost X17_raw_rmhost
    ## X1_raw_rmhost       0.8526070      0.7444287      0.8222234      0.8937260
    ## X4_raw_rmhost       0.8982738      0.7506635      0.8254953      0.9735796
    ## X8_raw_rmhost       0.8457272      0.7177244      0.8119455      0.8758203
    ## X10_raw_rmhost      0.7758957      0.6756243      0.8069033      0.8067703
    ## X11_raw_rmhost      1.0000000      0.7657269      0.7995143      0.8963971
    ## X12_raw_rmhost      0.7657269      1.0000000      0.7088499      0.7511500
    ## X14_raw_rmhost      0.7995143      0.7088499      1.0000000      0.8156735
    ## X17_raw_rmhost      0.8963971      0.7511500      0.8156735      1.0000000
    ## X21_raw_rmhost      0.8442174      0.6972224      0.7411435      0.8368365
    ## X23_raw_rmhost      0.9102129      0.7595069      0.8310519      0.9300897
    ## X24_raw_rmhost      0.9554863      0.7656755      0.7835535      0.8950684
    ## X25_raw_rmhost      0.8538017      0.7765075      0.8566512      0.8839933
    ##                X21_raw_rmhost X23_raw_rmhost X24_raw_rmhost X25_raw_rmhost
    ## X1_raw_rmhost       0.7822961      0.9328879      0.8202025      0.8214246
    ## X4_raw_rmhost       0.8572436      0.9286580      0.8994783      0.9035357
    ## X8_raw_rmhost       0.9029275      0.8499299      0.8144887      0.9022256
    ## X10_raw_rmhost      0.7623365      0.8417239      0.7613637      0.9224190
    ## X11_raw_rmhost      0.8442174      0.9102129      0.9554863      0.8538017
    ## X12_raw_rmhost      0.6972224      0.7595069      0.7656755      0.7765075
    ## X14_raw_rmhost      0.7411435      0.8310519      0.7835535      0.8566512
    ## X17_raw_rmhost      0.8368365      0.9300897      0.8950684      0.8839933
    ## X21_raw_rmhost      1.0000000      0.8193396      0.7935526      0.8034976
    ## X23_raw_rmhost      0.8193396      1.0000000      0.8862794      0.8744673
    ## X24_raw_rmhost      0.7935526      0.8862794      1.0000000      0.8633730
    ## X25_raw_rmhost      0.8034976      0.8744673      0.8633730      1.0000000

``` r
# beepr::beep(sound = "mario")
```

Conduct pairwise t-test on ordination

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:19:03 EDT"

``` r
## https://agroninfotech.blogspot.com/2020/06/one-way-repeated-measures-anova-in-r.html

summary(aov(ord1_df, formula = MDS1 ~ treatment + Error(panelist/treatment)))
```

    ## 
    ## Error: panelist
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Residuals  5 0.4276 0.08553               
    ## 
    ## Error: panelist:treatment
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## treatment  1 0.01099 0.01099   0.606  0.471
    ## Residuals  5 0.09060 0.01812

``` r
pairwise.t.test(x = ord1_df$MDS1, g = ord1_df$treatment, p.adjust.method = 'bonferroni')
```

    ## 
    ##  Pairwise comparisons using t tests with pooled SD 
    ## 
    ## data:  ord1_df$MDS1 and ord1_df$treatment 
    ## 
    ##          control
    ## stannous 0.66   
    ## 
    ## P value adjustment method: bonferroni

``` r
# beepr::beep(sound = "mario")
```

Construct lists of KOs for pathways and processes of interest

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:19:03 EDT"

``` r
## make lists of KOs to plot 

## use to output lists of ids for each set in the dataframe - using kegg here, but could use any annotation
  ## use complete kegg hierarchy file to ensure capture of duplicate identities
keggs <- khier %>%
  # filter(kegg2 == 'Cellular community - prokaryotes') %>%
  filter(kegg3 == 'Quorum sensing [PATH:ko02024]') #%>%
  # filter(Kegg != '00680')
# paste(shQuote(keggs$ko), collapse=", ") ## single quote
keggs <- khier %>%
  filter(kegg3 == 'Flagellar assembly [PATH:ko02040]')
cat(paste(shQuote(unique(keggs$ko), type="cmd"), collapse=", ")) ## double quote
```

    ## "K02386", "K02387", "K02388", "K02389", "K02390", "K02391", "K02392", "K02393", "K02394", "K02396", "K02397", "K02398", "K02399", "K02400", "K02401", "K02402", "K02403", "K02405", "K02406", "K02407", "K02408", "K02409", "K02410", "K02411", "K02412", "K02413", "K02414", "K02416", "K02417", "K02418", "K02419", "K02420", "K02421", "K02422", "K02423", "K02556", "K02557", "K13820", "K21217", "K21218"

``` r
keggs <- khier %>%
  filter(kegg3 == 'Bacterial motility proteins [BR:ko02035]')
cat(paste(shQuote(unique(keggs$ko), type="cmd"), collapse=", "))
```

    ## "K00575", "K02278", "K02279", "K02280", "K02281", "K02282", "K02283", "K02382", "K02383", "K02384", "K02385", "K02386", "K02387", "K02388", "K02389", "K02390", "K02391", "K02392", "K02393", "K02394", "K02395", "K02396", "K02397", "K02398", "K02399", "K02400", "K02401", "K02402", "K02403", "K02404", "K02405", "K02406", "K02407", "K02408", "K02409", "K02410", "K02411", "K02412", "K02413", "K02414", "K02415", "K02416", "K02417", "K02418", "K02419", "K02420", "K02421", "K02422", "K02423", "K02424", "K02425", "K02556", "K02557", "K02650", "K02651", "K02652", "K02653", "K02654", "K02655", "K02656", "K02657", "K02658", "K02659", "K02660", "K02661", "K02662", "K02663", "K02664", "K02665", "K02666", "K02667", "K02668", "K02669", "K02670", "K02671", "K02672", "K02673", "K02674", "K02675", "K02676", "K03406", "K03407", "K03408", "K03409", "K03410", "K03411", "K03412", "K03413", "K03414", "K03415", "K03516", "K03776", "K04562", "K05874", "K05875", "K05876", "K05877", "K06595", "K06596", "K06597", "K06598", "K06599", "K06600", "K06601", "K06602", "K06603", "K07324", "K07325", "K07327", "K07328", "K07329", "K07330", "K07331", "K07332", "K07333", "K07345", "K07346", "K07347", "K07348", "K07349", "K07350", "K07351", "K07352", "K07353", "K07354", "K07355", "K07356", "K07822", "K07991", "K10564", "K10565", "K11522", "K11523", "K11524", "K11525", "K11526", "K13626", "K13820", "K13924", "K18475", "K21217", "K21218", "K23985", "K23986"

``` r
lists <- list(
  `Biofilm formation - Escherichia coli`=c("K00688", "K00694", "K00703", "K00975", "K01991", "K02398", "K02402", "K02403", "K02405", "K02425", "K02777", "K03087", "K03563", "K03566", "K03567", "K04333", "K04334", "K04335", "K04336", "K04761", "K05851", "K06204", "K07173", "K07638", "K07648", "K07659", "K07676", "K07677", "K07678", "K07687", "K07689", "K07773", "K07781", "K07782", "K10914", "K11531", "K11931", "K11935", "K11936", "K11937", "K12687", "K14051", "K18502", "K18504", "K18509", "K18515", "K18516", "K18518", "K18521", "K18522", "K18523", "K18528", "K18968", "K21084", "K21085", "K21086", "K21087", "K21088", "K21089", "K21090", "K21091"),
  `Biofilm formation - Vibrio cholerae`=c("K00640", "K01791", "K01912", "K02405", "K02452", "K02453", "K02454", "K02455", "K02456", "K02457", "K02458", "K02459", "K02460", "K02461", "K02462", "K02463", "K02472", "K02777", "K02779", "K03087", "K03092", "K03557", "K03563", "K03606", "K03666", "K05851", "K05946", "K07173", "K07181", "K07678", "K07689", "K08604", "K08720", "K10909", "K10910", "K10911", "K10912", "K10913", "K10914", "K10915", "K10916", "K10917", "K10918", "K10919", "K10920", "K10921", "K10922", "K10923", "K10924", "K10925", "K10926", "K10927", "K10928", "K10929", "K10930", "K10931", "K10932", "K10933", "K10934", "K10935", "K10936", "K10937", "K10938", "K10939", "K10940", "K10941", "K10942", "K10943", "K10961", "K10962", "K10963", "K10964", "K10965", "K12276", "K13246", "K15851", "K16554", "K18515", "K20918", "K20919", "K20920", "K20921", "K20922", "K20945", "K20946", "K20947", "K20948", "K20949", "K20950", "K20951", "K20952", "K20953", "K20954", "K20955", "K20956", "K20957", "K20958", "K20959", "K20960", "K20961", "K20962", "K20963", "K20964", "K20965", "K20966", "K20988"),
  `Biofilm formation - Pseudomonas aeruginosa` = c("K01657", "K01658", "K01768", "K02398", "K02405", "K02657", "K02658", "K02659", "K02660", "K03563", "K03651", "K06596", "K06598", "K07678", "K07689", "K10914", "K10941", "K11444", "K11890", "K11891", "K11893", "K11895", "K11900", "K11901", "K11902", "K11903", "K11907", "K11912", "K11913", "K11915", "K12990", "K12992", "K13060", "K13061", "K13487", "K13488", "K13489", "K13490", "K13491", "K16011", "K17940", "K18000", "K18001", "K18002", "K18003", "K18099", "K18100", "K18101", "K18304", "K19291", "K19735", "K20257", "K20258", "K20259", "K20968", "K20969", "K20970", "K20971", "K20972", "K20973", "K20974", "K20975", "K20976", "K20977", "K20978", "K20987", "K20997", "K20998", "K20999", "K21000", "K21001", "K21002", "K21003", "K21004", "K21005", "K21006", "K21007", "K21008", "K21009", "K21010", "K21011", "K21012", "K21019", "K21020", "K21021", "K21022", "K21023", "K21024", "K21025", "K23127"),
  `Biofilm formation`=c("K00640", "K00688", "K00694", "K00703", "K00975", "K01657", "K01658", "K01768", "K01791", "K01912", "K01991", "K02398", "K02402", "K02403", "K02405", "K02425", "K02452", "K02453", "K02454", "K02455", "K02456", "K02457", "K02458", "K02459", "K02460", "K02461", "K02462", "K02463", "K02472", "K02657", "K02658", "K02659", "K02660", "K02777", "K02779", "K03087", "K03092", "K03557", "K03563", "K03566", "K03567", "K03606", "K03651", "K03666", "K04333", "K04334", "K04335", "K04336", "K04761", "K05851", "K05946", "K06204", "K06596", "K06598", "K07173", "K07181", "K07638", "K07648", "K07659", "K07676", "K07677", "K07678", "K07687", "K07689", "K07773", "K07781", "K07782", "K08604", "K08720", "K10909", "K10910", "K10911", "K10912", "K10913", "K10914", "K10915", "K10916", "K10917", "K10918", "K10919", "K10920", "K10921", "K10922", "K10923", "K10924", "K10925", "K10926", "K10927", "K10928", "K10929", "K10930", "K10931", "K10932", "K10933", "K10934", "K10935", "K10936", "K10937", "K10938", "K10939", "K10940", "K10941", "K10942", "K10943", "K10961", "K10962", "K10963", "K10964", "K10965", "K11444", "K11531", "K11890", "K11891", "K11893", "K11895", "K11900", "K11901", "K11902", "K11903", "K11907", "K11912", "K11913", "K11915", "K11931", "K11935", "K11936", "K11937", "K12276", "K12687", "K12990", "K12992", "K13060", "K13061", "K13246", "K13487", "K13488", "K13489", "K13490", "K13491", "K14051", "K15851", "K16011", "K16554", "K17940", "K18000", "K18001", "K18002", "K18003", "K18099", "K18100", "K18101", "K18304", "K18502", "K18504", "K18509", "K18515", "K18516", "K18518", "K18521", "K18522", "K18523", "K18528", "K18968", "K19291", "K19735", "K20257", "K20258", "K20259", "K20918", "K20919", "K20920", "K20921", "K20922", "K20945", "K20946", "K20947", "K20948", "K20949", "K20950", "K20951", "K20952", "K20953", "K20954", "K20955", "K20956", "K20957", "K20958", "K20959", "K20960", "K20961", "K20962", "K20963", "K20964", "K20965", "K20966", "K20968", "K20969", "K20970", "K20971", "K20972", "K20973", "K20974", "K20975", "K20976", "K20977", "K20978", "K20987", "K20988", "K20997", "K20998", "K20999", "K21000", "K21001", "K21002", "K21003", "K21004", "K21005", "K21006", "K21007", "K21008", "K21009", "K21010", "K21011", "K21012", "K21019", "K21020", "K21021", "K21022", "K21023", "K21024", "K21025", "K21084", "K21085", "K21086", "K21087", "K21088", "K21089", "K21090", "K21091", "K23127"),
  `Nitrogen metabolism`=c("K00260", "K00261", "K00262", "K00264", "K00265", "K00266", "K00284", "K00360", "K00362", "K00363", "K00366", "K00367", "K00368", "K00370", "K00371", "K00372", "K00374", "K00376", "K00459", "K00531", "K00926", "K01455", "K01501", "K01672", "K01673", "K01674", "K01725", "K01743", "K01915", "K01948", "K02305", "K02567", "K02568", "K02575", "K02586", "K02588", "K02591", "K03385", "K04561", "K05601", "K10534", "K10535", "K10944", "K10945", "K10946", "K15371", "K15576", "K15577", "K15578", "K15579", "K15864", "K15876", "K15877", "K17877", "K18245", "K18246", "K19823", "K20932", "K20933", "K20934", "K20935", "K22896", "K22897", "K22898", "K22899"),
  `Cell adhesion molecules`=c("K03160", "K03161", "K05412", "K05413", "K05689", "K05693", "K05695", "K05718", "K05719", "K06081", "K06087", "K06088", "K06089", "K06257", "K06449", "K06454", "K06456", "K06458", "K06459", "K06461", "K06464", "K06467", "K06470", "K06471", "K06474", "K06477", "K06478", "K06483", "K06485", "K06486", "K06487", "K06490", "K06491", "K06492", "K06494", "K06495", "K06496", "K06520", "K06523", "K06527", "K06531", "K06533", "K06538", "K06539", "K06544", "K06547", "K06548", "K06550", "K06567", "K06584", "K06585", "K06590", "K06591", "K06592", "K06708", "K06710", "K06713", "K06735", "K06736", "K06744", "K06745", "K06746", "K06747", "K06751", "K06752", "K06756", "K06757", "K06759", "K06760", "K06766", "K06770", "K06771", "K06775", "K06779", "K06780", "K06781", "K06785", "K06787", "K06791", "K06793", "K06796", "K06797", "K06809", "K06815", "K06816", "K07377", "K07378", "K07379", "K07380", "K07522", "K07523", "K10784", "K10785", "K16336", "K16337", "K16338", "K16350", "K16351", "K16359", "K16360", "K23268"),
  `Quorum sensing`=c("K00494", "K01114", "K01218", "K01318", "K01364", "K01399", "K01497", "K01580", "K01626", "K01635", "K01657", "K01658", "K01728", "K01897", "K01995", "K01996", "K01997", "K01998", "K01999", "K02031", "K02032", "K02033", "K02034", "K02035", "K02052", "K02053", "K02054", "K02055", "K02250", "K02251", "K02252", "K02253", "K02402", "K02403", "K02490", "K03070", "K03071", "K03073", "K03075", "K03076", "K03106", "K03110", "K03210", "K03217", "K03400", "K03666", "K06046", "K06352", "K06353", "K06354", "K06355", "K06356", "K06358", "K06359", "K06360", "K06361", "K06363", "K06364", "K06365", "K06366", "K06369", "K06375", "K06998", "K07173", "K07344", "K07645", "K07666", "K07667", "K07680", "K07691", "K07692", "K07699", "K07706", "K07707", "K07711", "K07715", "K07781", "K07782", "K07800", "K07813", "K08321", "K08605", "K08642", "K08777", "K09823", "K09936", "K10555", "K10556", "K10557", "K10558", "K10715", "K10823", "K10909", "K10910", "K10911", "K10912", "K10913", "K10914", "K10915", "K10916", "K10917", "K11006", "K11007", "K11031", "K11033", "K11034", "K11035", "K11036", "K11037", "K11039", "K11063", "K11216", "K11530", "K11531", "K11752", "K12257", "K12292", "K12293", "K12294", "K12295", "K12296", "K12415", "K12789", "K12990", "K13060", "K13061", "K13062", "K13063", "K13075", "K13815", "K13816", "K14051", "K14645", "K14982", "K14983", "K15580", "K15581", "K15582", "K15583", "K15654", "K15655", "K15656", "K15657", "K15850", "K15851", "K15852", "K15853", "K15854", "K16619", "K17940", "K18000", "K18001", "K18002", "K18003", "K18096", "K18098", "K18099", "K18100", "K18101", "K18139", "K18304", "K18306", "K18307", "K18315", "K18316", "K18317", "K18318", "K18319", "K19666", "K19731", "K19732", "K19733", "K19734", "K19735", "K20086", "K20087", "K20088", "K20089", "K20090", "K20248", "K20249", "K20250", "K20252", "K20253", "K20256", "K20257", "K20258", "K20259", "K20260", "K20261", "K20262", "K20263", "K20264", "K20265", "K20266", "K20267", "K20268", "K20269", "K20270", "K20271", "K20272", "K20273", "K20274", "K20275", "K20276", "K20277", "K20321", "K20322", "K20323", "K20324", "K20325", "K20326", "K20327", "K20328", "K20329", "K20330", "K20331", "K20332", "K20333", "K20334", "K20335", "K20336", "K20337", "K20338", "K20339", "K20340", "K20341", "K20342", "K20343", "K20344", "K20345", "K20373", "K20374", "K20375", "K20376", "K20377", "K20378", "K20379", "K20380", "K20381", "K20382", "K20383", "K20384", "K20385", "K20386", "K20387", "K20388", "K20389", "K20390", "K20391", "K20480", "K20481", "K20482", "K20483", "K20484", "K20485", "K20486", "K20487", "K20488", "K20489", "K20490", "K20491", "K20492", "K20494", "K20527", "K20528", "K20529", "K20530", "K20531", "K20532", "K20533", "K20539", "K20540", "K20552", "K20554", "K20555", "K22954", "K22955", "K22956", "K22957", "K22968", "K23133"),
  `Methane metabolism` = c("K00018", "K00024", "K00058", "K00093", "K00121", "K00122", "K00123", "K00124", "K00125", "K00126", "K00127", "K00148", "K00169", "K00170", "K00171", "K00172", "K00189", "K00192", "K00193", "K00194", "K00195", "K00196", "K00197", "K00198", "K00200", "K00201", "K00202", "K00203", "K00204", "K00205", "K00300", "K00317", "K00319", "K00320", "K00399", "K00400", "K00401", "K00402", "K00440", "K00441", "K00442", "K00443", "K00577", "K00578", "K00579", "K00580", "K00581", "K00582", "K00583", "K00584", "K00600", "K00625", "K00672", "K00830", "K00831", "K00850", "K00863", "K00918", "K00925", "K01007", "K01070", "K01079", "K01086", "K01499", "K01595", "K01622", "K01623", "K01624", "K01689", "K01834", "K01895", "K02203", "K02446", "K03388", "K03389", "K03390", "K03396", "K03421", "K03422", "K03532", "K03533", "K03841", "K04041", "K04480", "K05299", "K05884", "K05979", "K06034", "K06914", "K07072", "K07144", "K07811", "K07812", "K07821", "K08093", "K08094", "K08097", "K08264", "K08265", "K08685", "K08691", "K08692", "K09733", "K10713", "K10714", "K10944", "K10945", "K10946", "K10977", "K10978", "K11212", "K11260", "K11261", "K11529", "K11532", "K11645", "K11779", "K11780", "K11781", "K12234", "K13039", "K13788", "K13812", "K13831", "K13942", "K14028", "K14029", "K14067", "K14080", "K14081", "K14082", "K14083", "K14084", "K14126", "K14127", "K14128", "K14940", "K14941", "K15022", "K15228", "K15229", "K15633", "K15634", "K15635", "K16157", "K16158", "K16159", "K16160", "K16161", "K16162", "K16176", "K16177", "K16178", "K16179", "K16254", "K16255", "K16256", "K16257", "K16258", "K16259", "K16260", "K16305", "K16306", "K16370", "K16792", "K16793", "K17066", "K17067", "K17068", "K17100", "K18277", "K18933", "K19793", "K21071", "K22081", "K22082", "K22083", "K22084", "K22085", "K22086", "K22087", "K22305", "K22480", "K22481", "K22482", "K22515", "K22516", "K23995", "K24182"),
  `Denitrification` = c("K00370", "K00371", "K00374", "K02567", "K02568", "K00368", "K15864", "K04561", "K02305", "K00376"),
  `Dissimilatory nitrate reduction` = c("K00958", "K00394", "K00395", "K11180", "K11181"),
  `Assimilatory nitrate reduction` = c("K00367", "K10534", "K00372", "K00360", "K00366", "K17877"),
  `EPS` = c("K09688", "K09690", "K09689", "K09691", "K06041", "K07265", "K07266", "K00979", "K00694", "K19291", "K19292", "K16081", "K19296", "K01795", "K19294", "K01729", "K20541", "K20543", "K20542", "K05789", "K05790", "K10107", "K13009", "K13008", "K01991", "K01104", "K16692", "K02853", "K16693", "K16695", "K18799"),
  `EPS - ABC mechanisms` = c("K09688", "K09690", "K09689", "K09691", "K06041", "K07265", "K07266", "K00979"),
  `EPS - synthase mechanisms` = c("K00694", "K19291", "K19292", "K16081", "K19296", "K01795", "K19294", "K01729", "K20541", "K20543", "K20542"),
  `EPS - WZY mechanisms` = c("K05789", "K05790", "K10107", "K13009", "K13008", "K01991", "K01104", "K16692", "K02853", "K16693", "K16695", "K18799"),
  `N-glycosylation` = c("K07151", "K12666", "K12667", "K12668", "K12669", "K12670", "K00730", "K12691"),
  `Flagellar assembly` = c("K02386", "K02387", "K02388", "K02389", "K02390", "K02391", "K02392", "K02393", "K02394", "K02396", "K02397", "K02398", "K02399", "K02400", "K02401", "K02402", "K02403", "K02405", "K02406", "K02407", "K02408", "K02409", "K02410", "K02411", "K02412", "K02413", "K02414", "K02416", "K02417", "K02418", "K02419", "K02420", "K02421", "K02422", "K02423", "K02556", "K02557", "K13820", "K21217", "K21218"),
  `Bacterial chemotaxis` = c("K00575", "K02410", "K02416", "K02417", "K02556", "K02557", "K03406", "K03407", "K03408", "K03409", "K03410", "K03411", "K03412", "K03413", "K03414", "K03415", "K03776", "K05874", "K05875", "K05876", "K05877", "K10108", "K10439", "K10540", "K12368", "K13924"),
  `Bacterial motility proteins` = c("K00575", "K02278", "K02279", "K02280", "K02281", "K02282", "K02283", "K02382", "K02383", "K02384", "K02385", "K02386", "K02387", "K02388", "K02389", "K02390", "K02391", "K02392", "K02393", "K02394", "K02395", "K02396", "K02397", "K02398", "K02399", "K02400", "K02401", "K02402", "K02403", "K02404", "K02405", "K02406", "K02407", "K02408", "K02409", "K02410", "K02411", "K02412", "K02413", "K02414", "K02415", "K02416", "K02417", "K02418", "K02419", "K02420", "K02421", "K02422", "K02423", "K02424", "K02425", "K02556", "K02557", "K02650", "K02651", "K02652", "K02653", "K02654", "K02655", "K02656", "K02657", "K02658", "K02659", "K02660", "K02661", "K02662", "K02663", "K02664", "K02665", "K02666", "K02667", "K02668", "K02669", "K02670", "K02671", "K02672", "K02673", "K02674", "K02675", "K02676", "K03406", "K03407", "K03408", "K03409", "K03410", "K03411", "K03412", "K03413", "K03414", "K03415", "K03516", "K03776", "K04562", "K05874", "K05875", "K05876", "K05877", "K06595", "K06596", "K06597", "K06598", "K06599", "K06600", "K06601", "K06602", "K06603", "K07324", "K07325", "K07327", "K07328", "K07329", "K07330", "K07331", "K07332", "K07333", "K07345", "K07346", "K07347", "K07348", "K07349", "K07350", "K07351", "K07352", "K07353", "K07354", "K07355", "K07356", "K07822", "K07991", "K10564", "K10565", "K11522", "K11523", "K11524", "K11525", "K11526", "K13626", "K13820", "K13924", "K18475", "K21217", "K21218", "K23985", "K23986")
)

# beepr::beep(sound = "mario")
```

Using lists, construct bar plots for differential expression for each
category

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:19:03 EDT"

``` r
de <- read.delim("noiseq/nsb_ko_controlvsstannous_tmm_filtered.txt")
## ko_hier made in network rmd - need to move to one sheet
de_ann <- merge(de, ko_hier, by.x = "id", by.y = "ko", all.x = TRUE)
  de_ann$X <- NULL
  colnames(de_ann)
```

    ##  [1] "id"               "control_mean"     "stannous_mean"    "theta"           
    ##  [5] "prob"             "log2FC"           "P"                "kegg_gene"       
    ##  [9] "kegg1"            "kegg2"            "kegg3"            "kegg_description"

``` r
write_delim(x = de_ann, file = "tables/nsb_annotated_differential_expression.txt", delim = "\t", quote = NULL, col_names = TRUE)


## top 15 for each treatment by fold change
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  de_sub <- de_ann %>% filter(grepl(paste((lists[[i]]), collapse="|"), id)) %>%
              arrange(log2FC) %>%
              filter(row_number() > max(row_number()) - 15 | row_number() <= 15)
  p <- ggplot(de_sub, aes(x = log2FC, y = reorder(id, -log2FC), fill = -log10(P))) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient2(low="blue", mid = "gray88", high="red", midpoint = 1.3) +
    theme_classic() +
    labs(y = "KEGG Ortholog", x = expression("log"["2"]*" fold change"), title = paste0("<- higher in SnF2 | ",nm," | higher in MFP ->")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_legend_position(x = "topright")
  print(p)
  nm2 <- gsub(" ", "_", nm)
  ggsave(paste0("plots/nsb_pathway_",nm2,"_top15trmt.png"), dpi = 600)
  write.table(de_sub, paste0("tables/nsb_pathway_",nm2,"_top15trmt.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}
```

![](index_files/figure-gfm/de%20barplots-1.png)<!-- -->![](index_files/figure-gfm/de%20barplots-2.png)<!-- -->![](index_files/figure-gfm/de%20barplots-3.png)<!-- -->![](index_files/figure-gfm/de%20barplots-4.png)<!-- -->![](index_files/figure-gfm/de%20barplots-5.png)<!-- -->![](index_files/figure-gfm/de%20barplots-6.png)<!-- -->![](index_files/figure-gfm/de%20barplots-7.png)<!-- -->![](index_files/figure-gfm/de%20barplots-8.png)<!-- -->![](index_files/figure-gfm/de%20barplots-9.png)<!-- -->![](index_files/figure-gfm/de%20barplots-10.png)<!-- -->![](index_files/figure-gfm/de%20barplots-11.png)<!-- -->![](index_files/figure-gfm/de%20barplots-12.png)<!-- -->![](index_files/figure-gfm/de%20barplots-13.png)<!-- -->![](index_files/figure-gfm/de%20barplots-14.png)<!-- -->![](index_files/figure-gfm/de%20barplots-15.png)<!-- -->![](index_files/figure-gfm/de%20barplots-16.png)<!-- -->![](index_files/figure-gfm/de%20barplots-17.png)<!-- -->![](index_files/figure-gfm/de%20barplots-18.png)<!-- -->![](index_files/figure-gfm/de%20barplots-19.png)<!-- -->

``` r
## make table that includes all category annotations (many will be repeated)
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  nm2 <- gsub(" ", "_", nm)
  de_sub_all <- de %>% select(!X) %>% filter(grepl(paste((lists[[i]]), collapse="|"), id))
  # ko_hier_sub <- ko_hier %>% filter(grepl(paste((lists[[i]]), collapse="|"), ko))
  de_sub_all <- merge(de_sub_all, khier, by.x = "id", by.y = "ko", all.x = TRUE)
  write.table(de_sub_all, paste0("tables/nsb_pathway_",nm2,"_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}

## top 30 highest absolute fold changes
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  de_sub <- de_ann %>% filter(grepl(paste((lists[[i]]), collapse="|"), id)) %>%
              arrange(abs(log2FC)) %>%
              top_n(n = 30, wt = abs(log2FC))
  p <- ggplot(de_sub, aes(x = log2FC, y = reorder(id, -log2FC), fill = -log10(P))) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient2(low="blue", mid = "gray88", high="red", midpoint = 1.3) +
    theme_classic() +
    labs(y = "KEGG Ortholog", x = expression("log"["2"]*" fold change"), title = paste0("<- higher in SnF2 | ",nm," | higher in MFP ->")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_legend_position(x = "topright")
  print(p)
  nm2 <- gsub(" ", "_", nm)
  ggsave(paste0("plots/nsb_pathway_",nm2,"_top15trmt.png"), dpi = 600)
  write.table(de_sub, paste0("tables/nsb_pathway_",nm2,"_topFC.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}
```

![](index_files/figure-gfm/de%20barplots-20.png)<!-- -->![](index_files/figure-gfm/de%20barplots-21.png)<!-- -->![](index_files/figure-gfm/de%20barplots-22.png)<!-- -->![](index_files/figure-gfm/de%20barplots-23.png)<!-- -->![](index_files/figure-gfm/de%20barplots-24.png)<!-- -->![](index_files/figure-gfm/de%20barplots-25.png)<!-- -->![](index_files/figure-gfm/de%20barplots-26.png)<!-- -->![](index_files/figure-gfm/de%20barplots-27.png)<!-- -->![](index_files/figure-gfm/de%20barplots-28.png)<!-- -->![](index_files/figure-gfm/de%20barplots-29.png)<!-- -->![](index_files/figure-gfm/de%20barplots-30.png)<!-- -->![](index_files/figure-gfm/de%20barplots-31.png)<!-- -->![](index_files/figure-gfm/de%20barplots-32.png)<!-- -->![](index_files/figure-gfm/de%20barplots-33.png)<!-- -->![](index_files/figure-gfm/de%20barplots-34.png)<!-- -->![](index_files/figure-gfm/de%20barplots-35.png)<!-- -->![](index_files/figure-gfm/de%20barplots-36.png)<!-- -->![](index_files/figure-gfm/de%20barplots-37.png)<!-- -->![](index_files/figure-gfm/de%20barplots-38.png)<!-- -->

``` r
## all ko within the category
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  de_sub <- de_ann %>% filter(grepl(paste((lists[[i]]), collapse="|"), id))
  p <- ggplot(de_sub, aes(x = log2FC, y = reorder(id, -log2FC), fill = -log10(P))) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient2(low="blue", mid = "gray88", high="red", midpoint = 1.3) +
    theme_classic() +
    labs(y = "KEGG Ortholog", x = expression("log"["2"]*" fold change"), title = paste0("<- higher in SnF2 | ",nm," | higher in MFP ->")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_legend_position(x = "topright")
  print(p)
  nm2 <- gsub(" ", "_", nm)
  ggsave(paste0("plots/nsb_pathway_",nm2,"_all.png"), dpi = 600)
  write.table(de_sub, paste0("tables/nsb_pathway_",nm2,"_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}
```

![](index_files/figure-gfm/de%20barplots-39.png)<!-- -->![](index_files/figure-gfm/de%20barplots-40.png)<!-- -->![](index_files/figure-gfm/de%20barplots-41.png)<!-- -->![](index_files/figure-gfm/de%20barplots-42.png)<!-- -->![](index_files/figure-gfm/de%20barplots-43.png)<!-- -->![](index_files/figure-gfm/de%20barplots-44.png)<!-- -->![](index_files/figure-gfm/de%20barplots-45.png)<!-- -->![](index_files/figure-gfm/de%20barplots-46.png)<!-- -->![](index_files/figure-gfm/de%20barplots-47.png)<!-- -->![](index_files/figure-gfm/de%20barplots-48.png)<!-- -->![](index_files/figure-gfm/de%20barplots-49.png)<!-- -->![](index_files/figure-gfm/de%20barplots-50.png)<!-- -->![](index_files/figure-gfm/de%20barplots-51.png)<!-- -->![](index_files/figure-gfm/de%20barplots-52.png)<!-- -->![](index_files/figure-gfm/de%20barplots-53.png)<!-- -->![](index_files/figure-gfm/de%20barplots-54.png)<!-- -->![](index_files/figure-gfm/de%20barplots-55.png)<!-- -->![](index_files/figure-gfm/de%20barplots-56.png)<!-- -->![](index_files/figure-gfm/de%20barplots-57.png)<!-- -->

``` r
## top de by p-value
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  de_sub <- de_ann %>% filter(grepl(paste((lists[[i]]), collapse="|"), id)) %>%
              arrange(P) %>%
              top_n(n = -30, wt = P)
  p[[i]] <- ggplot(de_sub, aes(x = log2FC, y = reorder(id, -log2FC), fill = -log10(P))) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient2(low="blue", mid = "gray88", high="red", midpoint = 1.3, labels = c("1.0", "1.3", "2.0"), breaks=c(1,1.3,2)) +
    theme_classic() +
    labs(y = "KEGG Ortholog", x = expression("log"["2"]*" fold change"), 
         title = paste0("← higher in SnF2 | ",nm," | higher in MFP →")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10)) +
    theme_legend_position(x = "topright")
  print(p[[i]])
  nm2 <- gsub(" ", "_", nm)
  ggsave(paste0("plots/nsb_pathway_",nm2,"_low30pvalue_fl.png"), dpi = 600)
  ggsave(paste0("plots/nsb_pathway_",nm2,"_low30pvalue_fl.pdf"), dpi = 600)
  ggsave(paste0("plots/nsb_pathway_",nm2,"_low30pvalue_fl.tiff"), dpi = 600)
  write.table(de_sub, paste0("tables/nsb_pathway_",nm2,"_low30pvalue.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}
```

![](index_files/figure-gfm/de%20barplots-58.png)<!-- -->![](index_files/figure-gfm/de%20barplots-59.png)<!-- -->![](index_files/figure-gfm/de%20barplots-60.png)<!-- -->![](index_files/figure-gfm/de%20barplots-61.png)<!-- -->![](index_files/figure-gfm/de%20barplots-62.png)<!-- -->![](index_files/figure-gfm/de%20barplots-63.png)<!-- -->![](index_files/figure-gfm/de%20barplots-64.png)<!-- -->![](index_files/figure-gfm/de%20barplots-65.png)<!-- -->![](index_files/figure-gfm/de%20barplots-66.png)<!-- -->![](index_files/figure-gfm/de%20barplots-67.png)<!-- -->![](index_files/figure-gfm/de%20barplots-68.png)<!-- -->![](index_files/figure-gfm/de%20barplots-69.png)<!-- -->![](index_files/figure-gfm/de%20barplots-70.png)<!-- -->![](index_files/figure-gfm/de%20barplots-71.png)<!-- -->![](index_files/figure-gfm/de%20barplots-72.png)<!-- -->![](index_files/figure-gfm/de%20barplots-73.png)<!-- -->![](index_files/figure-gfm/de%20barplots-74.png)<!-- -->![](index_files/figure-gfm/de%20barplots-75.png)<!-- -->![](index_files/figure-gfm/de%20barplots-76.png)<!-- -->

``` r
# beepr::beep(sound = "mario")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:20:17 EDT"

``` r
p[[1]] <- p[[1]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[2]] <- p[[2]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[3]] <- p[[3]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[4]] <- p[[4]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())

plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], align = "vh", ncol = 2, rel_heights = c(2,2,0.75), labels = c("A", "B", "C", "D"))
```

![](index_files/figure-gfm/bar%20cowplot-1.png)<!-- -->

``` r
ggsave("plots/figure09_de_bars3.png", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars3.pdf", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars3.tiff", dpi = 600, height = 10, width = 10)

legend <- cowplot::get_legend(p[[4]] + theme(legend.position=c(0.3, 1), legend.direction = "horizontal"))
legend <- cowplot::get_legend(p[[4]] + theme(legend.position = "bottom", legend.direction = "horizontal"))

p[[1]] <- p[[1]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[2]] <- p[[2]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[3]] <- p[[3]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[4]] <- p[[4]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[5]] <- p[[5]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[6]] <- p[[6]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[7]] <- p[[7]] + theme(legend.position='none') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())

plot_grid(p[[1]], p[[2]], p[[3]], p[[6]], p[[4]], legend, align = "vh", ncol = 2, rel_heights = c(2,2,0.75), labels = c("A", "B", "C", "D", "E", NULL))
```

![](index_files/figure-gfm/bar%20cowplot-2.png)<!-- -->

``` r
ggsave("plots/figure09_de_bars.png", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars.pdf", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars.tiff", dpi = 600, height = 10, width = 10)

plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], legend, align = "vh", ncol = 2, rel_heights = c(2,2,0.75), labels = c("A", "B", "C", "D", NULL))
```

![](index_files/figure-gfm/bar%20cowplot-3.png)<!-- -->

``` r
ggsave("plots/figure09_de_bars2.png", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars2.pdf", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars2.tiff", dpi = 600, height = 10, width = 10)

p[[1]] <- p[[1]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[2]] <- p[[2]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[3]] <- p[[3]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[4]] <- p[[4]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[5]] <- p[[5]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[6]] <- p[[6]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())
p[[7]] <- p[[7]] + theme(legend.position='bottom') + theme(plot.title = element_text(size=10)) + theme(axis.title.y = element_blank())

plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], align = "vh", ncol = 2, rel_heights = c(2,2,0.75), labels = c("A", "B", "C", "D"))
```

![](index_files/figure-gfm/bar%20cowplot-4.png)<!-- -->

``` r
ggsave("plots/figure09_de_bars4.png", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars4.pdf", dpi = 600, height = 10, width = 10)
ggsave("plots/figure09_de_bars4.tiff", dpi = 600, height = 10, width = 10)

# beepr::beep(sound = "mario")
```

Plot mean normalized expression for 30 most significant KOs for each
category

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:20:40 EDT"

``` r
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  de_sub <- de_ann %>% filter(grepl(paste((lists[[i]]), collapse="|"), id)) %>%
              arrange(P) %>%
              top_n(n = -30, wt = P) %>%
    reshape::melt(measure.vars = 2:3)
  p <- ggplot(de_sub, aes(x = log(ceiling(value+1)), y = reorder(id, value), fill = variable)) +
    geom_bar(stat = 'identity', position=position_dodge(width=0.8), width = 0.8) +
    scale_fill_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous")) +
    theme_classic() +
    labs(y = "KEGG Ortholog", x = "Normalized Expression", title = nm) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  nm2 <- gsub(" ", "_", nm)
  ggsave(paste0("plots/ns_pathway_",nm2,"_top30sig_expression_ordered.png"), dpi = 600)
}
```

![](index_files/figure-gfm/mean%20expression-1.png)<!-- -->![](index_files/figure-gfm/mean%20expression-2.png)<!-- -->![](index_files/figure-gfm/mean%20expression-3.png)<!-- -->![](index_files/figure-gfm/mean%20expression-4.png)<!-- -->![](index_files/figure-gfm/mean%20expression-5.png)<!-- -->![](index_files/figure-gfm/mean%20expression-6.png)<!-- -->![](index_files/figure-gfm/mean%20expression-7.png)<!-- -->![](index_files/figure-gfm/mean%20expression-8.png)<!-- -->![](index_files/figure-gfm/mean%20expression-9.png)<!-- -->![](index_files/figure-gfm/mean%20expression-10.png)<!-- -->![](index_files/figure-gfm/mean%20expression-11.png)<!-- -->![](index_files/figure-gfm/mean%20expression-12.png)<!-- -->![](index_files/figure-gfm/mean%20expression-13.png)<!-- -->![](index_files/figure-gfm/mean%20expression-14.png)<!-- -->![](index_files/figure-gfm/mean%20expression-15.png)<!-- -->![](index_files/figure-gfm/mean%20expression-16.png)<!-- -->![](index_files/figure-gfm/mean%20expression-17.png)<!-- -->![](index_files/figure-gfm/mean%20expression-18.png)<!-- -->![](index_files/figure-gfm/mean%20expression-19.png)<!-- -->

``` r
for(i in 1:length(lists)){
  nm <- (names(lists)[i])
  de_sub <- de_ann %>% filter(grepl(paste((lists[[i]]), collapse="|"), id)) %>%
              arrange(log2FC) %>%
              top_n(n = 30, wt = abs(log2FC)) %>%
    reshape::melt(measure.vars = 2:3)
  p <- ggplot(de_sub, aes(x = log(ceiling(value+1)), y = reorder(id, value), fill = variable)) +
    geom_bar(stat = 'identity', position=position_dodge(width=0.8), width = 0.8) +
    scale_fill_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous")) +
    theme_classic() +
    labs(y = "KEGG Ortholog", x = "Normalized Expression", title = nm) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  nm2 <- gsub(" ", "_", nm)
  ggsave(paste0("plots/ns_pathway_",nm2,"_top30fc_expression_ordered.png"), dpi = 600)
}
```

![](index_files/figure-gfm/mean%20expression-20.png)<!-- -->![](index_files/figure-gfm/mean%20expression-21.png)<!-- -->![](index_files/figure-gfm/mean%20expression-22.png)<!-- -->![](index_files/figure-gfm/mean%20expression-23.png)<!-- -->![](index_files/figure-gfm/mean%20expression-24.png)<!-- -->![](index_files/figure-gfm/mean%20expression-25.png)<!-- -->![](index_files/figure-gfm/mean%20expression-26.png)<!-- -->![](index_files/figure-gfm/mean%20expression-27.png)<!-- -->![](index_files/figure-gfm/mean%20expression-28.png)<!-- -->![](index_files/figure-gfm/mean%20expression-29.png)<!-- -->![](index_files/figure-gfm/mean%20expression-30.png)<!-- -->![](index_files/figure-gfm/mean%20expression-31.png)<!-- -->![](index_files/figure-gfm/mean%20expression-32.png)<!-- -->![](index_files/figure-gfm/mean%20expression-33.png)<!-- -->![](index_files/figure-gfm/mean%20expression-34.png)<!-- -->![](index_files/figure-gfm/mean%20expression-35.png)<!-- -->![](index_files/figure-gfm/mean%20expression-36.png)<!-- -->![](index_files/figure-gfm/mean%20expression-37.png)<!-- -->![](index_files/figure-gfm/mean%20expression-38.png)<!-- -->

``` r
# beepr::beep(sound = "mario")
```

Prepare tables for GSEA analysis

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:21:08 EDT"

``` r
## filter counts
tmm_counts_filtered <- tmm_counts[rowSums(tmm_counts[,2:ncol(tmm_counts)]) > 1,]

## gene sets
## go_annotations.txt made with trinity script - extract_GO_assignments_from_Trinotate_xls.pl 
go_gsea <- read.table("data/go_annotations.txt", header=F, row.names=1,stringsAsFactors=F)
## count number of columns needed by counting commas in go term column for each row - find max and add 1 for total columns needed - 
## may need to run after filtering out low count transcripts
max(str_count(go_gsea$V2, ","))
```

    ## [1] 244

``` r
go_gsea %>% add_column(description = "na", .before = "V2") %>% 
  separate(V2, sep = ",", into = paste("go", 1:245, sep = "_")) %>%
  subset(row.names(.) %in% tmm_counts_filtered$X) -> go_gsea_out
write.table(go_gsea_out, "data/gsea_go.gmt", sep = "\t", col.names = FALSE, quote = FALSE, row.names = TRUE, na = "")

## expression counts table
entry_tmm_counts <- tmm_counts_filtered %>% 
  merge(., updated_annots, by.x = "X", by.y = "X.gene_id") %>%
  select(Entry, ends_with("_raw_rmhost")) %>%
  na_if("") %>% na.omit
entry_tmm_counts <- aggregate(entry_tmm_counts[,c(2:ncol(entry_tmm_counts))], list(entry_tmm_counts[,1]), FUN=sum)
entry_tmm_counts <- entry_tmm_counts %>%
  add_column(Description = "na", .after = "Group.1") %>%
  dplyr::rename(NAME=Group.1) %>% rbind(names(.), .)
write.table(entry_tmm_counts, "data/entry_tmm_counts.txt", sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, na = "")

## phenotypes
gsea3 <- tibble(V1 = paste(as.character(ncol(tmm_counts[2:ncol(tmm_counts)])), length(unique(myfactors$toothpaste)), "1", sep = " "))
gsea4 <- unique(myfactors$toothpaste)
gsea5 <- tibble(V1 = paste("#", gsea4[1], gsea4[2], sep = " "))
gsea6 <- myfactors %>% mutate(gsea = case_when(
  startsWith(toothpaste, "C") ~ 0,
  startsWith(toothpaste, "T") ~ 1
  )) %>% dplyr::select(gsea) %>% toString() %>% as.data.frame() %>%
  str_remove("c\\(") %>% str_remove("\\)") %>% str_remove_all(",") %>% data.frame(V1 = .)
gsea_pheno <- bind_rows(gsea3, gsea5, gsea6)
write.table(gsea_pheno, "data/gsea_pheno.cls", sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, na = "")

## chip - really annotations since not using microarray
## repeated from main analysis with modifications
annots <- read.delim("data/trinotate_annotation_report_1e3_rmhost_max1.xls", header = TRUE)
uniprot_map <- read.delim("data/uniprot-filtered-reviewed_yes_2022.tab", header = TRUE)
annots %>% separate(sprot_Top_BLASTX_hit, c("uniprot_x", NA, NA, "percent_id_x", "evalue_x", "recname_x", "taxonomy_x"), sep = "\\^", remove = FALSE) %>% 
  separate(sprot_Top_BLASTP_hit, c("uniprot_p", NA, NA, "percent_id_p", "evalue_p", "recname_p", "taxonomy_p"), sep = "\\^", remove = FALSE) %>% 
  separate(evalue_x, c(NA, "evalue_x"), sep = "\\:", remove = TRUE) %>%
  separate(evalue_p, c(NA, "evalue_p"), sep = "\\:", remove = TRUE) %>%
  separate(percent_id_x, c("percent_id_x", NA), sep = "\\%", remove = TRUE) %>%
  separate(percent_id_p, c("percent_id_p", NA), sep = "\\%", remove = TRUE) %>% 
  separate(Kegg, c(NA, "kegg1"), sep = ":", remove = FALSE, extra = "merge") %>%
  separate(kegg1, c("kegg1", "kegg2"), sep = ";", remove = TRUE) %>%
  arrange(X.gene_id, evalue_x) %>% group_by(X.gene_id) %>% 
  slice_head(n = 1) %>% ungroup() %>% 
  na_if(".") %>% mutate(., uniprot_id = coalesce(uniprot_x, uniprot_p)) %>%
  mutate(., taxa = coalesce(taxonomy_x, taxonomy_p)) %>%
  merge(., uniprot_map, by.x = "uniprot_id", by.y = "Entry.name", all.x = TRUE) -> updated_annots
gsea_chip <- updated_annots %>% 
  select(X.gene_id, uniprot_id, Protein.names) %>%
  subset(X.gene_id %in% row.names(go_gsea_out)) %>%
  dplyr::rename(`Probe Set ID`=X.gene_id, `Gene Symbol`=uniprot_id, `Gene Title`=Protein.names)
write.table(gsea_chip, "data/gsea_annots.chip", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, na = "")

# beepr::beep(sound = "mario")
```

After completing GSEA analysis, move output table to ./data and plot
heatmaps for 50 top features for each treatment

``` r
## after completing the analysis in gsea, output is imported and plotted

cat45 <- read.delim("data/entryid_for_heatmap_gsea.txt")
cat45bottom <- read.delim("data/entryid_for_heatmap_gsea_bottom.txt")
# cat45exp <- merge(cat45, entry_tmm_counts[-1,c(1,3:ncol(entry_tmm_counts))], by.x = "Uniprot.id", by.y = "NAME", all.x = TRUE)

hm_gct <- read.delim("data/heat_map_Top_50_Features.gct", skip = 2)
cat45_hm_gct <- merge(cat45, hm_gct[,c(1,3:ncol(hm_gct))], by.x = "Uniprot.id", by.y = "NAME", all.x = TRUE)
cat45_hm_gct_mat <- cat45_hm_gct %>% select(-Category, -Protein.name) %>% column_to_rownames(var = "Uniprot.id")

cat45bottom_hm_gct <- merge(cat45bottom, hm_gct[,c(1,3:ncol(hm_gct))], by.x = "Uniprot.id", by.y = "NAME", all.x = TRUE)
cat45bottom_hm_gct_mat <- cat45bottom_hm_gct %>% select(-Category, -Protein.name) %>% column_to_rownames(var = "Uniprot.id")

## back to pheatmap
category <- rbind(cat45, cat45bottom)
category <- category %>% column_to_rownames(var = "Uniprot.id") %>% select(-Protein.name)
cat_top <- cat45 %>% column_to_rownames(var = "Uniprot.id") %>% select(-Protein.name)
cat_bot <- cat45bottom %>% column_to_rownames(var = "Uniprot.id") %>% select(-Protein.name)
top_ordered <- cat45_hm_gct_mat[cat45$Uniprot.id, ]
bot_ordered <- cat45bottom_hm_gct_mat[cat45bottom$Uniprot.id, ]
annot_col <- myfactors %>% select(toothpaste)

newCols <- colorRampPalette(colorspace::diverge_hcl(length(unique(cat_top$Category)), palette = "Cork"))
mycolors <- newCols(length(unique(cat_top$Category)))
names(mycolors) <- unique(cat_top$Category)
mycolors <- list(Category = mycolors)

png(filename = "plots/figure11_pheatmap_top_labels.png", width = 10, height = 9, units = "in", res = 600)
# pdf(file = "plots/figure11_pheatmap_top_labels.pdf", width = 10, height = 9)
# tiff(filename = "plots/figure11_pheatmap_top_labels.tiff", width = 10, height = 9, units = "in", res = 600)
pheatmap(top_ordered, color = colorspace::diverge_hsv(n = 100, s = 2.5), cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
         annotation_row = cat_top, annotation_names_row = FALSE, annotation_colors = mycolors, 
         labels_row = paste(cat45$Uniprot.id, cat45$Protein.name, sep = ": "),
         cellwidth = 10, cellheight = 10, angle_col = 0, fontsize_col = 10, 
         # annotation_col = annot_col, annotation_names_col = FALSE,
         labels_col = c(expression(SnF[2]), "MFP"))
         # labels_col = c("p1", "p2", "p3", "p4", "p5", "p6", "p1", "p2", "p3", "p4", "p5", "p6"))
         # labels_col = c(expression(SnF[2]~"p1"), expression(SnF[2]~"p2"), expression(SnF[2]~"p3"), 
         #                expression(SnF[2]~"p4"), expression(SnF[2]~"p5"), expression(SnF[2]~"p6"),
         #                "CDC p1", "CDC p2", "CDC p3", "CDC p4", "CDC p5", "CDC p6"))
dev.off()
```

    ## quartz_off_screen 
    ##                 3

``` r
newCols <- colorRampPalette(colorspace::diverge_hcl(length(unique(cat_bot$Category)), palette = "Cork"))
mycolors <- newCols(length(unique(cat_bot$Category)))
names(mycolors) <- unique(cat_bot$Category)
mycolors <- list(Category = mycolors)

png(filename = "plots/figure12_pheatmap_bottom_labels.png", width = 10, height = 9, units = "in", res = 600)
# pdf(file = "plots/figure12_pheatmap_bottom_labels.pdf", width = 10, height = 9)
# tiff(filename = "plots/figure12_pheatmap_bottom_labels.tiff", width = 10, height = 9, units = "in", res = 600)
pheatmap(bot_ordered, color = colorspace::diverge_hsv(n = 100, s = 2.5), cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
         annotation_row = cat_bot, annotation_names_row = FALSE, annotation_colors = mycolors, 
         labels_row = paste(cat45bottom$Uniprot.id, cat45bottom$Protein.name, sep = ": "),
         cellwidth = 10, cellheight = 10, angle_col = 0, fontsize_col = 10, 
         labels_col = c(expression(SnF[2]), "MFP"))
         # labels_col = c("p1", "p2", "p3", "p4", "p5", "p6", "p1", "p2", "p3", "p4", "p5", "p6"))
         # labels_col = c(expression(SnF[2]~"p1"), expression(SnF[2]~"p2"), expression(SnF[2]~"p3"), 
         #                expression(SnF[2]~"p4"), expression(SnF[2]~"p5"), expression(SnF[2]~"p6"), 
         #                "CDC p1", "CDC p2", "CDC p3", "CDC p4", "CDC p5", "CDC p6"))
dev.off()
```

    ## quartz_off_screen 
    ##                 3

### Amplicon Analysis

dada2 analysis mainly follows tutorial from
<https://benjjneb.github.io/dada2/>

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:30 EDT"

``` r
library(dada2); packageVersion("dada2")
```

    ## [1] '1.24.0'

``` r
library(tidyr); packageVersion("tidyr")
```

    ## [1] '1.2.0'

``` r
library(ShortRead); packageVersion("ShortRead")
```

    ## [1] '1.54.0'

``` r
library(DECIPHER); packageVersion("DECIPHER")
```

    ## [1] '2.24.0'

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.40.0'

``` r
library(decontam); packageVersion("decontam")
```

    ## [1] '1.16.0'

``` r
library(vegan); packageVersion("vegan")
```

    ## [1] '2.6.2'

``` r
library(tidyverse); packageVersion("tidyverse")
```

    ## [1] '1.3.1'

``` r
library(wesanderson); packageVersion("wesanderson")
```

    ## [1] '0.3.6'

``` r
library(scales); packageVersion("scales")
```

    ## [1] '1.2.1'

``` r
library(plyr); packageVersion("plyr")
```

    ## [1] '1.8.7'

``` r
library(microbiome); packageVersion("microbiome")
```

    ## [1] '1.18.0'

``` r
library(reshape2); packageVersion("reshape2")
```

    ## [1] '1.4.4'

``` r
library(DESeq2); packageVersion("DESeq2")
```

    ## [1] '1.36.0'

``` r
library(cowplot); packageVersion("cowplot")
```

    ## [1] '1.1.1'

``` r
library(ranacapa); packageVersion("ranacapa")
```

    ## [1] '0.1.0'

``` r
library(data.table); packageVersion("data.table")
```

    ## [1] '1.14.2'

``` r
## prepared in Rmd for metatranscriptome - only needed if running independently
# dir.create("tables")
# dir.create("plots")
```

#### DADA2 pipeline

Read in data

``` r
# Sys.time()
# 
# path_16s <- "raw_reads"
# list.files(path_16s)
# 
# ## check pattern against actual reads
# fnFs <- sort(list.files(path_16s, pattern="_R1_", full.names = TRUE)) %>%
#   subset(grepl("MFP|SNF|Water", .))
# fnRs <- sort(list.files(path_16s, pattern="_R2_", full.names = TRUE)) %>%
#   subset(grepl("MFP|SNF|Water", .))
# 
# 
# sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# names(fnFs) <- sample.names
# names(fnRs) <- sample.names
```

Remove any remaining primers

``` r
# Sys.time()
# 
# ## used https://www.nemabiome.ca/dada2_workflow.html
# ## check that reads contain primers and count instances
# 
# ## count primers function
# count_primers <- function(primer, filename) {
#   num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
#   return(sum(num_hits > 0))
# }
# 
# ## primers used are 347F and 803R
# fwd_primer <- "GGAGGCAGCAGTRRGGAAT"
# rev_primer <- "CTACCRGGGTATCTAATCC"
# fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
# rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
# 
# 
# count_primers(fwd_primer, fnFs[[1]])
# count_primers(rev_primer, fnRs[[1]])
# count_primers(fwd_primer_rev, fnFs[[1]])
# count_primers(rev_primer_rev, fnRs[[1]])
# 
# ## prepare to trim reads with cutadapt (need to install program before running)
# cutadapt <- path.expand("~/anaconda3/bin/cutadapt")
# system2(cutadapt, args = "--version")
# 
# cut_dir <- file.path("cutadapt")
# if (!dir.exists(cut_dir)) dir.create(cut_dir)
# 
# fwd_cut <- file.path(cut_dir, basename(fnFs))
# rev_cut <- file.path(cut_dir, basename(fnRs))
# 
# names(fwd_cut) <- sample.names
# names(rev_cut) <- sample.names
# 
# cut_logs <- path.expand(file.path(cut_dir, paste0(sample.names, ".log")))
# 
# cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev,
#                    "-G", rev_primer, "-A", fwd_primer_rev,
#                    "-n", 2)
# 
# for (i in seq_along(fnFs)) {
#   system2(cutadapt,
#           args = c(cutadapt_args,
#                    "-o", fwd_cut[i], "-p", rev_cut[i],
#                    fnFs[i], fnRs[i]),
#           stdout = cut_logs[i])  
# }
```

Quality filter and trim reads

``` r
# Sys.time()
# 
# filt_dir <- file.path("cutadapt/filtered")
# if (!dir.exists(filt_dir)) dir.create(filt_dir)
# 
# fwd_filt <- file.path(filt_dir, basename(fnFs))
# rev_filt <- file.path(filt_dir, basename(fnRs))
# 
# names(fwd_filt) <- sample.names
# names(rev_filt) <- sample.names
# 
# length(fnFs)
# length(fnRs)
# 
# filtered_out <- filterAndTrim(fwd = fwd_cut, filt = fwd_filt, rev = rev_cut, filt.rev = rev_filt,
#                               truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# filtered_out
```

Denoise reads and merge pairs

``` r
# Sys.time()
# 
# err_fwd <- learnErrors(fwd_filt, multithread = TRUE)
# err_rev <- learnErrors(rev_filt, multithread = TRUE)
# plotErrors(err_fwd, nominalQ = TRUE)
# plotErrors(err_rev, nominalQ = TRUE)
# 
# dada_fwd <- dada(fwd_filt, err = err_fwd, multithread = TRUE)
# dada_rev <- dada(rev_filt, err = err_rev, multithread = TRUE)
# 
# mergers <- mergePairs(dadaF = dada_fwd, dadaR = dada_rev, derepF = fwd_filt, derepR = rev_filt,
#                       maxMismatch = 1, verbose=TRUE)
# 
# seqtab <- makeSequenceTable(mergers)
# dim(seqtab)
# 
# seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
# dim(seqtab_nochim)
# 
# table(nchar(getSequences(seqtab_nochim)))
# 
# getN <- function(x) sum(getUniques(x))
# 
# track <- cbind(filtered_out, sapply(dada_fwd, getN), sapply(dada_rev, getN), sapply(mergers, getN),
#                rowSums(seqtab_nochim))
# 
# colnames(track) <- c("raw", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chim")
# rownames(track) <- sample.names  
# 
# track
# write.table(track, file = "tables/track_amplicons.txt", quote = FALSE, sep = "\t", col.names=NA)
```

Assign taxa using the HOMD database  
Databases can be found at <https://homd.org>

``` r
# Sys.time()
# 
# taxa <- assignTaxonomy(seqtab_nochim, "databases/eHOMDv15.1_FL_Compilation_TS.fa", multithread=TRUE)
# taxa <- addSpecies(taxa, "databases/HOMD_16S_rRNA_RefSeq_V15.22.fasta") #2926
# taxa <- taxa[,-8]
```

For some of the ordination methods, a tree will be needed. This step is
very time consuming, so be sure to save the tree object to load in case
some steps need to be repeated later.

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:36 EDT"

``` r
# library(phangorn)
# 
# seqs <- getSequences(seqtab_nochim)
# names(seqs) <- seqs
# alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
# 
# phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
# dm <- dist.ml(phang.align)
# treeNJ <- NJ(dm)
# fit = pml(treeNJ, data=phang.align)
# 
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                       rearrangement = "stochastic", control = pml.control(trace = 0))
# 
# save(fitGTR, file = "data/fitGTR.RData")

## if repeating analysis and no need to remake tree, use -
load("data/fitGTR.RData")

# detach("package:phangorn", unload=TRUE)
```

Prepare metadata for analysis

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:36 EDT"

``` r
sampleid <- c("MFP1", "MFP2", "MFP3", "MFP4", "MFP5", "MFP6", "MFP7", "MFP8", "MFP9", "MFP10", "MFP11", "MFP12", "MFP13",
              "SNF1", "SNF2", "SNF3", "SNF4", "SNF5", "SNF6", "SNF7", "SNF8", "SNF9", "SNF10", "SNF11", "SNF12", "SNF13", "Water")
sampletype <- c("sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample",
              "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "sample", "control")              
source <- c("Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm",
                "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "Biofilm", "control")
toothpaste <- c("MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP", "MFP",
                "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "SnF2", "control")
treatment <- c("mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp", "mfp",
               "snf", "snf", "snf", "snf", "snf", "snf", "snf", "snf", "snf", "snf", "snf", "snf", "snf", "control")
study_design <- c("Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover",
                  "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover", "Crossover",
                  "Crossover", "Crossover", "Crossover", "Crossover", "control")
sex <- c("F", "F", "F", "F", "M", "F", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "F", "M", "F", "M", "M", "M", "M", "F", "F", "M", "control")
time <- c("-14", "-14", "-14", "-14", "-14", "-14", "-14", "-14", "-14", "-14", "-14", "-14", "-14", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "control")
panelist <- c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12", "p13",
              "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12", "p13", "control")

samdf <- data.frame(sampleid=sampleid, sampletype=sampletype, source=source, toothpaste=toothpaste, treatment=treatment,
                    study_design=study_design, sex=sex, time=time, panelist=panelist)
rownames(samdf) <- samdf$sampleid
```

#### Plots and stats with phyloseq

Make phyloseq object

``` r
# Sys.time()
# 
# ps_taxa <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
#                sample_data(samdf),
#                phy_tree(fitGTR$tree),
#                tax_table(taxa))
# ps_taxa
# 
# save(ps_taxa, file = "data/ps_taxa.RData")

## if repeating analysis and no need to remake tree, use -
load("data/ps_taxa.RData")
```

Remove contaminants using negative controls

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:36 EDT"

``` r
#identify contaminants
phyloseq_object <- ps_taxa
sample_data(phyloseq_object)$is.neg <- sample_data(phyloseq_object)$sampletype == "control"
contamdf.prev <- isContaminant(phyloseq_object, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ##  2303     1

``` r
#make dataframe
phyloseq_object.neg <- prune_samples(sample_data(phyloseq_object)$sampletype == "control", phyloseq_object)
phyloseq_object.neg.presence <- transform_sample_counts(phyloseq_object.neg, function(abund) 1*(abund>0))
phyloseq_object.pos <- prune_samples(sample_data(phyloseq_object)$sampletype == "control", phyloseq_object)
phyloseq_object.pos.presence <- transform_sample_counts(phyloseq_object.pos, function(abund) 1*(abund>0))

df.pres <- data.frame(prevalence.pos=taxa_sums(phyloseq_object.pos.presence),
                      prevalence.neg=taxa_sums(phyloseq_object.neg.presence),
                      contam.prev=contamdf.prev$contaminant)

grep("TRUE", df.pres$contam.prev)
```

    ## [1] 116

``` r
df.pres[grepl("TRUE",df.pres$contam.prev),] %>% rownames() -> BadTaxa

allTaxa <- taxa_names(ps_taxa)
allTaxa <- allTaxa[!(allTaxa %in% BadTaxa)]
ps <- prune_taxa(allTaxa, ps_taxa)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2303 taxa and 27 samples ]
    ## sample_data() Sample Data:       [ 27 samples by 9 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 2303 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 2303 tips and 2301 internal nodes ]

``` r
## remove chloroplasts and mitochondria
ps <- subset_taxa(ps, (Order!="Chloroplast"))
ps <- subset_taxa(ps, (Family!="Mitochondria"))

## remove controls
ps <- subset_samples(ps, sampletype != "control")
```

##### Ordination

Ordination for Bray-Curtis, UniFrac, and Euclidean
distance/dissimilarity matrices for PCoA and NMDS methods

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:37 EDT"

``` r
ps_prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ps_prop <- prune_samples(sample_sums(ps_prop) >= 0, ps_prop)
nmds_bray <- ordinate(ps_prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.1216233 
    ## Run 1 stress 0.1216233 
    ## ... New best solution
    ## ... Procrustes: rmse 5.626093e-05  max resid 0.0001675109 
    ## ... Similar to previous best
    ## Run 2 stress 0.1216233 
    ## ... New best solution
    ## ... Procrustes: rmse 5.91364e-05  max resid 0.0002301938 
    ## ... Similar to previous best
    ## Run 3 stress 0.1216233 
    ## ... Procrustes: rmse 6.031331e-05  max resid 0.0002336768 
    ## ... Similar to previous best
    ## Run 4 stress 0.1216233 
    ## ... New best solution
    ## ... Procrustes: rmse 1.933454e-05  max resid 6.400674e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.1216233 
    ## ... Procrustes: rmse 4.124228e-05  max resid 0.0001659204 
    ## ... Similar to previous best
    ## Run 6 stress 0.1216233 
    ## ... Procrustes: rmse 2.926069e-05  max resid 0.0001154564 
    ## ... Similar to previous best
    ## Run 7 stress 0.1216233 
    ## ... Procrustes: rmse 1.556173e-05  max resid 5.31405e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.1216233 
    ## ... Procrustes: rmse 2.765624e-05  max resid 0.0001033987 
    ## ... Similar to previous best
    ## Run 9 stress 0.1216234 
    ## ... Procrustes: rmse 0.0001375814  max resid 0.0005554136 
    ## ... Similar to previous best
    ## Run 10 stress 0.1641165 
    ## Run 11 stress 0.1216233 
    ## ... Procrustes: rmse 3.752749e-05  max resid 0.0001493728 
    ## ... Similar to previous best
    ## Run 12 stress 0.1216233 
    ## ... Procrustes: rmse 9.452091e-06  max resid 3.690964e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.1216233 
    ## ... New best solution
    ## ... Procrustes: rmse 1.341617e-06  max resid 4.140509e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.1216233 
    ## ... Procrustes: rmse 5.83248e-05  max resid 0.0002320699 
    ## ... Similar to previous best
    ## Run 15 stress 0.1216233 
    ## ... Procrustes: rmse 2.624369e-05  max resid 0.0001060169 
    ## ... Similar to previous best
    ## Run 16 stress 0.1216233 
    ## ... Procrustes: rmse 6.169051e-05  max resid 0.0002484098 
    ## ... Similar to previous best
    ## Run 17 stress 0.1216233 
    ## ... Procrustes: rmse 8.099291e-06  max resid 2.129092e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.1216233 
    ## ... Procrustes: rmse 5.658131e-05  max resid 0.0002283838 
    ## ... Similar to previous best
    ## Run 19 stress 0.1600309 
    ## Run 20 stress 0.1216233 
    ## ... Procrustes: rmse 8.912695e-05  max resid 0.0003598532 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
pcoa_bray <- ordinate(ps_prop, method="PCoA", distance="bray")
nmds_uni <- ordinate(ps_prop, method="NMDS", distance="unifrac")
```

    ## Run 0 stress 0.149848 
    ## Run 1 stress 0.149848 
    ## ... New best solution
    ## ... Procrustes: rmse 4.039733e-06  max resid 1.747753e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.2182688 
    ## Run 3 stress 0.149848 
    ## ... Procrustes: rmse 4.61255e-06  max resid 1.906993e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.1529609 
    ## Run 5 stress 0.149848 
    ## ... New best solution
    ## ... Procrustes: rmse 2.336447e-06  max resid 9.091041e-06 
    ## ... Similar to previous best
    ## Run 6 stress 0.149848 
    ## ... Procrustes: rmse 1.158076e-06  max resid 2.949575e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.2179913 
    ## Run 8 stress 0.2277548 
    ## Run 9 stress 0.1529609 
    ## Run 10 stress 0.2179913 
    ## Run 11 stress 0.1498481 
    ## ... Procrustes: rmse 1.214179e-05  max resid 3.380907e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.1529609 
    ## Run 13 stress 0.149848 
    ## ... Procrustes: rmse 2.525837e-06  max resid 6.380518e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.149848 
    ## ... Procrustes: rmse 6.875092e-06  max resid 2.547469e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.2174098 
    ## Run 16 stress 0.149848 
    ## ... Procrustes: rmse 7.6409e-06  max resid 2.917646e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.2097727 
    ## Run 18 stress 0.1529609 
    ## Run 19 stress 0.1529609 
    ## Run 20 stress 0.2243417 
    ## *** Solution reached

``` r
pcoa_uni <- ordinate(ps_prop, method="PCoA", distance="unifrac")
nmds_euc <- ordinate(ps_prop, method="NMDS", distance="euclidean")
```

    ## Run 0 stress 0.09783744 
    ## Run 1 stress 0.1006778 
    ## Run 2 stress 0.1152063 
    ## Run 3 stress 0.1220985 
    ## Run 4 stress 0.1159292 
    ## Run 5 stress 0.09867647 
    ## Run 6 stress 0.1091914 
    ## Run 7 stress 0.1052533 
    ## Run 8 stress 0.1171467 
    ## Run 9 stress 0.125427 
    ## Run 10 stress 0.1006778 
    ## Run 11 stress 0.09783747 
    ## ... Procrustes: rmse 8.534115e-05  max resid 0.0002952907 
    ## ... Similar to previous best
    ## Run 12 stress 0.09783747 
    ## ... Procrustes: rmse 9.523748e-05  max resid 0.000412945 
    ## ... Similar to previous best
    ## Run 13 stress 0.09783744 
    ## ... Procrustes: rmse 9.617684e-06  max resid 4.105059e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.1112083 
    ## Run 15 stress 0.1152829 
    ## Run 16 stress 0.1192242 
    ## Run 17 stress 0.1006779 
    ## Run 18 stress 0.1210003 
    ## Run 19 stress 0.1293761 
    ## Run 20 stress 0.1220985 
    ## *** Solution reached

``` r
pcoa_euc <- ordinate(ps_prop, method="PCoA", distance="euclidean")
```

##### Analysis of variance

Uses vegan adonis test

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:38 EDT"

``` r
## fit linear models for pc1 and pc2 for factors of interest
pcoa <- as.data.frame(pcoa_bray$vectors)
pcoa$merge <- rownames(pcoa)
pcoa_factors <- merge(pcoa, samdf, by.x = "merge", by.y = "sampleid")

PC1_tr <- lm(Axis.1 ~ treatment, data = pcoa_factors)
print("PC1 treatment")
```

    ## [1] "PC1 treatment"

``` r
anova(PC1_tr)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Axis.1
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## treatment  1 0.04738 0.047379  0.7069 0.4088
    ## Residuals 24 1.60865 0.067027

``` r
summary(PC1_tr)
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ treatment, data = pcoa_factors)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41455 -0.20987  0.03653  0.19881  0.42121 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   0.04269    0.07180   0.595    0.558
    ## treatmentsnf -0.08538    0.10155  -0.841    0.409
    ## 
    ## Residual standard error: 0.2589 on 24 degrees of freedom
    ## Multiple R-squared:  0.02861,    Adjusted R-squared:  -0.01186 
    ## F-statistic: 0.7069 on 1 and 24 DF,  p-value: 0.4088

``` r
PC1_sx <- lm(Axis.1 ~ sex, data = pcoa_factors)
print("PC1 sex")
```

    ## [1] "PC1 sex"

``` r
anova(PC1_sx)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Axis.1
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## sex        1 0.00031 0.000309  0.0045 0.9472
    ## Residuals 24 1.65572 0.068989

``` r
summary(PC1_sx)
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ sex, data = pcoa_factors)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36867 -0.24459  0.01565  0.22022  0.38629 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) -0.003190   0.070198  -0.045    0.964
    ## sexM         0.006912   0.103329   0.067    0.947
    ## 
    ## Residual standard error: 0.2627 on 24 degrees of freedom
    ## Multiple R-squared:  0.0001864,  Adjusted R-squared:  -0.04147 
    ## F-statistic: 0.004475 on 1 and 24 DF,  p-value: 0.9472

``` r
PC1_tm <- lm(Axis.1 ~ time, data = pcoa_factors)
print("PC1 time")
```

    ## [1] "PC1 time"

``` r
anova(PC1_tm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Axis.1
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## time       1 0.04738 0.047379  0.7069 0.4088
    ## Residuals 24 1.60865 0.067027

``` r
summary(PC1_tm)
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ time, data = pcoa_factors)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41455 -0.20987  0.03653  0.19881  0.42121 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  0.04269    0.07180   0.595    0.558
    ## time0       -0.08538    0.10155  -0.841    0.409
    ## 
    ## Residual standard error: 0.2589 on 24 degrees of freedom
    ## Multiple R-squared:  0.02861,    Adjusted R-squared:  -0.01186 
    ## F-statistic: 0.7069 on 1 and 24 DF,  p-value: 0.4088

``` r
PC1_pn <- lm(Axis.1 ~ panelist, data = pcoa_factors)
print("PC1 panelist")
```

    ## [1] "PC1 panelist"

``` r
anova(PC1_pn)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Axis.1
    ##           Df  Sum Sq  Mean Sq F value    Pr(>F)    
    ## panelist  12 1.51191 0.125993  11.365 5.283e-05 ***
    ## Residuals 13 0.14412 0.011086                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PC1_pn)
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ panelist, data = pcoa_factors)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.13607 -0.06106  0.00000  0.06106  0.13607 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  0.06848    0.07445   0.920  0.37444   
    ## panelistp10 -0.40189    0.10529  -3.817  0.00214 **
    ## panelistp11 -0.36909    0.10529  -3.505  0.00387 **
    ## panelistp12 -0.17325    0.10529  -1.645  0.12384   
    ## panelistp13  0.18807    0.10529   1.786  0.09739 . 
    ## panelistp2   0.10184    0.10529   0.967  0.35107   
    ## panelistp3  -0.34968    0.10529  -3.321  0.00552 **
    ## panelistp4   0.10509    0.10529   0.998  0.33643   
    ## panelistp5   0.31578    0.10529   2.999  0.01025 * 
    ## panelistp6   0.18340    0.10529   1.742  0.10513   
    ## panelistp7  -0.33680    0.10529  -3.199  0.00698 **
    ## panelistp8   0.08360    0.10529   0.794  0.44146   
    ## panelistp9  -0.23731    0.10529  -2.254  0.04210 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1053 on 13 degrees of freedom
    ## Multiple R-squared:  0.913,  Adjusted R-squared:  0.8326 
    ## F-statistic: 11.37 on 12 and 13 DF,  p-value: 5.283e-05

``` r
## adonis on each factor
ps_clean <- prune_taxa(!(taxa_names(ps) %in% "-1"), ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0.0, ps_clean)
ps_clean <- transform_sample_counts(ps_clean, function(x) x/sum(x))
ps_clean <- prune_samples(sample_sums(ps_clean) >= 0, ps_clean)
df <- data.frame(sample_data(ps_clean))
ps_clean_bray <- phyloseq::distance(ps_clean, method = "bray")
adonis_bt <- adonis2(ps_clean_bray ~ treatment, data = df)
adonis_bt
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_bray ~ treatment, data = df)
    ##           Df SumOfSqs      R2     F Pr(>F)
    ## treatment  1   0.1605 0.03202 0.794  0.589
    ## Residual  24   4.8529 0.96798             
    ## Total     25   5.0134 1.00000

``` r
adonis2(ps_clean_bray ~ sex, data = df)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_bray ~ sex, data = df)
    ##          Df SumOfSqs     R2      F Pr(>F)  
    ## sex       1   0.3389 0.0676 1.7401  0.093 .
    ## Residual 24   4.6745 0.9324                
    ## Total    25   5.0134 1.0000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis2(ps_clean_bray ~ time, data = df)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_bray ~ time, data = df)
    ##          Df SumOfSqs      R2     F Pr(>F)
    ## time      1   0.1605 0.03202 0.794  0.596
    ## Residual 24   4.8529 0.96798             
    ## Total    25   5.0134 1.00000

``` r
adonis_bp <- adonis2(ps_clean_bray ~ panelist, data = df)
adonis_bp
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_bray ~ panelist, data = df)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## panelist 12   4.3052 0.85873 6.5851  0.001 ***
    ## Residual 13   0.7083 0.14127                  
    ## Total    25   5.0134 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## you can also combine factors to look at interactions
adonis2(ps_clean_bray ~ treatment + sex, data = df)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_bray ~ treatment + sex, data = df)
    ##           Df SumOfSqs      R2      F Pr(>F)
    ## treatment  1   0.1605 0.03202 0.8180  0.571
    ## sex        1   0.3389 0.06760 1.7269  0.105
    ## Residual  23   4.5140 0.90037              
    ## Total     25   5.0134 1.00000

``` r
## adding euclidean and unifrac
ps_clean_euc <- phyloseq::distance(ps_clean, method = "euclidean")
ps_clean_uni <- phyloseq::distance(ps_clean, method = "unifrac")
adonis_et <- adonis2(ps_clean_euc ~ treatment, data = df)
adonis_et
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_euc ~ treatment, data = df)
    ##           Df SumOfSqs      R2      F Pr(>F)
    ## treatment  1  0.01933 0.01317 0.3202  0.954
    ## Residual  24  1.44855 0.98683              
    ## Total     25  1.46788 1.00000

``` r
adonis_ut <- adonis2(ps_clean_uni ~ treatment, data = df)
adonis_ut
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = ps_clean_uni ~ treatment, data = df)
    ##           Df SumOfSqs      R2      F Pr(>F)
    ## treatment  1  0.09248 0.03064 0.7586  0.895
    ## Residual  24  2.92595 0.96936              
    ## Total     25  3.01844 1.00000

``` r
## once you have your stats, you can insert them into the appropriate plot

## write tables for stats of interest
adonis2(ps_clean_bray ~ treatment, data = df) %>%
  write.table(., file = "tables/adonis_treatment_bray.txt", quote = FALSE, sep = "\t", col.names=NA, na = "")
```

##### Plot PCoA and NMDS

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:38 EDT"

``` r
## beta diversity plots using ordination output and anova output
bray <- plot_ordination(ps_prop, pcoa_bray, color="treatment", title="Bray PCoA") +
  theme_bw() + ggtitle("Beta Diversity: Bray-Curtis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  xlab(paste0("PC1: ", percent(pcoa_bray$values$Relative_eig[[1]], accuracy = 0.01))) + 
  ylab(paste0("PC2: ", percent(pcoa_bray$values$Relative_eig[[2]], accuracy = 0.01))) + 
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0,0,0,0))) +
  annotate("text", x = 0.01, y = 0.4, colour = "black", hjust = 0, size = 3, 
           label = paste0("adonis R2: ", round(adonis_bt$R2[1], digits = 3), "; p-value: ", round(adonis_bt$`Pr(>F)`[1], digits = 3)))
ggsave("plots/pcoa_treatment.pdf", height = 4.5, width = 5.5)
ggsave("plots/pcoa_treatment.jpeg", height = 4.5, width = 5.5, dpi = 500)
bray
```

![](index_files/figure-gfm/phyloseq%20beta%20diveristy%20plots-1.png)<!-- -->

``` r
braya <- plot_ordination(ps_prop, pcoa_bray, color="treatment", title="Bray PCoA") +
  geom_path(aes(group = panelist), color = "grey", arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  theme_bw() + ggtitle("Beta Diversity: Bray-Curtis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  xlab(paste0("PC1: ", percent(pcoa_bray$values$Relative_eig[[1]], accuracy = 0.01))) + 
  ylab(paste0("PC2: ", percent(pcoa_bray$values$Relative_eig[[2]], accuracy = 0.01))) + 
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0,0,0,0))) +
  annotate("text", x = 0.01, y = 0.4, colour = "black", hjust = 0, size = 3, 
           label = paste0("adonis R2: ", round(adonis_bt$R2[1], digits = 3), "; p-value: ", round(adonis_bt$`Pr(>F)`[1], digits = 3)))
ggsave("plots/pcoa_treatment_arrows.pdf", height = 4.5, width = 5.5)
ggsave("plots/pcoa_treatment_arrows.jpeg", height = 4.5, width = 5.5, dpi = 500)
braya
```

![](index_files/figure-gfm/phyloseq%20beta%20diveristy%20plots-2.png)<!-- -->

``` r
## make sure to check that you are adding the correct stats to the correct plot  
plot_ordination(ps_prop, nmds_bray, color="treatment", title="Bray NMDS") +
  theme_bw() + ggtitle("Beta Diversity: Bray-Curtis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.2, y = 0.76, colour = "darkgrey", 
           label = paste0("adonis R2: ", round(adonis_bt$R2[1], digits = 3), "; p-value: ", round(adonis_bt$`Pr(>F)`[1], digits = 3)))
```

![](index_files/figure-gfm/phyloseq%20beta%20diveristy%20plots-3.png)<!-- -->

``` r
ggsave("plots/nmds_treatment.pdf", height = 4.5, width = 5.5)
ggsave("plots/nmds_treatment.jpeg", height = 4.5, width = 5.5, dpi = 500)
plot_ordination(ps_prop, nmds_bray, color="treatment", title="Bray NMDS") +
  geom_path(aes(group = panelist), color = "grey", arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  theme_bw() + ggtitle("Beta Diversity: Bray-Curtis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.2, y = 0.75, colour = "darkgrey", 
           label = paste0("adonis R2: ", round(adonis_bt$R2[1], digits = 3), "; p-value: ", round(adonis_bt$`Pr(>F)`[1], digits = 3)))
```

![](index_files/figure-gfm/phyloseq%20beta%20diveristy%20plots-4.png)<!-- -->

``` r
ggsave("plots/nmds_treatment_arrows.pdf", height = 4.5, width = 5.5)
ggsave("plots/nmds_treatment_arrows.jpeg", height = 4.5, width = 5.5, dpi = 500)
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:44 EDT"

``` r
euc <- plot_ordination(ps_prop, pcoa_euc, color="treatment", title="Euclidean PCoA") +
  theme_bw() + ggtitle("Beta Diversity: Euclidean") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  xlab(paste0("PC1: ", percent(pcoa_euc$values$Relative_eig[[1]], accuracy = 0.01))) + 
  ylab(paste0("PC2: ", percent(pcoa_euc$values$Relative_eig[[2]], accuracy = 0.01))) + 
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0,0,0,0))) +
  annotate("text", x = 0.01, y = 0.19, colour = "black", hjust = 0, size = 3, 
           label = paste0("adonis R2: ", round(adonis_et$R2[1], digits = 3), "; p-value: ", round(adonis_et$`Pr(>F)`[1], digits = 3)))
ggsave("plots/pcoa_treatment_euc.pdf", height = 5.5, width = 6.5)
ggsave("plots/pcoa_treatment_euc.jpeg", height = 5.5, width = 6.5, dpi = 500)
euca <- plot_ordination(ps_prop, pcoa_euc, color="treatment", title="Euclidean PCoA") +
  geom_path(aes(group = panelist), color = "grey", arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  theme_bw() + ggtitle("Beta Diversity: Euclidean") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  xlab(paste0("PC1: ", percent(pcoa_euc$values$Relative_eig[[1]], accuracy = 0.01))) + 
  ylab(paste0("PC2: ", percent(pcoa_euc$values$Relative_eig[[2]], accuracy = 0.01))) + 
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0,0,0,0))) +
  annotate("text", x = 0.01, y = 0.19, colour = "black", hjust = 0, size = 3, 
           label = paste0("adonis R2: ", round(adonis_et$R2[1], digits = 3), "; p-value: ", round(adonis_et$`Pr(>F)`[1], digits = 3)))
ggsave("plots/pcoa_treatment_arrow_euc.pdf", height = 5.5, width = 6.5)
ggsave("plots/pcoa_treatment_arrow_euc.jpeg", height = 5.5, width = 6.5, dpi = 500)

## make sure to check that you are adding the correct stats to the correct plot  
plot_ordination(ps_prop, nmds_euc, color="treatment", title="Euclidean NMDS") +
  theme_bw() + ggtitle("Beta Diversity: Euclidean") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.19, y = 0.22, colour = "darkgrey", 
           label = paste0("adonis R2: ", round(adonis_et$R2[1], digits = 3), "; p-value: ", round(adonis_et$`Pr(>F)`[1], digits = 3)))
```

![](index_files/figure-gfm/plots%20for%20euc%20and%20uni-1.png)<!-- -->

``` r
ggsave("plots/nmds_treatment_euc.pdf", height = 5.5, width = 6.5)
ggsave("plots/nmds_treatment_euc.jpeg", height = 5.5, width = 6.5, dpi = 500)
plot_ordination(ps_prop, nmds_euc, color="treatment", title="Euclidean NMDS") +
  geom_path(aes(group = panelist), color = "grey", arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  theme_bw() + ggtitle("Beta Diversity: Euclidean") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.19, y = 0.22, colour = "darkgrey", 
           label = paste0("adonis R2: ", round(adonis_et$R2[1], digits = 3), "; p-value: ", round(adonis_et$`Pr(>F)`[1], digits = 3)))
```

![](index_files/figure-gfm/plots%20for%20euc%20and%20uni-2.png)<!-- -->

``` r
ggsave("plots/nmds_treatment_arrow_euc.pdf", height = 5.5, width = 6.5)
ggsave("plots/nmds_treatment_arrow_euc.jpeg", height = 5.5, width = 6.5, dpi = 500)

uni <- plot_ordination(ps_prop, pcoa_uni, color="treatment", title="Unifrac PCoA") +
  theme_bw() + ggtitle("Beta Diversity: Unifrac") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  xlab(paste0("PC1: ", percent(pcoa_uni$values$Relative_eig[[1]], accuracy = 0.01))) + 
  ylab(paste0("PC2: ", percent(pcoa_uni$values$Relative_eig[[2]], accuracy = 0.01))) + 
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0,0,0,0))) +
  annotate("text", x = 0.01, y = 0.19, colour = "black", hjust = 0, size = 3, 
           label = paste0("adonis R2: ", round(adonis_ut$R2[1], digits = 3), "; p-value: ", round(adonis_ut$`Pr(>F)`[1], digits = 3)))
ggsave("plots/pcoa_treatment_uni.pdf", height = 5.5, width = 6.5)
ggsave("plots/pcoa_treatment_uni.jpeg", height = 5.5, width = 6.5, dpi = 500)
unia <- plot_ordination(ps_prop, pcoa_uni, color="treatment", title="Unifrac PCoA") +
  geom_path(aes(group = panelist), color = "grey", arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  theme_bw() + ggtitle("Beta Diversity: Unifrac") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  xlab(paste0("PC1: ", percent(pcoa_uni$values$Relative_eig[[1]], accuracy = 0.01))) + 
  ylab(paste0("PC2: ", percent(pcoa_uni$values$Relative_eig[[2]], accuracy = 0.01))) + 
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "top", legend.title=element_blank(), legend.justification = "right",
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,-10,-10), legend.spacing.x = unit(0, "cm")) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0,0,0,0))) +
  annotate("text", x = 0.01, y = 0.19, colour = "black", hjust = 0, size = 3, 
           label = paste0("adonis R2: ", round(adonis_ut$R2[1], digits = 3), "; p-value: ", round(adonis_ut$`Pr(>F)`[1], digits = 3)))
ggsave("plots/pcoa_treatment_arrow_uni.pdf", height = 5.5, width = 6.5)
ggsave("plots/pcoa_treatment_arrow_uni.jpeg", height = 5.5, width = 6.5, dpi = 500)

## make sure to check that you are adding the correct stats to the correct plot  
plot_ordination(ps_prop, nmds_uni, color="treatment", title="Unifrac NMDS") +
  theme_bw() + ggtitle("Beta Diversity: Unifrac") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.15, y = 0.2, colour = "darkgrey", 
           label = paste0("adonis R2: ", round(adonis_ut$R2[1], digits = 3), "; p-value: ", round(adonis_ut$`Pr(>F)`[1], digits = 3)))
```

![](index_files/figure-gfm/plots%20for%20euc%20and%20uni-3.png)<!-- -->

``` r
ggsave("plots/nmds_treatment_uni.pdf", height = 5.5, width = 6.5)
ggsave("plots/nmds_treatment_uni.jpeg", height = 5.5, width = 6.5, dpi = 500)
plot_ordination(ps_prop, nmds_uni, color="treatment", title="Unifrac NMDS") +
  geom_path(aes(group = panelist), color = "grey", arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  theme_bw() + ggtitle("Beta Diversity: Unifrac") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous"),
                     labels = c("MFP", expression("SnF"["2"]))) +
  geom_hline(yintercept = 0, color = "lightgrey") +
  geom_vline(xintercept = 0, color = "lightgrey") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.15, y = 0.2, colour = "darkgrey", 
           label = paste0("adonis R2: ", round(adonis_ut$R2[1], digits = 3), "; p-value: ", round(adonis_ut$`Pr(>F)`[1], digits = 3)))
```

![](index_files/figure-gfm/plots%20for%20euc%20and%20uni-4.png)<!-- -->

``` r
ggsave("plots/nmds_treatment_arrow_uni.pdf", height = 5.5, width = 6.5)
ggsave("plots/nmds_treatment_arrow_uni.jpeg", height = 5.5, width = 6.5, dpi = 500)
```

##### Within sample diversity

Uses both phyloseq and vegan to calculate numerous alpha diversity
metrics and output to a table, then plot Observed species, Shannon
index, and Simpson index

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:49 EDT"

``` r
## alpha diversity in vegan - add to sample data
samdf_wstats <- data.frame(sample_data(ps_clean))
asv_table <- ps_clean@otu_table
samdf_wstats$shannon_vegan <- vegan::diversity(asv_table, index="shannon")
samdf_wstats$simpson_vegan <- vegan::diversity(asv_table, index="simpson")
samdf_wstats$invsimpson_vegan <- vegan::diversity(asv_table, index="inv")
## add alpha diversity from phyloseq to sample data
## these save as df within the column and include standard error columns
samdf_wstats[,ncol(samdf_wstats)+1] <- estimate_richness(ps, measures="Observed")
samdf_wstats[,ncol(samdf_wstats)+1] <- estimate_richness(ps, measures="ACE")
samdf_wstats[,ncol(samdf_wstats)+1] <- estimate_richness(ps, measures="Chao1")
samdf_wstats[,ncol(samdf_wstats)+1] <- estimate_richness(ps, measures="Fisher")
samdf_wstats[,ncol(samdf_wstats)+1] <- estimate_richness(ps, measures="Shannon")
samdf_wstats[,ncol(samdf_wstats)+1] <- estimate_richness(ps, measures="Simpson")

## add alpha diversity metrics to the phyloseq sample data table (no need to do both this and above)
ps_alpha <- ps
sample_data(ps_alpha)$shannon_physeq <- estimate_richness(ps_alpha, measures="Shannon")
sample_data(ps_alpha)$simpson_physeq <- estimate_richness(ps_alpha, measures="Simpson")
sample_data(ps_alpha)$observed_physeq <- estimate_richness(ps_alpha, measures="Observed")

## test for normalcy - microbiome data is rarely normal
shapiro.test(samdf_wstats$shannon_vegan)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$shannon_vegan
    ## W = 0.9672, p-value = 0.5522

``` r
shapiro.test(samdf_wstats$simpson_vegan)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$simpson_vegan
    ## W = 0.87441, p-value = 0.004397

``` r
shapiro.test(samdf_wstats$invsimpson_vegan)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$invsimpson_vegan
    ## W = 0.94364, p-value = 0.164

``` r
shapiro.test(samdf_wstats$Observed)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$Observed
    ## W = 0.92095, p-value = 0.04733

``` r
shapiro.test(samdf_wstats$ACE)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$ACE
    ## W = 0.92086, p-value = 0.0471

``` r
shapiro.test(samdf_wstats$Chao1)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$Chao1
    ## W = 0.92105, p-value = 0.04757

``` r
shapiro.test(samdf_wstats$Fisher)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  samdf_wstats$Fisher
    ## W = 0.91625, p-value = 0.03676

``` r
## calculate anova from vegan alpha diversity stats (best for normal)
## in this case, inverted simpson was the closest to normal
inv_tp <- aov(invsimpson_vegan ~ treatment, data = samdf_wstats)
summary(inv_tp)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    1    3.2   3.159   0.136  0.715
    ## Residuals   24  556.0  23.166

``` r
## tukey's test on anova for pairwise comparisons
TukeyHSD(inv_tp)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = invsimpson_vegan ~ treatment, data = samdf_wstats)
    ## 
    ## $treatment
    ##               diff       lwr      upr     p adj
    ## snf-mfp -0.6971869 -4.593519 3.199146 0.7151411

``` r
## kruskal wallace test on vegan alpha diversity stats (best for non-normal)
kruskal.test(shannon_vegan ~ treatment, data = samdf_wstats)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  shannon_vegan by treatment
    ## Kruskal-Wallis chi-squared = 0.00065746, df = 1, p-value = 0.9795

``` r
kruskal.test(shannon_vegan ~ sex, data = samdf_wstats)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  shannon_vegan by sex
    ## Kruskal-Wallis chi-squared = 1.0582, df = 1, p-value = 0.3036

``` r
kruskal.test(shannon_vegan ~ panelist, data = samdf_wstats)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  shannon_vegan by panelist
    ## Kruskal-Wallis chi-squared = 23.026, df = 12, p-value = 0.02751

``` r
## wilcoxon rank sum test is used for within groups
pairwise.wilcox.test(samdf_wstats$shannon_vegan, samdf_wstats$treatment, p.adjust.method = "fdr")
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  samdf_wstats$shannon_vegan and samdf_wstats$treatment 
    ## 
    ##     mfp
    ## snf 1  
    ## 
    ## P value adjustment method: fdr

``` r
pairwise.wilcox.test(samdf_wstats$simpson_vegan, samdf_wstats$treatment, p.adjust.method = "fdr")
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  samdf_wstats$simpson_vegan and samdf_wstats$treatment 
    ## 
    ##     mfp 
    ## snf 0.96
    ## 
    ## P value adjustment method: fdr

``` r
pairwise.wilcox.test(samdf_wstats$Observed, samdf_wstats$treatment, p.adjust.method = "fdr")
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  samdf_wstats$Observed and samdf_wstats$treatment 
    ## 
    ##     mfp 
    ## snf 0.61
    ## 
    ## P value adjustment method: fdr

``` r
## alpha diversity (does not change with annotation method)
plot_richness(ps, x="treatment", measures=c("Shannon", "Simpson", "Observed"), color="sex") +
  theme_bw() + ggtitle("Alpha Diversity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3, type = "continuous"))
```

![](index_files/figure-gfm/phyloseq%20alpha%20diversity%20plot-1.png)<!-- -->

``` r
ggsave("plots/alpha-div.pdf")
ggsave("plots/alpha-div.jpeg", dpi = 500)

write.table(samdf_wstats, "tables/samples_wstats.txt", quote = FALSE, sep = "\t", row.names = FALSE)

## remake alpha as boxplots
## add stats to plots
plota <- samdf_wstats %>% select(sampleid, treatment, Observed, Shannon, Simpson)
plota %>%
  pivot_longer(cols = 3:5, names_to = "stat", values_to = "value") -> plota
alpha <- ggplot(plota, aes(x = treatment, y = value, color = treatment)) +
  geom_jitter(alpha = 0.5, size = 0.3) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Measure", color = "") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2, type = "continuous")) +
  scale_x_discrete(labels=c("mfp" = "MFP", "snf" = expression("SnF"["2"]))) +
  theme(legend.position="none", axis.title.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ stat, scales = "free_y")
ggsave("plots/alpha_boxplots_legend.jpeg", height = 5, width = 6.5, dpi = 1000)

# wes_palette("Darjeeling1", n = 2, type = "continuous")[]
```

##### Taxonomy

Plot taxa by all levels from Phyla down to ASV

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:50 EDT"

``` r
## ggplot of top 10 by taxa level grouped by factor and showing percent abundance
ps_melt <- psmelt(ps)

phylum <- ps_melt %>% select(Sample, Phylum, Abundance) %>% dcast(Phylum ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,2:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- phylum %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)
unique(ps_melt$Phylum)
```

    ## [1] "Firmicutes"               "Proteobacteria"          
    ## [3] "Bacteroidetes"            "Actinobacteria"          
    ## [5] "Fusobacteria"             "Absconditabacteria (SR1)"
    ## [7] "Synergistetes"            "Saccharibacteria (TM7)"

``` r
p <- ggplot(filt, aes(x = treatment, y = percentage, fill = Phylum)) +
  theme_bw() +
  geom_bar(position="fill", stat="summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  scale_x_discrete(labels=c("mfp" = "MFP", "snf" = expression("SnF"["2"]))) +
  labs(y = "Percent abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
p
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-1.png)<!-- -->

``` r
# ggsave("plots/bar_phylum.pdf")
# ggsave("plots/bar_phylum.jpeg", dpi = 500, height = 4, width = 5)

class <- ps_melt %>% select(Sample, Class, Abundance) %>% dcast(Class ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,2:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- class %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

ggplot(filt, aes(x = treatment, y = percentage, fill = Class)) +
  theme_bw() +
  geom_bar(stat = "summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(x = "Toothpaste", y = "Percent abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-2.png)<!-- -->

``` r
# ggsave("plots/bar_class.pdf")
# ggsave("plots/bar_class.jpeg", dpi = 500, height = 4, width = 5)

order <- ps_melt %>% select(Sample, Order, Abundance) %>% dcast(Order ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,2:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- order %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

ggplot(filt, aes(x = treatment, y = percentage, fill = Order)) +
  theme_bw() +
  geom_bar(stat = "summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(x = "Toothpaste", y = "Percent abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-3.png)<!-- -->

``` r
# ggsave("plots/bar_order.pdf")
# ggsave("plots/bar_order.jpeg", dpi = 500, height = 4, width = 5)

family <- ps_melt %>% select(Sample, Family, Abundance) %>% dcast(Family ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,2:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- family %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

ggplot(filt, aes(x = treatment, y = percentage, fill = Family)) +
  theme_bw() +
  geom_bar(stat = "summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(x = "Toothpaste", y = "Percent abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-4.png)<!-- -->

``` r
# ggsave("plots/bar_family.pdf")
# ggsave("plots/bar_family.jpeg", dpi = 500, height = 4, width = 5)

genus <- ps_melt %>% select(Sample, Genus, Abundance) %>% dcast(Genus ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,2:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- genus %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

g <- ggplot(filt, aes(x = treatment, y = percentage, fill = Genus)) +
  theme_bw() +
  geom_bar(position="fill", stat="summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  scale_x_discrete(labels=c("mfp" = "MFP", "snf" = expression("SnF"["2"]))) +
  labs(y = "Percent abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
g
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-5.png)<!-- -->

``` r
# ggsave("plots/bar_genus.pdf")
# ggsave("plots/bar_genus.jpeg", dpi = 500, height = 4, width = 5)

genus <- ps_melt %>% select(Sample, Genus, Abundance) %>% dcast(Genus ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.)
write.table(genus, "tables/genera.txt", quote = FALSE, row.names = FALSE, sep = "\t")

## for top 10 asvs, add on a simple asv name to replace nucleotides first
ps2 <- ps
taxa_names(ps2) <- paste("ASV", 1:ntaxa(ps2), sep = "")
ps2_melt <- psmelt(ps2)
asv <- ps_melt %>% select(Sample, OTU, Abundance) %>% dcast(OTU ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,2:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- asv %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

ggplot(filt, aes(x = treatment, y = percentage, fill = OTU)) +
  theme_bw() +
  geom_bar(stat = "summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(x = "Toothpaste", y = "Percent abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-6.png)<!-- -->

``` r
# ggsave("plots/bar_asv.pdf")
# ggsave("plots/bar_asv.jpeg", dpi = 500, height = 4, width = 5)

## by species
ps2_melt$Species <- ifelse(is.na(ps2_melt$Species), ps2_melt$OTU, ps2_melt$Species)
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

gs <- ggplot(filt, aes(x = treatment, y = percentage, fill = paste(Genus, Species))) +
  theme_bw() +
  geom_bar(position="fill", stat="summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", y = "Percent abundance") +
  scale_x_discrete(labels=c("mfp" = "MFP", "snf" = expression("SnF"["2"]))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(labels = scales::percent)
gs
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-7.png)<!-- -->

``` r
# ggsave("plots/bar_species.pdf")
# ggsave("plots/bar_species.jpeg", dpi = 500, height = 4.7, width = 5)

# species <- dcast(ps2_melt[,c(2,20,21,3)], Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.)
# write.table(genus, "tables/species.txt", quote = FALSE, row.names = FALSE, sep = "\t")

## by species
ps2_melt$Species <- ifelse(is.na(ps2_melt$Species), ps2_melt$OTU, ps2_melt$Species)
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 25, rowSums(.[,3:ncol(.)])) %>% melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(percentage = value/sum(value)*100)

gs20 <- ggplot(filt, aes(x = treatment, y = percentage, fill = paste(Genus, Species))) +
  theme_bw() +
  geom_bar(position="fill", stat="summary") +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 25, type = "continuous")) +
  labs(fill = "Genus species", y = "Percent abundance") +
  scale_x_discrete(labels=c("mfp" = "MFP", "snf" = expression("SnF"["2"]))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.key.size = unit(0.25, 'cm'))
gs20
```

![](index_files/figure-gfm/phyloseq%20taxonomy%20plots-8.png)<!-- -->

``` r
ggsave("plots/bar_species_25.pdf")
ggsave("plots/bar_species_25.jpeg", dpi = 500, height = 4.7, width = 8)
```

Additionally, plot specific taxonomic groups of interest

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:30:58 EDT"

``` r
ps2_melt$Species <- ifelse(is.na(ps2_melt$Species), ps2_melt$OTU, ps2_melt$Species)
ps2_melt <- separate(ps2_melt, Species, into = "Species", remove = TRUE, sep = "_")
ps2_melt <- separate(ps2_melt, Genus, into = "Genus", remove = TRUE, sep = "\\[")

## Prevotella
species <- subset(ps2_melt, Genus=="Prevotella")
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>%
  melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(cpm = value/sum(value)*1000000)

ggplot(filt, aes(x = treatment, y = log(cpm), fill = paste(Genus, Species))) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", x = "Toothpaste", y = "Relative abundance (cpm)") +
  labs(title = "Top 10 Prevotella (of 236 ASVs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![](index_files/figure-gfm/individual%20taxa%20plots-1.png)<!-- -->

``` r
ggsave("plots/box_prevotella.pdf")
ggsave("plots/box_prevotella.jpeg", dpi = 500, height = 4.5, width = 5.25)

## Porphyromonas
species <- subset(ps2_melt, Genus=="Porphyromonas")
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>%
  melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(cpm = value/sum(value)*1000000)

ggplot(filt, aes(x = treatment, y = log(cpm), fill = paste(Genus, Species))) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", x = "Toothpaste", y = "Relative abundance (cpm)") +
  labs(title = "Top 10 Porphyromonas (of 73 ASVs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![](index_files/figure-gfm/individual%20taxa%20plots-2.png)<!-- -->

``` r
ggsave("plots/box_porphyromonas.pdf")
ggsave("plots/box_porphyromonas.jpeg", dpi = 500, height = 4.5, width = 5.25)

## Porphyromonas gingivalis
species <- subset(ps2_melt, Genus=="Porphyromonas" & Species=="gingivalis")
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>%
  melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(cpm = value/sum(value)*1000000)

ggplot(filt, aes(x = treatment, y = log(cpm), fill = paste(Genus, Species))) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", x = "Toothpaste", y = "Relative abundance (cpm)") +
  labs(title = "Porphyromonas gingivalis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![](index_files/figure-gfm/individual%20taxa%20plots-3.png)<!-- -->

``` r
ggsave("plots/box_pgingivalis.pdf")
ggsave("plots/box_pgingivalis.jpeg", dpi = 500, height = 4.5, width = 5.25)

## Firmicutes
species <- subset(ps2_melt, Phylum=="Firmicutes")
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>%
  melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(cpm = value/sum(value)*1000000)

ggplot(filt, aes(x = treatment, y = log(cpm), fill = paste(Genus, Species))) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", x = "Toothpaste", y = "Relative abundance (cpm)") +
  labs(title = "Top 10 Firmicutes (of 434 ASVs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![](index_files/figure-gfm/individual%20taxa%20plots-4.png)<!-- -->

``` r
ggsave("plots/box_firmicutes.pdf")
ggsave("plots/box_firmicutes.jpeg", dpi = 500, height = 4.5, width = 5.25)

## Bacteroidetes
species <- subset(ps2_melt, Phylum=="Bacteroidetes")
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>%
  melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(cpm = value/sum(value)*1000000)

ggplot(filt, aes(x = treatment, y = log(cpm), fill = paste(Genus, Species))) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", x = "Toothpaste", y = "Relative abundance (cpm)") +
  labs(title = "Top 10 Bacteroidetes (of 442 ASVs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![](index_files/figure-gfm/individual%20taxa%20plots-5.png)<!-- -->

``` r
ggsave("plots/box_bacteroidetes.pdf")
ggsave("plots/box_bacteroidetes.jpeg", dpi = 500, height = 4.5, width = 5.25)

## Spirochaetota - none identified in data

## Actinobacteria
species <- subset(ps2_melt, Phylum=="Actinobacteria")
species <- ps_melt %>% select(Sample, Genus, Species, Abundance) %>% dcast(Genus + Species ~ Sample, value.var = 'Abundance', fun.aggregate = sum, na.rm = TRUE) %>% na.omit(.) %>%
  top_n(., 10, rowSums(.[,3:ncol(.)])) %>%
  melt(.) %>% merge(., samdf, by.x = "variable", by.y = "sampleid")
filt <- species %>% group_by(treatment) %>% dplyr::mutate(cpm = value/sum(value)*1000000)

ggplot(filt, aes(x = treatment, y = log(cpm), fill = paste(Genus, Species))) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 10, type = "continuous")) +
  labs(fill = "Genus species", x = "Toothpaste", y = "Relative abundance (cpm)") +
  labs(title = "Top 10 Actinobacteria (of 149 ASVs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![](index_files/figure-gfm/individual%20taxa%20plots-6.png)<!-- -->

``` r
ggsave("plots/box_actinobacteria.pdf")
ggsave("plots/box_actinobacteria.jpeg", dpi = 500, height = 4.5, width = 5.25)
```

##### Rarefaction curve

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:31:07 EDT"

``` r
rarecurve(as.data.frame(ps@otu_table), step = 100, sampleid, xlab = "Sample Size", ylab = "Species",
          label = TRUE)
```

![](index_files/figure-gfm/rarefaction%20curve-1.png)<!-- -->

``` r
## wrapper to start with phyloseq, rarefaction with vegan rarecurve(), and plot with ggplot2
p <- ggrare(ps, step = 100, color = "treatment", label = "sampleid", se = TRUE)
```

    ## rarefying sample MFP1
    ## rarefying sample MFP10
    ## rarefying sample MFP11
    ## rarefying sample MFP12
    ## rarefying sample MFP13
    ## rarefying sample MFP2
    ## rarefying sample MFP3
    ## rarefying sample MFP4
    ## rarefying sample MFP5
    ## rarefying sample MFP6
    ## rarefying sample MFP7
    ## rarefying sample MFP8
    ## rarefying sample MFP9
    ## rarefying sample SNF1
    ## rarefying sample SNF10
    ## rarefying sample SNF11
    ## rarefying sample SNF12
    ## rarefying sample SNF13
    ## rarefying sample SNF2
    ## rarefying sample SNF3
    ## rarefying sample SNF4
    ## rarefying sample SNF5
    ## rarefying sample SNF6
    ## rarefying sample SNF7
    ## rarefying sample SNF8
    ## rarefying sample SNF9

![](index_files/figure-gfm/rarefaction%20curve-2.png)<!-- -->

``` r
p + coord_cartesian(xlim = c(0,10000))
```

![](index_files/figure-gfm/rarefaction%20curve-3.png)<!-- -->

##### Differential abundance

DE calculated using deseq2 and plotted

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:34:22 EDT"

``` r
## repeated from above - no need to run twice
# ps2 <- ps
# taxa_names(ps2) <- paste("ASV", 1:ntaxa(ps2), sep = "")

sample_data(ps2)$treatment <- relevel(factor(sample_data(ps2)$treatment), ref = "mfp")
dds <- phyloseq_to_deseq2(ps2, ~ treatment)
dds <- DESeq2::estimateSizeFactors(dds, type = 'poscounts')
dds <- DESeq(dds, fitType = "local")
results(dds)
```

    ## log2 fold change (MLE): treatment snf vs mfp 
    ## Wald test p-value: treatment snf vs mfp 
    ## DataFrame with 2233 rows and 6 columns
    ##          baseMean log2FoldChange     lfcSE       stat    pvalue      padj
    ##         <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
    ## ASV1     46053.77      0.0857540  0.468595  0.1830025  0.854796  0.926454
    ## ASV2     14161.61      0.2849254  0.642664  0.4433506  0.657512  0.926454
    ## ASV3      5155.23     -0.3574503  1.317759 -0.2712561  0.786194  0.926454
    ## ASV4      4713.39     -0.0263356  1.122141 -0.0234691  0.981276  0.987989
    ## ASV5      3969.29     -0.7953041  1.166567 -0.6817477  0.495399  0.926454
    ## ...           ...            ...       ...        ...       ...       ...
    ## ASV2229 0.0258934       0.199162  0.840082   0.237074  0.812599  0.926454
    ## ASV2230 0.0258934       0.199162  0.840082   0.237074  0.812599  0.926454
    ## ASV2231 0.0258934       0.199162  0.840082   0.237074  0.812599  0.926454
    ## ASV2232 0.0258934       0.199162  0.840082   0.237074  0.812599  0.926454
    ## ASV2233 0.0334575       0.199162  0.895640   0.222368  0.824027  0.926454

``` r
## McMurdie and Holmes (2014) Waste Not, Want Not: Why Rarefying Microbiome Data is Inadmissible. PLoS Computational Biology in press
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:34:29 EDT"

``` r
alpha_cutoff <- 0.1 ## switched to 1 for plots labeled diff-abund_phylum_*_all

resultsNames(dds)
```

    ## [1] "Intercept"            "treatment_snf_vs_mfp"

``` r
res <- results(dds, name = "treatment_snf_vs_mfp", cooksCutoff = FALSE)
sigtab <- res[which(res$pvalue < alpha_cutoff), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps2)[rownames(sigtab), ], "matrix"))
head(sigtab)
```

    ##          baseMean log2FoldChange     lfcSE      stat     pvalue      padj
    ## ASV11  2075.45068       1.331704 0.6979907  1.907911 0.05640268 0.9264542
    ## ASV13  1972.30410      -1.847196 1.0143215 -1.821115 0.06858939 0.9264542
    ## ASV46    81.33976      -3.227817 1.8503352 -1.744450 0.08108069 0.9264542
    ## ASV74   298.18241       1.301787 0.7626082  1.707019 0.08781844 0.9264542
    ## ASV100   11.98080      -3.337108 1.7744995 -1.880591 0.06002755 0.9264542
    ## ASV108   14.40669      -5.916225 2.9366921 -2.014588 0.04394782 0.9264542
    ##         Kingdom                   Phylum                          Class
    ## ASV11  Bacteria               Firmicutes                  Negativicutes
    ## ASV13  Bacteria           Proteobacteria             Betaproteobacteria
    ## ASV46  Bacteria Absconditabacteria (SR1) Absconditabacteria (SR1) [C-1]
    ## ASV74  Bacteria               Firmicutes                        Bacilli
    ## ASV100 Bacteria               Firmicutes                  Negativicutes
    ## ASV108 Bacteria           Proteobacteria            Gammaproteobacteria
    ##                                 Order                         Family
    ## ASV11                  Veillonellales                Veillonellaceae
    ## ASV13                    Neisseriales                  Neisseriaceae
    ## ASV46  Absconditabacteria (SR1) [O-1] Absconditabacteria (SR1) [F-1]
    ## ASV74                 Lactobacillales               Streptococcaceae
    ## ASV100                 Veillonellales                Veillonellaceae
    ## ASV108                 Pasteurellales                Pasteurellaceae
    ##                                 Genus    Species
    ## ASV11                     Veillonella    parvula
    ## ASV13                       Neisseria flavescens
    ## ASV46  Absconditabacteria (SR1) [G-1]       <NA>
    ## ASV74                   Streptococcus       <NA>
    ## ASV100                    Veillonella       <NA>
    ## ASV108                    Haemophilus       <NA>

``` r
write.table(res, "tables/diffabund.txt", quote = FALSE, sep = "\t")
res_taxa <- merge(res, ps2@tax_table, by = 'row.names')
write.table(res_taxa, "tables/diffabund_wtaxa.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:34:29 EDT"

``` r
# Phylum order
x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x)) %>% sort(., TRUE)
sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x)) %>% sort(., TRUE)
sigtab$Genus <- factor(as.character(sigtab$Genus), levels=names(x))
sigtab$shape <- ifelse(sigtab$pvalue<=0.05, "sig", "nonsig")
write.table(sigtab, "tables/sigtable_0.1.txt", quote = FALSE, sep = "\t")

levels(sigtab$Genus) <- sub(" ", "\n", levels(sigtab$Genus))

da <- ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Phylum, shape = shape)) +
  geom_point(size=2) +
  scale_shape_manual(values = c(16,8), guide = "none") +
  theme_bw() +
  geom_vline(xintercept = 0, color = "grey") +
  labs(title = expression("higher in MFP  |  higher in SnF"["2"]*""), x = expression("log"["2"]*" fold change"), y = NULL) +
  labs(x = expression("log"["2"]*" fold change"), y = NULL) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=10), axis.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # geom_hline(yintercept = 0, color = "darkgrey") +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 9, type = "continuous"))
da
```

![](index_files/figure-gfm/diff%20abund%20plots-1.png)<!-- -->

``` r
# ggsave("plots/diff-abund_phylum.jpeg", dpi = 500, width = 6, height = 3, units = "in")
# ggsave("plots/diff-abund_phylum.pdf", width = 6, height = 3, units = "in")
ggsave("plots/diff-abund_phylum.jpeg", dpi = 500)
```

##### Stats and plots for manuscript

``` r
lp <- get_legend(p)
lg <- get_legend(g)
lgs <- get_legend(gs20)
lda <- get_legend(da)

plot_grid(p, g, gs20, da, labels = "AUTO", align = "vh")
```

![](index_files/figure-gfm/for%20manuscript%20taxa/da-1.png)<!-- -->

``` r
ggsave("plots/taxa-da_25.jpeg", width = 11, height = 8, units = "in", dpi = 500)

cp <- plot_grid(p + theme(legend.position="none"), NULL,
          g + theme(legend.position="none"), NULL,
          gs20 + theme(legend.position="none"), NULL,
          da + theme(legend.position="none"), NULL,
          labels = c("A", "", "B", "", "C", "", "D", ""),
          ncol = 4, rel_widths = c(3,2,3,2)) +
  draw_grob(lp, 0.14, 0.5, 0.5, 0.5, hjust = 0, vjust = 0) +
  draw_grob(lg, 0.615, 0.5, 0.5, 0.5, hjust = 0, vjust = 0) +
  draw_grob(lgs, 0.15, 0.05, 0.5, 0.5, hjust = 0, vjust = 0) +
  draw_grob(lda, 0.64, 0.05, 0.5, 0.5, hjust = 0, vjust = 0)
save_plot("plots/taxa-da_aligned_25.jpeg", cp, base_width = 11, base_height = 8)

plot_grid(p + theme(legend.position="none"), NULL,
          g + theme(legend.position="none"), NULL,
          gs20 + theme(legend.position="none"), NULL,
          da + theme(legend.position="none"), NULL,
          labels = c("A", "", "B", "", "C", "", "D", ""),
          ncol = 4, rel_widths = c(3,2,3,2)) +
  draw_grob(lp, 0.14, 0.5, 0.5, 0.5, hjust = 0, vjust = 0) +
  draw_grob(lg, 0.615, 0.5, 0.5, 0.5, hjust = 0, vjust = 0) +
  draw_grob(lgs, 0.15, 0.05, 0.5, 0.5, hjust = 0, vjust = 0) +
  draw_grob(lda, 0.64, 0.05, 0.5, 0.5, hjust = 0, vjust = 0)
```

![](index_files/figure-gfm/for%20manuscript%20taxa/da-2.png)<!-- -->

``` r
ggsave("plots/figure05_taxa-da_aligned_gg_notitle_25.jpeg", width = 11, height = 8, units = "in", dpi = 500, bg = "white")
ggsave("plots/figure05_taxa-da_aligned_gg_notitle_25.pdf", width = 11, height = 8, units = "in", dpi = 500, bg = "white")
ggsave("plots/figure05_taxa-da_aligned_gg_notitle_25.tiff", width = 11, height = 8, units = "in", dpi = 500, bg = "white")
```

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:34:45 EDT"

``` r
kruskal.test(Observed ~ treatment, data = samdf_wstats)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Observed by treatment
    ## Kruskal-Wallis chi-squared = 0.28994, df = 1, p-value = 0.5903

``` r
kruskal.test(Shannon ~ treatment, data = samdf_wstats)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Shannon by treatment
    ## Kruskal-Wallis chi-squared = 0.00065746, df = 1, p-value = 0.9795

``` r
kruskal.test(Simpson ~ treatment, data = samdf_wstats)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Simpson by treatment
    ## Kruskal-Wallis chi-squared = 0.0059172, df = 1, p-value = 0.9387

``` r
ddply(samdf_wstats, ~treatment, summarise, mean = mean(Observed), sd = sd(Observed))
```

    ##   treatment     mean       sd
    ## 1       mfp 220.3077 72.02243
    ## 2       snf 214.5385 78.78940

``` r
ddply(samdf_wstats, ~treatment, summarise, mean = mean(Shannon), sd = sd(Shannon))
```

    ##   treatment     mean        sd
    ## 1       mfp 2.934808 0.5422159
    ## 2       snf 3.001520 0.4248979

``` r
ddply(samdf_wstats, ~treatment, summarise, mean = mean(Simpson), sd = sd(Simpson))
```

    ##   treatment      mean         sd
    ## 1       mfp 0.8444823 0.10438684
    ## 2       snf 0.8645097 0.06352908

``` r
summarize_phyloseq(subset_samples(ps2, treatment=="mfp"))
```

    ## [[1]]
    ## [1] "1] Min. number of reads = 83199"
    ## 
    ## [[2]]
    ## [1] "2] Max. number of reads = 222780"
    ## 
    ## [[3]]
    ## [1] "3] Total number of reads = 1919347"
    ## 
    ## [[4]]
    ## [1] "4] Average number of reads = 147642.076923077"
    ## 
    ## [[5]]
    ## [1] "5] Median number of reads = 143018"
    ## 
    ## [[6]]
    ## [1] "7] Sparsity = 0.901340039271074"
    ## 
    ## [[7]]
    ## [1] "6] Any OTU sum to 1 or less? YES"
    ## 
    ## [[8]]
    ## [1] "8] Number of singletons = 801"
    ## 
    ## [[9]]
    ## [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)0.626959247648903"
    ## 
    ## [[10]]
    ## [1] "10] Number of sample variables are: 9"
    ## 
    ## [[11]]
    ## [1] "sampleid"     "sampletype"   "source"       "toothpaste"   "treatment"   
    ## [6] "study_design" "sex"          "time"         "panelist"

``` r
summarize_phyloseq(subset_samples(ps2, treatment=="snf"))
```

    ## [[1]]
    ## [1] "1] Min. number of reads = 118808"
    ## 
    ## [[2]]
    ## [1] "2] Max. number of reads = 204407"
    ## 
    ## [[3]]
    ## [1] "3] Total number of reads = 2040168"
    ## 
    ## [[4]]
    ## [1] "4] Average number of reads = 156936"
    ## 
    ## [[5]]
    ## [1] "5] Median number of reads = 157989"
    ## 
    ## [[6]]
    ## [1] "7] Sparsity = 0.903923662544352"
    ## 
    ## [[7]]
    ## [1] "6] Any OTU sum to 1 or less? YES"
    ## 
    ## [[8]]
    ## [1] "8] Number of singletons = 875"
    ## 
    ## [[9]]
    ## [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)0.850873264666368"
    ## 
    ## [[10]]
    ## [1] "10] Number of sample variables are: 9"
    ## 
    ## [[11]]
    ## [1] "sampleid"     "sampletype"   "source"       "toothpaste"   "treatment"   
    ## [6] "study_design" "sex"          "time"         "panelist"

``` r
summarize_phyloseq(ps2)
```

    ## [[1]]
    ## [1] "1] Min. number of reads = 83199"
    ## 
    ## [[2]]
    ## [1] "2] Max. number of reads = 222780"
    ## 
    ## [[3]]
    ## [1] "3] Total number of reads = 3959515"
    ## 
    ## [[4]]
    ## [1] "4] Average number of reads = 152289.038461538"
    ## 
    ## [[5]]
    ## [1] "5] Median number of reads = 146668.5"
    ## 
    ## [[6]]
    ## [1] "7] Sparsity = 0.902631850907713"
    ## 
    ## [[7]]
    ## [1] "6] Any OTU sum to 1 or less? YES"
    ## 
    ## [[8]]
    ## [1] "8] Number of singletons = 36"
    ## 
    ## [[9]]
    ## [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)1.47783251231527"
    ## 
    ## [[10]]
    ## [1] "10] Number of sample variables are: 9"
    ## 
    ## [[11]]
    ## [1] "sampleid"     "sampletype"   "source"       "toothpaste"   "treatment"   
    ## [6] "study_design" "sex"          "time"         "panelist"

``` r
plot_grid(alpha, bray, uni, euc, labels = "AUTO")
```

![](index_files/figure-gfm/additional%20for%20manuscript-1.png)<!-- -->

``` r
ggsave("plots/alpha-beta.jpeg", width = 10, height = 10, units = "in", dpi = 500)

plot_grid(alpha, braya, unia, euca, labels = "AUTO")
```

![](index_files/figure-gfm/additional%20for%20manuscript-2.png)<!-- -->

``` r
ggsave("plots/figure06_alpha-beta_arrows.jpeg", width = 10, height = 10, units = "in", dpi = 500)
ggsave("plots/figure06_alpha-beta_arrows.pdf", width = 10, height = 10, units = "in", dpi = 500)
ggsave("plots/figure06_alpha-beta_arrows.tiff", width = 10, height = 10, units = "in", dpi = 500)
```

``` r
# View(ps2@tax_table)
table_s1 <- merge(t(ps2@otu_table), ps2@tax_table, by.x = 'row.names', by.y = 'row.names')
write.table(table_s1, "tables/table_s1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

### Output citations for packages loaded in analysis

``` r
Sys.time()
```

    ## [1] "2022-10-03 18:34:51 EDT"

``` r
citations <- function(includeURL = TRUE, includeRStudio = TRUE) {
  if(includeRStudio == TRUE) {
    ref.rstudio <- rstudioapi::versionInfo()$citation
    if(includeURL == FALSE) {
      ref.rstudio$url <- NULL;
    }
    print(ref.rstudio, style = 'text')
    cat('\n')
  }

  cit.list <- c('base', names(sessionInfo()$otherPkgs))
  for(i in 1:length(cit.list)) {
    ref <- citation(cit.list[i])
    if(includeURL == FALSE) {
      ref$url <- NULL;
    }
    print(ref, style = 'text')
    cat('\n')
  }
}

citations()
```

    ## RStudio Team (2022). _RStudio: Integrated Development Environment for
    ## R_. RStudio, PBC, Boston, MA. <http://www.rstudio.com/>.
    ## 
    ## R Core Team (2022). _R: A Language and Environment for Statistical
    ## Computing_. R Foundation for Statistical Computing, Vienna, Austria.
    ## <https://www.R-project.org/>.

    ## Warning in citation(cit.list[i]): no date field in DESCRIPTION file of package
    ## 'ranacapa'

    ## Kandlikar G (2022). _ranacapa: Utility Functions and 'shiny' App for
    ## Simple Environmental DNA Visualizations and Analyses_. R package
    ## version 0.1.0, <https://github.com/gauravsk/ranacapa>.
    ## 
    ## Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change
    ## and dispersion for RNA-seq data with DESeq2." _Genome Biology_, *15*,
    ## 550. doi:10.1186/s13059-014-0550-8
    ## <https://doi.org/10.1186/s13059-014-0550-8>.
    ## 
    ## Wickham H (2007). "Reshaping Data with the reshape Package." _Journal
    ## of Statistical Software_, *21*(12), 1-20.
    ## <http://www.jstatsoft.org/v21/i12/>.
    ## 
    ## Lahti L, Shetty S (2012-2019). "microbiome R package."
    ## 
    ## Wickham H (2011). "The Split-Apply-Combine Strategy for Data Analysis."
    ## _Journal of Statistical Software_, *40*(1), 1-29.
    ## <https://www.jstatsoft.org/v40/i01/>.
    ## 
    ## Wickham H, Seidel D (2022). _scales: Scale Functions for
    ## Visualization_. R package version 1.2.1,
    ## <https://CRAN.R-project.org/package=scales>.
    ## 
    ## Davis NM, Proctor D, Holmes SP, Relman DA, Callahan BJ (2017). "Simple
    ## statistical identification and removal of contaminant sequences in
    ## marker-gene and metagenomics data." _bioRxiv_, 221499.
    ## doi:10.1101/221499 <https://doi.org/10.1101/221499>.
    ## 
    ## McMurdie PJ, Holmes S (2013). "phyloseq: An R package for reproducible
    ## interactive analysis and graphics of microbiome census data." _PLoS
    ## ONE_, *8*(4), e61217.
    ## <http://dx.plos.org/10.1371/journal.pone.0061217>.
    ## 
    ## Wright ES (2016). "Using DECIPHER v2.0 to Analyze Big Biological
    ## Sequence Data in R." _The R Journal_, *8*(1), 352-359.
    ## 
    ## Müller K, Wickham H, James DA, Falcon S (2022). _RSQLite: SQLite
    ## Interface for R_. R package version 2.2.14,
    ## <https://CRAN.R-project.org/package=RSQLite>.
    ## 
    ## Morgan M, Anders S, Lawrence M, Aboyoun P, Pagès H, Gentleman R (2009).
    ## "ShortRead: a Bioconductor package for input, quality assessment and
    ## exploration of high-throughput sequence data." _Bioinformatics_, *25*,
    ## 2607-2608. doi:10.1093/bioinformatics/btp450
    ## <https://doi.org/10.1093/bioinformatics/btp450>,
    ## <http://dx.doi.org10.1093/bioinformatics/btp450>.
    ## 
    ## Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan
    ## M, Carey V (2013). "Software for Computing and Annotating Genomic
    ## Ranges." _PLoS Computational Biology_, *9*.
    ## doi:10.1371/journal.pcbi.1003118
    ## <https://doi.org/10.1371/journal.pcbi.1003118>,
    ## <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
    ## 
    ## Morgan M, Obenchain V, Hester J, Pagès H (2022). _SummarizedExperiment:
    ## SummarizedExperiment container_. R package version 1.26.1,
    ## <https://bioconductor.org/packages/SummarizedExperiment>.
    ## 
    ## Ahlmann-Eltze C, Hickey P, Pagès H (2022). _MatrixGenerics: S4 Generic
    ## Summary Statistic Functions that Operate on Matrix-Like Objects_. R
    ## package version 1.8.1,
    ## <https://bioconductor.org/packages/MatrixGenerics>.
    ## 
    ## Bengtsson H (2022). _matrixStats: Functions that Apply to Rows and
    ## Columns of Matrices (and to Vectors)_. R package version 0.62.0,
    ## <https://CRAN.R-project.org/package=matrixStats>.
    ## 
    ## Morgan M, Pagès H, Obenchain V, Hayden N (2022). _Rsamtools: Binary
    ## alignment (BAM), FASTA, variant call (BCF), and tabix file import_. R
    ## package version 2.12.0, <https://bioconductor.org/packages/Rsamtools>.
    ## 
    ## Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan
    ## M, Carey V (2013). "Software for Computing and Annotating Genomic
    ## Ranges." _PLoS Computational Biology_, *9*.
    ## doi:10.1371/journal.pcbi.1003118
    ## <https://doi.org/10.1371/journal.pcbi.1003118>,
    ## <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
    ## 
    ## Pagès H, Aboyoun P, Gentleman R, DebRoy S (2022). _Biostrings:
    ## Efficient manipulation of biological strings_. R package version
    ## 2.64.0, <https://bioconductor.org/packages/Biostrings>.
    ## 
    ## Arora S, Morgan M, Carlson M, Pagès H (2022). _GenomeInfoDb: Utilities
    ## for manipulating chromosome names, including modifying them to follow a
    ## particular naming style_. R package version 1.32.2,
    ## <https://bioconductor.org/packages/GenomeInfoDb>.
    ## 
    ## Pagès H, Aboyoun P (2022). _XVector: Foundation of external vector
    ## representation and manipulation in Bioconductor_. R package version
    ## 0.36.0, <https://bioconductor.org/packages/XVector>.
    ## 
    ## Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan
    ## M, Carey V (2013). "Software for Computing and Annotating Genomic
    ## Ranges." _PLoS Computational Biology_, *9*.
    ## doi:10.1371/journal.pcbi.1003118
    ## <https://doi.org/10.1371/journal.pcbi.1003118>,
    ## <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
    ## 
    ## Pagès H, Lawrence M, Aboyoun P (2022). _S4Vectors: Foundation of
    ## vector-like and list-like containers in Bioconductor_. R package
    ## version 0.34.0, <https://bioconductor.org/packages/S4Vectors>.
    ## 
    ## Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022).
    ## _BiocParallel: Bioconductor facilities for parallel evaluation_. R
    ## package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.
    ## 
    ## Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP
    ## (2016). "DADA2: High-resolution sample inference from Illumina amplicon
    ## data." _Nature Methods_, *13*, 581-583. doi:10.1038/nmeth.3869
    ## <https://doi.org/10.1038/nmeth.3869>.
    ## 
    ## Eddelbuettel D, François R (2011). "Rcpp: Seamless R and C++
    ## Integration." _Journal of Statistical Software_, *40*(8), 1-18.
    ## doi:10.18637/jss.v040.i08 <https://doi.org/10.18637/jss.v040.i08>.
    ## 
    ## Eddelbuettel D (2013). _Seamless R and C++ Integration with Rcpp_.
    ## Springer, New York. doi:10.1007/978-1-4614-6868-4
    ## <https://doi.org/10.1007/978-1-4614-6868-4>, ISBN 978-1-4614-6867-7.
    ## 
    ## Eddelbuettel D, Balamuta JJ (2018). "Extending extitR with extitC++: A
    ## Brief Introduction to extitRcpp." _The American Statistician_, *72*(1),
    ## 28-36. doi:10.1080/00031305.2017.1375990
    ## <https://doi.org/10.1080/00031305.2017.1375990>.

    ## Warning in citation(cit.list[i]): no date field in DESCRIPTION file of package
    ## 'EnhancedVolcano'

    ## Blighe K, Rana S, Lewis M (2022). _EnhancedVolcano: Publication-ready
    ## volcano plots with enhanced colouring and labeling_. R package version
    ## 1.13.2, <https://github.com/kevinblighe/EnhancedVolcano>.
    ## 
    ## Slowikowski K (2021). _ggrepel: Automatically Position Non-Overlapping
    ## Text Labels with 'ggplot2'_. R package version 0.9.1,
    ## <https://CRAN.R-project.org/package=ggrepel>.
    ## 
    ## Kolde R (2019). _pheatmap: Pretty Heatmaps_. R package version 1.0.12,
    ## <https://CRAN.R-project.org/package=pheatmap>.
    ## 
    ## Hamilton NE, Ferry M (2018). "ggtern: Ternary Diagrams Using ggplot2."
    ## _Journal of Statistical Software, Code Snippets_, *87*(3), 1-17.
    ## doi:10.18637/jss.v087.c03 <https://doi.org/10.18637/jss.v087.c03>.
    ## 
    ## Wei T, Simko V (2021). _R package 'corrplot': Visualization of a
    ## Correlation Matrix_. (Version 0.92),
    ## <https://github.com/taiyun/corrplot>.
    ## 
    ## Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P,
    ## O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward
    ## M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M, Durand S,
    ## Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M,
    ## Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter
    ## Braak C, Weedon J (2022). _vegan: Community Ecology Package_. R package
    ## version 2.6-2, <https://CRAN.R-project.org/package=vegan>.
    ## 
    ## Sarkar D (2008). _Lattice: Multivariate Data Visualization with R_.
    ## Springer, New York. ISBN 978-0-387-75968-5,
    ## <http://lmdvr.r-forge.r-project.org>.
    ## 
    ## Simpson G (2022). _permute: Functions for Generating Restricted
    ## Permutations of Data_. R package version 0.9-7,
    ## <https://CRAN.R-project.org/package=permute>.
    ## 
    ## Ram K, Wickham H (2018). _wesanderson: A Wes Anderson Palette
    ## Generator_. R package version 0.3.6,
    ## <https://CRAN.R-project.org/package=wesanderson>.
    ## 
    ## Wilke C (2020). _cowplot: Streamlined Plot Theme and Plot Annotations
    ## for 'ggplot2'_. R package version 1.1.1,
    ## <https://CRAN.R-project.org/package=cowplot>.
    ## 
    ## Bhatnagar S (2021). _manhattanly: Interactive Q-Q and Manhattan Plots
    ## Using 'plotly.js'_. R package version 0.3.0,
    ## <https://CRAN.R-project.org/package=manhattanly>.
    ## 
    ## Tarazona S, Garcia-Alcalde F, Dopazo J, Ferrer A, Conesa A (2011).
    ## "Differential expression in RNA-seq: a matter of depth." _Genome
    ## Research_, *21*(12), 4436.
    ## 
    ## Tarazona S, Furio-Tari P, Turra D, Pietro AD, Nueda MJ, Ferrer A,
    ## Conesa A (2015). "Data quality aware analysis of differential
    ## expression in RNA-seq with NOISeq R/Bioc package." _Nucleic Acids
    ## Research_, *43*(21), e140.
    ## 
    ## Bates D, Maechler M, Jagan M (2022). _Matrix: Sparse and Dense Matrix
    ## Classes and Methods_. R package version 1.4-1,
    ## <https://CRAN.R-project.org/package=Matrix>.
    ## 
    ## Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo
    ## HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry
    ## RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Ole's AK, Pag`es H,
    ## Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M (2015).
    ## "Orchestrating high-throughput genomic analysis with Bioconductor."
    ## _Nature Methods_, *12*(2), 115-121.
    ## <http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
    ## 
    ## Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M.,
    ## Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo,
    ## R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M,
    ## MacDonald, J., Obenchain, V., Ole's, K. A, Pag`es, H., Reyes, A.,
    ## Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan, M.
    ## (2015). "Orchestrating high-throughput genomic analysis with
    ## Bioconductor." _Nature Methods_, *12*(2), 115-121.
    ## <http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
    ## 
    ## Stubben C (2016). _trinotateR: Trinotate annotation report summaries_.
    ## R package version 1.0.
    ## 
    ## Dowle M, Srinivasan A (2021). _data.table: Extension of `data.frame`_.
    ## R package version 1.14.2,
    ## <https://CRAN.R-project.org/package=data.table>.
    ## 
    ## Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor
    ## package for differential expression analysis of digital gene expression
    ## data." _Bioinformatics_, *26*(1), 139-140.
    ## doi:10.1093/bioinformatics/btp616
    ## <https://doi.org/10.1093/bioinformatics/btp616>.
    ## 
    ## McCarthy DJ, Chen Y, Smyth GK (2012). "Differential expression analysis
    ## of multifactor RNA-Seq experiments with respect to biological
    ## variation." _Nucleic Acids Research_, *40*(10), 4288-4297.
    ## doi:10.1093/nar/gks042 <https://doi.org/10.1093/nar/gks042>.
    ## 
    ## Chen Y, Lun AAT, Smyth GK (2016). "From reads to genes to pathways:
    ## differential expression analysis of RNA-Seq experiments using Rsubread
    ## and the edgeR quasi-likelihood pipeline." _F1000Research_, *5*, 1438.
    ## doi:10.12688/f1000research.8987.2
    ## <https://doi.org/10.12688/f1000research.8987.2>.
    ## 
    ## Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
    ## "limma powers differential expression analyses for RNA-sequencing and
    ## microarray studies." _Nucleic Acids Research_, *43*(7), e47.
    ## doi:10.1093/nar/gkv007 <https://doi.org/10.1093/nar/gkv007>.
    ## 
    ## Wickham H (2021). _forcats: Tools for Working with Categorical
    ## Variables (Factors)_. R package version 0.5.1,
    ## <https://CRAN.R-project.org/package=forcats>.
    ## 
    ## Wickham H (2019). _stringr: Simple, Consistent Wrappers for Common
    ## String Operations_. R package version 1.4.0,
    ## <https://CRAN.R-project.org/package=stringr>.
    ## 
    ## Wickham H, François R, Henry L, Müller K (2022). _dplyr: A Grammar of
    ## Data Manipulation_. R package version 1.0.9,
    ## <https://CRAN.R-project.org/package=dplyr>.
    ## 
    ## Henry L, Wickham H (2020). _purrr: Functional Programming Tools_. R
    ## package version 0.3.4, <https://CRAN.R-project.org/package=purrr>.
    ## 
    ## Wickham H, Hester J, Bryan J (2022). _readr: Read Rectangular Text
    ## Data_. R package version 2.1.2,
    ## <https://CRAN.R-project.org/package=readr>.
    ## 
    ## Müller K, Wickham H (2022). _tibble: Simple Data Frames_. R package
    ## version 3.1.8, <https://CRAN.R-project.org/package=tibble>.
    ## 
    ## Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_.
    ## Springer-Verlag New York. ISBN 978-3-319-24277-4,
    ## <https://ggplot2.tidyverse.org>.
    ## 
    ## Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
    ## Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E,
    ## Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi
    ## K, Vaughan D, Wilke C, Woo K, Yutani H (2019). "Welcome to the
    ## tidyverse." _Journal of Open Source Software_, *4*(43), 1686.
    ## doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
    ## 
    ## Wickham H, Girlich M (2022). _tidyr: Tidy Messy Data_. R package
    ## version 1.2.0, <https://CRAN.R-project.org/package=tidyr>.

``` r
beepr::beep(sound = "mario")
```

##### END
