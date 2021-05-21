###The Function is this one:
library(tidyverse)

read_mot_output <- function(.file){
  #two first lines raw output:
  #                             V1      V2                  V3               V4                  V5
  #  1                       173-201       +            0.305227 Sigma24_spacer17   ecoli_sodc_22_F_-
  #  2 ATGAAACGTTTTAGTCTGGCTATTCTGGC 4.58445   ecoli_sodc_22_F_-  
  .file <- .file[file.info(.file)$size != 0]
  .f <- lapply(.file, function(r){
    .sites_raw <- read.table(r, fill = T, header= F)
    .even_indexes <- seq(2,nrow(.sites_raw), 2)
    .odd_indexes <- seq(1,nrow(.sites_raw)-1, 2)
    ind <- length(.sites_raw[.even_indexes, ]) + 3
    #.odd <- .sites_raw[.odd_indexes,][1:5]
    #.odd <- .odd %>% select(V1, V2, V3, V4, V5) 
    .sites_motevo <- cbind(.sites_raw[.odd_indexes,][1:5], .sites_raw[.even_indexes, ])[,1:8]
    
    colnames(.sites_motevo) <- c('pos', 'strand', 'post_prob', 'tf', 'name', 'seq', 'WM_score', 'seq_name')
    .sites_motevo$WM_score <- as.numeric(as.character(.sites_motevo$WM_score))
    .sites_motevo$post_prob <- as.numeric(as.character(.sites_motevo$post_prob))
    .sites_motevo$pos <- as.character(.sites_motevo$pos)
    .sites_motevo$seq_name <- gsub('_Alon_22_F_', '', .sites_motevo$seq_name)
    .sites_motevo$seq_name <- gsub('_22_F_', '', .sites_motevo$seq_name)
    .sites_motevo$seq_name <- gsub('ecoli_', '', .sites_motevo$seq_name)
    .sites_motevo$seq <- as.character(.sites_motevo$seq)
    .sites_motevo$name <- as.character(.sites_motevo$name)
    .sites_motevo$tf <- as.character(.sites_motevo$tf)
    .sites_motevo$strand <- as.character(.sites_motevo$strand)
    
    return(.sites_motevo)
  })
  return(do.call(rbind, .f))
}

fls <- list.files(path="/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/constitexpr_noisefloor_20200915_motevo_001",pattern = '^sites', ignore.case = T, full.names = T)
mot <- read_mot_output(fls)
mot <- as.tibble(mot)

# Tidying motevo results and filtering for sigma70 factor only and relevant sequences
#pos   strand post_prob tf               name  seq                    WM_score seq_name
prom_list <- c("hi1","hi2","hi3","med1","med2","med3","rpmB","rpsB","rrnB","rplN")
motSigma <- mot %>% 
  filter(seq_name %in% prom_list) %>% 
  filter(grepl("Sigma70",tf)) %>%
  filter(strand=="+") %>% 
  separate(pos,into=c("beg","end"),convert=TRUE,remove=FALSE) %>% 
  extract(tf, "spacer", "Sigma70_spacer([0-9]{2})") %>%
  mutate(spacer=as.double(spacer)) %>% 
  select(-c(name,strand,post_prob)) %>% 
  select(c("seq_name","pos","beg","end","WM_score","spacer","seq"))

motSigma

# Importing data from Dorde jupyternotebook and tidy data to compare with motevo results
pathToData <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/constitexpr_noisefloor_20201013_dorde_001/20201016_results.csv"
dordeResults <- readr::read_delim(pathToData,",")
dordeResults <- dordeResults %>% 
  rename(id=X1) %>% 
  separate(seq_name,sep="_",c("specie","seq_name")) %>%
  mutate(beg=str_length(sequence)-(p_all+spacer_all+12)+1,
         end=str_length(sequence)-p_all+1,
         pos=paste(as.character(beg),as.character(end),sep="-"),
         seq=substr(sequence,beg,end)) %>% 
  rename(WM_score=E_all,spacer=spacer_all) %>% 
  select(-c(id,E_tot,E_max,p_max,spacer_max,specie,sequence,energy_entropy,p_all)) %>% 
  select(c("seq_name","pos","beg","end","WM_score","spacer","seq"))
  
  
