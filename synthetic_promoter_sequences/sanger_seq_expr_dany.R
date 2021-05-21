library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1")

setwd("~/Documents/Biozentrum/People/Tobias")
setwd("/Volumes/BC2_VNG/julou/Documents/Biozentrum/People/Tobias/")
myseqs <- read.table("data_matrix_sanger.txt", header=TRUE) %>%
  mutate(seqID=substring(as.character(seqID), 2)) %>%
  separate(seqID, c("plate", "well", "lib"), sep="_") %>%
  mutate(plate=substring(plate, 2))

myexprs <- read.table("med_var_evolved_ex_varcorr")
# fix p16ex problem
spl <- strsplit(as.character(myexprs[[1]]), "_")
for (ispl in which(sapply(spl, length) == 3))
  spl[[ispl]] <- spl[[ispl]][c(1, 3)]
myexprs[[1]] <- sapply(spl, function(.x) paste(.x, collapse="_"))

myexprs <- setNames(myexprs, 
                    c("cloneID", "mean_log_original", "var_log_original", "mean_log", "var_log", "rho_log", "mean_abs_orig", 
                      "var_abs_orig", "mean_abs", "var_abs", "rho_abs", "number_of_datapoints", "var_log_corrected")) %>%
  separate(cloneID, c("plate", "well"), sep="_") %>%
  extract(plate, "plate", "p(\\d+)ex")

selections <- c("35"="medium", "36"="high", "37"="medium", 
                "38"="high", "39"="medium", "40"="high",
                "65"="medium", "66"="high", "67"="medium",
                "68"="high", "69"="medium", "70"="high")
rounds <- c("35"=3, "36"=3, "37"=3, 
            "38"=3, "39"=3, "40"=3,
            "65"=5, "66"=5, "67"=5,
            "68"=5, "69"=5, "70"=5)
mydata <- full_join(myseqs, myexprs, by=c("plate", "well")) %>%
  mutate(plate=factor(plate), well=factor(well), lib=factor(lib)) %>%
  mutate(selection=factor(selections[lib]), round=factor(rounds[lib]))

qplot(mean_log, E_f_max_all, col=selection, data=mydata %>% mutate(round=paste("round", round)) %>% na.omit, facets=round~.) +
  geom_point(col="black", size=8, shape=43, data=filter(mydata, plate==15, well=="C03") %>% mutate(round=paste("round", round)))
# qplot(mean_log, E_f_max_all, col=selection, data=mydata %>% mutate(round=paste("round", round)) %>% filter(plate!=16) %>% na.omit, facets=round~.)
# qplot(mean_log, E_f_max_all, col=selection, data=mydata %>% mutate(lib=paste("lib", lib)) %>% na.omit, facets=~lib)

filter(mydata, plate==15, well=="C03")
mean( filter(mydata, round==5, selection=="medium")$mean_log)
mean( filter(mydata, round==5, selection=="medium")$E_f_max_all)


tmp <- group_by(mydata, plate) %>%
  do(layout=(function(.df) {
    .l <- arrange(.df, plate, well) %>%
      .$lib %>% as.character %>% as.numeric %>%
      matrix(8, 12, byrow=TRUE)
#     browser()
    print(.l)
    return(.l)
  })(.) )


# WITH EXTENDED SEQS ####

mydata_ext <- read.table("data_matrix_sanger_extseqs.txt", header=TRUE, sep="\t") %>%
  extract(cloneId, c("plate", "well"), "p(\\d+).*_(.*)") %>%
  left_join(select(mydata, plate, well, lib)) %>%
  extract(round, c("selection", "round"), "([a-z]+)(\\d)") %>%
  mutate(plate=factor(plate), well=factor(well), selection=factor(selection), round=factor(round))

mydata_ext %>% 
  filter(lib %in% c(65, 66), as.numeric(as.character(plate))<=3, pos_10_box!="down") %>% 
  mutate(round=paste("round", round)) %>% na.omit %>% 
  ggplot(aes(mean_log, E_f_max_all, col=selection)) +
  facet_grid(round~.) +
  geom_point(col="black", size=8, shape=43, data=filter(mydata, plate==15, well=="C03") %>% mutate(round=paste("round", round))) +
  geom_vline(aes(col="med"), xintercept=mean( filter(mydata_ext, pos_10_box!="down", round==5, selection=="med")$mean_log)) +
  geom_vline(aes(col="hi"), xintercept=mean( filter(mydata_ext, pos_10_box!="down", round==5, selection=="hi", mean_log>8.5)$mean_log)) +
  geom_text(aes(label=paste(plate, well, sep='_')))
# qplot(mean_log, E_f_max_all, col=selection, data=mydata %>% mutate(round=paste("round", round)) %>% filter(plate!=16) %>% na.omit, facets=round~.)
# qplot(mean_log, E_f_max_all, col=selection, data=mydata %>% mutate(lib=paste("lib", lib)) %>% na.omit, facets=~lib)

filter(mydata_ext, plate==15, well=="C03")
mean( filter(mydata_ext, pos_10_box!="down", round==5, selection=="med")$mean_log)
mean( filter(mydata_ext, pos_10_box!="down", round==5, selection=="med")$E_f_max_all)



