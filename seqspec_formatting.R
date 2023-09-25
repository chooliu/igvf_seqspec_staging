# small R script showing how onlists were obtained =============================


library(tidyverse)
library(tools)

sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# attached base packages:
#   [1] tools     stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4     readr_2.1.0     tidyr_1.1.4     tibble_3.1.6
# [8] ggplot2_3.3.5   tidyverse_1.3.1
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.7       cellranger_1.1.0 pillar_1.6.4     compiler_4.0.2   dbplyr_2.1.1     bit_4.0.4
# [7] jsonlite_1.7.2   lubridate_1.8.0  lifecycle_1.0.1  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.12
# [13] reprex_2.0.1     cli_3.1.0        rstudioapi_0.13  DBI_1.1.1        parallel_4.0.2   haven_2.4.3
# [19] xml2_1.3.2       withr_2.4.2      httr_1.4.2       fs_1.5.0         generics_0.1.1   vctrs_0.3.8
# [25] hms_1.1.1        bit64_4.0.5      grid_4.0.2       tidyselect_1.1.1 glue_1.5.0       R6_2.5.1
# [31] fansi_0.5.0      readxl_1.3.1     vroom_1.5.6      tzdb_0.2.0       modelr_0.1.8     magrittr_2.0.1
# [37] backports_1.3.0  scales_1.1.1     ellipsis_0.3.2   rvest_1.0.2      assertthat_0.2.1 colorspace_2.0-2
# [43] utf8_1.2.2       stringi_1.7.5    munsell_0.5.0    broom_0.7.10     crayon_1.4.2





# i5 & i7 ======================================================================


# valid i5 and i7 barcodes that each identify a 384-well plate
# note that a paired i5 & i7 are reverse complements of one another
# source: supplementary files in Heffel, et al (2022) && Liu, et al (2023) preprints

valid_i5i7 <-
  read_tsv("~/Desktop/valid_i5i7.cpp", col_names = c("name")) %>%
  separate(name, into = c("indexname", "i7", "i5")) %>%
  arrange(indexname)

write_lines(valid_i5i7$i7, "i7_onlist.txt")
write_lines(valid_i5i7$i5, "i5_onlist.txt")



md5sum("i5_onlist.txt")  # 4f15e5aea0ea639b94a3ede217701ad3
md5sum("i7_onlist.txt") # d21cac93bee474f5113020c0d8f71c3b





# cell/well barcodes ===========================================================

# 8bp cell barcodes: or more technically, well-IDs
# these identify the well of the 384-well plate in the assay
# if desired, can be used to check for technical artifacts by position
set1 <- read_lines("A02a_cellbarcodes_subset1.fa")
set2 <- read_lines("A02a_cellbarcodes_subset2.fa")

# format valid cell barcodes into single sequence
# arranged in order A01, A02, ..., P23, P24
# source: https://github.com/chooliu/snmCTseq_Pipeline --> Scripts/
tibble(seq = set1[c(F, T)] %>% gsub("\\^", "", .),
       well = set1[c(T, F)] %>% gsub(">", "", .)) %>%
bind_rows(
  tibble(seq = set2[c(F, T)] %>% gsub("\\^", "", .),
       well = set2[c(T, F)] %>% gsub(">", "", .))) %>%
  mutate(row = str_sub(well, 0, 1),
         col = str_sub(well, 2, 3) %>% str_pad(., width = 2, side = "left", pad = "0")) %>%
  arrange(row, col) %>%
  .$seq %>% unlist %>%
  write_lines("cb_onlist.txt")




md5sum("cb_onlist.txt") # 10a4583b09d8c97bd5b1c6fd30be3aa0




