library(tidyverse)

valid <- read_tsv("~/Desktop/valid_i5i7.cpp", col_names = c("name")) %>%
  separate(name, into = c("indexname", "i7", "i5")) %>%
  arrange(indexname) %>%
  transmute(i7 = paste0(i7, " # ", indexname),
            i5 = paste0(i5, " # ", indexname))


write_lines(valid$i7, "i7_onlist.txt")
write_lines(valid$i5, "i5_onlist.txt")

system("head i7_onlist.txt")


tools::md5sum("i5_onlist.txt")  # a0aa9b1754c4b7950816dc212a081d14
tools::md5sum("i7_onlist.txt") # 7bf9d57e65e6f74cfc6f4b10f4b679fc


set1 <- read_lines("A02a_cellbarcodes_subset1.fa")
set2 <- read_lines("A02a_cellbarcodes_subset2.fa")


tibble(seq = set1[c(F, T)] %>% gsub("\\^", "", .),
       well = set1[c(T, F)] %>% gsub(">", "", .)) %>%
bind_rows(
  tibble(seq = set2[c(F, T)] %>% gsub("\\^", "", .),
       well = set2[c(T, F)] %>% gsub(">", "", .))) %>%
  mutate(row = str_sub(well, 0, 1),
         col = str_sub(well, 2, 3) %>% str_pad(., width = 2, side = "left", pad = "0")) %>%
  arrange(row, col) %>%
  transmute(out = paste0(seq, " # ", row, col)) %>%
  write_lines("cb_onlist.txt")


sl("NNNNNNNNNNNNNNNNN")

tools::md5sum("cb_onlist.txt")  # df0aa7d170263a7435cd1074ac69b89a


seqspec init -n "snmCT-seq" -m 2 -o "snmCTseq.yaml" \
"(((barcode:8,linker:9,methyl:134)R1.fastq.gz,
(linker:15,methyl:136)R2.fastq.gz,
(index7:10)I1.fastq.gz,
(index5:10)I2.fastq.gz)methyl,
((barcode:8,linker:9,methyl:134)R1.fastq.gz,(linker:15,methyl:136)R2.fastq.gz,(index7:10)I1.fastq.gz,(index5:10)I2.fastq.gz)rna)"


# seqspec init -n "snmCT-seq" -m 2 -o "snmCTseq.yaml" \
# "(((illumina_p5:29)illumina_p5,
# (truseq_read1:33)truseq_read1,
# (barcode:10)I1.fastq.gz,
# (barcode:8,linker:9,methyl:134)R1.fastq.gz,
# (methyl:141,linker:10)R2.fastq.gz,
# (truseq_read2:27)truseq_read2_trunc,
# (illumina_index7:10)I2.fastq.gz,
# (illumina_p7:24)illumina_p7)methyl,
# ((illumina_p5:29)illumina_p5,
# (truseq_read1:33)truseq_read1,
# (barcode:10)I1.fastq.gz,
# (barcode:8,linker:9,cdna:134)R1.fastq.gz,
# (cdna:141,linker:10)R2.fastq.gz,
# (truseq_read2:27)truseq_read2_trunc,
# (illumina_index7:10)I2.fastq.gz,
# (illumina_p7:24)illumina_p7)rna)"


seqspec print multiome.yaml 


str_length("AATGATACGGCGACCACCGAGATCTACAC")
str_length("ACACTCTTTCCCTACACGACGCTCT")
str_length("ATCTCGTATGCCGTCTTCTGCTTG") # 24 P7

seqspec init -n "snmCT-seq" -m 2 -o "snmCTseq.yaml" \
"((illumina_p5:29,
truseq_read1_trunc:25,
(i5_barcode:10)I2.fastq.gz,
(barcode:8,linker:9,methyl:134)R1.fastq.gz,
(methyl:141,linker:10)R2.fastq.gz,
truseq_read2_trunc:25,
(i7_barcode:10)I1.fastq.gz,
illumina_p7:24)methyl,

(illumina_p5:29,
truseq_read1:33,
(i5_barcode:10)I2.fastq.gz,
(barcode:8,linker:9,cdna:134)R1.fastq.gz,
(cdna:141,linker:10)R2.fastq.gz,
truseq_read2:27,
(i7_barcode:10)I1.fastq.gz,
illumina_p7:24)rna)"



str_length("GATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ATCTCGTATGCCGTCTTCTGCTTG p7



AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNN


!Assay
seqspec_version: 0.0.0
assay: 'snmCT-seq'
sequencer: 'Illumina NextSeq500'
name: snmCT-seq
doi: 'https://doi.org/10.1016/j.xgen.2022.100107'
publication_date: '09 March 2022'
description: 'Single-nucleus methylcytosine and transcriptome sequencing (snmCT-seq). Specification assumes 2x151bp paired-end reads with 10bp dual indices based on snmC-seq3 and SMART-seq2. Bisulfite-converted gDNA reads (methyl) and cDNA reads (RNA) are separated in silico.'
modalities:
  - methyl
  - rna