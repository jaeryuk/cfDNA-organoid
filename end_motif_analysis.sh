##################################################
#Extract end sequences of cfDNA from BAM files
#Bash script
#################################################

end_motif() {
id="$1"

##path
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

##NFR##
##Extract 5 prime 20 bp sequence from 1st paired-reads
samtools view -f 67 -e '(tlen <=118) & (tlen >= -118)' ${input_directory}/${id}.bam \
| cut -f10 | grep -o ^.................... \
> ${output_directory}/${id}.NFR.1st.fasta

##Extract 3 prime 20 bp fasta sequence from 2nd paired-reads
samtools view -f 131 -e '(tlen <=118) & (tlen >= -118)' ${input_directory}/${id}.bam \
| cut -f10 | grep -o ....................$ \
> ${output_directory}/${id}.NFR.2nd.fasta

##Merge end sequences from 1st and 2nd reads
paste ${output_directory}/${id}.NFR.1st.fasta ${output_directory}/${id}.NFR.2nd.fasta \
| sed 's/\t//g' | awk 'length($1)==40' > ${output_directory}/${id}.NFR.fasta

##NBR##
##Extract 5 prime 20 bp sequence from 1st paired-reads
samtools view -f 67 -e '(tlen >=119) | (tlen <= -119)' ${input_directory}/${id}.bam \
| cut -f10 | grep -o ^.................... \
> ${output_directory}/${id}.NBR.1st.fasta

##Extract 3 prime 20 bp fasta sequence from 2nd paired-reads
samtools view -f 131 -e '(tlen >=119) | (tlen <= -119)' ${input_directory}/${id}.bam \
| cut -f10 | grep -o ....................$ \
> ${output_directory}/${id}.NBR.2nd.fasta

##Merge end sequences from 1st and 2nd reads.
paste ${output_directory}/${id}.NBR.1st.fasta ${output_directory}/${id}.NBR.2nd.fasta \
| sed 's/\t//g' | awk 'length($1)==40' > ${output_directory}/${id}.NBR.fasta

##Randomly extract 5,000,000 reads
shuf -n 5000000 ${output_directory}/${id}.NFR.fasta > ${output_directory}/${id}.NFR.5M.fasta
shuf -n 5000000 ${output_directory}/${id}.NBR.fasta >${output_directory}/${id}.NBR.5M.fasta
}
export -f end_motif

##path
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

##run gnu parallel
ls ${input_directory}/*.bam | xargs -n 1 basename |
sed 's/.bam//' | parallel --jobs ${jobs} --verbose --resume-failed --joblog ${output_directory}/end_motif.log end_motif

##################################################
#Plot end motifs wih DNA sequence logo
#R script
#################################################

library(ggseqlogo)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

#Input path
input_directory <- "/path/to/input_directory/"

#Load end sequence file
junction_motif_df <-
fread(paste0(input_directory,"ID.short.5M.fasta"), header=FALSE)

#Adjust parameter
junction_motif <- junction_motif_df$V1
number <- paste0("n = ", length(junction_motif))

grob <- grobTree(textGrob(number, x=0.04,  y=0.95, hjust=0,
  gp=gpar(col="black", fontsize=20)))

##Plot sequence logo
p1<-
junction_motif %>%
ggseqlogo(method = 'bits' , col_scheme='nucleotide', seq_type='dna')

p1 +
scale_x_continuous(breaks= seq(1,40,by=1),
                   labels=c("-20", "-19", "-18", "-17", "-16", "-15", "-14", "-13", "-12", "-11",
                           "-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1",
                           "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
     axis.text.y = element_text(size=20),
     axis.title.y = element_text(size=30)) +
ylim(c(0,0.10)) +
annotation_custom(grob)
