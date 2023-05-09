##################################################
#RepeatMasker region annotation
#Bash script
#################################################

repeat_annotation() {
id="$1"

hg38_repeatmasker_bed="RepeatMasker.hg38.bed"
reference_size="/Homo_sapiens_assembly38.chrom.sizes"
reference_fasta="Homo_sapiens_assembly38.fasta"
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

##NFR##
#NFR reads extraction and convert to genomic region in bed format
samtools view -b -e ' ( (tlen <= 118) & (tlen >= -118) )' -f66 ${input_directory}/${id}.bam | bedtools bamtobed -i stdin | cut -f1,2,3,6 > ${output_directory}/${id}.NFR.bed

##extract positive strands and extract 5 prime end positon
cat ${output_directory}/${id}.NFR.bed | awk '$4 == "+" {print $1 "\t" $2 "\t" $2+1}' > ${output_directory}/${id}.NFR.pos_str.bed

##extract negative strands and extract 5 prime end positon
cat ${output_directory}/${id}.short.bed | awk '$4 == "-" {print $1 "\t" $3 "\t" $3+1}' > ${output_directory}/${id}.NFR.neg.bed

##concatenate end positions from both strands
cat ${output_directory}/${id}.NFR.pos_str.bed ${output_directory}/${id}.NFR.neg.bed | sort -k1,1 -k2,2n > ${output_directory}/${id}.NFR.both_str.bed

##Annotate RepeatMasker region
bedtools intersect -a ${output_directory}/${id}.NFR.both_str.bed -b ${hg38_repeatmasker_bed} -wa -wb \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
> ${output_directory}/${id}.NFR.both_str.inRepeat.txt

##Annotate regions not in RepeatMasker
bedtools intersect -v -a ${output_directory}/${id}.NFR.both_str.bed -b ${hg38_repeatmasker_bed} \
> ${output_directory}/${id}.NFR.both_str.outRepeat.txt

##Concatenate two files
cat ${output_directory}/${id}.NFR.both_str.txt ${output_directory}/${id}.NFR.both_str.outRepeat.txt \
> ${output_directory}/${id}.NFR.both_str.RepeatAnno.txt

##Count reads of each annotated repeat region
cat ${output_directory}/${id}.NFR.both_str.RepeatAnno.txt |  awk -v c=8 'BEGIN{FS=OFS="\t"} {for(i=NF+1; i<=c; i++) $i="No_repeat"} 1' | cut -f5,6,8 | sort | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  | sort -rnk4 > ${output_directory}/${id}.NFR.summary.txt


##NBR##
#NFR reads extraction and convert to genomic region in bed format
samtools view -b -e ' ( (tlen <= 118) & (tlen >= -118) )' -f66 ${input_directory}/${id}.bam | bedtools bamtobed -i stdin | cut -f1,2,3,6 > ${output_directory}/${id}.NBR.bed

##extract positive strands and extract 5 prime end positon
cat ${output_directory}/${id}.NBR.bed | awk '$4 == "+" {print $1 "\t" $2 "\t" $2+1}' > ${output_directory}/${id}.NBR.pos_str.bed

##extract negative strands and extract 5 prime end positon
cat ${output_directory}/${id}.short.bed | awk '$4 == "-" {print $1 "\t" $3 "\t" $3+1}' > ${output_directory}/${id}.NBR.neg.bed

##concatenate end positions from both strands
cat ${output_directory}/${id}.NBR.pos_str.bed ${output_directory}/${id}.NBR.neg.bed | sort -k1,1 -k2,2n > ${output_directory}/${id}.NBR.both_str.bed

##Annotate RepeatMasker region
bedtools intersect -a ${output_directory}/${id}.NBR.both_str.bed -b ${hg38_repeatmasker_bed} -wa -wb \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
> ${output_directory}/${id}.NBR.both_str.inRepeat.txt

##Annotate regions not in RepeatMasker
bedtools intersect -v -a ${output_directory}/${id}.NBR.both_str.bed -b ${hg38_repeatmasker_bed} \
> ${output_directory}/${id}.NBR.both_str.outRepeat.txt

##Concatenate two files
cat ${output_directory}/${id}.NBR.both_str.txt ${output_directory}/${id}.NBR.both_str.outRepeat.txt \
> ${output_directory}/${id}.NBR.both_str.RepeatAnno.txt

##Count reads of each annotated repeat region
cat ${output_directory}/${id}.NBR.both_str.RepeatAnno.txt |  awk -v c=8 'BEGIN{FS=OFS="\t"} {for(i=NF+1; i<=c; i++) $i="No_repeat"} 1' | cut -f5,6,8 | sort | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  | sort -rnk4 > ${output_directory}/${id}.NBR.summary.txt

}
export -f repeat_annotation

##path
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

##run gnu parallel
ls ${input_directory}/*.bam | xargs -n 1 basename |
sed 's/.bam//' | parallel --jobs ${jobs} --verbose --resume-failed --joblog ${output_directory}/bamCoverage_parallel.log bamCoverage_parallel
