##################################################
#Inser size and read legth extraction from BAM file
#Bash script
#################################################

CollectInsertSizeMetrics() {
id="$1"

#path
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

echo ${id}

#insert size collection
samtools view -h -b -e 'rname != "chrM"' ${input_directory}/${id}.bam \
| picard CollectInsertSizeMetrics \
     -I /dev/stdin \
     -O ${output_directory}/${id}.insert_size_metrics.txt \
     -H ${output_directory}/${id}.insert_size_histogram.pdf \
     --MINIMUM_PCT 0.5 \
     -W 1000

#insert size file processing
 cat ${output_directory}/${id}.insert_size_metrics.txt \
 | tail -n +11 \
 > ${output_directory}/${id}.insert_size_hist.txt

#size < 150 read legth collect
samtools view -f 66 -e 'rname != "chrM"' ${input_directory}/${id}.bam \
| awk '{if (length($10) <= 150) print length($10)}' \
| sort -n | uniq -c | awk '{print $2,"\t",$1}' \
> ${output_directory}/${id}.size.below150.txt

}
export -f CollectInsertSizeMetrics

#path
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

jobs="100%"

#run gnu parallel
ls ${input_fastq_directory}/*.fastq.gz | xargs -n 1 basename |
sed 's/.fastq.gz//' | parallel --jobs ${jobs} --verbose --resume-failed --joblog ${output_directory}/CollectInsertSizeMetrics.log CollectInsertSizeMetrics


##################################################
#Insert size and read legth data processing
#python3 script
#################################################

import pandas as pd
from pathlib import Path

##Insert size data processing
insert_size = '/path/to/insert_size'
source_files = sorted(Path(insert_size).glob('*.insert_size_hist.txt'))

dataframes = []
for file in source_files:
    df = pd.read_csv(file, sep="\t")
    sample_ID = file.name.replace(".insert_size_hist.txt","")
    df['Sample_ID'] = sample_ID
    dataframes.append(df)

insert_size1 = pd.concat(dataframes)

##Read length data processing
read_length = '/path/to/read_length'
source_files = sorted(Path(read_length).glob('*.size.below150.txt'))

dataframes = []
for file in source_files:
    colnames=['insert_size', 'All_Reads.fr_count']
    df = pd.read_csv(file, sep="\t", names=colnames, header=None)
    df = df[["insert_size", "All_Reads.fr_count"]]
    sample_ID = file.name.replace(".size.below150.txt","")
    df['Sample_ID'] = sample_ID
    dataframes.append(df)

insert_size2 = pd.concat(dataframes)

##Concatenate insert size data and read length data
insert_size = pd.concat([insert_size1, insert_size2], axis=0)

##save
insert_size.to_csv("insert_size.read_length.below150.txt", sep="\t")
