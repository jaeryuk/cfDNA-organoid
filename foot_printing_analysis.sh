##################################################
#convert bam file to bigwigfile using deeptools
#Bash script
#################################################

bamCoverage_parallel() {
id="$1"

  input_directory="/path/to/input_directory"
  output_directory="/path/to/output_directory"

  blacklist="hg38-blacklist.v3.bed.gz"

bamCoverage --bam ${bam_directory}/${id}.bam \
    -o ${output_directory}/${id}.bw \
    --binSize 10 \
    --normalizeUsing BPM \
    --centerReads \
    --ignoreForNormalization chrX \
    --extendReads \
    --smoothLength 30 \
    --blackListFileName ${blacklist} \
    --minMappingQuality 30
  }
export -f bamCoverage_parallel

##path
input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

##run gnu parallel
ls ${input_directory}/*.bam | xargs -n 1 basename |
sed 's/.bam//' | parallel --jobs ${jobs} --verbose --resume-failed --joblog ${output_directory}/bamCoverage_parallel.log bamCoverage_parallel


##########################################################
#compute and plot BPM-normalized coverage around DNA binding proteins using deepTools
#Bash script
#########################################################

input_directory="/path/to/input_directory"
output_directory="/path/to/output_directory"

refGene="hg38.refGene.gtf.gz"
super_enhancer="/hg38_super_enhancer.liftover.bed"
TFBS="encRegTfbsClustered.hg38.bed"
blacklist="hg38-blacklist.v3.bed.gz"

for bw in ${input_directory}/*.bw
do


##Transcription units
  computeMatrix scale-regions \
      -b 2000 -a 2000 \
      --regionBodyLength 2000 \
      -R ${refGene} \
      -S ${input_directory}/*.bw \
      --skipZeros \
      -o ${input_directory}/TSS_TES.gz \
      --outFileSortedRegions ${input_directory}/TSS_TES.bed.gz \
      -bl ${blacklist}

  plotProfile -m ${input_directory}/TSS_TES.gz \
       -out ${input_directory}/TSS_TES.png \
       --perGroup \
       --colors red green blue \
       --legendLocation upper-right \
       --yAxisLabel "BPM-normalized coverage" \
       --regionsLabel " "

##Transcription start sites
  computeMatrix reference-point \
        --referencePoint TSS \
        -b 1000 -a 1000 \
        -R ${refGene} \
        -S ${input_directory}/*.bw \
        --skipZeros \
        -o ${input_directory}/TSS.gz \
        --outFileSortedRegions ${input_directory}/TSS.bed.gz \
        -bl ${blacklist}

  plotProfile -m ${input_directory}/TSS.gz \
         -out ${input_directory}/TSS.png \
         --perGroup \
         --colors red green blue \
         --legendLocation upper-right \
         --yAxisLabel "BPM-normalized coverage" \
         --regionsLabel " "

##Transcription factor binding sites
  computeMatrix reference-point \
        --referencePoint center \
        -b 1000 -a 1000 \
        -R ${TFBS} \
        -S ${input_directory}/*.bw \
        --skipZeros \
        -o ${input_directory}/TFBS.gz \
        --outFileSortedRegions ${input_directory}/TFBS.bed.gz \
        -bl ${blacklist}

  plotProfile -m ${input_directory}/TFBS.gz \
         -out ${input_directory}/TFBS.png \
         --perGroup \
         --colors red green blue \
         --legendLocation upper-right \
         --yAxisLabel "BPM-normalized coverage" \
         --regionsLabel " " \
         --refPointLabel "Center of TFBS"

##super enhancer regions
  computeMatrix scale-regions \
         -b 2000 -a 2000 \
         --regionBodyLength 2000 \
         -R ${super_enhancer} \
         -S ${input_directory}/*.bw \
         --skipZeros \
         -o ${input_directory}/super_enhancer.gz \
         --outFileSortedRegions ${input_directory}/super_enhancer.bed.gz \
         -bl ${blacklist}

  plotProfile -m ${input_directory}/super_enhancer.gz \
          -out ${input_directory}/super_enhancer.png \
          --perGroup \
          --colors red green blue \
          --legendLocation upper-right \
          --yAxisLabel "BPM-normalized coverage" \
          --startLabel "Start" \
          --endLabel "End" \
          --regionsLabel " "
done
