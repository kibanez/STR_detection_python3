# Bash script to generate pileup plots - important to see the reads that EH takes into account when calling the expansion
# Note this requires running python3 (from January'2020 python2.7 will stop existing)
# https://github.com/Illumina/GraphAlignmentViewer

# Requirements for runnning this script:
# Python3
# Matplotlib -- pip3.7 install matplotlib
# Pysam -- pip3.7 install pysam
# Numpy -- pip3.7 install numpy
# PyYAML (if visualizing output from versions older than v3) -- pip3.7 install pyyaml

# Path to the binary python file
GRAPH_SCRIPT=/Users/kibanez/git/GraphAlignmentViewer/GraphAlignmentViewer.py

# Fasta file with .fai index for reference sequence. If not provided, flanks are set to 'N'. Default: None--reference_fasta REFERENCE_FASTA
#REFERENCE_FASTA=/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa
REFERENCE_FASTA=/Volumes/KIKU_STRs/data/reference_GEL_GRCh38/genome.fa

# Path to variant catalog used to run EH v3
VARIANT_CATALOG_V3=/Users/kibanez/Documents/STRs/specs/EHv3.2.2/batch_march2020/GRCh38/variant_catalog_GRCh38_2march2020.json

# Input folder where the json/vcf/log or bam files are
INPUT_FOLDER=/Volumes/KIKU_STRs/data/gnomAD/EHv322/data/

# Output folder where we want the plots
OUTPUT_FOLDER=/Users/kibanez/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/pileup_gnomAD/figures

# Parameters passing by console
if [ $# -ne 1 ]
  then
    echo "A list with platekeys together with LOCUS_ID is required (separated by commas)"
fi

FILE_INPUT=$1

module load python/3.6.5
# Load the virtual environment for dependencies
source /genomes/scratch/kgarikano/GEL_STR/Visualisation/GraphAlignmentViewer/venv/bin/activate


mkdir -p ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}


while read line; do

   IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
    ID_NAME='EH_'${NAMES[0]}
    LOCUS_ID=${NAMES[1]}

    #INPUT_BAM=${INPUT_FOLDER}${ID_NAME}'_realigned.bam'
    #INPUT_VCF=${INPUT_FOLDER}${ID_NAME}'.vcf'
    INPUT_BAM=${INPUT_FOLDER}${ID_NAME}'.expansion_hunter3_realigned.bam'
    INPUT_VCF=${INPUT_FOLDER}${ID_NAME}'.expansion_hunter3.vcf'

    OUTPUT_FILE_NAME=${ID_NAME}'_'${LOCUS_ID}

    python3 ${GRAPH_SCRIPT} \
    --variant_catalog ${VARIANT_CATALOG_V3} \
    --read_align ${INPUT_BAM} \
    --output_prefix ${OUTPUT_FILE_NAME} \
    --output_dir ${OUTPUT_FOLDER} \
    --dpi 600 \
    --title_prefix ${ID_NAME}'_EHv3_pileup' \
    --reference_fasta ${REFERENCE_FASTA} \
    --locus_id ${LOCUS_ID} \
    --file_format v3

                                                                                                                                                                             68,0-1        97%
done <${FILE_INPUT}