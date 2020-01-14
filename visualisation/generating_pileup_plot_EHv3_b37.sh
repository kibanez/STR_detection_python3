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
GRAPH_SCRIPT=/genomes/scratch/kgarikano/GEL_STR/Visualisation/GraphAlignmentViewer/GraphAlignmentViewer.py

# Fasta file with .fai index for reference sequence. If not provided, flanks are set to 'N'. Default: None--reference_fasta REFERENCE_FASTA
REFERENCE_FASTA=/genomes/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa

# Path to variant catalog used to run EH v3
VARIANT_CATALOG_V3=/genomes/scratch/kgarikano/GEL_STR/specs/EHv3.1.2/GRCh37/variant_catalog_GRCh37_12oct2019.json

# Input folder where the json/vcf/log or bam files are
INPUT_FOLDER=/genomes/scratch/kgarikano/GEL_STR/Validation_golden_table/output_EHv3.1.2/

# Output folder where we want the plots
#OUTPUT_FOLDER=/genomes/scratch/kgarikano/GEL_STR/Visualisation/paper_positive/EHv3
OUTPUT_FOLDER=/genomes/scratch/kgarikano/GEL_STR/Visualisation/population/AFR

# Parameters passing by console
if [ $# -ne 1 ]
  then
    echo "A list together with platekeys and LOCUS_ID is required, separated by commas"
fi

list_ids=$1

module load python/3.6.5
# Load the virtual environment for dependencies
source /genomes/scratch/kgarikano/GEL_STR/Visualisation/GraphAlignmentViewer/venv/bin/activate


mkdir -p ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}

cat ${list_ids} | while read line; do

   IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
    ID_NAME='EH_'${NAMES[0]}
    LOCUS_ID=${NAMES[1]}

    INPUT_BAM=${INPUT_FOLDER}${ID_NAME}'_realigned.bam'
    INPUT_VCF=${INPUT_FOLDER}${ID_NAME}'.vcf'
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

done