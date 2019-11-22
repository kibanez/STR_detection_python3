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
REFERENCE_FASTA=/Users/kibanez/Documents/human_reference/b38/genome.fa
# Path to variant catalog used to run EH v3
VARIANT_CATALOG_V3=/Users/kibanez/Documents/STRs/specs/EHv3.1.2/GRCh38/variant_catalog_GRCh38_12oct2019.json
# Input folder where the json/vcf/log or bam files are
INPUT_FOLDER=/Users/kibanez/Documents/STRs/ANALYSIS/random_pileups
# Output folder where we want the plots
OUTPUT_FOLDER=/Users/kibanez/Documents/STRs/ANALYSIS/random_pileups/plots
# List of the platekeys to generate pileups
list_ids=/Users/kibanez/Downloads/borratzeko/list_platekeys_poster.txt


# Load the virtual environment for dependencies
source /Users/kibanez/git/GraphAlignmentViewer/venv/bin/activate

mkdir -p ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}


cat ${list_ids} | while read line; do

    IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line" 
    ID_NAME=${NAMES[0]}'.str_eh_v3'

    INPUT_BAM=${INPUT_FOLDER}${ID_NAME}'_realigned.bam'
	  INPUT_VCF=${INPUT_FOLDER}${ID_NAME}'.vcf'


    python3 ${GRAPH_SCRIPT} \
    --variant_catalog ${VARIANT_CATALOG_V3} \
    --read_align ${INPUT_BAM} \
    --vcf ${INPUT_VCF} \
    --reference_fasta ${REFERENCE_FASTA} \
    --file_format v3

done
