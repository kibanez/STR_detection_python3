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
GRAPH_SCRIPT=/genomes/scratch/kgarikano/GEL_STR/STR_Visualization/GraphAlignmentViewer/GraphAlignmentViewer.py

# Fasta file with .fai index for reference sequence. If not provided, flanks are set to 'N'. Default: None--reference_fasta REFERENCE_FASTA
REFERENCE_FASTA=/genomes/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa

# Path to variant catalog used to run EH v3
VARIANT_CATALOG_V3=/genomes/scratch/kgarikano/GEL_STR/specs/EHv3.1.2/GRCh37/variant_catalog_GRCh37_12oct2019.json

# Input folder where the json/vcf/log or bam files are
INPUT_FOLDER=/home/dpasko/STRs/research_113K/EH_output_v3.1.2_October2019/

# Output folder where we want the plots
OUTPUT_FOLDER=/genomes/scratch/kgarikano/GEL_STR/STR_Visualization/random_plots/

# List of the platekeys to generate pileups
list_ids=/genomes/scratch/kgarikano/GEL_STR/STR_Visualization/list_platekey_EHv3_GRCh37.txt

# Locus id
LOCUS_ID=ATXN1

module load python/3.6.5
# Load the virtual environment for dependencies
source /genomes/scratch/kgarikano/GEL_STR/STR_Visualization/GraphAlignmentViewer/venv/bin/activate


mkdir -p ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}


cat ${list_ids} | while read line; do

   IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
    ID_NAME='EH_'${NAMES[0]}

    INPUT_BAM=${INPUT_FOLDER}${ID_NAME}'_realigned.bam'
    INPUT_VCF=${INPUT_FOLDER}${ID_NAME}'.vcf'


    python3 ${GRAPH_SCRIPT} \
    --variant_catalog ${VARIANT_CATALOG_V3} \
    --read_align ${INPUT_BAM} \
    --output_prefix ${ID_NAME} \
    --output_dir ${OUTPUT_FOLDER} \
    --dpi 600 \
    --title_prefix ${ID_NAME}'_EHv3_pileup' \
    --reference_fasta ${REFERENCE_FASTA} \
    --locus_id ${LOCUS_ID} \
    --file_format v3

done