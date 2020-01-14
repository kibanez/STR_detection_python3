# Bash script to generate pileup plots - important to see the reads that EH takes into account when calling the expansion
# Note this is a script that runs the always used pileup-generator. I call it VINTAGE, since this has been replaced by https://github.com/Illumina/GraphAlignmentViewer
# BUT I have experienced some bugs when running this on various EH-v2.5.5 plots

GRAPH_SCRIPT=/genomes/scratch/kgarikano/GEL_STR/Visualisation/scripts/STR_pileup.py
INPUT_FOLDER=/home/dpasko/STRs/research_80K/EH_output_v2.5.5_July2019/
OUTPUT_FOLDER=/genomes/scratch/kgarikano/GEL_STR/Visualisation/random_plots/

# Parameters passing by console
if [ $# -ne 1 ]
  then
    echo "A list containing a PLATEKEY together with the LOCUS_ID is required, separated by commas"
fi

# Change to output folder to create pileup plots there
cd ${OUTPUT_FOLDER}

FILE_INPUT=$1

while read line; do

   IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
   PLATEKEY='EH_'${NAMES[0]}
   LOCUS_ID=${NAMES[1]}

   JSON_FILE=${INPUT_FOLDER}${PLATEKEY}'.json'
   LOG_FILE=${INPUT_FOLDER}${PLATEKEY}'_alignments_relevant_reads.log'


   echo ${GRAPH_SCRIPT} \
   --json ${JSON_FILE}  \
   --read_align ${LOG_FILE} \
   --color \
   --repeat_id ${LOCUS_ID}


   ${GRAPH_SCRIPT} \
   --json ${JSON_FILE}  \
   --read_align ${LOG_FILE} \
   --color \
   --repeat_id ${LOCUS_ID}

done <${FILE_INPUT}
