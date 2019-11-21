# Calculate the number of `repeatNumber` estimated by EHdn-v0.8.0

# Output folders for EHdn GRCh38
OUTPUT_PATH_b38=/genomes/scratch/kgarikano/GEL_STR/EHdn/output_EHdnv0.8.0/numberRepeats/GRCh38

# Input folders for EHdn v0.8.0 for GRCh38
INPUT_PATH_b38=/genomes/scratch/kgarikano/GEL_STR/EHdn/output_EHdnv0.8.0/GRCh38


if [ $# -lt 1 ]
  then
    echo "A file containing a list of json files is required (<BAM_NAME>,<BAM_GENDER>,<BAM_PATH>\n)"
    exit
fi

# mkdir output directories
mkdir -p ${OUTPUT_PATH_b38}

input_list=$1

cat $input_list | while read line; do

    IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
    JSON_FILE=${NAMES[0]}

    #?| Let's take the ID, being the last character of the string JSON_FILE
    PLATEKEY=$(echo $JSON_FILE | cut -d '/' -f 9-)

    echo ${PLATEKEY}
    OUTPUT_WGS=${OUTPUT_PATH_b38}${PLATEKEY}'number.repeatUnit.txt'

    grep RepeatUnit ${JSON_FILE} > ${OUTPUT_WGS}

done
