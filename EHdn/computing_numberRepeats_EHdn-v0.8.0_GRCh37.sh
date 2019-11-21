# Calculate the number of `repeatNumber` estimated by EHdn-v0.8.0

# Output folders for EHdn GRCh37
OUTPUT_PATH_b37=/genomes/scratch/kgarikano/GEL_STR/EHdn/output_EHdnv0.8.0/numberRepeats/GRCh37

# Input folders for EHdn v0.8.0 for GRCh37
INPUT_PATH_b37=/genomes/scratch/kgarikano/GEL_STR/EHdn/output_EHdnv0.8.0/GRCh37


if [ $# -lt 1 ]
  then
    echo "A file containing a list of json files is required (<BAM_NAME>,<BAM_GENDER>,<BAM_PATH>\n)"
    exit
fi

# mkdir output directories
mkdir -p ${OUTPUT_PATH_b37}

input_list=$1

cat $input_list | while read line; do

    IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
    JSON_FILE=${NAMES[0]}

    #?| Let's take the ID, being the last character of the string JSON_FILE
    PLATEKEY=$(echo $JSON_FILE | cut -d '/' -f 9-)

    echo ${PLATEKEY}
    OUTPUT_WGS=${OUTPUT_PATH_b37}${PLATEKEY}'number.repeatUnit.txt'

    grep RepeatUnit ${JSON_FILE} > ${OUTPUT_WGS}

done
