
if [[ $# -eq 0 ]] ; then
    echo 'Two arguments are required: 1) the VCF file to be annotated and 2) the output directory to save the corresponding annotated csv file'
    exit 1
fi


# Input parameters
INPUT_VCF=$1
OUTPUT_DIR=$2


OUTPUT_FILE=${OUTPUT_DIR}/"${INPUT_VCF##*/}"_annovar.csv

cd /home/kgarikano/GEL_STR/annovar/

perl /home/kgarikano/GEL_STR/annovar/table_annovar.pl \
-vcfinput ${INPUT_VCF} \
humandb/ \
-buildver hg19 \
-out ${OUTPUT_FILE} \
-remove \
-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
-operation gx,r,f,f,f \
-nastring . \
-polish \