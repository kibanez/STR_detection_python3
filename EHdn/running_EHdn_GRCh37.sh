# Bash script to run ExpansionHunter denovo (EHdn) across all genomes aligned with GRCh37 human genome assembly

eh_binary='/genomes/scratch/kgarikano/GEL_STR/EHdn/sw/ExpansionHunterDenovo-v0.8.6-linux_x86_64/bin/ExpansionHunterDenovo-v0.8.6'
reference_fasta='/genomes/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
output_folder='/genomes/scratch/kgarikano/GEL_STR/EHdn/output_EHdnv0.8.0/GRCh37'

# Create dir if does not exist
mkdir -p $output_folder

export LC_ALL=C; unset LANGUAGE

# Parameters passing by console

if [ $# -ne 1 ]
  then
    echo "A file with a list of platekeys/genomes and their path (separated by comma) is required"
fi

input_list=$1

cat $input_list | while read line; do

        IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
        lp_id=${NAMES[0]}
        path_to_bam=${NAMES[1]}
        #path_to_bam=${path_to_bam}'/Assembly/'${lp_id}'.bam'

        output_json=${output_folder}'/'${lp_id}'_EHdeNovo'

    echo ${eh_binary} profile --reads ${path_to_bam} --reference ${reference_fasta} --output-prefix ${output_json} --min-anchor-mapq 50 --max-irr-mapq 60
    ${eh_binary} profile --reads ${path_to_bam} --reference ${reference_fasta} --output-prefix ${output_json} --min-anchor-mapq 50 --max-irr-mapq 60

done
