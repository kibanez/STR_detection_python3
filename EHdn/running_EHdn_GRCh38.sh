# Bash script to run ExpansionHunter denovo (EHdn) across all genomes aligned with GRCh37 human genome assembly

eh_binary='/genomes/scratch/kgarikano/GEL_STR/sw/latest_EHdn'
reference_fasta='/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa'
output_folder='/genomes/scratch/kgarikano/GEL_STR/EHdn/output_EHdnv0.8.0/GRCh38'

# Create dir if does not exist
mkdir -p $output_folder

export LC_ALL=C; unset LANGUAGE

input_list=''

cat $input_list | while read line; do

        IFS=’,’ read -ra NAMES <<< "$line"
        lp_id=${NAMES[0]}
        path_to_bam=${NAMES[1]}

        output_json=$output_folder$lp_id'_EHdeNovo'

    echo $eh_binary profile --reads $path_to_bam --reference $reference_fasta --output-prefix $output_json --min-anchor-mapq 50 --max-irr-mapq 60
    $eh_binary profile --reads $path_to_bam --reference $reference_fasta --output-prefix $output_json --min-anchor-mapq 50 --max-irr-mapq 60
done

