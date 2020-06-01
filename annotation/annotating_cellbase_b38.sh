module load cellbase

if [[ $# -eq 0 ]] ; then
    echo 'Two arguments are required: 1) the VCF file to be annotated and 2) the output directory to save the corresponding annotated json file'
    exit 1
fi


# Input parameters
INPUT_VCF=$1
OUTPUT_DIR=$2


OUTPUT_FILE=${OUTPUT_DIR}/"${INPUT_VCF##*/}".json.gz

cellbase.sh variant-annotation \
--input-file ${INPUT_VCF} \
--num-threads 8 \
--assembly GRCh38 \
--output ${OUTPUT_FILE} \
--species hsapiens \
--local \
--exclude expression,geneDisease,drugInteraction,conservation,functionalScore \
--custom-file /genomes/resources/genomeref/data/human/GRCh38_annotation/agg_db/GEL_GL_7273/GEL_GL_6628.duprem.sites.annot.subset.atomic.left.split.vcf.gz,/home/kgarikano/annotation/external_datasets/HGMD_PRO_2018/HGMD_PRO_2018.1_hg38_id_2_info.duprem.atomic.left.split_chrNomenclatureChanged.vcf.gz,/home/kgarikano/annotation/external_datasets/dbSNP/All_20161122.vcf,/home/kgarikano/annotation/external_datasets/gnomad/gnomad.genomes.r2.0.1.sites.coding.autosomes.picard.hg38.final.atomic.left.split.dedupped.vcf \
--custom-file-id GEL.GL.6628,HGMD_2018,dbSNP_human_9606_b149,gnomad_170228 \
--custom-file-fields AF,AF_RD,AF_Cancer,AC,AC_RD,AC_Cancer,AN,AN_RD,AN_Cancer,GN:CLASS,MUT,GENE,DNA,PROT,DB,tu,HGMD_IDs:RS,dbSNPBuildID,GENEINFO:AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_Female,AF_Male,AN_AFR,AN_AMR,AN_ASJ,AN_EAS,AN_FIN,AN_NFE,Hom_AFR,Hom_AMR,Hom_ASJ,Hom_EAS,Hom_FIN,Hom_NFE,CSQ
