module load bio/cellbase/v4.7.1

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
--assembly GRCh37 \
--output ${OUTPUT_FILE} \
--species hsapiens \
--local \
--exclude expression,geneDisease,drugInteraction,conservation,functionalScore \
--custom-file /genomes/resources/genomeref/data/human/GRCh37_annotation/agg4_db/GEL_GL_5523/GEL_GL_5277.sites.annot.subset.duprem.atomic.left.split.vcf.gz,/genomes/resources/genomeref/data/human/GRCh37_annotation/GEL_Platypus/GEL_Platypus_RD_1777/GEL_Platypus_RD_1777.sites.vcf.gz,/genomes/resources/genomeref/data/human/GRCh37_annotation/HGMD/2017.1/HGMD_PRO_2017.1_hg19_id_2_info.duprem.atomic.left.split.vcf.gz,/home/kgarikano/annotation/external_datasets/gnomad/gnomad.genomes.r2.0.1.sites.coding.autosomes.vcf \
--custom-file-id GEL.GL.5277,GEL.Platypus.RD.1777,HGMD.2015.4,gnomad \
--custom-file-fields AF,AC,AN,GN,RD_AC,RD_AN,CG_AC,CG_AN,V1_AC,V1_AN,V2_AC,V2_AN,125_AC,125_AN,150_AC,150_AN:AF,AC,AN,N_HET,N_HOM:CLASS,MUT,PHEN,HGMD_IDs:AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_Female,AF_Male,AN_AFR,AN_AMR,AN_ASJ,AN_EAS,AN_FIN,AN_NFE,Hom_AFR,Hom_AMR,Hom_ASJ,Hom_EAS,Hom_FIN,Hom_NFE,CSQ

