"""
Created on 17/09/20

@author: kristina ibanez-garikano

Script that converts an annotated VCF file into a csv/tsv file in order to open in Excel
"""
import logging
import optparse
import os
import sys
import vcf


class OptionParser(optparse.OptionParser):
    def check_required(self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)

        if atrib is None:
            return False
        else:
            return True


def print_tables(hash_table, f_output):
    """
    Function that prints two tsv files: one containing all the EH VCF, already enriched somehow
    and the other table, containing only the STR loci that are spotted by having more repetitions that theoretically
    and in practice have been seen

    :param hash_table: it contains the frequencies for each locus
    :param f_output: file in where the function will write the info in hash_table
    :return:
    """

    l_fields = ['chr', 'start', 'end', 'allele', 'gene', 'ref', 'alt', 'Repeat_Motif',
                'num_samples', 'AF', 'list_samples']

    l_chr = set([item[0] for item in hash_table.keys()])

    chr_specials = []
    if 'X' in l_chr:
        l_chr.remove('X')
        chr_specials.append('X')
    elif 'Y' in l_chr:
        l_chr.remove('Y')
        chr_specials.append('Y')

    l_chr = sorted(l_chr)
    l_chr = [str(i) for i in l_chr]

    for i in chr_specials:
        l_chr.append(i)

    fo = open(f_output, 'w')
    fo.write('\t'.join(l_fields) + '\n')
    for key in sorted(sorted(hash_table.keys(), key=itemgetter(1)), key=lambda x: l_chr.index(x[0])):
        fo.write('\t'.join(map(lambda field: hash_table[key].get(field, '.'), l_fields)) + '\n')
    fo.close()


def read_vcf(input_vcf, logger):
    hash_table = {}
    vcf_reader = vcf.Reader(filename=input_vcf)

    for i, r in enumerate(vcf_reader):
        hash_fields = dict(r.INFO)
        hash_fields.update(dict(zip(r.samples[0].data._fields, r.samples[0].data)))

        chrom = r.CHROM
        pos = r.POS
        alt = r.ALT
        l_samples = len(r.samples)

        hash_variant = {}
        hash_variant['QUAL'] = str(hash_fields.get('QUAL', '.'))
        hash_variant['FILTER'] = str(hash_fields.get('FILTER', '.'))
        hash_variant['Func.refGene'] = str(hash_fields.get('Func.refGene', '.'))
        hash_variant['Gene.refGene'] = str(hash_fields.get('Gene.refGene', '.'))
        hash_variant['GeneDetail.refGene'] = str(hash_fields.get('GeneDetail.refGene', '.'))
        hash_variant['ExonicFunc.refGene'] = str(hash_fields.get('ExonicFunc.refGene', '.'))
        hash_variant['AAChange.refGene'] = str(hash_fields.get('AAChange.refGene', '.'))
        hash_variant['cytoBand'] = str(hash_fields.get('cytoBand', '.'))
        hash_variant['ExAC_ALL'] = str(hash_fields.get('ExAC_ALL', '.'))
        hash_variant['ExAC_AFR'] = str(hash_fields.get('ExAC_AFR', '.'))
        hash_variant['ExAC_AMR'] = str(hash_fields.get('ExAC_AMR', '.'))
        hash_variant['ExAC_EAS'] = str(hash_fields.get('ExAC_EAS', '.'))
        hash_variant['ExAC_FIN'] = str(hash_fields.get('ExAC_FIN', '.'))
        hash_variant['ExAC_NFE'] = str(hash_fields.get('ExAC_NFE', '.'))
        hash_variant['ExAC_OTH'] = str(hash_fields.get('ExAC_OTH', '.'))
        hash_variant['ExAC_SAS'] = str(hash_fields.get('ExAC_SAS', '.'))
        hash_variant['avsnp147'] = str(hash_fields.get('avsnp147', '.'))
        hash_variant['SIFT_score'] = str(hash_fields.get('SIFT_score', '.'))
        hash_variant['SIFT_pred'] = str(hash_fields.get('SIFT_pred', '.'))
        hash_variant['Polyphen2_HDIV_score'] = str(hash_fields.get('Polyphen2_HDIV_score', '.'))
        hash_variant['Polyphen2_HDIV_pred'] = str(hash_fields.get('Polyphen2_HDIV_pred', '.'))
        hash_variant['Polyphen2_HVAR_score'] = str(hash_fields.get('Polyphen2_HVAR_score', '.'))
        hash_variant['Polyphen2_HVAR_pred'] = str(hash_fields.get('Polyphen2_HVAR_pred', '.'))
        hash_variant['LRT_score'] = str(hash_fields.get('LRT_score', '.'))
        hash_variant['LRT_pred'] = str(hash_fields.get('LRT_pred', '.'))
        hash_variant['MutationTaster_score'] = str(hash_fields.get('MutationTaster_score', '.'))
        hash_variant['MutationTaster_pred'] = str(hash_fields.get('MutationTaster_pred', '.'))
        hash_variant['MutationAssessor_score'] = str(hash_fields.get('MutationAssessor_score', '.'))
        hash_variant['MutationAssessor_pred'] = str(hash_fields.get('MutationAssessor_pred', '.'))
        hash_variant['FATHMM_score'] = str(hash_fields.get('FATHMM_score', '.'))
        hash_variant['FATHMM_pred'] = str(hash_fields.get('FATHMM_pred', '.'))
        hash_variant['PROVEAN_score'] = str(hash_fields.get('PROVEAN_score', '.'))
        hash_variant['PROVEAN_pred'] = str(hash_fields.get('PROVEAN_pred', '.'))
        hash_variant['VEST3_score'] = str(hash_fields.get('VEST3_score', '.'))
        hash_variant['CADD_raw'] = str(hash_fields.get('CADD_raw', '.'))
        hash_variant['CADD_phred'] = str(hash_fields.get('CADD_phred', '.'))
        hash_variant['DANN_score'] = str(hash_fields.get('DANN_score', '.'))
        hash_variant['fathmm-MKL_coding_score'] = str(hash_fields.get('fathmm-MKL_coding_score', '.'))
        hash_variant['fathmm-MKL_coding_pred'] = str(hash_fields.get('fathmm-MKL_coding_pred', '.'))
        hash_variant['MetaSVM_score'] = str(hash_fields.get('MetaSVM_score', '.'))
        hash_variant['MetaSVM_pred'] = str(hash_fields.get('MetaSVM_pred', '.'))
        hash_variant['MetaLR_score'] = str(hash_fields.get('MetaLR_score', '.'))
        hash_variant['MetaLR_pred'] = str(hash_fields.get('MetaLR_pred', '.'))
        hash_variant['integrated_fitCons_score'] = str(hash_fields.get('integrated_fitCons_score', '.'))
        hash_variant['integrated_confidence_value'] = str(hash_fields.get('integrated_confidence_value', '.'))
        hash_variant['GERP++_RS'] = str(hash_fields.get('GERP++_RS', '.'))
        hash_variant['phyloP7way_vertebrate'] = str(hash_fields.get('phyloP7way_vertebrate', '.'))
        hash_variant['phyloP20way_mammalian'] = str(hash_fields.get('phyloP20way_mammalian', '.'))
        hash_variant['phastCons7way_vertebrate'] = str(hash_fields.get('phastCons7way_vertebrate', '.'))
        hash_variant['phastCons20way_mammalian'] = str(hash_fields.get('phastCons20way_mammalian', '.'))
        hash_variant['SiPhy_29way_logOdds'] = str(hash_fields.get('SiPhy_29way_logOdds', '.'))

        for vcf in l_samples:
            hash_variant[vcf] = "."

        hash_table[(chrom, pos, alt)] = hash_variant



def run(argv=None):

    if argv is None: argv = sys.argv

    parser = OptionParser(add_help_option=True, description="")
    parser.add_option("--v", default=None,
                      help="Absolute path to annotated VCF file", dest="input_vcf")
    parser.add_option("--o", default=None,
                      help="Output file name",
                      dest="f_output")
    parser.add_option("--O", default=None,
                      help="The output directory in which the tabulated data will be written",
                      dest="d_output")
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)

    if not parser.check_required("--v"):
        raise IOError('Absolute path to annotated VCF file is missing')

    if not parser.check_required("--o"):
        raise IOError('Output file name is missing')

    if not parser.check_required("--O"):
        raise IOError('The output directory in which the tabulated data will be written is missing')



    input_vcf = options.input_vcf

    if not os.path.isfile(input_vcf):
        raise IOError('The path to the annotated VCF file %s does not exist') \
              % input_vcf

    # Output folder in which the tabulated output will be written
    output_folder = options.d_output

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_csv = options.f_output

    if output_csv is None:
        raise IOError('The name for the tabulated CSV file is missing %s') % output_csv

    output_file = os.path.join(output_folder, output_csv)

    # Configure logger
    formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)
    logger = logging.getLogger("preprocess")
    logger.setLevel(logging.INFO)
    logger.addHandler(console)

    hash_table = read_vcf(input_vcf, logger)
    print_tables(hash_table, output_file)


if __name__ == '__main__':
    run()