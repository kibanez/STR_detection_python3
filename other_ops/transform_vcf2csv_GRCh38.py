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
from operator import itemgetter


class OptionParser(optparse.OptionParser):
    def check_required(self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)

        if atrib is None:
            return False
        else:
            return True


def print_tables(hash_table, f_output, l_samples):
    """
    Function that prints csv file from annotated VCF file

    :param hash_table: it contains the whole VCF file records
    :param f_output: csv file to write to
    :param l_samples: list of sample IDs
    :return:
    """

    l_fields = ['chr', 'pos', 'ref', 'alt', 'QUAL', 'FILTER',
                'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
                'cytoBand', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH',
                'ExAC_SAS',
                'avsnp147', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred',
                'Polyphen2_HVAR_score',
                'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred',
                'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score',
                'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score',
                'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred',
                'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate',
                'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds']
    l_fields = l_fields + l_samples
    
    l_chr = set([item[0] for item in hash_table.keys()])

    fo = open(f_output, 'w')
    fo.write(','.join(l_fields) + '\n')
    for key in sorted(hash_table.keys(), key=itemgetter(1)):
        fo.write(','.join(map(lambda field: hash_table[key].get(field, '.'), l_fields)) + '\n')
    fo.close()


def read_annovar_vcf(input_vcf):
    """
    Function that reads annotated VCF file with PyVCF
    :param input_vcf: annotated VCF with annovar
    :return: hash table with all records (defined as tuple of (chr, pos, alt)) with all corresponding info, genomic
    and functional
    """
    hash_table = {}
    vcf_reader = vcf.Reader(filename=input_vcf)

    for i, r in enumerate(vcf_reader):
        hash_variant = {}

        hash_fields = dict(r.INFO)
        hash_fields.update(dict(zip(r.samples[0].data._fields, r.samples[0].data)))

        chrom = r.CHROM
        pos = str(r.POS)
        ref = str(r.REF)
        alt = str(r.ALT[0])
        l_samples = len(r.samples)

        if r.FILTER == []:
            hash_variant['FILTER'] = "PASS"
        else:
            hash_variant['FILTER'] = str(r.FILTER)

        hash_variant['QUAL'] = str(r.QUAL)

        hash_variant['chr'] = chrom.strip()
        hash_variant['pos'] = pos.strip()
        hash_variant['ref'] = ref.strip()
        hash_variant['alt'] = alt.strip()
        hash_variant['Func.refGene'] = str(hash_fields.get('Func.refGene', '.')[0])
        hash_variant['Gene.refGene'] = str(hash_fields.get('Gene.refGene', '.')[0])
        hash_variant['GeneDetail.refGene'] = str(hash_fields.get('GeneDetail.refGene', '.')[0])
        hash_variant['ExonicFunc.refGene'] = str(hash_fields.get('ExonicFunc.refGene', '.')[0])
        hash_variant['AAChange.refGene'] = str(hash_fields.get('AAChange.refGene', '.')[0])
        hash_variant['cytoBand'] = str(hash_fields.get('cytoBand', '.')[0])
        hash_variant['ExAC_ALL'] = str(hash_fields.get('ExAC_ALL', '.'))
        hash_variant['ExAC_AFR'] = str(hash_fields.get('ExAC_AFR', '.'))
        hash_variant['ExAC_AMR'] = str(hash_fields.get('ExAC_AMR', '.'))
        hash_variant['ExAC_EAS'] = str(hash_fields.get('ExAC_EAS', '.'))
        hash_variant['ExAC_FIN'] = str(hash_fields.get('ExAC_FIN', '.'))
        hash_variant['ExAC_NFE'] = str(hash_fields.get('ExAC_NFE', '.'))
        hash_variant['ExAC_OTH'] = str(hash_fields.get('ExAC_OTH', '.'))
        hash_variant['ExAC_SAS'] = str(hash_fields.get('ExAC_SAS', '.'))
        hash_variant['avsnp147'] = str(hash_fields.get('avsnp147', '.')[0])
        hash_variant['SIFT_score'] = str(hash_fields.get('SIFT_score', '.')[0])
        hash_variant['SIFT_pred'] = str(hash_fields.get('SIFT_pred', '.')[0])
        hash_variant['Polyphen2_HDIV_score'] = str(hash_fields.get('Polyphen2_HDIV_score', '.')[0])
        hash_variant['Polyphen2_HDIV_pred'] = str(hash_fields.get('Polyphen2_HDIV_pred', '.')[0])
        hash_variant['Polyphen2_HVAR_score'] = str(hash_fields.get('Polyphen2_HVAR_score', '.')[0])
        hash_variant['Polyphen2_HVAR_pred'] = str(hash_fields.get('Polyphen2_HVAR_pred', '.')[0])
        hash_variant['LRT_score'] = str(hash_fields.get('LRT_score', '.')[0])
        hash_variant['LRT_pred'] = str(hash_fields.get('LRT_pred', '.')[0])
        hash_variant['MutationTaster_score'] = str(hash_fields.get('MutationTaster_score', '.')[0])
        hash_variant['MutationTaster_pred'] = str(hash_fields.get('MutationTaster_pred', '.')[0])
        hash_variant['MutationAssessor_score'] = str(hash_fields.get('MutationAssessor_score', '.')[0])
        hash_variant['MutationAssessor_pred'] = str(hash_fields.get('MutationAssessor_pred', '.')[0])
        hash_variant['FATHMM_score'] = str(hash_fields.get('FATHMM_score', '.')[0])
        hash_variant['FATHMM_pred'] = str(hash_fields.get('FATHMM_pred', '.')[0])
        hash_variant['PROVEAN_score'] = str(hash_fields.get('PROVEAN_score', '.')[0])
        hash_variant['PROVEAN_pred'] = str(hash_fields.get('PROVEAN_pred', '.')[0])
        hash_variant['VEST3_score'] = str(hash_fields.get('VEST3_score', '.')[0])
        hash_variant['CADD_raw'] = str(hash_fields.get('CADD_raw', '.')[0])
        hash_variant['CADD_phred'] = str(hash_fields.get('CADD_phred', '.')[0])
        hash_variant['DANN_score'] = str(hash_fields.get('DANN_score', '.')[0])
        hash_variant['fathmm-MKL_coding_score'] = str(hash_fields.get('fathmm-MKL_coding_score', '.')[0])
        hash_variant['fathmm-MKL_coding_pred'] = str(hash_fields.get('fathmm-MKL_coding_pred', '.')[0])
        hash_variant['MetaSVM_score'] = str(hash_fields.get('MetaSVM_score', '.')[0])
        hash_variant['MetaSVM_pred'] = str(hash_fields.get('MetaSVM_pred', '.')[0])
        hash_variant['MetaLR_score'] = str(hash_fields.get('MetaLR_score', '.')[0])
        hash_variant['MetaLR_pred'] = str(hash_fields.get('MetaLR_pred', '.')[0])
        hash_variant['integrated_fitCons_score'] = str(hash_fields.get('integrated_fitCons_score', '.')[0])
        hash_variant['integrated_confidence_value'] = str(hash_fields.get('integrated_confidence_value', '.')[0])
        hash_variant['GERP++_RS'] = str(hash_fields.get('GERP++_RS', '.')[0])
        hash_variant['phyloP7way_vertebrate'] = str(hash_fields.get('phyloP7way_vertebrate', '.')[0])
        hash_variant['phyloP20way_mammalian'] = str(hash_fields.get('phyloP20way_mammalian', '.')[0])
        hash_variant['phastCons7way_vertebrate'] = str(hash_fields.get('phastCons7way_vertebrate', '.')[0])
        hash_variant['phastCons20way_mammalian'] = str(hash_fields.get('phastCons20way_mammalian', '.')[0])
        hash_variant['SiPhy_29way_logOdds'] = str(hash_fields.get('SiPhy_29way_logOdds', '.')[0])

        l_samples = r.samples[::]
        l_sample_ids = []
        for sample in l_samples:
            sample_id = sample.sample
            sample_gt = sample.data.GT
            hash_variant[sample_id] = sample_gt
            l_sample_ids.append(sample_id)

        hash_table[(chrom, pos, alt)] = hash_variant

    return hash_table, l_sample_ids


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

    hash_table, l_samples = read_annovar_vcf(input_vcf)
    print_tables(hash_table, output_file, l_samples)


if __name__ == '__main__':
    run()