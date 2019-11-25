"""
Created on 20/08/2019

@author: kristina ibanez-garikano

Function that takes all VCF generated by EH v3.0.0, merging them and computing AF for each locus
"""

import copy
import logging
import optparse
import os
import re
import sys
import vcf
from configparser import ConfigParser
from operator import itemgetter
from os import path as osp


localModulesBase = osp.dirname(osp.realpath(__file__))
modulesRelDirs = ["../modules/"]

for moduleRelDir in modulesRelDirs:
    sys.path.insert(0, osp.join(localModulesBase, moduleRelDir))


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
                'num_samples', 'AF', 'LC', 'list_samples']

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


def merging_vcf(l_vcf, path_vcf, logger):
    """
    Function that receives a list of individual VCF files to merge them all

    The main objective here is to extract the GT frequencies for each STR loci (being trinucleotide regions)
    It returns a VCF file containing the AF for each STR detected in the whole cohort
    chr   pos ref alt frequency

    :param l_vcf: list of VCF files generated by EH
    :param path_vcf: folder in where l_vcf are located
    :param logger:
    :return:
    """

    hash_table = {}

    total_samples = len(l_vcf)

    for vcf_input in l_vcf:

        name_vcf = vcf_input

        vcf_input = os.path.join(path_vcf, vcf_input)

        if not os.path.isfile(vcf_input):
            raise IOError("The VCF file does not exist %s " % vcf_input)

        logger.info("Analysing the STR regions within %s" % vcf_input)

        # Sometimes the VCF generated by EH is empty. We need to check whether the VCF is empty or not
        if os.path.getsize(vcf_input) == 0:
            print(vcf_input + '\n')
            continue

        vcf_reader = vcf.Reader(filename=vcf_input)

        for i, r in enumerate(vcf_reader):
            hash_fields = dict(r.INFO)
            hash_fields.update(dict(zip(r.samples[0].data._fields, r.samples[0].data)))

            # We want to take the max value from CI field
            ci_field = hash_fields.get('REPCI', '0').split('/')
            max_ci_allele1 = ci_field[0].split('-')[1]
            if len(ci_field) > 1:
                max_ci_allele2 = ci_field[1].split('-')[1]

            # We only analyse STR markers >=3 length
            if len(str(hash_fields.get('RU', 0))) >= 3:
                pos = str(r.INFO['END'] - 1)
                gene = str(hash_fields.get('REPID', ""))

                # we retrieve all the info contained in the INFO fields
                hash_variant = {}
                hash_variant['Repeat_Motif'] = str(hash_fields.get('RU', 0))
                hash_variant['Reference_length_bp'] = str(hash_fields.get('RL', 0))
                hash_variant['gene'] = gene
                hash_variant['gt'] = str(hash_fields.get('GT', ""))
                hash_variant['chr'] = r.CHROM
                hash_variant['start'] = str(r.POS)
                hash_variant['end'] = str(r.INFO['END'])
                hash_variant['ref'] = r.REF

                hash_variant['alt_size'] = str(max_ci_allele1)
                hash_variant['ref_size'] = str(max_ci_allele2)
                hash_variant['num_samples'] = '1'
                hash_variant['list_samples'] = name_vcf
                hash_variant['LC'] = str(hash_fields.get('LC', 0))

                if hash_variant.get('gt') == '0/0':
                    allele = max_ci_allele1
                    hash_variant['allele'] = allele

                    if (r.CHROM, pos, gene, allele) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele)]['num_samples'] = str(int(hash_table.get((r.CHROM, pos, gene, allele))['num_samples']) + 2)
                        hash_table[(r.CHROM, pos, gene, allele)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele))['list_samples'] + ';' + \
                            name_vcf + \
                            '_x2'

                    else:
                        # we add the new alternative allele info
                        hash_variant['num_samples'] = '2'
                        # we update the number of the genome + '_x2)
                        hash_variant['list_samples'] = name_vcf + '_x2'
                        hash_table[(r.CHROM, pos, gene, allele)] = hash_variant

                elif hash_variant.get('gt') == '1/1':
                    allele = max_ci_allele1
                    hash_variant['allele'] = allele

                    if (r.CHROM, pos, gene, allele) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele)]['num_samples'] = \
                            str(int(hash_table.get((r.CHROM, pos, gene, allele))['num_samples']) + 2)
                        hash_table[(r.CHROM, pos, gene, allele)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele))['list_samples'] + \
                            ';' + \
                            name_vcf + \
                            '_x2'
                    else:
                        # we add the new alternative allele info
                        hash_variant['num_samples'] = '2'
                        # we update the number of the genome + '_x2)
                        hash_variant['list_samples'] = name_vcf + '_x2'
                        hash_table[(r.CHROM, pos, gene, allele)] = hash_variant

                elif hash_variant.get('gt') == '1' or hash_variant.get('gt') == '0':
                    allele = max_ci_allele1
                    hash_variant['allele'] = allele

                    if (r.CHROM, pos, gene, allele) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele)]['num_samples'] = str(
                            int(hash_table.get((r.CHROM, pos, gene, allele))['num_samples']) + 1)
                        hash_table[(r.CHROM, pos, gene, allele)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele))['list_samples'] \
                            + ';' + name_vcf
                    else:
                        # we add the new alternative allele info
                        hash_variant['num_samples'] = '1'
                        hash_variant['list_samples'] = name_vcf
                        hash_table[(r.CHROM, pos, gene, allele)] = hash_variant

                elif hash_variant.get('gt') == '0/1':
                    allele_ref = max_ci_allele1
                    allele_alt = max_ci_allele2

                    if (r.CHROM, pos, gene, allele_ref) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele_ref)]['num_samples'] = str(
                            int(hash_table.get((r.CHROM, pos, gene, allele_ref))['num_samples']) + 1)
                        hash_table[(r.CHROM, pos, gene, allele_ref)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele_ref))['list_samples'] + ';' + name_vcf
                    else:
                        # we specify the allele
                        hash_variant['allele'] = allele_ref
                        # we add the new alternative allele info
                        hash_table[(r.CHROM, pos, gene, allele_ref)] = hash_variant

                    if (r.CHROM, pos, gene, allele_alt) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele_alt)]['num_samples'] = str(
                            int(hash_table.get((r.CHROM, pos, gene, allele_alt))['num_samples']) + 1)
                        hash_table[(r.CHROM, pos, gene, allele_alt)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele_alt))['list_samples'] + ';' + name_vcf
                    else:
                        hash_variant_alt = copy.deepcopy(hash_variant)
                        hash_variant_alt['allele'] = allele_alt
                        # we add the new alternative allele info
                        hash_table[(r.CHROM, pos, gene, allele_alt)] = hash_variant_alt

                elif hash_variant.get('gt') == '1/2':
                    allele_alt1 = max_ci_allele1
                    allele_alt2 = max_ci_allele2

                    if (r.CHROM, pos, gene, allele_alt1) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele_alt1)]['num_samples'] = str(
                            int(hash_table.get((r.CHROM, pos, gene, allele_alt1))['num_samples']) + 1)
                        hash_table[(r.CHROM, pos, gene, allele_alt1)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele_alt1))['list_samples'] + ';' + name_vcf
                    else:
                        # we specify the allele
                        hash_variant['allele'] = allele_alt1
                        # we add the new alternative allele info
                        hash_table[(r.CHROM, pos, gene, allele_alt1)] = hash_variant

                    if (r.CHROM, pos, gene, allele_alt2) in hash_table:
                        # we add info to an existing key - repeat size allele
                        hash_table[(r.CHROM, pos, gene, allele_alt2)]['num_samples'] = str(
                            int(hash_table.get((r.CHROM, pos, gene, allele_alt2))['num_samples']) + 1)
                        hash_table[(r.CHROM, pos, gene, allele_alt2)]['list_samples'] = \
                            hash_table.get((r.CHROM, pos, gene, allele_alt2))['list_samples'] + ';' + name_vcf
                    else:
                        # we specify the allele
                        hash_variant_alt = copy.deepcopy(hash_variant)
                        hash_variant_alt['allele'] = allele_alt2
                        # we add the new alternative allele info
                        hash_table[(r.CHROM, pos, gene, allele_alt2)] = hash_variant_alt

    for key, value in iter(hash_table.items()):
        # we calculate the AF for each STR site (for each alternate allele)
        af = float(hash_table.get(key)['num_samples']) / float(total_samples)
        hash_table[key]['AF'] = str(af)

    logger.info("Computing the allele frequencies after having digested all the individual VCF files")

    return hash_table


def read_cfg_file(cfg_filename):
    '''
    Function that reads the configuration file which includes the paths to all the fundamental input

    :param cfg_filename: configuration file
    :return: hash_table containing info in the config file
    '''
    fi = open(cfg_filename, 'r')

    config = ConfigParser.ConfigParser()
    config.readfp(fi)

    hash_cfg = {}

    for field in config.options('GENERAL'):
        hash_cfg[field] = config.get('GENERAL', field)

    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT', field)

    for field in config.options('REFERENCE'):
        hash_cfg[field] = config.get('REFERENCE', field)

    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE', field)

    fi.close()

    return hash_cfg


def run(argv=None):

    if argv is None: argv = sys.argv

    parser = OptionParser(add_help_option=True, description="")
    parser.add_option("--s", default=None,
                      help="The path in which the resulting EH VCF files are", dest="f_samples")
    parser.add_option("--o", default=None,
                      help="The output VCF name in which the merged VCF will be write",
                      dest="f_output")
    parser.add_option("--O", default=None,
                      help="The output directory in which the merged VCF will be write",
                      dest="d_output")
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)

    if not parser.check_required("--s"):
        raise IOError('The path to the VCF files generated by running EH is missing')

    if not parser.check_required("--o"):
        raise IOError('The merged VCF file is missing')

    if not parser.check_required("--O"):
        raise IOError('The output directory in which the merged VCF files will be write is missing')

    if options.f_samples != None:

        path_samples = options.f_samples

        if not os.path.exists(path_samples):
            raise IOError('The path to the VCF files generated by running EH %s does not exist') \
                  % path_samples

        # Output folder in which the output of EH will be saved
        output_folder = options.d_output

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        merged_vcf = options.f_output

        if merged_vcf is None:
            raise IOError('The name for the merged VCF is missing %s') % merged_vcf

        output_file = os.path.join(output_folder, merged_vcf)


        # Configure logger
        formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        console.setLevel(logging.INFO)
        logger = logging.getLogger("preprocess")
        logger.setLevel(logging.INFO)
        logger.addHandler(console)

        logger.info("The process of merging all the individual VCF files has started")

        # Here we retrieve all the VCF files generated by running EH algorithm in a cohort
        # EH for each sample generates a *.json, *.vcf and *.log file
        # we are interested in the VCF files

        l_vcf = [f for f in os.listdir(path_samples) if f.endswith('.vcf')]

        l_samples = []

        for vcf in l_vcf:
            sample_name = re.sub('^EH_', '', vcf)
            sample_name = re.sub('.vcf$', '', sample_name)
            l_samples.append(sample_name)

        hash_table = merging_vcf(l_vcf, path_samples, logger)

        print_tables(hash_table, output_file)

        logger.info('The merged VCF files has been annotated and enriched - %s' % output_file)


if __name__ == '__main__':
    run()

