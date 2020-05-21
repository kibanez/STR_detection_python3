import json
import sys
from operator import itemgetter

if len(sys.argv) < 2:
    sys.exit('ERROR: Usage: %s <json> <outFile>' % sys.argv[0])


def print_tables(hash_table, f_output):
    """
    Function that prints the meat of hash_table

    :param hash_table: it contains the frequencies and enrichment from cellbase annotation
    :param f_output: file in where the function will write the info in hash_table
    :return:
    """

    l_fields = ['chr', 'start', 'ref', 'alt', 'type', 'rsid', 'consequenceType',
                'GEL.GL.5277', 'GEL.Platypus.RD.1777',
                'GNOMAD_GENOMES_ALL', 'GNOMAD_GENOMES_AFR', 'GNOMAD_GENOMES_AMR', 'GNOMAD_GENOMES_EAS',
                'GNOMAD_GENOMES_NFE', 'GNOMAD_GENOMES_FIN', 'GNOMAD_GENOMES_FEMALE', 'GNOMAD_GENOMES_MALE',
                '1kG_phase3_ALL', '1kG_phase3_SAS', '1kG_phase3_AFR', '1kG_phase3_EUR', '1kG_phase3_AMR',
                '1kG_phase3_EAS', 'UK10K_ALL', 'UK10K_TWINSUK_NODUP', 'UK10K_TWINSUK', 'ALSPAC',
                'HGMD_version', 'HGMD_ID', 'HGMD_PHEN', "HGMD_CLASS", "HGMD_MUT"]

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
    fo.write(','.join(l_fields) + '\n')
    for key in sorted(sorted(hash_table.keys(), key=itemgetter(1)), key=lambda x: l_chr.index(x[0])):
        fo.write(','.join(map(lambda field: hash_table[key].get(field, '.'), l_fields)) + '\n')
    fo.close()


# Read input JSON
annotFile = sys.argv[1]
outFile = sys.argv[2]

# Parse JSON
hash_table = {}
with open(annotFile) as f:
    for line in f:
        hash_variant = {}

        annot = json.loads(line)
        cellbaseFrq = ""
        noStudy = ""
        outLine = annot["annotation"]["chromosome"] + "\t" + str(annot["annotation"]["start"]) + "\t" + \
                  annot["annotation"]["reference"] + "\t" + annot["annotation"]["alternate"]

        # Principal keys
        chr = str(annot['chromosome'])
        start = str(annot['start'])
        ref = str(annot['reference'])
        alt = str(annot['alternate'])

        # Include single information
        hash_variant['chr'] = chr
        hash_variant['start'] = start
        hash_variant['ref'] = ref
        hash_variant['alt'] = alt
        hash_variant['type'] = annot.get("type", ".")
        hash_variant['rsid'] = annot.get('id', ".")

        if "additionalAttributes" in annot["annotation"]:
            if "GEL.GL.5277" in annot["annotation"]["additionalAttributes"]:
                hash_variant['GEL.GL.5277'] = \
                    annot["annotation"]["additionalAttributes"]["GEL.GL.5277"]["attribute"]["AF"]
            if "GEL.Platypus.RD.1777" in annot["annotation"]["additionalAttributes"]:
                hash_variant['GEL.Platypus.RD.1777'] = \
                    annot["annotation"]["additionalAttributes"]["GEL.Platypus.RD.1777"]['attribute']['AF']
        else:
            hash_variant['GEL.GL.5277'] = '.'
            hash_variant['GEL.Platypus.RD.1777'] = '.'

        # Extract CellBase annotation for GNOMAD_EXOMES-ALL
        if "populationFrequencies" in annot["annotation"] and len(annot["annotation"]["populationFrequencies"]) > 0:
            for pf in annot["annotation"]["populationFrequencies"]:
                if pf["study"] == "GNOMAD_GENOMES":
                    if pf["population"] == "ALL":
                        hash_variant['GNOMAD_GENOMES_ALL'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "AFR":
                        hash_variant['GNOMAD_GENOMES_AFR'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "AMR":
                        hash_variant['GNOMAD_GENOMES_AMR'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "EAS":
                        hash_variant['GNOMAD_GENOMES_EAS'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "NFE":
                        hash_variant['GNOMAD_GENOMES_NFE'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "FIN":
                        hash_variant['GNOMAD_GENOMES_FIN'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "FEMALE":
                        hash_variant['GNOMAD_GENOMES_FEMALE'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "MALE":
                        hash_variant['GNOMAD_GENOMES_MALE'] = str(pf["altAlleleFreq"])
                elif pf["study"] == "1kG_phase3":
                    if pf["population"] == "ALL":
                        hash_variant['1kG_phase3_ALL'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "SAS":
                        hash_variant['1kG_phase3_SAS'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "AFR":
                        hash_variant['1kG_phase3_AFR'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "EUR":
                        hash_variant['1kG_phase3_EUR'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "AMR":
                        hash_variant['1kG_phase3_AMR'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "EAS":
                        hash_variant['1kG_phase3_EAS'] = str(pf["altAlleleFreq"])
                elif pf["study"] == "UK10K":
                    if pf["population"] == "ALL":
                        hash_variant['UK10K_ALL'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "TWINSUK_NODUP":
                        hash_variant['UK10K_TWINSUK_NODUP'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "TWINSUK":
                        hash_variant['UK10K_TWINSUK'] = str(pf["altAlleleFreq"])
                    elif pf["population"] == "ALSPAC":
                        hash_variant['ALSPAC'] = str(pf["altAlleleFreq"])

        if "additionalAttributes" in annot["annotation"]:
            if "HGMD.2015.4" in annot["annotation"]["additionalAttributes"]:
                hash_variant["HGMD_version"] = "HGMD.2015.4"
                hash_variant["HGMD_ID"] = \
                    annot["annotation"]["additionalAttributes"]["HGMD.2015.4"]["attribute"]["HGMD_IDs"]
                hash_variant["HGMD_PHEN"] = \
                    annot["annotation"]["additionalAttributes"]["HGMD.2015.4"]["attribute"]["PHEN"]
                hash_variant["HGMD_CLASS"] = \
                    annot["annotation"]["additionalAttributes"]["HGMD.2015.4"]["attribute"]["CLASS"]
                hash_variant["HGMD_MUT"] = \
                    annot["annotation"]["additionalAttributes"]["HGMD.2015.4"]["attribute"]["MUT"]
            else:
                hash_variant["HGMD_version"] = "HGMD.2015.4"
                hash_variant["HGMD_ID"] = '.'
                hash_variant["HGMD_PHEN"] = '.'
                hash_variant["HGMD_CLASS"] = '.'
                hash_variant["HGMD_MUT"] = '.'

        hash_table[(chr, start, ref, alt)] = hash_variant

# Writing down hash_table into a file
    print_tables(hash_table, outFile)

