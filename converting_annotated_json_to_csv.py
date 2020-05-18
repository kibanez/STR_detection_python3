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

    l_fields = ['rsid', 'chr', 'ref', 'alt', '1kG_phase3_ALL']

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


# Read input JSON
annotFile = sys.argv[1]
outFile = open(sys.argv[2], 'w')

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
        rsid = str(annot['id'])
        chr = str(annot['chromosome'])
        ref = str(annot['reference'])
        alt = str(annot['alternate'])

        # Extract CellBase annotation for GNOMAD_EXOMES-ALL
        if "populationFrequencies" in annot["annotation"] and len(annot["annotation"]["populationFrequencies"]) > 0:
            for pf in annot["annotation"]["populationFrequencies"]:
                if pf["study"] == "1kG_phase3" and pf["population"] == "ALL":
                    hash_variant['1kG_phase3_ALL'] = str(pf["altAlleleFreq"])
                else:
                    noStudy = "NO_STUDY"
            if hash_variant['1kG_phase3_ALL'] == "" and noStudy != "":
                hash_variant['1kG_phase3_ALL'] = noStudy
#        else:
#            outLine += "\t" + "NOT_ANNOTATED"
        hash_table[rsid, chr, ref, alt] = hash_variant

# Writing down hash_table into a file
    print_tables(hash_table, outFile)

