import sys
import json

if len(sys.argv) < 2:
    sys.exit('ERROR: Usage: %s <json> <outFile>' % sys.argv[0])

# Read input JSON
annotFile = sys.argv[1]

# Open outfile
outFile = open(sys.argv[2], 'w')

# Parse JSON
with open(annotFile) as f:
    for line in f:
        annot = json.loads(line)
        cellbaseFrq = ""
        noStudy = ""
        outLine = annot["annotation"]["chromosome"] + "\t" + str(annot["annotation"]["start"]) + "\t" + \
                  annot["annotation"]["reference"] + "\t" + annot["annotation"]["alternate"]

        # Extract CellBase annotation for GNOMAD_EXOMES-ALL
        if "populationFrequencies" in annot["annotation"] and len(annot["annotation"]["populationFrequencies"]) > 0:
            for pf in annot["annotation"]["populationFrequencies"]:
                if pf["study"] == "1kG_phase3" and pf["population"] == "ALL":
                    cellbaseFrq = str(pf["altAlleleFreq"])
                else:
                    noStudy = "NO_STUDY"
            if cellbaseFrq == "" and noStudy != "":
                cellbaseFrq = noStudy
        else:
            outLine += "\t" + "NOT_ANNOTATED"
        outLine += "\t" + cellbaseFrq
        # # Extract annotated frequencies from GNOMAD
        # for stu in annot["studies"]:
        #     for fi in stu["files"]:
        #         if "AF" in fi["attributes"]:
        #             outLine += "\t" + fi["attributes"]["AF"]
        outFile.write(outLine + "\n")
outFile.close()

