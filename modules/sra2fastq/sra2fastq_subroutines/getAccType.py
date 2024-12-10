import re

# define a function to get the SRA/ERA/DRA type from the accession number
def getAccType(acc):
    if re.match(r'^(SRX|ERX|DRX)', acc):
        return "ByExp"
    elif re.match(r'^(SRS|ERS|DRS)', acc):
        return "BySample"
    elif re.match(r'^(SRP|ERP|DRP)', acc):
        return "ByStudy"
    else:
        return "ByRun"
