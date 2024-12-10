import re

# define a function to check if the given accession number is valid
def isValidAcc(acc: str):
    return bool(re.match(r'^(SR|ER|DR)', acc))