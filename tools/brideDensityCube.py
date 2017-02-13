import sys


def lineToTab(line):
	line=line.strip()
	while "  " in line:
		line = line.replace("  ", " ")
	return line.split(" ")

fileNameIn = sys.argv[1]
fileNameOut = fileNameIn.replace(".cube", "_bride.cube")

with open(fileNameIn, 'r') as fileIn:
	with open(fileNameOut, 'w') as fileOut:
		for line in fileIn:
			tab=lineToTab(line)
			if(len(tab)!=1):
				fileOut.write(line)
			else:
				value = float(tab[0])
				value = min(value, 2)
				value = max(value, 0)
				fileOut.write(str(value)+"\n")
