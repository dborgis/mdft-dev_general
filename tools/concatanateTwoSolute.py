import sys


def lineToTab(line):
	line = line.strip()
	while("  " in line):
		line = line.replace("  ", " ")
	return line.split()


if len(sys.argv)<=2:
	sys.exit( "ERROR: Missing parameters\n\nusage: python " + sys.argv[0] + " firstSoluteFile secondSoluteFile concatanateSoluteFile\n" )
	


fileNameOne = sys.argv[1]
fileNameTwo = sys.argv[2]
outFileName = sys.argv[3]


with open(fileNameOne, 'r') as fileInOne:
	linesInFirstFile = fileInOne.readlines()
with open(fileNameTwo, 'r') as fileInTwo:
	linesInSecondFile = fileInTwo.readlines()
	
nAtomsInFirstFile = int(linesInFirstFile[1].strip().split()[0])
nAtomsInSecondFile = int(linesInSecondFile[1].strip().split()[0])

with open(outFileName, 'w') as fileOut:
	
	fileOut.write(linesInFirstFile[0])  # comment line
	fileOut.write(str(nAtomsInFirstFile+nAtomsInSecondFile)+"\n")   # number of particules
	fileOut.write(linesInFirstFile[2])   # headers
	
	for line in linesInFirstFile[3:]:
		fileOut.write(line)

	for line in linesInSecondFile[3:]:
		tab = lineToTab(line)
		tab[0] = str(int(tab[0])+nAtomsInFirstFile)
		fileOut.write(" ".join(tab)+"\n")


