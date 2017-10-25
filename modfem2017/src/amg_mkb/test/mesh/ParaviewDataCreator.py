from sys import argv

def handleLevel(coordinates, fileName, line):
    dataFile = open(fileName,'w+')
    dataFile.write('x coord, y coord, z coord\n')
    numbers = line.split()
    for number in numbers:
	dataFile.write(",".join(coordinates[number]) + '\n')
    dataFile.close()
	    

script, filename, workDir = argv

f = open(filename)

coordinates = {}
for line in f:
    if(line.startswith('LEVEL')):
	break
    else:
	numbers = line.split()
	coordinates[numbers[0]] = (numbers[1],numbers[2],numbers[3])

dataFile = open(workDir + '/data' + '0.csv','w+')
dataFile.write('x coord, y coord, z coord\n')
for point in coordinates:
    dataFile.write(",".join(coordinates[point]) + '\n')

fileCnt = 1
for line in f:
    if(line.startswith('LEVEL')):
	fileCnt+=1
    else:
	handleLevel(coordinates,workDir + '/data' + str(fileCnt) + '.csv', line)

print 'END'
