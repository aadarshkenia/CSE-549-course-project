from Bio import SeqIO
from datetime import datetime
input_file = 'SuperRead/superReadSequences.fasta'

# This function constructs the kmer-read matrix.
def generateMatrix(kmerToReadMap,superReads):
	matrix = []
	sortedKeys = kmerToReadMap.keys()
	print ('Number of Keys: ' + str(len(sortedKeys)))
	#sortedKeys.sort()
	for key in sortedKeys:
		row = [0] * len(superReads)
		for value in kmerToReadMap[key]:
			row[value] = 1
		matrix.append(row)
	return matrix

# This function constructs unique-kmer(key) - read(value) pair from superreads.
def prepareMatrix(superReads, kmerSize=3):

	kmerToReadMap = {}
	i = 0
	for read in superReads:
		if i%10000 is 0:
			print str(i) + ' reads have been processed'
		rowNum = read[0]
		temp = read[1]
		for k in range(0,len(read[1])+1-kmerSize):
			key = temp[k:k+kmerSize]
			try:
				value = kmerToReadMap.pop(key)
				value.append(rowNum)
			except KeyError:
				value = [rowNum]
			kmerToReadMap[key] = value
		i += 1

	return kmerToReadMap


def sortKmerReadMap(kmerToReadMap):

	kmerToReadList = [];
	print 'entering sortKmerReadMap function' + str(datetime.now())
	for key in kmerToReadMap.keys():
		kmerToReadList.append([key, kmerToReadMap[key], len(kmerToReadMap[key])])
	kmerToReadListSortedDes= sorted(kmerToReadList, key=lambda x: x[2], reverse=True)
	print "sortKmerRead list created" + str(datetime.now())
	newKmerToReadMap = {}
	for i in range(0, 10):
		print str(kmerToReadListSortedDes[i][2])

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
superReads = []
kmerSize = 101
print 'Kmer size:' + str(kmerSize)
i = 0
for fasta in fasta_sequences:
	if(len(str(fasta.seq)) > kmerSize):
		name, sequence = fasta.id, str(fasta.seq)
		superReads.append([i, sequence])
		i = i + 1;
#print superReads
print 'Number of reads : ' + str(len(superReads))
print('Going to prepare matrix')

count = len(superReads[0][1])


kmerToReadMap = prepareMatrix(superReads,kmerSize)

print 'going to generate Matrix' + str(datetime.now())

#Sort the kmer-read pairs so that we consider only the most occuring kmers.

sortKmerReadMap(kmerToReadMap)

'''
#print kmerToReadMap
print('Going to generate matrix')
matrix = generateMatrix(kmerToReadMap,superReads)
f = open('Matrix.txt','w')
f.write(matrix)
'''
