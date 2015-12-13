from Bio import SeqIO
from datetime import datetime
input_file = 'SuperRead/superReadSequences.fasta'

# This function constructs the kmer-super read matrix.
def generateMatrix(kmerToReadMap,superReads):
	matrix = []
	sortedKeys = kmerToReadMap.keys()
	print ('Number of Keys: ' + str(len(sortedKeys)))
	#sortedKeys.sort()
	i = 0
	for key in sortedKeys:
		if i%100 is 0:
			print str(i) + ' keys processed'
		row = [0] * len(superReads)
		for value in kmerToReadMap[key]:
			row[value] = 1
		matrix.append(row)
		i += 1
	return matrix

# This function constructs unique-kmer(key) - read(value) pair from superreads.
def prepareMatrix(superReads, kmerSize):

	kmerToReadMap = {}
	i = 0
	for read in superReads:
		if i%10000 is 0:
			print str(i) + ' reads have been processed'
		rowNum = read[0]
		temp = read[1]
		for k in range(0,len(read[1])+1-kmerSize):
			key = temp[k:k+kmerSize]
			if key in kmerToReadMap:
				kmerToReadMap[key] = kmerToReadMap[key] + [rowNum]
			else:
				kmerToReadMap[key] = [rowNum]
			#try:
			#	value = kmerToReadMap.pop(key)
			#	value.append(rowNum)
			#except KeyError:
			#	value = [rowNum]
			#kmerToReadMap[key] = value
		i += 1

	return kmerToReadMap


def mapKmerToLmer(kmerToReadMap, lmerSize):
	print "inside mapKmerToLmer"
	lmerToReadMap = {}
	i = 0
	for key in kmerToReadMap.keys():
		if i%100000 is 0:
			print str(i) + ' k-mers have been processed'
		flag = False
		for k in range(0, len(key) + 1 - lmerSize):
			temp = key[k : k+lmerSize]
			if flag:
				if(lmer > temp):
					lmer = temp
			else:
				flag = True
				lmer = temp

		if lmer in lmerToReadMap:
			lmerToReadMap[lmer] = lmerToReadMap[lmer] + kmerToReadMap[key]
			#value = lmerToReadMap.pop(lmer)
			#value.append(kmerToReadMap[key])
			#lmerToReadMap[lmer] = value
		else:
			lmerToReadMap[lmer] = kmerToReadMap[key]
		i += 1
	return lmerToReadMap


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

print "unique kmer count: " + str(len(kmerToReadMap.keys()))


#Sort the kmer-read pairs so that we consider only the most occuring kmers.
#sortKmerReadMap(kmerToReadMap)

lmerToReadMap = mapKmerToLmer(kmerToReadMap, 7)
print "First unique lmer count: " + str(len(lmerToReadMap.keys()))

#lmerToReadMap = mapKmerToLmer(lmerToReadMap, (kmerSize)/4)
#print "second unique lmer count: " + str(len(lmerToReadMap.keys()))

#print lmerToReadMap
print "matrix size will be: " + str(len(lmerToReadMap.keys()) * len(superReads))


print('Going to generate matrix')
#print("samle key:"  + str(lmerToReadMap.keys()[0]) + "\t sample value:" + str(lmerToReadMap.values()[0]))
matrix = generateMatrix(lmerToReadMap,superReads)

#f = open('Matrix.txt','w')
#f.write(matrix)
