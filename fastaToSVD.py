from Bio import SeqIO
input_file = 'SuperRead/superReadSequences.fasta'

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
	
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
superReads = []
i = 0
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
        superReads.append([i, sequence])
	i += 1
#print superReads
print 'Number of reads : ' + str(len(superReads))
print('Going to prepare matrix')


count = len(superReads[0][1]) - 5000
print 'Kmer size is: '+str(count);

kmerToReadMap = prepareMatrix(superReads,len(superReads[0][1]) - 3000)

raw_input('Going to generate Matrix')
#print kmerToReadMap
print('Going to generate matrix')
matrix = generateMatrix(kmerToReadMap,superReads)
f = open('Matrix.txt','w')
f.write(matrix)

