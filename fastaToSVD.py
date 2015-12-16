from Bio import SeqIO
from datetime import datetime
import numpy as np
from copy import copy, deepcopy
from scipy import spatial

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


def computeSVD(X, featureCount):

	print "computing SVD"
	P, D, QT = np.linalg.svd(X, full_matrices=False)
	P1 = copy(P[:, 0:featureCount])
	D1 = copy(D[0:featureCount])
	QT1 = copy(QT[0:featureCount, :])
	superReadVectorMatrix = np.transpose(QT1)
	return superReadVectorMatrix

def computeAdjacencyList(vectorMatrix):

	print "computing adjacenylist"

	adjacencyList = {}
	for i in range(0, len(vectorMatrix)):
		if i%100 is 0:
			print str(i) + ' k-mers have been processed'
		for j in range(0, len(vectorMatrix)):
			if(i != j):
				result = 1 - spatial.distance.cosine(vectorMatrix[i], vectorMatrix[j])
				#print result
				if(result > 0.6):
					if i in adjacencyList:
						value = adjacencyList.pop(i)
						value.append(j)
						adjacencyList[i] = value
					else:
						adjacencyList[i] = [j]

	return adjacencyList

def clusterReads(vectorMatrix, featureCount):
	cluster = {}
	for i in range(0, len(vectorMatrix)):
		if i%100 is 0:
			print str(i) + ' reads have been clustered'
		max = -1
		bucketIndex = -1
		for j in range(0, featureCount):
			bucket = np.zeros((1, featureCount))
			bucket[0][j] = 1
			result = 1 - spatial.distance.cosine(bucket, vectorMatrix[i])
			if(result > max):
				max = result
				bucketIndex = j
		if bucketIndex in cluster:
			cluster[bucketIndex] = cluster[bucketIndex] + [i]
		else:
			cluster[bucketIndex] = [i]
		#print "bucket of " + str(i) + "  is :" + str(bucketIndex)
	return cluster




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


#print lmerToReadMap
print "matrix size will be: " + str(len(lmerToReadMap.keys()) * len(superReads))


print('Going to generate matrix')

matrix = generateMatrix(lmerToReadMap,superReads)


featureCount = 50
superReadVectorMatrix = computeSVD(matrix, featureCount)

cluster = clusterReads(superReadVectorMatrix, featureCount)
#print cluster
filename = "output-" + str(featureCount) + ".txt"
target = open(filename, 'w')
for i in cluster:
  target.write("\n >>> \n %s " % cluster[i])
target.close()
#print superReadVectorMatrix
#adjacencyList = computeAdjacencyList(superReadVectorMatrix)
#print QT1
#print vectorMatrix

#print adjacencyList
