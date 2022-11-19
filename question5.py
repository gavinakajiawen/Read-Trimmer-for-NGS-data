#!/usr/bin/env python3

nucleotides = "ATGACGTGCGCGTCGCTCTGA"
quality = "GHIJKKLLLOOOPPPQQQRRR"

def get_ASCII_number(letter, phred_version = 64):
	''' This function takes a letter and transforms it into an phred score'''
	answer = ord(letter) - phred_version
	
	return answer


ASCII_scores = []

for i in quality:
	number = get_ASCII_number(i, 64)
	ASCII_scores.append(number)

nuc_list=[]
for i in nucleotides:
	nuc_list.append(i.split())

#print(ASCII_scores)

enum = []
#for i in enumerate(nucleotides):
#	print(i, end = "")
	
# put them in a dictionary
moving_window = 3
threshold = 10
# make the key three letters

# values is quality

for i in range(len(nuc_list)):
	print(nucleotides[i:i+3])
	#print(nuc_list[i])
''' 
new_list =[] 
while value < 10:
	
	append(new_list)
	
	discard number 
'''
