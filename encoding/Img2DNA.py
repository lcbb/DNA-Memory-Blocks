#!/usr/bin/env python3
from sys import argv
import math
import random
from PIL import Image
import imagehash

# Relevant constants (must match those used in DNA2Img.py)
SEED_MIN = 2
SEED_MAX = 2**31 - 2
SEED_NDIGITS = math.ceil(math.log(SEED_MAX)/math.log(3))
RNG_MULTIPLIER = 5
RNG_OFFSET = 0
RNG_MODULO = 2**31-1
NUCLEOTIDES = 'GATC'
RUN_LENGTH_MAX = 15
RUN_LENGTH_NDIGITS = 2
def make_rng(multiplier, offset, modulo, seed = 1):
	x = seed
	def rng():
		nonlocal x
		x = (multiplier * x + offset) % modulo
		return x
	return rng


def Int2Ternary(val, num_digits):

	assert(0 <= val < 3**num_digits)

  	# Converts the number to a base 3 integer.
	num_base3 = ''
	for i in range(num_digits-1, -1, -1):
		digit_val = 3**i
		num_base3 += str(val // digit_val)
		val %= digit_val
	return num_base3
def Int2Quaternary(val, num_digits):

	assert(0 <= val < 4**num_digits)

	# Converts the number of a base 4 integer.
	num_base4 = ''
	for i in range(num_digits-1, -1, -1):
		digit_val = 4**i
		num_base4 += str(val // digit_val)
		val %= digit_val
	return num_base4

def Quaternary2Nt(quaternary, rng = None):
	seq = ''
	for digit in quaternary:
		rand = rng() if rng is not None else 0
		offset = (int(digit) + rand) % 4
		seq += NUCLEOTIDES[offset]
	return seq

def Ternary2Nt(ternary, prefix = '', rng = None):
	if prefix == '':
		seq = random.choice('ATCG')
	else:
		seq = prefix

	for digit in ternary:
		rand = rng() if rng is not None else 0
		offset = (NUCLEOTIDES.index(seq[-1]) + 1 + int(digit) + rand) % 4
		seq += NUCLEOTIDES[offset]

	return seq[len(prefix):]



def Int2Nt(num, num_digits, prefix = 'G'):
	""" Converts the given number to a sequence of nucleotides using the Huffman
	encoding scheme to encode a base-3 integer, in which each digit of the integer
	is stored by the transition graph described in the huffman_encoding dict.
	The number must satisfy 0 <= num < 3^num_digits. """

	ternary = Int2Ternary(num, num_digits)
	return Ternary2Nt(ternary, prefix)

def Color2Nt(color, prefix = ''):
	""" Note: not used in current encoding. """
	return 'ACTG'[2*color + random.randint(0,1)]

def Img2RLE(img):
	""" Returns the run-length encoding scheme for the given image,
	in the form of a list of 2-tuples (color, run-length). """
	imgdata = [v%2 for v in img.getdata()] # for some reason, BW images still have pixel vals to 255...
	idx = 0
	rledata = []
	while idx < len(imgdata):
		cur_color = imgdata[idx]
	
		try:  run_length = imgdata.index(1-cur_color, idx) - idx
		except:  run_length = len(imgdata) - idx
		run_length = min(run_length, RUN_LENGTH_MAX)
	
		rledata.append((cur_color, run_length))
		idx += run_length

	return rledata

def RLE2Nt(rledata):
	""" Returns an encoding of the RLE in nucleotides. """

	# Get random number generator
	seed = random.randint(SEED_MIN, SEED_MAX)
	rng = make_rng(RNG_MULTIPLIER, RNG_OFFSET, RNG_MODULO, seed = seed)

	# Compute quaternary representation of RLE
	quaternary = ''
	store_color = True
	for color, run_length in rledata:
	  if store_color:  quaternary += Int2Quaternary(color, 1)
	  quaternary += Int2Quaternary(run_length, RUN_LENGTH_NDIGITS)
	  store_color = run_length==RUN_LENGTH_MAX

	# Turn run-length encoding into nucleotide sequence
	imgseq = Int2Nt(seed, SEED_NDIGITS, prefix='G')
	imgseq += Quaternary2Nt(quaternary, rng)
	return seed, imgseq

def DNAImageHash(img):
	""" Returns the 'hash' of the input image (a PIL.Image object)
	as two DNA oligo sequences. These sequences are computed from 
	the most significant and least significant halves of the image's 'whash'. """

	whash = str(imagehash.whash(img))
	assert(len(whash) == 16)

	pri5=Int2Nt(int(whash[:8], 16), 24, prefix='G')
	pri3=Int2Nt(int(whash[8:], 16), 24, prefix='G')
	primers=[pri5,pri3]
	return(primers)

def Img2Nt(img):
	# New run-length encoding scheme:
	#   1. Start at first pixel.
	#   2. Store color of current pixel.
	#   3. Traversing in a column-first order, find run length and store x=min(run length, MAX).
	#   4. Advance x pixels.
	#   5. If x==MAX, goto 2. Otherwise, goto 3.
	width,height = img.size

	rledata = Img2RLE(img)
	seed, imgseq = RLE2Nt(rledata)
	
	# Get primer sequences
	master_forward='AAATTTGAATTCGTCGTCGTCCCCTCAAACT'
	master_reverse='GCTGAAAAGGTGGCATCAATCTGCAGGCAAAC'
	whash_forward, whash_reverse = DNAImageHash(img)
	
	# Get complete sequence
	finalSeq = master_forward + whash_forward + imgseq + whash_reverse + master_reverse

	return seed, finalSeq

def checkSeq(sequence):
	import re
	to_check = [
		('GAATTC', 1),
		('CTGCAG', 1),
		('GGGG', 0),
		('AAAAAA', 0),
		('TTTTTT', 0),
		('CCCC', 1)
	]
	for seq, count in to_check:
		res = len(re.findall('(?={})'.format(seq), sequence))
		if res != count:
			return False
	return True
	

fIn = argv[1]
img=Image.open(fIn).convert(mode='1').transpose(method=Image.TRANSPOSE)

max_tries = 5000

# Check sequence for common problems
seed, finalSeq = Img2Nt(img)
i = 1
while not checkSeq(finalSeq):
	seed, finalSeq = Img2Nt(img)
	i += 1
	if i >= max_tries:
		print('WARNING: Could not find a satisfactory sequence after {} attempts'.format(i))
		break

print('>' + fIn[:fIn.rfind('.')] + '_' + str(seed))
print(finalSeq)
print()
