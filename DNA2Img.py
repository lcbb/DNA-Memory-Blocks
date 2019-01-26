import sys
import math
from PIL import Image

if len(sys.argv) != 4:
  print("Usage: python3 DNA2Img.py <seq> <imgwidth> <imgheight>")
  sys.exit(0)

seq = sys.argv[1]
width = int(sys.argv[2])
height = int(sys.argv[3])

# Relevant constants (must match those used in Img2DNA.py)
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


# Convert dimers to digits
dimer_to_digit = {
  'GA': 0, 'GT': 1, 'GC': 2, 'GG': -1,
  'AT': 0, 'AC': 1, 'AG': 2, 'AA': -1,
  'TC': 0, 'TG': 1, 'TA': 2, 'TT': -1,
  'CG': 0, 'CA': 1, 'CT': 2, 'CC': -1
}

# Get primer sequences
master_forward='AAATTTGAATTCGTCGTCGTCCCCTCAAACT'
master_reverse='GCTGAAAAGGTGGCATCAATCTGCAGGCAAAC'
#assert(seq[:len(master_forward)] == master_forward)
#assert(seq[-len(master_reverse):] == master_reverse)

# Find location of the master primers, which may be truncated to no more than 10nt in length
master_forward_trunc = master_forward[-10:]
master_reverse_trunc = master_reverse[:10]

assert(master_forward_trunc in seq)
assert(master_reverse_trunc in seq)
master_forward_pos = seq.index(master_forward_trunc)
master_reverse_pos = seq.rindex(master_reverse_trunc)

print("Master forward primer found at position", master_forward_pos)
print("Master reverse primer found at position", master_reverse_pos)

# Extract relevant portion of the sequence
master_forward_seq = seq[:master_forward_pos + len(master_forward_trunc)]
master_reverse_seq = seq[master_reverse_pos:]
hash_forward_seq = 'G' + seq[len(master_forward_seq):len(master_forward_seq)+24]
hash_reverse_seq = 'G' + seq[-24+master_reverse_pos:master_reverse_pos]
prefix_len = len(master_forward_seq) + len(hash_forward_seq) - 1
suffix_len = len(master_reverse_seq) + len(hash_reverse_seq) - 1
seedseq = 'G' + seq[prefix_len:prefix_len+SEED_NDIGITS]
imgseq = seq[prefix_len+SEED_NDIGITS:-suffix_len]

print(master_forward_seq)
print(master_reverse_seq)
print(hash_forward_seq)
print(hash_reverse_seq)

# Interpret seed nucleotides
seed_digits = [dimer_to_digit[seedseq[i:i+2]] for i in range(len(seedseq)-1)]
seed = sum([3**(SEED_NDIGITS-i-1) * digit for i,digit in enumerate(seed_digits)])
rng = make_rng(RNG_MULTIPLIER, RNG_OFFSET, RNG_MODULO, seed=seed)

#print(seed)

# Convert the nucleotide sequence to quaternary digits
img_digits = [(NUCLEOTIDES.index(nt)+4-rng())%4 for nt in imgseq]
#img_digits = [dimer_to_digit[imgseq[i:i+2]] for i in range(len(imgseq)-1)]
#img_digits = [(d+8-rng())%4 for i,d in enumerate(img_digits)] # account for random offsets

# Decode the run-length encoding
imgdata = []
idx = 1
cur_color = img_digits[0]
while idx < len(img_digits):
  quaternary_str = ''.join([str(v) for v in img_digits[idx:idx+RUN_LENGTH_NDIGITS]])
  try:  run_length = int(quaternary_str, 4)
  except:  run_length = 0

  if run_length == 0:
    print("WARNING: Sequencing error, found invalid run length")
    run_length = 1

  imgdata.extend([cur_color]*run_length)

  idx += RUN_LENGTH_NDIGITS

  if run_length == RUN_LENGTH_MAX and idx < len(img_digits):
    cur_color = img_digits[idx]
    idx += 1
  else:
    cur_color = 1-cur_color
imgdata = imgdata[:width*height]

# Place image data into image
img = Image.new(size = (width, height), mode = '1')
img.putdata(imgdata)
img = img.transpose(Image.TRANSPOSE)
img.show()

# Verify hashes
