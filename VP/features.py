#!usr/bin/python

import os
import sys
import click
import pickle
import argparse

from quantiprot.utils.mapping import simplify
from quantiprot.utils.io import load_fasta_file
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.utils.sequence import SequenceSet, subset, columns

def save_to_disc(fname, item, **kwargs):
	""" A function to save any object to disc
		using pickle.

		Parameters
		----------
		fname: str
			Name of the file to save to.
		
		item: any object
			The item which needs to be saved
	"""

	# Open pickle file and dump
	out = open(fname, 'w')
	pickle.dump(item, out)
	out.close()

def write_to_log_file(log, fname='features.log'):
	""" A function to write errors to a log file.

		Parameters
		----------
		log: list
			The log of errors.

		fname: str
			Name of the log file to write to.
	"""


	# Open log file
	f = open(fname, 'w')

	# Write to log file
	for a in log:
		a += '\n'
		f.write(a)

	f.close()

def encode_sequences_from_file(fname, dirname):
	""" Function to get the encodings for all
		sequences from a fasta file which is aligned.

		This encoding method uses AAIndex to get the
		encodings. The amino acid index used for this
		currently is 'JOND920101' - Relative frequency 
		of occurrence (Jones et al., 1992).

		Parameters
		----------
		fname: str
			Name of the file from which to read sequences.

		dirname: str
			Name of the file in which 'fname' is present.

		Returns
		-------
		encoded: list
			A list of the encoded sequences.

		log: list
			A log containing any errors encountered.	
	"""

	# Define the path to the file
	fpath = dirname + '/' + fname

	# Read file as a set of Quantiprot sequences
	sequences = load_fasta_file(fpath)

	# Create a Feature object
	aa2freq_map = get_aaindex_file("JOND920101")
	aa2freq_map.mapping['-'] = 0.0
	aa2freq_map.mapping['X'] = 0.0
	freq_feat = Feature(aa2freq_map)
	
	# Encode sequences using relative frequency
	encoded, log = [], []
	for seq in sequences:
		try:
			f = freq_feat(seq)
			encoded.append(f.data)
		except:
			# For logging errors
			e, message, _tb = sys.exc_info()
			log.append(fname + ': line ' + str(e) + ': ' + str(message))
			continue

	return encoded, log

def pad_sequences(seqs):
	lengths = [len(x) for x in seqs]
	m = max(lengths)

	for i in range(len(seqs)):
		while (len(seqs[i]) < max):
			seqs[i].append(0.)

	return seqs

def generate_sequences_from_dir(dir):

	files = os.listdir(dir)

	# Encode all sequences in all files
	# in the given directory
	full_log, enc_seqs = [], []
	for fname in files:
		enc, log = encode_sequences_from_file(fname)
		enc_seqs.append(enc)
		full_log += log

	# Log if log is not empty
	if full_log:
		write_to_log_file(full_log)

	# Pad the sequences
	seqs = pad_sequences(enc_seqs)

	save_to_disc(seqs)

def main():
	get_features()

if __name__ == "__main__":
	main()