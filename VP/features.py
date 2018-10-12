#!usr/bin/python

import os
import sys
import click
import pickle
import argparse
import numpy as np

from keras.preprocessing.sequence import pad_sequences

from quantiprot.utils.mapping import simplify
from quantiprot.utils.io import load_fasta_file
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.utils.sequence import SequenceSet, subset, columns

# Retrieve ids set
with open('data/virus_ids.pkl') as f:
	ids_set = pickle.load(f)

def remove_seq_id(iden):
	""" Remove the Seq-ID from the identifier
		of any sequence.

		For more info about Seq-ID, refer to 
		virus_doc.md in VirulencePredictor/docs.

		Parameters
		----------
		iden: str
			The identifier to remove Seq-ID from.

		Returns
		-------
		iden: str
			Identifier without Seq-ID in it.
	"""
	iden = iden.split('|')[1]
	iden = iden.split('_')
	iden.pop(1)
	iden.pop(2)
	iden = '_'.join(iden)
	
	return iden

def get_feature_map(index='JOND920101'):
	""" To get the feature mapping object 
		using the amino acid index given. 

		The mapping is created using AAindex.
		'-' is mapped to 0.0.

		Parameters
		----------
		index: str
			Index of the amino acid.

		Returns
		-------
		feat_map: Feature obj
			A feature object which can transform
			any sequence to a sequence of numbers.
	"""
	
	# Create a Feature object
	aaindex_map = get_aaindex_file(index)
	aaindex_map.mapping['-'] = 0.0
	feat_map = Feature(aaindex_map)
	
	return feat_map

def pad_encoding(enc, pad_len):
	""" A function to pad all the values in a 
		dictionary.

		Parameters
		----------
		enc: dict
			A dictionary which contains keys as 
			identifiers and values as numerical
			sequences.

		pad_len: int
			Length to which to pad or truncate to.

		Returns
		-------
		enc: dict
			The padded dictionary
	"""
	for key in enc.keys():
		val = enc[key]
		val = np.reshape(val, (1, len(val)))
		val = pad_sequences(val, maxlen=pad_len, 
							dtype='float32', padding='pre', truncating='pre', value=0.0)
		enc[key] = val[0]
		
	return enc

def encoded_seq_from_file(fname, dirname):
	""" Function to encode the sequences from
		a file using AAindex. The encoded sequences
		are then padded to maximum length.

		Parameters
		----------
		fname: str
			Name of file to read from.

		dirname: str
			Name of directory the file is present in.

		Returns
		-------
		enc: dict
			The encodings for all the ids.
	"""

	# Load the fasta file
	f = load_fasta_file(dirname + '/' + fname)
	feat_map = get_feature_map()
	
	# Get the sequences in a dataset
	dataset = []
	for i in range(len(f)):
		dataset.append(f[i])
		
	# Create a dictionary with keys as identifiers
	# and their values as the data.
	enc = {}
	for seq in dataset:
		seq_id = remove_seq_id(seq.identifier)
		enc[seq_id] = feat_map(seq).data
	
	# Pad all sequences to maximum value in the
	# dataset.
	maxlen = max([len(val) for val in enc.values()])
	enc = pad_encoding(enc, maxlen)

	# Check if all values have lengths
	# equal to maxlen.
	for val in enc.values():
		assert len(val) == maxlen
	
	return enc

def pad_dict(enc, max_len):
	""" Pad a dictionary so that it resembles a
		complete dictionary with all the ids present
		for the viruses and then pad those values to 
		the max_len.

		Parameters
		----------
		enc: dict
			Dictionary of identifiers: encodings

		max_len: int
			The length to which to pad or truncate the 
			encodings.

		Returns
		-------
		enc: dict
			A padded dictionary of encodings
	"""
	
	# Add all keys which are not already present in 
	# the dictionary from the set of ids for the viruses.
	fill_len = len(enc.values()[0])

	fill_dict = {}
	for k in ids_set - set(enc.keys()):
		fill_dict[k] = np.zeros((fill_len,))
		
	# Add the created dictionary to the pre-existing dictionary
	enc.update(fill_dict)
	
	# Pad all values in the dictionary to the maximum length.
	for k in enc.keys():
		val = enc[k]
		val = np.reshape(val, (1, len(val)))
		val = pad_sequences(val, maxlen=max_len,
						dtype='float32', padding='pre', truncating='pre', value=0.0)
		enc[k] = val[0]

	return enc

def recombine(encs):
	""" A function to recombine the sequences
		in the segments with the original viruses.
		This function uses the identifiers to relocate
		the sequences to their original viruses.
		
		Parameters
		----------
		encs: list of dicts
			A list containing the encodings for all the 
			different segments.
			
		Returns
		-------
		features_dict: dict
			A dictionary with keys as the identifiers and
			the values as the concatenated list of embeddings.
	"""
	# Some dimensions to check with later
	enc_dim = len(encs[0].values()[0])
	total_segments = len(encs)
	
	# Initialise the dictionary
	features_dict = {}
	for k in ids_set:
		features_dict[k] = []
		
	# Add the correct sequence to its parent virus
	for enc in encs:
		for k in enc.keys():
			features_dict[k].append(enc[k])
		
	# Check if initial shape is maintained
	for val in features_dict.values():
		assert len(val) == total_segments
			
	# Concatenate each value in the dict to form
	# one large array
	for k in features_dict.keys():
		features_dict[k] = np.concatenate(features_dict[k], axis=0)
		
	# Check if the reshaping happened properly
	for val in features_dict.values():
		assert len(val) == enc_dim * total_segments
			
	return features_dict

def test_enc_shape(encs, 
					total_segments, 
					total_entities, 
					global_max_len):

	""" A function to test if the encodings 
		are the correct shape or not.

		Parameters
		----------
		encs: list of dicts
			The encodings.

		total_segments: int
			Total number of segments (13 for Influenza).

		total_entities: int
			Total number of virus entities (including variations
			on the same one).

		global_max_len: int
			Length of the largest sequence.
	"""

	try:
		assert len(encs) == total_segments
		for enc in encs:
			assert len(enc) == total_entities
			for val in enc.values():
				assert len(val) == global_max_len
			keys = set(enc.keys())
			assert not ids_set - keys
	except:
		err = "Shape of encoding is not proper."
		raise ValueError(err)

	print "Shape of encodings is correct."

def extract_features():
	""" Function to extract features from the 
		preprocessed files.
	"""

	# Get encodings from AAindex
	print ""
	encs = []
	print "Encoding using AAindex..."
	with click.progressbar(os.listdir('data/prep')) as bar:
		for fname in bar:
			encs.append(encoded_seq_from_file(fname, 'data/prep'))

	# Length of the longest sequence
	global_max_len = max([len(enc.values()[0]) for enc in encs])

	# Pad all encodings
	print "Padding encodings..."
	with click.progressbar(range(len(encs))) as bar:
		for i in bar:
			encs[i] = pad_dict(encs[i], max_len=global_max_len)

	# Test that the shape of the encodings is uniform
	test_enc_shape(encs, 13, 215, 772)

	# Recombine
	print "Recombining features..."
	features = recombine(encs)

	# Save to disc
	print "Saving to disc..."
	f = open('data/features_long.pkl', 'w')
	pickle.dump(features, f)
	f.close()

	print "Done."

def main():
	extract_features()

if __name__ == "__main__":
	main()