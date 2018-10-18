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
	ids_set_virus = pickle.load(f)

with open('data/mouse_ids.pkl') as f:
	ids_set_mouse = pickle.load(f)

# DONE
def remove_seq_id_virus(iden):
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

# DONE
def remove_seq_id_mouse(iden):
	""" Remove the Seq-ID from the identifier
		of any sequence.

		Parameters
		----------
		iden: str
			The identifier to remove Seq-ID from.

		Returns
		-------
		iden: str
			Identifier without Seq-ID in it.
	"""

	iden = iden.split('_')
	if len(iden) == 3:
		iden.pop()
	elif len(iden) == 2:
		l = list(iden[1])
		for i in range(5):
			l.pop()
		iden[1] = ''.join(l)
	iden = '_'.join(iden)

	return iden

# DONE
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

# DONE
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

# DONE
def encoded_seq_from_file(fname, dirname, particle):
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
	if particle == 'virus':
		for seq in dataset:
			seq_id = remove_seq_id_virus(seq.identifier)
			enc[seq_id] = feat_map(seq).data
	elif particle == 'mouse':
		for seq in dataset:
			seq_id = remove_seq_id_mouse(seq.identifier)
			if seq_id not in ids_set_mouse:
				print seq.identifier
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


# DONE
def pad_dict(enc, max_len, particle):
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

		particle: str
			'virus' or 'mouse'.

		Returns
		-------
		enc: dict
			A padded dictionary of encodings
	"""
	
	# Add all keys which are not already present in 
	# the dictionary from the set of ids for the viruses.
	fill_len = len(enc.values()[0])

	fill_dict = {}
	if particle == 'virus':
		for k in ids_set_virus - set(enc.keys()):
			fill_dict[k] = np.zeros((fill_len,))
	elif particle == 'mouse':
		for k in ids_set_mouse - set(enc.keys()):
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

# DONE
def recombine(encs, particle):
	""" A function to recombine the sequences
		in the segments with the original viruses.
		This function uses the identifiers to relocate
		the sequences to their original viruses.
		
		Parameters
		----------
		encs: list of dicts
			A list containing the encodings for all the 
			different segments.

		particle: str
			'virus' or 'mouse'.
			
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
	if particle == 'virus':
		for k in ids_set_virus:
			features_dict[k] = []
	elif particle == 'mouse':
		for k in ids_set_mouse:
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

# DONE
def test_enc_shape(encs, 
					total_segments, 
					total_entities,
					global_max_len,
					particle):

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

		particle: str
			'virus' or 'mouse'.
	"""
	return 

# DONE
def extract_features(particle):
	""" Function to extract features from the 
		preprocessed files.

		Parameters
		----------
		particle: str
			Can be 'virus' or 'mouse' depending on what
			we want to extract features from.
	"""

	# Get encodings from AAindex
	encs = []
	print "Encoding using AAindex..."

	if particle == 'virus':
		read_dir = 'data/virus'
	else:
		read_dir = 'data/mouse/dir_002'

	with click.progressbar(os.listdir(read_dir)) as bar:
		for fname in bar:
			encs.append(encoded_seq_from_file(fname, read_dir, particle))

	# Length of the longest sequence
	global_max_len = max([len(enc.values()[0]) for enc in encs])

	# Pad all encodings
	print "Padding encodings..."
	with click.progressbar(range(len(encs))) as bar:
		for i in bar:
			encs[i] = pad_dict(encs[i], global_max_len, particle)

	# Test that the shape of the encodings is uniform
	if particle == 'virus':
		test_enc_shape(encs, 13, 215, global_max_len, particle)
	elif particle == 'mouse':
		test_enc_shape(encs, 1000, 12, global_max_len, particle)

	# Recombine
	print "Recombining features..."
	features = recombine(encs, particle)

	print "Done."

	return features