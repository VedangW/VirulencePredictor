#!usr/bin/python

import os
import sys
import click
import pickle
import argparse
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from config import DATA_REPO
from keras.preprocessing.sequence import pad_sequences

from quantiprot.utils.io import load_fasta_file
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.utils.sequence import SequenceSet

# Retrieve ids set
with open(DATA_REPO + '/orders') as f:
	virus_order, mouse_order, inv_virus_order, inv_mouse_order = pickle.load(f)
	ids_set_virus = set(virus_order.keys())
	ids_set_mouse = set(mouse_order.keys())

# Retrive viruses_dict
with open(DATA_REPO + '/viruses_dict') as f:
	viruses_dict = pickle.load(f)

def _change_format_virus(iden):
	""" Remove the Seq-ID from the identifier
		of any sequence.

		For more info about Seq-ID, refer to 
		virus_doc.md in VirulencePredictor/docs.
	"""
	iden = iden.split('|')[1]
	iden = iden.split('_')
	strain = iden[0]

	return viruses_dict[strain]

def _change_format_mouse(iden):
	""" Remove the Seq-ID from the identifier
		of any sequence.
	"""

	iden = iden.split('_')

	if len(iden) == 3:
		iden.pop()
	elif len(iden) == 2:
		l = list(iden[1])
		for i in range(5):
			l.pop()
		iden[1] = ''.join(l)

	iden.pop(0)
	iden = iden[0]

	return iden

def _get_feature_map(index='JOND920101'):
	""" To get the feature mapping object 
		using the amino acid index given. 

		The mapping is created using AAindex.
		'-' is mapped to 0.0.
	"""
	
	# Create a Feature object
	aaindex_map = get_aaindex_file(index)
	aaindex_map.mapping['-'] = 0.0
	feat_map = Feature(aaindex_map)
	
	return feat_map

def _pad_encoding(enc, pad_len):
	""" A function to pad all the values in a 
		dictionary.
	"""
	for key in enc.keys():
		val = enc[key]
		val = np.reshape(val, (1, len(val)))
		val = pad_sequences(val, 
							maxlen=pad_len, 
							dtype='float32', 
							padding='pre', 
							truncating='pre', 
							value=0.0)
		enc[key] = val[0]
		
	return enc

def encoded_seq_from_file(fname, dirname, particle, index):
	""" Function to encode the sequences from
		a file using AAindex. The encoded sequences
		are then padded to maximum length.
	"""

	# Load the fasta file
	f = load_fasta_file(dirname + '/' + fname)
	feat_map = _get_feature_map(index)
	
	# Get the sequences in a dataset
	dataset = []
	for i in range(len(f)):
		dataset.append(f[i])
		
	# Create a dictionary with keys as identifiers
	# and their values as the data.
	enc = {}
	if particle == 'virus':
		for seq in dataset:
			seq_id = _change_format_virus(seq.identifier)
			enc[seq_id] = feat_map(seq).data
	elif particle == 'mouse':
		for seq in dataset:
			seq_id = _change_format_mouse(seq.identifier)
			if seq_id not in ids_set_mouse:
				print seq.identifier, seq_id
			enc[seq_id] = feat_map(seq).data
	
	# Pad all sequences to maximum value in the
	# dataset.
	maxlen = max([len(val) for val in enc.values()])
	enc = _pad_encoding(enc, maxlen)

	# Check if all values have lengths
	# equal to maxlen.
	for val in enc.values():
		assert len(val) == maxlen
	
	return enc

def pad_dict(enc, max_len, particle):
	""" Pad a dictionary so that it resembles a
		complete dictionary with all the ids present
		for the viruses and then pad those values to 
		the max_len.
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

def recombine(encs, particle):
	""" A function to recombine the sequences
		in the segments with the original viruses.
		This function uses the identifiers to relocate
		the sequences to their original viruses.
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

	print len(features_dict)

	# Add the correct sequence to its parent virus
	for enc in encs:
		for k in enc.keys():
			try:
				features_dict[k].append(enc[k])
			except Exception as e:
				print "Couldn't preprocess", e.message

	# Check if initial shape is maintained
	for key in features_dict.keys():
		try:
			assert len(features_dict[key]) == total_segments
		except:
			print key
			print len(features_dict[key]), total_segments

	return features_dict

def extract_features(particle, index):
	""" Function to extract features from the 
		preprocessed files.
	"""

	# Get encodings from AAindex
	encs = []
	print "Encoding using AAindex..."

	if particle == 'virus':
		read_dir = 'data/virus'
	else:
		read_dir = 'data/mouse'

	with click.progressbar(os.listdir(read_dir)) as bar:
		for fname in bar:
			encs.append(encoded_seq_from_file(fname, 
				read_dir, particle, index))

	# Length of the longest sequence
	global_max_len = max([len(enc.values()[0]) for enc in encs])

	# Pad all encodings
	print "Padding encodings..."
	with click.progressbar(range(len(encs))) as bar:
		for i in bar:
			encs[i] = pad_dict(encs[i], 
							global_max_len, 
							particle)

	# Recombine
	print "Recombining features..."
	features = recombine(encs, particle)

	return features