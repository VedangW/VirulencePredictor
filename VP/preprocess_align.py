#!usr/bin/python

import os
import click
import argparse
import numpy as np
import pandas as pd

from config import AMINO_ACIDS
from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.sequence import SequenceSet, subset, columns

# File IO Utility functions

def get_fname(fname):
	""" Get a new file name for any file
		by just adding '_prep' to it.
	"""

	fname = fname.split('.')
	fname = fname[0]
	fname += '_prep.fasta'
	return fname

def save_as_fasta(df, fname='prep.fasta'):
	""" Save a series of sequences as 
		a fasta file.
	"""
	
	# Gets the sequences in fasta format
	seqs = []
	for index, row in df.iterrows():
		seq = '>' + index + '\n' + row.str.cat(sep='') + '\n'
		seqs.append(seq)

	# Write to file
	with open(fname, 'w') as f:
		for seq in seqs:
			f.write(seq)

# Preprocessing functions

def check_anomaly(df):
	""" Checks if there is an anomaly
		anywhere in the sequences.
	"""

	to_change = {}
	for col in df.columns:
		# Check if any anomaly characters exist
		diff = set(df[col].unique()) - AMINO_ACIDS
		if diff:
			to_change[col] = list(diff)
			
	return to_change

def replace_with_X(df, change_cols):
	""" Replaces any anomaly present with
		an X in the dataframe.
		
		Changes are to be calculated in
		check_anomaly().
	"""

	for k in change_cols.keys():
		to_change = change_cols[k]
		for anomaly in to_change:
			df[k].replace(anomaly, 'X', inplace=True)

	return df

def replace_X_with_mode(df):
	""" Replaces 'X' in each column with the
		mode of that column.
		
		If 'X' is the mode of some column, throws
		an error. Also if 'X' is present in some 
		column after preprocessing, throws an error.
	"""
	
	# Check which columns have anomalies
	cols = check_anomaly(df).keys()
	
	# Replace 'X' with mode for each column
	for col in cols:
		mode = df[col].mode()[0]
		if mode == 'X':
			# If mode of column is 'X' then drop it.
			df.drop([col], axis=1, inplace=True)
		else:
			df[col].replace('X', mode, inplace=True)
			
	# Check that no anomalies are left
	try:
		assert not check_anomaly(df)
	except:
		err = "Anomalies still present in df after preprocessing.\n" + \
				"anomalies = " + str(check_anomaly(df))
		raise ValueError(err)
		
	return df

def get_df_from_file(fname):
	""" Loads data from a file to a dataframe.
	
		The identifiers are set as a primary
		key of the dataframe.
		
		The sequence entries are a1, a2,... and so on.
	"""

	# Load the file
	f = load_fasta_file(fname)
	
	# Get the identifiers and the sequences
	names, dataset = [], []
	for i in range(len(f)):
		dataset.append(f[i].data)
		names.append(f[i].identifier)
	
	# Generate a header for the dataframe
	headers = ['a' + str(i + 1) for i in range(np.shape(dataset)[1])]

	# Generate dataframe
	df = pd.DataFrame(dataset, columns=headers)
	df['names'] = names
	df = df.set_index('names')
	
	return df

def preprocess_alignments(dirname, ow, particle):
	""" Function to preprocess all files in a directory.
	"""

	# Create directories if not already existing
	if particle == 'virus': 
		# Check directory
		if not os.path.exists('data/virus'):
			# If doesn't exist, make a new one
			os.makedirs('data/virus')
		elif ow == True:
			# Overwrite all pre-existing files
			pass
		else:
			# Warn the user and stop the program
			raise ValueError("Directory 'virus' already exists.")
	else:
		# Check directory
		if not os.path.exists('data/mouse'):
			# If doesn't exist, make a new one
			os.makedirs('data/mouse')
		elif ow == True:
			# Overwrite all pre-existing files
			pass
		else:
			# Warn the user and stop the program
			raise ValueError("Directory 'mouse' already exists.")

	# Call on functions to preprocess the files
	print "Preprocessing files..."
	with click.progressbar(os.listdir(dirname)) as bar:
		for fname in bar:
			df = get_df_from_file(dirname + '/' + fname)
			anomalies = check_anomaly(df)
			if anomalies:
				df = replace_with_X(df, anomalies)
				df = replace_X_with_mode(df)

			# Save file as fasta
			fname = get_fname(fname)
			if particle == 'virus':
				save_as_fasta(df, 'data/virus/' + fname)
			else:
				save_as_fasta(df, 'data/mouse/' + fname)