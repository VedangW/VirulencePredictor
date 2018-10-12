#!usr/bin/python

import os
import click
import argparse
import numpy as np
import pandas as pd

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.sequence import SequenceSet, subset, columns

# The set of amino acid single letter notations
amino_acids = frozenset(['G', 'P', 'A', 'V', 'L', 
						 'I', 'M', 'C', 'F', 'Y', 
						 'W', 'H', 'K', 'R', 'Q', 
						 'N', 'E', 'D', 'S', 'T', 
						 '-'])

def get_fname_virus(fname):
	""" Get a new file name by removing
		the 'aligned' and adding the 'prep'.
	"""

	fname = fname.split('_')
	fname.pop()
	fname = ''.join(fname)
	fname += '_prep'

	return fname

def check_anomaly(df):
	""" Checks if there is an anomaly
		anywhere in the sequences.
	"""

	to_change = {}
	for col in df.columns:
		# Check if any anomaly characters exist
		diff = set(df[col].unique()) - amino_acids
		if diff:
			to_change[col] = list(diff)
			
	return to_change 

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
	f = open(fname, 'w')
	
	for seq in seqs:
		f.write(seq)

	f.close()

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
			# Throw error if 'X' is the mode of the column
			err = "Mode of column " + col + " is 'X'" 
			raise ValueError(err)
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

def preprocess_alignments(dirname, ow):
	""" Function to preprocess all files in a directory.
	"""

	# Check directory
	if not os.path.exists('data/prep'):
		# If doesn't exist, make a new one
		os.makedirs('data/prep')
	elif ow == True:
		# Overwrite all pre-existing files
		pass
	else:
		# Warn the user and stop the program
		raise ValueError("Directory 'prep' already exists.")
		
	label = "Preprocessing files..."
	with click.progressbar(os.listdir(dirname),
		label=label) as bar:
		for fname in bar:
			df = get_df_from_file(dirname + '/' + fname)
			anomalies = check_anomaly(df)
			if anomalies:
				df = replace_with_X(df, anomalies)
				df = replace_X_with_mode(df)
			fname = get_fname_virus(fname)
			save_as_fasta(df, 'data/prep/' + fname)

	print ("Done.")

def main():
	# Parser arguments
	parser = argparse.ArgumentParser(
		description='Preprocesses the aligned segments.')

	parser.add_argument('-d', '--aligned_dir', nargs='?', default='aligned', 
		type=str, help='Name of directory containing alignments.')

	parser.add_argument('-ow', '--overwrite', nargs='?', default=True, 
		type=bool, help='Set True to overwrite.')

	args = parser.parse_args()
	preprocess_alignments(args.aligned_dir, args.overwrite)

if __name__ == "__main__":
	main()