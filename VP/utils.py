#!usr/bin/python

import pickle
import argparse
import numpy as np
import pandas as pd

def remove_zeros(X):
	""" Remove columns where all values are zeros.
		This function is not being used right now.
	"""
	X_t = X.T
	where_zeros = np.where(~X_t.any(axis=1))[0]
	X_t = np.delete(X_t, where_zeros, axis=0)
	return X_t.T

def pad_with_zero_vectors(X, threshold_len):
	""" Pad a 2D array with 1D zero arrays to
		achieve a threshold len.
		This function is not being used right now. 
	"""
	_n_rows, n_cols = np.shape(X)
	new_shape = ((threshold_len, n_cols))
	X_changed = np.zeros(new_shape)
	X_changed[:X.shape[0],:X.shape[1]] = X

	return X_changed

def change_bracket(x):
	""" A utility function to change brackets to slash
		Eg. 'A/Netherlands/602/2009(H1N1)' is changed
		to 'A/Netherlands/602/2009/H1N1'.

		Primarily used in change_format().

		Parameters
		----------
		x: str
			String to process.
	"""

	if x.endswith(')'):
		x = x.split('(')
		l = x[-1]
		x.pop()
		l = l.split(')')[0]
		x = '('.join(x)
		x += '/' + l
	return x

def change_format(x):
	""" A utitlity function to add '/1' to all viruses
		without more than one version and remove brackets
		around influenza version.
		Eg. 'A/Netherlands/602(H1N2)/3' is changed to
		'A/Netherlands/602/H1N2/3' while
		'A/Netherlands/602(H1N2)' is changed to
		'A/Netherlands/602/H1N2/1'

		Primarily used for virus_ids in feature extraction
		of influenza.

		Parameters
		---------
		x: str
			String to process.
	"""
	if x.endswith(')'):
		x += '/1'
	x = x.split('/')
	l = x[-1]
	x.pop()
	x = '/'.join(x)
	x = change_bracket(x)
	x += '/' + l
	return x

def create_ids(v, m):
	""" Generate and store the ids from the graph
		dataset. 

		Parameters
		----------
		v: bool
			If virus_ids are to be generated.

		m: bool
			If mouse_ids are to be generated.
	"""
	df = pd.read_csv('data/merged.csv')

	# Generate virus_ids
	if v:
		influenza = df['Influenza_virus_name'].unique()
		virus_ids = set(influenza)

		with open('data/virus_ids.pkl', 'w') as f:
			pickle.dump(virus_ids, f)

	# Generate mouse_ids
	if m:
		mouse = df['Host_strain'].unique()
		mouse_ids = set(mouse)

		with open('data/mouse_ids.pkl', 'w') as f:
			pickle.dump(mouse_ids, f)

def main():
	# Parser arguments
	parser = argparse.ArgumentParser(
		description='To generate virus and mouse ids.')

	parser.add_argument('-v', '--virus', nargs='?', default=True, 
		type=bool, help='Set True to generate virus_ids.')
	parser.add_argument('-m', '--mouse', nargs='?', default=True, 
		type=bool, help='Set True to generate mouse_ids.')

	args = parser.parse_args()
	create_ids(args.virus, args.mouse)

if __name__ == "__main__":
	main() 