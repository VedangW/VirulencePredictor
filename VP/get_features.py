#!usr/bin/python

""" Extracts features from virus or mouse particles
	as needed.

	Uses files: preprocess_align, features, reduce_dimension
	in that order.
"""
import os
import pickle
import argparse
import numpy as np
import features as feat
import reduce_dimension as red
import preprocess_align as prep

import warnings
warnings.filterwarnings("ignore")

from sklearn.decomposition import PCA
from utils import pad_with_zero_vectors
from keras.preprocessing.sequence import pad_sequences

def extract_reduced_features_from_alignments(aligned_dir, overwrite, particle, preprocess=True):
	""" Calls functions from different files to extract
		features from alignments.

		Following is done:
		1. Alignments are preprocessed (replace rogue values with 'X').
		2. Alignments are converted to encodings.
		3. Encodings are reduced in size using AE or VAE.

		Parameters
		----------
		aligned_dir: str
			Directory where the alignments are present.

		overwrite: bool
			Whether to overwrite preprocessed files.

		particle: str
			'virus' or 'mouse'

		preprocess: bool
			True if you want to preprocess.
	"""
	if preprocess == True:
		print ""
		print "Preprocessing alignments ==>"
		print ""
		prep.preprocess_alignments(aligned_dir, overwrite, particle)

	if particle == 'virus':
		print ""
		print "Extracting features ==>"
		print ""
		features = feat.extract_features(particle, "")
		print ""
		print "Reducing dimension ==>"
		print ""
		red.reduce_dimension_ae(features, particle, 0)

	elif particle == 'mouse':
		dir_path = '/home/vedang/Documents/mouse_prep/'
		mouse_dirs = os.listdir(dir_path)

		reduced = []
		for i in range(len(mouse_dirs)):
			print ""
			print "Extracting features from batch", i, "==>"
			print ""
			features = feat.extract_features(particle, dir_path + mouse_dirs[i])
			mouse_order = features.keys()
			print ""
			print "Reducing dimension ==>"
			print ""
			reduced.append(red.reduce_dimension_pca(features, particle, i))

		reduced = np.concatenate(reduced, axis=1)

		print ""
		if len(reduced) < 100:
			print "Length is not 100."
			print "Shape of array is =", reduced.shape
			print "Padding the sequence..."
			reduced = pad_sequences(reduced, maxlen=100, dtype='float32', 
				padding='post', truncating='post', value=0.0)

			print "Shape of array after padding =", reduced.shape
		elif len(reduced) > 100:
			print "Length is greated than 100."
			print "Shape of array is =", reduced.shape
			n_particles = reduced.shape[0]

			print "Reducing length by PCA..."
			reduced = pad_with_zero_vectors(reduced, 100)
			pca = PCA(n_components=100)
			reduced = pca.fit_transform(reduced)
			reduced = reduced[:n_particles, :]
			del pca
			print "Final shape of array =", reduced.shape

		print "Saving to file..."
		fname = 'data/mouse_enc_pca.pkl'
		with open(fname, 'w') as f:
			pickle.dump(reduced, f)

		for i in range(len(mouse_order)):
			mouse_order[i] = [i, mouse_order[i]]

		with open('data/mouse_order.pkl', 'w') as f:
			pickle.dump(mouse_order, f)

		print ""
		print "Done."

def main():
	# Parser arguments
	parser = argparse.ArgumentParser(
		description='Preprocesses the aligned segments.')

	parser.add_argument('-d', '--aligned_dir', nargs='?', default='data/aligned', 
		type=str, help='Name of directory containing alignments.')
	parser.add_argument('-ow', '--overwrite', nargs='?', default=True, 
		type=bool, help='Set True to overwrite.')
	parser.add_argument('-p', '--particle', nargs='?', default='virus',
		type=str, help="'mouse' or 'virus'")
	parser.add_argument('-pr', '--preprocess', nargs='?', default=True,
		type=bool, help="Set True to preprocess")

	args = parser.parse_args()

	# Check if the wrong flags have not been given
	if args.particle == 'virus' or args.particle == 'mouse':
		pass
	else:
		raise ValueError("Particle flag can only be 'virus' or 'mouse'.")

	extract_reduced_features_from_alignments(args.aligned_dir, 
		args.overwrite, args.particle, False)

if __name__ == "__main__":
	main()