#!usr/bin/python

import os
import pickle
import argparse
import numpy as np

from features import *
from reduce_dimension import *
from preprocess_align import *

class EmbeddingsByAAindex():

	def __init__(self, particle,
				index, preprocess=True,
				reducer='ipca',
				overwrite=False, save=True):
		""" Class to convert proteomes to AAindex embeddings. """

		self.particle = particle
		self.index = index
		self.preprocess = preprocess
		self.reducer = reducer
		self.overwrite = overwrite
		self.save = save

	def __call__(self, aligned_dir):
		""" Convert directory containing 
			alignments to features. """

		self.aligned_dir = aligned_dir
		print

		if self.preprocess:
			print "Preprocessing data ->"
			preprocess_alignments(self.aligned_dir, 
									self.overwrite, 
									self.particle)
			print "Done."
			print
			
		self._extract_and_reduce()

	def _extract_and_reduce(self,):
		# Extract features from preprocessed data
		print "Extracting features ->"
		features = extract_features(self.particle, 
									self.index)
		print "Done."

		print 

		# Reduce dimension by IPCA
		print "Reducing dimension ->"
		reduce_dimension(features, self.particle, 
			self.index, self.reducer, self.save)
		print "Done."


# Parser arguments
parser = argparse.ArgumentParser(
	description='Preprocesses the aligned segments.')

parser.add_argument('-d', '--aligned_dir', nargs='?', 
					default='data/aligned', type=str, 
					help="Directory of alignments.")

parser.add_argument('-p', '--particle', nargs='?', 
					choices=['virus', 'mouse'], type=str, 
					help="'mouse' or 'virus'")

parser.add_argument('-i', '--index', nargs='?', 
					default='JOND920101', type=str, 
					help="Index from AAindex to use")

parser.add_argument('-rm', '--reduction_method', nargs='?', 
					choices=['ipca', 'pca', 't-sne'], 
					default='ipca', type=str, 
					help="Reduction method to use.")

parser.add_argument('-s', '--save', 
					help="Set true to save features",
					action='store_true')

parser.add_argument('-pr', '--preprocess',
					help="Set True to preprocess", 
					action='store_true')

parser.add_argument('-ow', '--overwrite',
					help="Set True to overwrite", 
					action='store_true')

args = parser.parse_args()

emb_aaindex = EmbeddingsByAAindex(particle=args.particle,
								index=args.index,
								preprocess=args.preprocess,
								overwrite=args.overwrite,
								reducer=args.reduction_method,
								save=args.save)

emb_aaindex(args.aligned_dir)