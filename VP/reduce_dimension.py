#!usr/bin/python

import pickle
import numpy as np

from config import THRESHOLD
from sklearn.decomposition import PCA, IncrementalPCA

# Length of reduced data
threshold_len = THRESHOLD

def remove_zeros(X):
	""" Remove columns where all values are zeros.
		This function is not being used right now.
	"""
	X_t = X.T
	where_zeros = np.where(~X_t.any(axis=1))[0]
	X_t = np.delete(X_t, where_zeros, axis=0)
	return X_t.T

def reduce_dimension(features, particle, index,
					method='ipca', save=True):
	""" Function to reduce the dimension of 
		features using Incremental PCA.

		Use when number of alignment files is
		too large.
	"""

	# Get names of the mice/viruses
	keys = np.array(features.keys())

	# Concatenate each set of feature vectors
	X = features.values()
	for i in range(len(X)):
		X[i] = np.concatenate(X[i], axis=0)
	X = np.vstack(X)

	# Remove zeros
	print "Removing zeros..."
	X = remove_zeros(X)
	print "Shape of X after removing zeros =", X.shape

	# Reduce dimension using IPCA
	if method == 'ipca':
		reducer = IncrementalPCA(n_components=threshold_len, 
			batch_size=30)
	elif method == 'pca':
		reducer = PCA(n_components=threshold_len)

	# Reduce the features
	X_reduced = reducer.fit_transform(X)

	# Check if any of the features is all zeros
	for i in range(len(X_reduced)):
		if np.array_equal(X_reduced[i], np.zeros(threshold_len,)):
			err = "All zeros in sample for " + keys[i]
			raise ValueError(err)

	# Concatenate keys and reduced features
	keys = keys.reshape((keys.shape[0], 1))
	X = np.concatenate([keys, X_reduced], axis=1)

	print "Shape of features after reduction =", X.shape

	# Save to disc
	if save:
		print "Saving data..."
		if particle == 'virus':
			fname = 'data/virus_' + index
			with open(fname, 'w') as f:
				pickle.dump(X, f)
	 	elif particle == 'mouse':
			fname = 'data/mouse_' + index
			with open(fname, 'w') as f:
				pickle.dump(X, f)

	return X