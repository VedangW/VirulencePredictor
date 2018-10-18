#!usr/bin/python

import pickle
import numpy as np

from keras.models import Model
from keras.layers import Input, Dense
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

threshold_len = 100

def remove_zeros(X):
	X_t = X.T
	where_zeros = np.where(~X_t.any(axis=1))[0]
	X_t = np.delete(X_t, where_zeros, axis=0)
	return X_t.T

def pad_with_zero_vectors(X):
	_n_rows, n_cols = np.shape(X)
	new_shape = ((threshold_len, n_cols))
	X_changed = np.zeros(new_shape)
	X_changed[:X.shape[0],:X.shape[1]] = X

	return X_changed

def build_and_train_AE(X, encoding_dim=100,
						input_len=10036,
						num_epochs=2):

	""" Function which builds a deep autoencoder 
		to reduce the dimensions of the training data.

		Parameters
		----------
		X: list of lists or similar
			Contains all the vectors to reduce in dimension.

		input_len: int
			Length of input vector.

		encoding_dim: int
			Dimension of the latent layer.

		num_epochs: int
			Number of epochs to run on model.

		Returns
		-------
		X_enc: list of lists or similar
			Dimensionally reduced X.
	"""
	print "Reducing dimension."
	print "Input length is", input_len
	print "Output length will be", encoding_dim
	
	print ""
	print "Building and training autoencoder..."

	n_rows, n_cols = np.shape(X)
	print n_rows, n_cols

	if input_len > threshold_len:
		print "Input length too high, removing zero columns..."
		X = remove_zeros(X)
		print "Applying PCA..."
		print "Padding..."
		X = pad_with_zero_vectors(X)
		print "Transforming..."
		pca = PCA(n_components=threshold_len)
		X = pca.fit_transform(X)

	print "Feature sized reduced by PCA."
	print "New shape =", np.shape(X)
	input_len = np.shape(X)[1]

	# Input layer
	input_layer = Input(shape=(input_len,))

	# Encoder layers
	encoded1 = Dense(512, activation='relu')(input_layer)
	encoded2 = Dense(256, activation='relu')(encoded1)
	encoded3 = Dense(128, activation='relu')(encoded2)
	latent = Dense(encoding_dim, activation='relu')(encoded3)

	# Decoder layers
	decoded1 = Dense(128, activation='relu')(latent)
	decoded2 = Dense(256, activation='relu')(decoded1)
	decoded3 = Dense(512, activation='relu')(decoded2)
	reconstruction = Dense(input_len, activation='sigmoid')(decoded3)

	# Define the whole autoencoder model
	autoencoder = Model(inputs=input_layer, outputs=reconstruction)
	autoencoder.compile(optimizer='adadelta', loss ='binary_crossentropy')

	# Fit the values in the AE
	autoencoder.fit(X, X, epochs=num_epochs,
		batch_size=2, shuffle=False,
		validation_split=0.2)

	# Define the encoder model
	encoder = Model(inputs=input_layer, outputs=latent)
	encoded_input = Input(shape =(encoding_dim,))

	# Predict the training set on it again
	X_enc = encoder.predict(X)

	print ""
	print "Dimension has been reduced."
	# Check shape of the encoded model
	print "Shape of data now =", np.shape(X_enc)
	print ""

	return X_enc

def reduce_dimension(features, particle):
	""" Reduce the dimension of the original
		set of vectors using one of AE or VAE.

		VAE has not been implemented yet.

		Parameters
		----------
		particle: str
			'virus' or 'mouse'.
	"""
	# Create a training set
	X = np.vstack(features.values())
	input_len = len(X[0])

	# # Reduce dimension using AE
	X_enc = build_and_train_AE(X, input_len=input_len)

	print "Saving file to disc..."
	# Save to disc
	if particle == 'virus':
		with open('data/features_virus_reduced_ae.pkl', 'w') as f:
			pickle.dump(X_enc, f)
	elif particle == 'mouse':
		with open('data/features_mouse_reduced_ae.pkl', 'w') as f:
			pickle.dump(X_enc, f)

	print "Done."