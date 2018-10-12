#!usr/bin/python

import pickle
import numpy as np

from keras.models import Model
from keras.layers import Input, Dense

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
	print ""
	print "Reducing dimension."
	print "Input length is", input_len
	print "Output length will be", encoding_dim
	
	print ""
	print "Building and training autoencoder..."

	# Input layer
	input_layer = Input(shape=(input_len,))

	# Encoder layers
	encoded1 = Dense(4096, activation='relu')(input_layer)
	encoded2 = Dense(2048, activation='relu')(encoded1)
	encoded3 = Dense(1024, activation='relu')(encoded2)
	encoded4 = Dense(512, activation='relu')(encoded3)
	encoded5 = Dense(256, activation='relu')(encoded4)
	encoded6 = Dense(128, activation='relu')(encoded5)
	latent = Dense(encoding_dim, activation='relu')(encoded6)

	# Decoder layers
	decoded1 = Dense(128, activation='relu')(latent)
	decoded2 = Dense(256, activation='relu')(decoded1)
	decoded3 = Dense(512, activation='relu')(decoded2)
	decoded4 = Dense(1024, activation='relu')(decoded3)
	decoded5 = Dense(2048, activation='relu')(decoded4)
	decoded6 = Dense(4096, activation='relu')(decoded5)
	reconstruction = Dense(10036, activation='sigmoid')(decoded6)

	# Define the whole autoencoder model
	autoencoder = Model(inputs=input_layer, outputs=reconstruction)
	autoencoder.compile(optimizer='adadelta', loss ='binary_crossentropy')

	# Fit the values in the AE
	autoencoder.fit(X, X, nb_epoch=num_epochs, 
		batch_size=32, shuffle=False, 
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

def reduce_dimension():
	""" Reduce the dimension of the original
		set of vectors using one of AE or VAE.

		VAE has not been implemented yet.
	"""

	# Load from file with full vectors
	with open('data/features_full.pkl') as f:
		features = pickle.load(f)

	# Create a training set
	X = np.vstack(features.values())

	# Reduce dimension using AE
	X_enc = build_and_train_AE(X)

	print "Saving file to disc..."
	# Save to disc
	with open('data/features_reduced_ae.pkl', 'w') as f:
		pickle.dump(X_enc, f)

	print "Done."

def main():
	reduce_dimension()

if __name__ == "__main__":
	main()