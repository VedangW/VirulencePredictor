#!usr/bin/python

import pickle
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from utils import remove_zeros

from keras.models import Model
from keras.layers import Input, Dense
from sklearn.decomposition import PCA

threshold_len = 9

def reduce_dimension_pca(features, particle, i):
	X = features.values()
	for i in range(len(X)):
		X[i] = np.concatenate(X[i], axis=0)
	X = np.vstack(X)

	print "Shape of features before reduction =", X.shape

	if X.shape[1] > 1000000:
		print "Length too high, removing zeros..."
		X = remove_zeros(X)
		print "Shape of X after removing zeros =", X.shape

	print "Applying PCA..."
	pca = PCA(n_components=threshold_len)
	X = pca.fit_transform(X)

	del pca

	print "Shape of features after reduction =", X.shape
	print "Done."
	return X

def reduce_dimension_ae(features, particle, i):
	vals = features.values()
	X = np.vstack([item for sublist in vals for item in sublist])

	num_epochs = 100

	print "Training Intermediate autoencoder..."

	input_len = np.shape(X)[1]
	inter_dim = 100

	# Input layer
	input_layer = Input(shape=(input_len,))

	# Encoder layers
	encoded1 = Dense(512, activation='relu')(input_layer)
	encoded2 = Dense(256, activation='relu')(encoded1)
	encoded3 = Dense(128, activation='relu')(encoded2)
	latent = Dense(inter_dim, activation='relu')(encoded3)

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
	encoded_input = Input(shape =(inter_dim,))

	print "Encoding batches..."
	final_dict = {}
	virus_order = []
	viruses = features.keys()
	for i in range(len(viruses)):
		batch = np.vstack(features[viruses[i]])
		batch_enc = encoder.predict(batch)
		final_dict[viruses[i]] = np.concatenate(batch_enc, axis=0)
		virus_order.append([i, viruses[i]])

	X_final = np.vstack(final_dict.values())

	"Training final autoencoder..."

	input_len = np.shape(X_final)[1]
	final_dim = 100

	# Input layer
	input_layer = Input(shape=(input_len,))

	# Encoder layers
	encoded1 = Dense(512, activation='relu')(input_layer)
	encoded2 = Dense(256, activation='relu')(encoded1)
	encoded3 = Dense(128, activation='relu')(encoded2)
	latent = Dense(final_dim, activation='relu')(encoded3)

	# Decoder layers
	decoded1 = Dense(128, activation='relu')(latent)
	decoded2 = Dense(256, activation='relu')(decoded1)
	decoded3 = Dense(512, activation='relu')(decoded2)
	reconstruction = Dense(input_len, activation='sigmoid')(decoded3)

	# Define the whole autoencoder model
	autoencoder = Model(inputs=input_layer, outputs=reconstruction)
	autoencoder.compile(optimizer='adadelta', loss ='binary_crossentropy')

	# Fit the values in the AE
	autoencoder.fit(X_final, X_final, epochs=num_epochs,
		batch_size=2, shuffle=False,
		validation_split=0.1)

	# Define the encoder model
	encoder = Model(inputs=input_layer, outputs=latent)
	encoded_input = Input(shape =(final_dim,))

	# Predict the training set on it again
	X_final_enc = encoder.predict(X_final)

	print "Shape of final encoded feature vector =", np.shape(X_final_enc)

	if particle == 'virus':
		fname = 'data/virus_enc_ae.pkl'
		with open(fname, 'w') as f:
			pickle.dump(X_final_enc, f)
		with open('data/virus_order.pkl', 'w') as f:
			pickle.dump(virus_order, f)
 	elif particle == 'mouse':
		fname = 'data/mouse_enc_ae' + i + '.pkl'
		with open(fname, 'w') as f:
			pickle.dump(X_final_enc, f)