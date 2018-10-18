#!usr/bin/python

""" Extracts features from virus or mouse particles
	as needed.

	Uses files: preprocess_align, features, reduce_dimension
	in that order.
"""
import argparse
import features as feat
import reduce_dimension as red
import preprocess_align as prep

def extract_reduced_features_from_alignments(aligned_dir, overwrite, particle):
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
	"""
	print ""
	print "Preprocessing alignments ==>"
	print ""
	prep.preprocess_alignments(aligned_dir, overwrite, particle)
	print ""
	print "Extracting features ==>"
	print ""
	features = feat.extract_features(particle)
	print ""
	print "Reducing dimension ==>"
	print ""
	red.reduce_dimension(features, particle)

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

	args = parser.parse_args()

	# Check if the wrong flags have not been given
	if args.particle == 'virus' or args.particle == 'mouse':
		pass
	else:
		raise ValueError("Particle flag can only be 'virus' or 'mouse'.")

	extract_reduced_features_from_alignments(args.aligned_dir, 
		args.overwrite, args.particle)

if __name__ == "__main__":
	main()