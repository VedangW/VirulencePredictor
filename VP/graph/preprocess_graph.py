#!usr/bin/python

import os
import pickle
import argparse
import numpy as np
import pandas as pd

# Disable warnings
import warnings
warnings.filterwarnings("ignore")

def save_to_csv(df, fname='collected.csv', delimiter=','):
	""" A function which saves a dataframe to a csv
		file.

		Parameters
		----------
		df: pd.DataFrame
			The dataframe to save.

		fname: str
			The name of the csv file to store in.

		delimiter: str
			The character which seperates the column entries
			in the csv file.
	"""
	df.to_csv(fname, sep=delimiter)

def LD50parse(ld50):
	""" A function to return the log to the base 10
		of the ld50 value. This function is used with
		df.apply in this context.

		For example: '10E5.5' is returned as 5.5.

		Parameters
		----------
		ld50: str
			An ld50 value written as a string.

		Returns
		-------
		log ld50: float
			The log to the base 10 of the ld50 input.
	"""
	return float(ld50.split('E')[1])

def preprocess_df(df, out_file, delim):
	""" A function to perform DataFrame operations on the df passed
		to it. Eventually the dataframe is saved as a CSV file.
		
		The function removes redundancies related to virus-mouse pairs
		and the concentration (LD50).Extra details such as Author name,
		Notes on LD50 and so on are dropped from the df.The LD50 
		value is converted to ceil(log LD50)).

		Parameters
		----------
		df: pd.DataFrame
			A dataframe read from the source text file.

		out_file: str
			Name of csv file to save in.

		delim: str
			Delimiter for csv file.
	"""

	print ("Preprocessing dataframe...")

	# The column 'virus_mouse' is to remove redundancies where the mouse
	# and the virus are same. 'vmc' is for removing redundancies where
	# virus, mouse and concentration are the same.
	print ("Removing redundancies...")
	df['virus_mouse'] = df['Influenza_virus_name'] + df['Host_strain']
	df['vmc'] = df['virus_mouse'] + df['LD50']
	df.drop_duplicates(subset='vmc', keep="first", inplace=True)

	# Drop these columns at this point because they are useless and have a lot of 
	df.drop(['LD50_unit', 'Notes on LD50', 'vmc'], axis=1, inplace=True)

	print ("Dropping NaNs...")
	# Remove wherever LD50 is 'NA' or NaN or None
	df = df[df['LD50'] != 'NA']
	df.dropna(inplace=True)

	print ("Transforming LD50 values...")
	# Take the log of the LD50 column
	df['LD50'] = df['LD50'].apply(lambda x: LD50parse(x))
	df['LD50'] = df['LD50'].apply(lambda x: float(x))

	# Get all instances of virus_mouse which repeat
	mult = []
	for v in df['virus_mouse'].unique():
		samps = df[df['virus_mouse'] == v]
		if len(samps) > 1:
			mult.append(v)

	# Create a dictionary for all values of 'virus_mouse'
	# which occur multiple times. The value against them in
	# the dictionary will be the mean of all the corresponding
	# values of LD50 for them.
	mult_dict = {}
	for v in mult:
		try:
			mult_dict[v] = df[df['virus_mouse'] == v]['LD50'].mean()
		except: 
			print df[df['virus_mouse'] == v]['LD50']

	for v in mult:
		df['LD50'].loc[df['virus_mouse'] == v] = mult_dict[v]

	# Create a vmc column again and drop according to redundancies in that
	df['vmc'] = df['virus_mouse'] + df['LD50'].apply(lambda x: str(x))
	df.drop_duplicates(subset='vmc', keep='first', inplace=True)

	print ("Dropping unwanted columns...")
	# Drop other redundant columns
	df.drop(['Reference', 'virus_mouse', 'vmc'], axis=1, inplace=True)

	# Apply the ceiling function to all values of LD50
	df['LD50'] = df['LD50'].apply(lambda x: int(x) + 1)

	print ("Saving dataframe to csv...")

	# Save dataframe to a csv file
	save_to_csv(df=df, fname=out_file, delimiter=delim)

def preprocess(fname, out_file, delim):
	""" A function that reads data from the source file
		and stores it in a dataframe which is passed to preprocess_df()
		for further preprocessing.

		Parameters
		----------
		fname: str
			Name of file to read.

		out_file: str
			Name of file to store data into.

		delim: str
			Delimiter for csv file.
	""" 

	print ("Reading data into dataframe...")
	# Read txt file for data
	f = open(fname, 'r')
	lines = f.readlines()

	# Strip and split data into a list
	data = []
	for i in range(len(lines)):
		data.append(lines[i].strip().split('\t'))

	# Define a dataframe
	cols = data.pop(0)
	df = pd.DataFrame(data, columns=cols)

	preprocess_df(df, out_file, delim)

def main():
	# Parser arguments
	parser = argparse.ArgumentParser(
		description='Preprocesses the data containing host-pathogen interaction.')

	parser.add_argument('--delim', nargs='?', default=',', 
		type=str, help='Delimiter for the CSV file.')

	parser.add_argument('--data', nargs='?', default='flu-infections-for-recomb.txt', 
		type=str, help='Name of data file.')

	parser.add_argument('--out', nargs='?', default='collected.csv', 
		type=str, help='Name of output file.')

	parser.add_argument('--ow', nargs='?', default=False,
		type=bool, help='Set to true to overwrite out file if already exists.')

	args = parser.parse_args()

	in_file = 'data/' + args.data
	out_file = 'data/' + args.out

	# Check if mentioned file exists in the directory.
	if not os.path.isfile(in_file):
		raise ValueError('Mentioned data file is not present.')

	# Check if output file already exists.
	if os.path.isfile(out_file) and not args.ow :
		raise ValueError('A file by the same name already exists.')

	preprocess(in_file, out_file, args.delim)

	print ("Done.")

if __name__ == "__main__":
	main()