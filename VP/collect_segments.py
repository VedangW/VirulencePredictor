#!usr/bin/python

import os
import click
import argparse

def get_seg_num(name, default):
	""" Get segment names according to swapping rules. """ 

	# Create a dictionary for swapping
	swap_dict = {'Seg1p1': 'Seg2p1', 
	'Seg2p1': 'Seg1p1', 
	'Seg6p1': 'Seg6p2',
	'Seg6p2': 'Seg6p1'}

	# Keys in the dictionary
	keys = list(swap_dict.keys())

	# Swap seg_num
	if name == 'B' and default in keys:
		return swap_dict[default]

	return default

def store_segments(f, fname, delimiter):
	""" Function to write the segments into the respective file name. """ 

	lines = f.readlines()

	# Can be 'A', 'B', 'maA', 'rA'
	name = fname.split('.')[1]

	for i in range(len(lines)):
		if not i % 2:
			default = lines[i].strip().split('_')[-1]
			seg_num = get_seg_num(name, default)
			seg_name = lines[i].strip()
		elif i % 2:
			f_seg = open('Segments/' + seg_num, 'a')
			write_str = seg_name + delimiter + lines[i]
			f_seg.write(write_str)
			f_seg.close()

def collect_segments(delimiter):
	""" The caller function which iterates through all files. """

	new_dir = 'Segments'
	if not os.path.exists(new_dir):
		os.makedirs(new_dir)

	with click.progressbar(os.listdir('proteomes'), label='Processing files') as bar:
		for fname in bar:
			f = open('proteomes/' + fname, 'r')
			store_segments(f, fname, delimiter)
			f.close()

def main():
	# Parser argument for delimiter
	parser = argparse.ArgumentParser(description='Collects segments from proteomes')
	parser.add_argument('--delim', nargs='?', const='\n', default='\n', 
		type=str, help='Delimiter for the segment files')
	args = parser.parse_args()

	collect_segments(args.delim)

if __name__ == "__main__":
	main()