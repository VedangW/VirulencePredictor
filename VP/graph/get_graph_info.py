#!usr/bin/python

import pickle
import numpy as np
import pandas as pd

def generate_graph_info():
	""" A function to retrieve relevant information from
		the graph to pass to GCMC.
	"""
	# Read full merged dataset
	df = pd.read_csv('data/merged.csv')

	# List of unique mice and viruses
	mice = df['Host_strain'].unique()
	viruses = df['Influenza_virus_name'].unique()

	# Total unique mice and viruses
	num_mice = len(mice)
	num_viruses = len(viruses)

	virus_dict = {}
	mouse_dict = {}

	count = 0
	# Add viruses
	for v in viruses:
		virus_dict[v] = count
		count += 1
		
	count = 0
	# Add mice
	for m in mice:
		mouse_dict[m] = count
		count += 1

	v_nodes = df['Host_strain'].apply(lambda x: mouse_dict[x])
	u_nodes = df['Influenza_virus_name'].apply(lambda x: virus_dict[x])

	u_nodes = np.array(u_nodes.tolist())
	v_nodes = np.array(v_nodes.tolist())

	# LD50 values as edge weights
	ld50 = np.array(df['LD50'].tolist())

	# Store all together and save
	graph_info = (num_viruses, num_mice, u_nodes, v_nodes, ld50)

	with open('data/graph_info.pkl', 'w') as f:
		pickle.dump(graph_info, f)

def main():
	generate_graph_info()

if __name__ == "__main__":
	main()