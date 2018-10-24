## Storing of Data

Here is the way that the data is stored for the project. The data is stored at ```VP/data```.

Following are the details of the files:

#### Files with data to use

- ```graph_info.pkl```: Read using pickle. The object retrieved, say ```g``` will be a tuple of 5 elements. ```g[0]``` is the number of viruses, ```g[1]``` is the number of mice. ```g[2], g[3] and g[4]``` make up the flu dataset, that is, they represent ```Influenza_virus_name```, ```Host_strain``` and ```LD50``` respectively. However, the influenza viruses are encoded as integers instead of their names. So for example ```g[4][20]``` is the LD50 value for virus ```g[2][20]``` and mouse ```g[3][20]```.

- ```virus_enc_ae.pkl```: Read using pickle. Has a numpy array of features for the viruses of shape (202, 100) (202 virus genomes and 100 features for each).

- ```mouse_enc_pca.pkl```: Read using pickle. Has a numpy array of features for the mice of shape (12, 100) (12 mouse genomes and 100 features for each).

#### Host-Pathogen-LD50 dataset:

- ```merged.csv```: Contains the final flu dataset after merging all versions of the available viruses. Rows are of the form - ```Host_strain``` - ```Influenza_virus_name``` - ```LD50```.

#### Files which describe the data

- ```mouse_ids.pkl```: Contains the names of all the mice in the dataset.

- ```virus_ids.pkl```: Contains the names of all the viruses in the dataset.

- ```mouse_order.pkl```: Contains the mapping between the mouse IDs and the numbers which are used to represent them in ```graph_info.pkl```.

- ```virus_order.pkl```: Contains the mapping between the virus IDs and the numbers which are used to represent them in ```graph_info.pkl```.