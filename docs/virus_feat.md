## Preprocessing for Viruses

Initially, the data for viruses is given as a set of fasta files (let the number of these be ```num_virus```, each containing a set of protein sequences, representing the individual segments of that virus (```PB2```, ```PB1```, ```PB1-F2```, ```PA```, ```HA``` and so on). Following are the steps taken to convert this data into features for each of these virus files. It should be noted that the final number of rows in the feature matrix will be exactly equal to the number of files present in the input directory.

#### Alignment of segments

Each virus has in it a certain number of segments and this can be upto 13. These segments are constant throughout all proteins in the dataset. The first step is to collect these segments into 13 segment files, each of which ideally should have sequences from all the viruses (but don't and this problem is tackled later). This is performed by the ```collect_segments.py``` module.

As these collected segments may have a terminal codon at the end (```'*'```), this is removed first before alignment using the module ```remove_terminals.py```.

These segments are then individually aligned using the alignment tool [```MUSCLE```](https://www.ebi.ac.uk/Tools/msa/muscle/).

#### Removing anomalies

There might be entries in the sequences which don't correspond to any of the [amino acid characters](http://www.cryst.bbk.ac.uk/education/AminoAcid/the_twenty.html). All of these are first converted to an ```X``` in the sequences. The two special cases, ```B``` and ```Z``` are also replaced. This is done by ```preprocess_align```. The characters ```X``` are then replaced with the mode of the column they exist in. Should the mode itself be ```X```, the column is dropped altogether.

#### Encoding of sequences

This is done using [AAindex](https://www.ncbi.nlm.nih.gov/pubmed/9847231), the mapping index being ```JOND920101``` (Relative frequency of occurrence (Jones et al., 1992)). The mapping for ```'-'``` is given as ```0.0```. The sequences are then padded so that they reach the same length, say ```unified_len```, using the value ```0.0``` for padding. 

As mentioned before, the number of sequences for each segment in the dataset is not the same at this point of time. To make them equal, sequences of zeros of length ```unified_len``` are inserted in them. All of this is achieved by ```features.py```.
Let the number of sequences in each segment at this point be ```segment_len```.

#### Dimension reduction and recombination

The dimensions of these sequences are then reduced as follows:

1. All the sequences from all segments are stacked together to form an array of sequences of length ```13 * segment_len```. This array is used to train an autoencoder, with the hidden length ```inter_dim```.
2. To encode all sequences in the array, batches of length ```13``` are created, each batch corresponding to the sequences from one particular virus. This is then passed through the encoder network to encode them. All of these vectors in the batch are then concatenated to form a single vector. Thus, at the end of this step, there would be ```num_virus``` number of vectors.
3. Another autoencoder is used to reduce the length of these vectors to ```100``` (an arbitrary small number).

Thus, the virus features are generated. 