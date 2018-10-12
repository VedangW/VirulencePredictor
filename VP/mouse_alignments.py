# ALIGNMENT

def alignProt(faa_in,faa_aln_out):
    
    # INPUT
    #   faa_in: the path for faa (protein) file
    #   faa_aln_out: the path for the file containing the aligned proteins
    # OUTPUT
    #   File faa_aln_out containing the aligned proteins
    
    from Bio import SeqIO
    from Bio.Align.Applications import MuscleCommandline

    faarecords = list(SeqIO.parse(faa_in, "fasta"))
    fsize = os.path.getsize(faa_in)

    if len(faarecords) == 0:
        open(faa_aln_out,'a').close()
    else:
        muscle_exe = os.path.join(os.getcwd(),'muscle')
        if fsize < 1000000:
            cline = MuscleCommandline(muscle_exe,input=faa_in,out=faa_aln_out)
        else:
            cline = MuscleCommandline(muscle_exe,input=faa_in,out=faa_aln_out,\
            maxiters=1,diags=True,sv=True)
        print(cline)
        cline()

import os

# Directory for input and output
dir_in = os.path.join(parent_dir, "RECOMB2019", "per-mouse-protein")
dir_out = os.path.join(parent_dir, "RECOMB2019", "per-mouse-protein-aligned4")
if not os.path.exists(dir_out): os.makedirs(dir_out)

files = [re.sub('.fasta','',x) for x in os.listdir(dir_in)]

#for f in files[347:348]:
for f in files[15835:15836]:
    
    faa_in = os.path.join(dir_in,f+'.fasta')
    faa_aln_out = os.path.join(dir_out,f+'.fasta')
    alignProt(faa_in,faa_aln_out)