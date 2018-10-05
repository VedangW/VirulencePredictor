#!usr/bin/python

import os
import sys
import click
import argparse

def write_to_log_file(log, fname='terminals.log'):
    """ A function to write errors to a log file.

        Parameters
        ----------
        log: list
            The log of errors.

        fname: str
            Name of the log file to write to.
    """


    # Open log file
    f = open(fname, 'w')

    # Write to log file
    for a in log:
        a += '\n'
        f.write(a)

    f.close()

def remove_terminals(fname, dirname):
    """ A function which removes the terminal symbol, '*'
        from a Segment file (fasta format).
        
        The symbol needs to be present at the end of sequences
        for this to happen.
        
        All errors are stored in a list called 'log' and returned.
        
        Parameters
        ----------
        fname: str
            A string which is the path of the original file.

        dirname: str
            Name of folder in which the file exists.
            
        Returns
        -------
        log: list
            A log containing the line numbers and the errors
            associated, should they arise.
    """
    
    # Read the lines from fname
    fpath = dirname + '/' + fname
    f = open(fpath, 'r')
    lines = f.readlines()
    f.close()
    
    # new_ is the new set of lines.
    new_ = []

    log = []
    for i in range(len(lines)):
        try:
            # Even line numbers contain the identifier.
            if i % 2 == 0:
                new_.append(lines[i])
            # Odd line numbers contain the sequences so remove the *.
            else:
                l = lines[i].strip().split('*')[0] + '\r\n'
                new_.append(l)
        except:
            # Add error messages to log
            e, message, _tb = sys.exc_info()
            log.append(fname + ': line ' + str(i) + ': ' + str(e) + ': ' + str(message))
            continue
    
    # Open a new file to write into
    new_fname = fname + '_new'
    new_fpath = dirname + '/' + new_fname
    f = open(new_fpath, 'w')
    
    # Write the new lines to the new file.
    for n in new_:
        f.write(n)

    f.close()
    return log

def remove_terminals_from_dir(dir):
    """ Remove the terminals from all files in a directory.

        Errors are stored in a log file, 'terminals.log'.

        Parameters
        ----------
        dir: str
            Name of the directory from which to read files.
    """

    files = os.listdir(dir)

    full_log = []
    label = 'Processing files...'
    with click.progressbar(files, label=label) as bar:
        for fname in bar:
            log = remove_terminals(fname, dir)
            full_log += log

    # If the log is not empty, create a log file
    if full_log: 
        write_to_log_file(full_log)

def main():
    # Parser argument for delimiter
    parser = argparse.ArgumentParser(description='Removes terminals (*) from all files in a directory')
    parser.add_argument('-d', '--dir', nargs='?', default='segments', 
        type=str, help='Name of directory to read files from.')
    args = parser.parse_args()

    remove_terminals_from_dir(args.dir)

if __name__ == "__main__":
    main()