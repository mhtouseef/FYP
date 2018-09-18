#get all the files present in that directory
import glob
list_of_files = glob.glob('./*.txt')           # create the list of file
for file_name in list_of_files:
    seq=file_name
    print(seq)
    #fasta_reader(seq)    #function call to remove header