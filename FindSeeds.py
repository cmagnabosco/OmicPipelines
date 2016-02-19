from Bio import SeqIO
import os
from optparse import OptionParser

############################################  
## PARSE COMMAND LINE ARGUMENTs
############################################  

def parse_args():
    parser = OptionParser(description = 'Generate a seed file for TargetedAssembly.py', epilog="Example: python FindSeeds.py --ref_file seed.fa --ref_type nucl --fwd metagenome.fa --outFile myseeds.fa [options]")
    parser.set_defaults(n_threads='1',evalue='1e-5')
## Required options
    parser.add_option('--ref_file', type = 'string',dest='ref_file',help='FASTA (nucleotide or protein) file of reference sequences')
    parser.add_option('--ref_type',type='choice',choices=['nucl','prot'], dest='ref_type',help='Indicate whether reference file is a nucleotide (nucl) or protein (prot) file.')
    parser.add_option('--fwd',type='string', dest='fwd', help='FASTA/Q file of the metagenome. For paired end reads, -fwd will indicate the file with forward reads and -rev indicate the reverse reads')
    parser.add_option('--outFile', type='string',dest='outFile',help='Name of the output file that will contain the seed sequences (FASTA format).')
## Other options
    parser.add_option('--rev', dest='rev',help = 'FASTA/Q file that contains the reverse reads of paired end sequencing')
    parser.add_option('--n_threads',type='string',dest='n_threads',help="Number of threads to use [1]")
    parser.add_option('--evalue',type='string',dest='evalue',help='evalue cut-off to use if --search Y [1e-5]')
## Parse args
    (options, args) = parser.parse_args()
## Check for required args
    if not options.ref_file:
        parser.error('--ref_file not given. Input a seed file')
    if not options.ref_type:
        parser.error('--ref_type not given. Input either nucl or prot.')
    if not options.fwd:
        parser.error('--fwd not given. Input a reads file') 
    if not options.outFile:
        parser.error('--outFile not given.')
    if os.stat(options.ref_file).st_size ==0:
        parser.error('--ref_file empty. Check file for sequences')
    if os.stat(options.fwd).st_size == 0:
        parser.error('--fwd file is empty. Check file for sequences')
    with open(options.ref_file, 'r') as f:
        if f.readline()[0]!='>':
            parser.error('--ref_file is not FASTA')
    return options 


options = parse_args()

if options.rev:
    print 'Mapping %s and %s to %s file %s' % (options.fwd,options.rev,options.ref_type,options.ref_file)
    paired = 1
else:
    print 'Mapping %s to %s file %s' % (options.fwd,options.ref_type,options.ref_file)
    paired = 0


#######################
## Build blast database
#######################

def makeSeedDB(ref,dbtype):
    myDB = ref
    if '.' in myDB:
        myDB=myDB[:myDB.rfind('.')]
    command = 'makeblastdb -in %s -input_type fasta -dbtype %s -out %s' % (ref,dbtype,myDB)
    print 'Building Blast Database: %s' % command
    os.system(command)
    print 'Building Blast Database: COMPLETE'
    return myDB

####################
## Search 4 Seeds
###################

def blast(myDB):
    if options.ref_type == 'prot':
        blast = 'blastx'
    elif options.ref_type == 'nucl':
        blast = 'blastn'
    input = open(options.fwd,'r')
    print input
    for l in input:
        if l[0] == '>':
            blastfile = options.fwd
        else:
            blastfile = 'FLAGFLAG23592835.fa'
        break
    input.close()
    if 'FLAGFLAG23592835.fa' == blastfile:
        temp = open(options.fwd,'r')
        output = open(blastfile,'w')
        for rec in SeqIO.parse(temp,'fastq'):
            output.write('>'+rec.description+'\n'+str(rec.seq)+'\n')
        temp.close()
        output.close()
    blastOut = 'seeds_%s.%s.txt' % (myDB,blast)
    command = blast + ' -db %s -query %s -evalue %s -outfmt 6 -num_threads %s -out %s' % (myDB, blastfile, options.evalue, options.n_threads, blastOut)
    print 'Executing Blast: %s' % command
    os.system(command)
    print 'Executing Blast: COMPLETE'
    try:
        os.remove('FLAGFLAG23592835.fa')
        print '%s complete' % blast
    except:
        print '%s complete' % blast
    return blastOut

def blast2fasta(bf):
    input = open(bf,'r')
    h=[]
    for l in input:
        t=l.rstrip().split('\t')
        h.append(t[0])
    input.close()
    fasta = options.fwd
    ff = open(fasta,'r')
    for l in ff:
        if l[0] == '>':
            seqType = 'fasta'
        else:
            seqType = 'fastq'
        break
    ff.close()
    ff = open(fasta,'r')
    output = open(options.outFile,'w')
    for rec in SeqIO.parse(ff,seqType):
        if rec.id in set(h):
            output.write('>'+rec.id+'\n'+str(rec.seq)+'\n')
    output.close()
    ff.close()


if __name__ == '__main__':
    db = makeSeedDB(options.ref_file,options.ref_type)
    blastOutputFile = blast(db) 
    if os.stat(blastOutputFile).st_size > 0:
        blast2fasta(blastOutputFile)
    else:
        print 'No significant hits found. Try increasing the evalue'
