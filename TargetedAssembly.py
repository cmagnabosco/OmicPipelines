################################################ 
## Cara Magnabosco
#### Pipeline for running PRICE
################################################


from Bio import SeqIO
import os
from optparse import OptionParser
import networkx as nx
import sys


############################################  
## PARSE COMMAND LINE ARGUMENTs
############################################  

def parse_args():
    parser = OptionParser(description = 'A pipeline for performing targeted assembly.', epilog="Example:\n python TargetedAssembly.py --seed_file seed.fa --fwd metagenome.fa --insert_size=N --output_file=output.fasta [options]")
    parser.set_defaults(n_threads='1')#,evalue='1e-5',search='N')
## Required options
    parser.add_option('--seed_file', type = 'string',dest='seed_file',help='FASTA (nucleotide or protein) file of sequences to target')
#    parser.add_option('--seed_type',type='choice',choices=['nucl','prot'], dest='seed_type',help='Indicate whether mapping file is a nucleotide (nucl) or protein (prot) file.')
    parser.add_option('--fwd',type='string', dest='fwd', help='FASTA/Q file of the metagenome. For paired end reads, -fwd will indicate the file with forward reads and -rev indicate the reverse reads')
    parser.add_option('--insert_size', type = 'int', dest='insert_size',help='Insert size (Required). If reads are single end (--fwd only), use the max sequence length')
    parser.add_option('--output_file',type = 'string', dest = 'priceOut', help = 'Name for output file. Must end in ".fasta"')
## Other options
    parser.add_option('--rev', dest='rev',help = 'FASTA/Q file that contains the reverse reads of paired end sequencing')
    parser.add_option('--n_threads',type='string',dest='n_threads',help="Number of threads to use [1]")
    parser.add_option('--mate',type='int',help = 'Designate --mate=1 if --fwd and --rev are mate-paired (pointing away from one another) rather than paired-end (pointing towards one another)') 
#    parser.add_option('--evalue',type='string',dest='evalue',help='evalue cut-off to use if --search Y [1e-5]')
#    parser.add_option('--search',choices=('Y','N'), dest = 'search', help='Controls whether or not to use the target_file as the seeds (--search N; DEFAULT) or as a database to search for seeds within the metagenomic data (--search Y) using BLAST')
## Parse args
    (options, args) = parser.parse_args()
## Check for required args
    if not options.seed_file: 
        parser.error('--seed_file not given. Input a seed file')
#    if not options.seed_type:
#        parser.error('--seed_type not given. Input either nucl or prot.')
    if not options.fwd:
        parser.error('--fwd not given. Input a reads file') 
    if not options.insert_size:
        parser.error('--insert_size not provided. Please provide')
    if not options.priceOut:
        parser.error('--output_file not provided. Please provide')
    if os.stat(options.seed_file).st_size ==0:
        parser.error('--seed_file empty. Check file for sequences')    
    if os.stat(options.fwd).st_size ==0:
        parser.error('--fwd empty. Check file for sequences')
    if options.rev:
        if os.stat(options.rev).st_size == 0:
            parser.error('--rev empty. Check file for sequences')
    if options.priceOut[options.priceOut.rfind('.')+1:] != 'fasta':
        parser.error('--output_file not in correct format. See python TargetedAssembly.py -h for help')
    return options 




def PricePaired(readType):
    priceCommand = 'PriceTI -%s %s %s %d 99 -icf %s 1 1 5 -nc 30 -dbmax 72 -mol %d -mpi 99 -MPI 90 -target 95 2 1 1 -a 16 -o %s' % (readType, options.fwd, options.rev, options.insert_size, options.seed_file,options.insert_size/10, options.priceOut)
    os.system(priceCommand)


def PriceSingle():
    priceCommand = 'PriceTI -spf %s %d %d -icf %s 1 1 5 -nc 30 -dbmax 72 -mol %d -mpi 99 -MPI 90 -target 95 2 1 1 -a 16 -o %s' % (options.fwd, (options.insert_size/2+options.insert_size/4),options.insert_size, options.seed_file, options.insert_size/10, options.priceOut)
    os.system(priceCommand)



def run_assembly_pipeline():
## Execute PRICE PROTOCOL depending on whether paired end or not
    if paired == 1:
        if options.mate: # If reads are mate-paired require mpp
            PricePaired('mpp')
        else:
            PricePaired('fpp')
    if paired == 0:
        PriceSingle()


#def makeBlastdb(file):
#    db = file[:file.rfind('.')]
#    os.system('makeblastdb -in %s -input_type fasta -dbtype nucl -out %s' % (file,db))
#    return db

#def blastn(file, db):
#    out = file[:file.rfind('.')]+'.txt'
#    os.system('blastn -db %s -query %s -evalue 1e-30 -outfmt 6 -out %s' % (db,file,out))
#    return out

#def findOverlaps(file,fasta):
#    eDges = []
#    input = open(file,'r')
#    for l in input:
#        t = l.rstrip().split('\t')
#        if float(t[2])>99:
#            if t[0] != t[1]:
#                eDges.append((t[0],t[1]))
#        G=nx.Graph()
#        G.add_edges_from(eDges)
#        print nx.number_connected_components(G)
#        for i,g in enumerate(nx.connected_components(G)):
#            # g is a set of headers that were connected                       
#            print i,g
#            input = open(fasta,'r')
#            output = open('Clustered_Seqs_'+str(i+1)+'.fasta','w')
#            for rec in SeqIO.parse(input,'fasta'):
#                print rec.id
#                if rec.id in g:
#                    print rec.id
#                    SeqIO.write(rec,output,'fasta')
#            output.close()
#            input.close()
#        input = open(fasta,'r')
#        header=[]
#        for rec in SeqIO.parse(input,'fasta'):
#            header.append(rec.id)
#            input.close()
#            missed = list(set(header)-set(g))
#            input = open(fasta,'r')
#            output = open('Seqs_not_clustered.faa','w')
#            for rec in SeqIO.parse(input,'fasta'):
#                if rec.id in set(missed):
#                    SeqIO.write(rec,output,'fasta')
#            output.close()
#            input.close()

#def merge_assemblies(tmpF):
#    myDB = makeBlastdb(tmpF)
#    blastOutput = blastn(tmpF,myDB)
#    findOverlaps(blastOutput,tmpF)


options = parse_args()
if options.rev:
    print 'Targeted Assembly of %s and %s based on %s' % (options.fwd,options.rev,options.seed_file)
    paired = 1
else:
    print 'Targeted Assembly of %s based on %s' % (options.fwd,options.seed_file)
    paired = 0



if __name__ == '__main__':
    run_assembly_pipeline()
    myMerge = options.priceOut.replace('.fasta','.cycle30.fasta')
    print myMerge
    import MergeAssemblies
    MergeAssemblies.merge_assemblies(myMerge)

