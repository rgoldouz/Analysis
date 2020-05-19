# MakeGFAL-COPYTxtFiles 

# Get grid proxy first
#voms-proxy-init --voms cms

## Use this python file to make any gfal-copy file you want
# Options
# File location
# Txt File name
import argparse , os

# set up an argument parser                                                                                                                                                                         
parser = argparse.ArgumentParser()

parser.add_argument('--loc', dest = 'LOC', default= 'srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/' )
parser.add_argument('--name', dest='NAME', default='2016_DYM10to50_ext0')# used to name txt file and eos directory 
parser.add_argument('--dest', dest='EOSDEST', default= '/eos/cms/store/user/asparker/TopLFV/' ) # '/eos/user/a/asparker/TopLFV/')# eos directory to move files to
parser.add_argument('--v', dest='VERBOSE', default=True)
parser.add_argument('--cp', dest='CPTOEOS', default=True) # Move these files to my EOS space

ARGS = parser.parse_args()


# Check that files exist and you can contact their server location using gfal-ls
#print "Test gfal-ls"

testlstxt = "gfal-ls " + ARGS.LOC

os.system( testlstxt )






# Loop over gfal-ls output and add ARGS.LOC to file name so they are in correct format for gfal-copy
#gfal-copy -f --from-file files.txt file://$PWD 
# where files.txt is a file where every line is a source in srm url syntax


# create list of files in SRM URL syntax
txtfileList = []

for fn in filter(None,os.popen( testlstxt  ).read().split('\n')):
#    #print "FILE X : {}".format(fn)
    txtfileList.append( ARGS.LOC + fn )
    if ARGS.VERBOSE :  print "Adding to txt file :  {}".format( ARGS.LOC + fn  ) 

#if ARGS.VERBOSE :
#    print "txtfileList"
#    print txtfileList
 

#'\n'

# Write the list to a text file
# use 1 txt file for  every X files since it tends to fail often
nFilespertxt = 50

txtfiles = []
print "printing txt file names"
for i in xrange(0, int(len(txtfileList)/ nFilespertxt)+1   ):
    txtfiles.append( ARGS.NAME +'_' + str(i) )
    print ARGS.NAME +'_' + str(i) 

# split the txtfileList ito the same number to pieces
chunks = [ txtfileList[nFilespertxt*i:nFilespertxt*(i+1)] for i in range(len(txtfileList)/nFilespertxt + 1)]

# make a list of gfalcopy commands
gfcoms = []

print "Writing contents to txt files"
# open each txt file and write X files or less to the txt file
for it, t in enumerate(txtfiles) :
    gfcomt = "gfal-copy -f --from-file  " + t + ".txt  " + ARGS.EOSDEST+ARGS.NAME + "/"
    gfcoms.append(gfcomt)
    with open(  t + ".txt", "w") as output:
        for f in chunks[it] :
            output.write(str(f) + '\n')

print "Complete list of gfal-copy commands for this sample :"
print gfcoms

#print"gfal-ls output written to file {}".format( ARGS.NAME + '.txt'  )
#myfile.close()




if ARGS.CPTOEOS :
    # Move these files to my EOS space
    #gfal-copy -f --from-file files.txt file://$PWD                                                                                             # where files.txt is a file where every line is a source in srm url syntax 
    lsnew = "ls "+  ARGS.EOSDEST+ARGS.NAME                                      
    isreal = os.system( lsnew  )  
    print  "ls "+ ARGS.EOSDEST+ARGS.NAME 
    #print isreal
    
    if isreal > 0 :  
        print "create eos dir because it doesn't exist: {}".format(ARGS.EOSDEST+ARGS.NAME)
        os.system( "mkdir "+ ARGS.EOSDEST+ARGS.NAME )

        # Copy files to their new home
        for g in gfcoms :
            os.system( g )
        print "Done moving all files from {} to {}".format ( ARGS.LOC , ARGS.EOSDEST+ARGS.NAME  )
    elif isreal == 0 :
        print "eos dir exists, skipping..."
        # Instead we need to check which files have already copied and copy the remainng ones

        # create list of files that have already been copied
        copied = []
        for fn in filter(None,os.popen( lsnew  ).read().split('\n')):
            copied.append( fn )                                                                                                                                                                   

        print "Last file successfully copied was : {}".format( copied[-1] )

        # determine which of split txt files this file is located inside
        for ic,c in enumerate(chunks) :
            if str(copied[-1]) in str(c) :
                #print c
                #print gfcoms[ic-1]
                #print gfcoms[ic]
                print "The following commands must be run for this sample to complete copying :"
                for com in gfcoms[ic:]:
                    print com
                    os.system( com )

print "Finished copying sample!"



