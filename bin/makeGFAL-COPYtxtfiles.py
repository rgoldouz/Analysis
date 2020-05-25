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
 

# Write the list to a text file
# use 1 txt file for  every X files since it tends to fail often
nFilespertxt = 50

txtfiles = []
print "printing txt file names"
for i in xrange(0, int(len(txtfileList)/ nFilespertxt)+1   ):
    txtfiles.append( ARGS.NAME +'_' + str(i) )
    print ARGS.NAME +'_' + str(i) 

# split the txtfileList into the same number to pieces
chunks = [ txtfileList[nFilespertxt*i:nFilespertxt*(i+1)] for i in range(len(txtfileList)/nFilespertxt + 1)]


## Now decide where to copy the files before creating the copy commands


##Check CERNbox EOS for existing files
os.system( "mkdir "+ "/eos/user/a/asparker/TopLFV/"+ ARGS.NAME )
lsnew1 = "ls " +'/eos/user/a/asparker/TopLFV/'+ ARGS.NAME
copied1 = []


for fn in filter(None,os.popen( lsnew1  ).read().split('\n')):
    copied1.append( fn ) 
if len(copied1) > 0 :
    copiedfile =  copied1[-1]
    print "Last file successfully copied was : {}".format( copied1[-1] )
    ARGS.EOSDEST = '/eos/user/a/asparker/TopLFV/'
    print ARGS.EOSDEST
else:
    print "None of the files were copied to CERNbox yet, Copy to CMS EOS"

    # Check CMS EOS for existing files
    os.system( "mkdir "+ ARGS.EOSDEST+ARGS.NAME )    
    ## CMS EOS
    lsnew = "ls "+  ARGS.EOSDEST+ARGS.NAME                                      
    copied = []

    print "Check which files were already copied to CMS EOS"

    for fn in filter(None,os.popen( lsnew  ).read().split('\n')):
        copied.append( fn ) 
    if len(copied) > 0 :
        copiedfile =  copied[-1]
        print "Last file successfully copied was : {}".format( copied[-1] )
        print ARGS.EOSDEST
    else:
        #print "None of the files were copied to CMS EOS yet, Copy entire  sample to CMS EOS"
        copiedfile = 'Def_not_a_file_name_in_real_list'




# make a list of gfalcopy commands now that we have decided where to copy to
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

remainingCommands = []
# Set remaining commands 
if copiedfile == 'Def_not_a_file_name_in_real_list':
    remainingCommands  = gfcoms




if ARGS.CPTOEOS :

    # Copy everything
    if remainingCommands == gfcoms :
        print "Copy everything..."
        print "The following commands must be run for this sample to complete copying :" 
        for com in gfcoms:
            print com
            #print "If this looks right then turn on the os.system( com ) commands and run again"
            os.system( com ) #test without this turned on
    # Only some of the commands must be redone
    else :

        print " determine which of split txt files the last copied file is located inside"
        for ic,c in enumerate(chunks) :
            if str(copiedfile) in str(c) :
                print copiedfile
                print c
                print c[-1]
 
                if str(copiedfile) in str(c[-1]):
                    print "Last file copied was the last one we needed, No need to rerun copy commands"
                    remainingCommands = []
                else :
                    print "The following commands must be run for this sample to complete copying :"
                    for com in gfcoms[ic:]  :
                        print com
                        remainingCommands.append(com)
                        #print "If this looks right then turn on the os.system( com ) commands and run again"
                        os.system( com ) #test without this turned on
                
            
#    # If files are on CERNbox(newdir == True )   AND not all are done copying                                                                                                                                 
#    #then copy entire sample to CMS EOS (Because my CERNbox is full)  
#    if (newdir == True and len(remainingCommands) > 0  and remainingCommands != gfcoms   ):
#        print "If files are on CERNbox  AND not all are done copying then copy entire sample to CMS EOS (Because my CERNbox is full)  "
#        #os.system( "mkdir "+ ARGS.EOSDEST+ARGS.NAME )
#        print "Created directory : "
#        print ARGS.EOSDEST+ARGS.NAME
#        for com in remainingCommands:
#            print com
#            print "If this looks right then turn on the os.system( com ) commands and run again"
#            #os.system( com ) #test without this turned on  

print "Finished copying sample!"



