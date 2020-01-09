#!/usr/bin/env python
#$ -S /usr/bin/python
"""
Run (almost) any command in parallel across a cluster (using SGE or SLURM) or on a single server. The only requirement is that the input file can be broken into pieces (using regex or line counts). Many bioinformatics file types (fasta, fastq, gbk) are implicitly understood. The resulting output fragments will simply be concatenated.

If the command to run is one of the currently recognized commands: (blastall, rpsblast, softbery, meta_rna) you can simply give the command as you would run it locally. EG:
    batch_launcher.py -- blastall -i some_file.fasta -d some_db -p blastx -o some_file.vs.db.blastx

If the script does not recognize your command, you can tell it where to find the input and output files by indicating with the -i and -o flags, where in your command these files are named:
    batch_launcher.py -i -i -o -o -- blastall -i some_file.fasta -d some_db -p blastx -o some_file.vs.db.blastx

The -i and -o flags (in addition to accepting flags) will also take a number indicating which element of the command has the file name. The executable or script is element number 0 of the command. ("-i 2 -o 8" would have worked in the above example)

The -i flag can be omitted or set to '%stdin' to read from standard input. (You must supply the input type and the chunk size!)

The -o flag can be ommitted to print to stdout (if --wait is used) or set to one of the following:
    %stdout (must use --wait)
    %e.xxx (replace the extension of the input file with '.xxx')
    %a.xxx (add '.xxx' to the input file name)
    %E.xxx and %A.xxx (same as above, but strip path from inputfile, and write in ./)
    %s/foo/bar/ and %S/foo/bar/ (apply a regular expression sub to input file(name))

Multiple -o flags can be given if multiple output files are created.

Additionally, if any output file is a directory, the final output will be a directory containing output directories created by each sub-task. There will be no merging done. 

The script will try to recognize the record type of the input file based on its extension. If it can't you will have to either supply a regex pattern to split the file on (-p "^>" for fasta, -p "^LOCUS" for gbk, etc) or manually set the file type (eg: -T fasta).
"""
version="batch_launcher 1.0.0"
SGE='sge'
SLURM='slurm'
LOCAL='local'
AUTO='auto'
queueBinaryMap={SGE:'qsub',
                SLURM:'sbatch'}

import argparse
import sys, re, logging, os, tempfile, subprocess, shutil, traceback, datetime, shlex, time, fileinput, io, threading, queue
from math import ceil

###########
# edl.util (selected functions)
#import re, logging, sys, os
#from edl.expressions import accessionRE
#logger=logging.getLogger(__name__)
##########
def countBasesInFasta(fastaFile):
    """
    Given a fasta file, return a dict where the number of records and the total number of bases are given by 'records' and 'bases' respectively.
    """
    recordRE=re.compile(r'^>')
    whiteSpaceRE=re.compile(r'\s+')
    totalBases=0
    totalSeqs=0
    with open(fastaFile) as f:
        for line in f:
            if recordRE.match(line):
                totalSeqs+=1
                continue
            totalBases+=len(whiteSpaceRE.sub('',line))

urlRE=re.compile(r'[a-z]+\:\/\/')
def openInputFile(infile, *args):
    """
    return input stream. Allow for text, gzipped text, or standard input if None given.
    """
    if infile is None:
        logging.info("Reading input from STDIN")
        return sys.stdin

    if isinstance(infile, str):
        if urlRE.match(infile):
            import urllib2
            return urllib2.urlopen(infile)
        if len(infile)>3 and infile[-3:]=='.gz':
            import gzip
            return gzip.GzipFile(infile,'rb')
        elif len(infile)>4 and infile[-4:]=='.bz2':
            import bz2
            return bz2.BZ2File(infile,'rb')
        else:
            return open(infile,'rt')
    else:
        return infile

############
# edl.batch
#import tempfile, re, os, logging, sys
#import numpy as np
#from edl.util import openInputFile
############

def checkTmpDir(tmpDir,jobName):
    """
    Make sure tmp dir is empty.
    Create it if necessary.
    If name is None, create dir name based on job name
    """
    if tmpDir is None or tmpDir=="":
        tmpDir=tempfile.mkdtemp(suffix=jobName,dir=".")
    else:
        if os.path.exists(tmpDir):
            raise Exception("Temporary directory already exists! (%s)" % (tmpDir))
        else:
            os.makedirs(tmpDir)

    logging.debug("Created temporary directory: %s" % tmpDir)
    return tmpDir

def create_parent_dir(path):
    """ create the parent dir for path if it doesn't already exist """
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

def getFileType(options, infile):
    """
    figure out how to split up the input file and return a FileType object with the necessary info.

    The options paramter is assumed to be an object with the following values:
        .infileType (None or a string indicating the input file type)
        .pattern (None or a regex pattern to find a record boundary)
        .numLines (None or an integer number of lines per record)
    if all of the above are None, the extension of the name in infile is used to guess
    """
    if options.pattern is not None:
        if options.numLines is not None:
            logging.error("Cannot use BOTH pattern and number of lines to split records!")
            raise Exception("Conflicting options: pattern and numLines")


    # if there is a file type, use that (And override if other options present)
    if options.infileType is not None:
        fileType = fileTypeMap[options.infileType]
        if fileType.sepRE is not None:
            logging.info('Using user chosen "%s" pattern to split files: /%s/' % (options.infileType, fileType.sepRE.pattern))
        else:
            logging.info('Using user chosen "%s" number of lines (%d) to split files' % (options.infileType, fileType.numLines))
    else:
        fileType=None

    if options.pattern is not None:
        sepRE = re.compile(options.pattern)
        logging.info('Using user supplied pattern to split files: /%s/' % (options.pattern))
        if fileType is None:
            fileType=FragmentableFileType(name=options.pattern,sepRE=sepRE)
        else:
            fileType.sepRE=sepRE
            fileType.numLines=None
    elif options.numLines is not None:
        logging.info('Using user supplied number of lines to split files: /%s/' % (options.numLines))
        if fileType is None:
            fileType=FragmentableFileType(str(options.numLines), numLines=options.numLines)
        else:
            logging.info("Overriding record line count (%d)!" % (options.numLines))
            fileType.sepRE=None
            fileType.numLines=options.numLines
    elif fileType is None:
        # try to guess from the extension
        if infile is None:
            sys.exit("We cannot infer the record separator from the filename when using standar input. Please specify the input file type with -T or a pattern to match the first line with -p")
        fileType = getTypeFromFileName(infile)
        if fileType.sepRE is not None:
            logging.info('File looks like %s, using pattern to split files: /%s/' % (fileType.name, fileType.sepRE.pattern))
        else:
            logging.info('File looks like %s, using number of lines (%d) to split files' % (fileType.name, fileType.numLines))
    return fileType

def getTypeFromFileName(filename):
    """
    find the file's extension in the fileExtensionMap
    return the separator RE for that file type

    If file not recognized, inform user how to set the pattern manually
    """
    ext = os.path.splitext(filename)[1]
    fileType = fileExtensionMap.get(ext,None)
    if fileType is not None:
        return fileType
    else:
        logging.warning("""Cannot guess type of %s!

The input file does not have a recognized extension:
%s

You must manually set the file type with '-T TYPE' where TYPE is in:
%s
OR set the record separator pattern with '-p PATTERN' where PATTERN is a regular expression. E.G. '^>" for fasta, or '^LOCUS' for gbk.
""" % (filename, list(fileExtensionMap.keys()), list(fileExtensionMap.keys())))
        sys.exit(65)

def getSizePerChunk(infile, splits, fileType, splitOnSize=False):
    """
    Get total size of all records and return target size for each chunk to end up with number of chunks specified by 'splits'
    """
    if infile is None:
        raise Exception("We cannot determine chunk size from STDIN!")

    if splitOnSize:
        # get a custom function that returns the size of this type of record
        recordSizer=fileType.sizer
    else:
        # just return 1 for each record
        recordSizer=recordCounter

    # loop through records
    inhandle = openInputFile(infile)
    totalSize = 0
    for record in fileType.recordStreamer(inhandle):
        totalSize+=recordSizer(record)
    inhandle.close()

    return calculateChunkSize(totalSize,splits)

def calculateChunkSize(size,splits):
    """
    how big should the fragments be?
    """
    chunk = int(ceil(size/float(splits)))
    logging.info("Setting chunk to: %d=ceil(%d/%d)" % (chunk,size,splits))
    return chunk

# Simply count records
recordCounter=lambda x: 1

def defaultRecordSizer(recordLines):
    """
    return the total number of characters in the record text
    """
    size=0
    for line in recordLines:
        size+=len(line)
    return size

def fastaRecordSizer(recordLines):
    """
    Returns the number of charcters in every line excluding:
        the first (header) line
        whitespace at the start and end of lines
    """
    size=0
    for i in range(1,len(recordLines)):
        size+=len(recordLines[i].strip())
    return size

def fastqRecordSizer(recordLines):
    """
    Returns the number of charaters in the lines between the sequence header (@) and the quality header (@) excluding whitespace at the start and end of lines
    """
    size=0
    for i in range(1,len(recordLines)):
        line=recordLines[i]
        if len(line)>0 and line[0]=='+':
            return size
        size+=len(line.strip())

    logging.warning("Did not find quality line in record: %s" % recordLines)
    return size

def gbRecordSizer(recordLines):
    """
    uses biopython to parse record
    """
    from Bio import SeqIO
    record = SeqIO.read(recordLines,'genbank')
    return len(record)

def fragmentInput(infile, options, frag_base):
    """
    Wraps the following methods into one:
        getFileType(options, infile)
        getSizePerChunk(infile, options.splits, fileType, splitOnSize=options.splitOnSize)
        fragmentInputBySize(infile, tmpdir, options.chunk, fileType, frag_base,   splitOnSize=options.splitOnSize, suffix=)
    """
    fileType=getFileType(options, infile)
    # if we had to get from file name, save:
    if options.infileType is None and options.pattern is None and options.numLines is None:
        options.infileType = fileType.name

    dont_fragment=False
    if options.chunk is None:
        if options.splits is None:
            options.splits=DEFAULT_SPLITS
        if options.splits==1:
            # Don't split if one fragment requested
            dont_fragment=True
        options.chunk=getSizePerChunk(infile, options.splits, fileType, splitOnSize=options.splitOnSize)

    # Are we going to save fragments for later reuse?
    if options.keepFragments:
        options.frag_dir = get_resuable_fragment_dir(infile, options.chunk,
                                                     options.splitOnSize)
        if options.fragBase is None:
            options.fragBase = 'chunk'
        suffix = os.path.splitext(infile)[1]
        if options.fragSuff is None:
            options.fragSuff = suffix
        num = count_existing_fragments(options.frag_dir, options.fragBase,
                                       options.fragSuff, infile)
        if num > 0:
            # fragments already exist, let's continue
            return num
    else:
        options.frag_dir = options.tmpDir
        if options.fragBase is None:
            options.fragBase = "file.fragment"
        suffix = '.in'
        if options.fragSuff is None:
            options.fragSuff = suffix

    # fragment input file
    if dont_fragment:
        # Just symlink to tmp_dir
        logging.debug("Since only one fragment requested, entire source file will be linked to tmp_dir")
        link_source_into_tmp(infile, options.frag_dir, options.fragBase,
                             options.fragSuff)
        options.keepFragments = False
        return 1
    else:
        return fragmentInputBySize(infile, 
                                   options.frag_dir, 
                                   options.chunk, 
                                   fileType, 
                                   options.fragBase, 
                                   splitOnSize=options.splitOnSize,
                                   suffix=options.fragSuff)

def get_resuable_fragment_dir(infile, chunk, splitOnSize):
    """
    ge consitent place to put re-usable file fragments. Something like:
        # {dirname(infile)}/.batch_launcher/{basename(infile)}/chunk_size
    """
    fdir = os.path.join(os.path.dirname(infile),
                        '.batch_launcher',
                        os.path.splitext(os.path.basename(infile))[0],
                        'chunk_{}{}'.format(chunk, 
                                            'size' if splitOnSize else ''))
    if not os.path.exists(fdir):
        os.makedirs(fdir)
    elif not os.path.isdir:
        raise Exception("Desired fragment dir ({}) alredy exists as a file!  "
                        "Rename your input file or omit the -K "
                        "option").format(fdir)
    return fdir


def link_source_into_tmp(infile, tmpdir, frag_base,
                          suffix='.in'):
    """
    create symlink in tmp dir that points to input file
    """
    logging.debug("Linking input input: %r" % ({'infile':infile,
                                                'tmpDir':tmpdir,
                                                'base':frag_base}))
    tmp_file_name=getFragmentPath(tmpdir, frag_base, 1, suffix)
    os.link(infile, tmp_file_name)
    return 1

def count_existing_fragments(folder, frag_base, suffix, input_file):
    """
    Return count of fragments newer than input, return -1
    if input is newer than any
    """
    input_file_time = os.path.getmtime(input_file)
    num = 0
    count = 0
    logging.info("Looking for existing fragments in %s", folder)
    while True:
        num += 1
        tmp_file_name=getFragmentPath(folder, frag_base, num, suffix)
        logging.debug("checkking %s", tmp_file_name)
        if os.path.exists(tmp_file_name):
            count += 1
            if os.path.getmtime(tmp_file_name) < input_file_time:
                # if any file is older, return -1
                return -1
        else:
            # at first missing file, return count
            return count
        
def fragmentInputBySize(infile, tmpdir, chunk, fileType, frag_base, splitOnSize=True, suffix='.in'):
    """
    Break up input into files of size chunk in tmpdir. Return number of fragments.
    """
    logging.debug("Fragmenting input: %r" %
                  ({'infile':infile,'tmpDir':tmpdir,'chunk':chunk,'base':frag_base}))
    inhandle = openInputFile(infile)
    num = fragmentInputStreamBySize(inhandle, tmpdir, chunk, fileType,
                                    frag_base, splitOnSize=splitOnSize, suffix=suffix)
    if infile is not None:
        inhandle.close()
    return num

def fragmentInputStreamBySize(inhandle, tmpdir, chunk, fileType, frag_base, splitOnSize=True, suffix='.in'):
    if splitOnSize:
        # get a custom function that returns the size of this type of record
        recordSizer=fileType.sizer
    else:
        # just return 1 for each record
        recordSizer=lambda x: 1

    count=0
    num=1
    tmp_file_name=getFragmentPath(tmpdir, frag_base, num, suffix)
    #logging.debug('Writing fragment (%d,%d,%d): %s' % (chunk,count,num,tmpFileName))
    tmpFile = open(tmp_file_name, 'w')
    for record in fileType.recordStreamer(inhandle):
        recordSize=recordSizer(record)
        count+=recordSize

        if count>chunk:
            # close previous chunk and open new one
            tmpFile.close
            num+=1
            tmp_file_name=getFragmentPath(tmpdir, frag_base, num, suffix)
            #logging.debug('Writing fragment (%d,%d,%d): %s' % (chunk,count,num,tmpFileName))
            tmpFile = open(tmp_file_name, 'w')
            count=recordSize

        # write record
        tmpFile.writelines(record)

    tmpFile.close()

    return num

def getFragmentPath(directory, base, index, suffix='.in'):
    return "%s%s%s" % (directory, os.sep, getFragmentName(base,index,suffix=suffix))

def getFragmentName(base, index, suffix='.in'):
    return "%s%05d%s" % (base, index, suffix)

def getFragmentPrefix(base,index):
    return getFragmentName(base,index,suffix='.pre')

def formatCommand(command):
    """
    given a list of command elements, print string approximating what you'd type at a shell prompt
    """
    cmdstr=""
    logging.debug(repr(command))
    for arg in command:
        if " " in arg:
            cmdstr=cmdstr+" \""+arg+"\""
        else:
            cmdstr=cmdstr+" "+arg
    return cmdstr

def regexRecordGenerator(fileType, stream):
    """
    Using the sepRE setting in fileType, break the input stream into records
    """
    lastRecord=[]
    for line in stream:
        if fileType.sepRE.match(line):
            if len(lastRecord)>0:
                yield lastRecord
                del lastRecord[:]
        lastRecord.append(line)

    if len(lastRecord)>0:
        yield lastRecord

def linedRecordGenerator(fileType, stream):
    """
    Using the numLines setting in fileType, break the input stream into records
    """
    lastRecord=[]
    for index,line in enumerate(stream):
        if index%fileType.numLines==0:
            if len(lastRecord)>0:
                yield lastRecord
                del lastRecord[:]
        lastRecord.append(line)

    if len(lastRecord)>0:
        yield lastRecord

def addFragmentingOptions(parser, defaults={"splits":400}):
    option_group = parser.add_argument_group('Fragmenting Options')

    option_group.add_argument("-L", "--recordLines", metavar="NUMLINES",
                              dest='numLines', default=None, type=int,
                       help="Number of lines per record")
    option_group.add_argument("-P", "--pattern", metavar="PATTERN", dest='pattern', default=None,
                       help="Regular expression to split records")
    choices = list(fileTypeMap.keys())
    option_group.add_argument("-T","--infileType", dest='infileType', default=None,
                              help='Type of input file. Otherwise, choose by '
                              'extension. Known types are: ' + 
                              ','.join(choices))
    option_group.add_argument("-C", "--chunkSize", type=int, dest='chunk',  metavar="FRAG_SIZE",
                      help="The number of records per fragment. Overrides NUM_FRAGS")
    default=defaults.get("splits",None)
    option_group.add_argument("-N", "--numChunks", dest='splits', type=int, metavar="NUM_FRAGS", default=default,
                      help="The number of fragments to create (defaults to %s)"
                             % (default,))
    option_group.add_argument("-s", "--splitOnSize",  default=False, action='store_true',
                      help="create chunks based on record size, not number of records. For known sequence types (fasta, fastq, gb), sequence length is used, otherwize the full size of the record text is counted")
    option_group.add_argument("-K", "--keepFragments", default=False,
                            action='store_true',
                            help=("keep fragments for later re-use (useful for"
                                  "HMMs)"))
#################
# Classes
#################
class FragmentableFileType:
    def __init__(self, name, sepRE=None, numLines=None, sizer=None):
        self.name=name
        self.sepRE=sepRE
        if sepRE is not None:
            self._recordStreamer=regexRecordGenerator
        self.numLines=numLines
        if numLines is not None:
            self._recordStreamer=linedRecordGenerator
        if sizer==None:
            self.sizer=defaultRecordSizer
        else:
            self.sizer=sizer

    def recordStreamer(self, stream):
        return self._recordStreamer(self, stream)

##################
# Constants from edl.batch
##################
FASTA=FragmentableFileType('fasta',sizer=fastaRecordSizer,sepRE=re.compile(r'^>(\S+)'))
FASTQ=FragmentableFileType('fastq',sizer=fastqRecordSizer,numLines=4)
GENBANK=FragmentableFileType('gb',sizer=gbRecordSizer,sepRE=re.compile(r'^LOCUS'))
TABLE=FragmentableFileType('table',sepRE=re.compile(r'^'))
HMM=FragmentableFileType('hmm',sepRE=re.compile(r'^HMMER'))
fileTypeMap={FASTA.name:FASTA,
             FASTQ.name:FASTQ,
             GENBANK.name:GENBANK,
             TABLE.name:TABLE,
             HMM.name:HMM,
            }
fileExtensionMap={'.fa':FASTA,
                  '.fna':FASTA,
                  '.faa':FASTA,
                  '.ffn':FASTA,
                  '.fasta':FASTA,
                  '.fastq':FASTQ,
                  '.gb':GENBANK,
                  '.gbk':GENBANK,
                  '.gbff':GENBANK,
                  '.gpff':GENBANK,
                  '.tab':TABLE,
                  '.tsv':TABLE,
                  '.csv':TABLE,
                  '.m8':TABLE,
                  '.hmm':HMM,
                 }

def main():
    ## set up CLI
    usage = "usage: %prog [options] -- [command to run on cluster]"
    description = """
Given an input file with multiple records, an output file, and a command: process the records in parallel across a cluster or multiprocesser server.
    """
    parser = argparse.ArgumentParser(usage, description=description)
    parser.add_argument("--version", dest="version", default=False,
                    action="store_true", help="print version")
    parser.add_argument("-A", "--about",
                  action="store_true", dest="about", default=False,
                  help="Print description")

    # logging
    option_group = parser.add_argument_group('Logging Options')
    option_group.add_argument("-l", "--logToFile", default=False, action="store_true",
                      help="print logging messages to file")
    option_group.add_argument("-v", "--verbose",
                      action="count", dest="verbose", default=1,
                      help="Print log messages. Use twice for debugging")
    option_group.add_argument("-q", '--quiet', dest='verbose',
                      action="store_const", const=0,
                      help="Suppress warnings. Only print fatal messages")
    option_group.add_argument("-V","--loglevel", type=int, dest="verbose",
                     help="Shortcut to set verbosity directly")
    # options for processing command
    option_group = parser.add_argument_group('Command Parsing Options')
    option_group.add_argument("-G", "--ignore", default=False, action='store_true',
    help="Do not perform checks based on the program name (in the command). Batch_launcher will need you to indicate input/output/threads with -i,-o, and -t, and it will skip any program specific processing (e.g. db checks for blast)")
    option_group.add_argument("-i", "--inputFlag", metavar="FLAG",
                      help="Indicate where in command to find the input file to fragment. The value can be an option flag (eg '-i' or '-fasta') or a number indicating a positional argument. Only needed if program in command is not recognized.")
    option_group.add_argument("-o", "--outputFlags", action='append',metavar='FLAG',
                      help="Indicate where in command (or not) to find the "
                           "output file. The value can be an option flag or "
                           "positional argument (as with -i, above), a file "
                           "name (if the command always produces the same "
                           "output file), '%%e.EXT' if the output file "
                           "is the input file with a new extension (e.g. "
                           "%%e.faa), '%%a.ext' if the output file is the "
                           "input fie with an additional extension, or "
                           "'%%s/foo/bar/' if the output file is the input "
                           "file with a regexp substitution. Multiple values "
                           "are permitted to account for multiple outputs.")
    option_group.add_argument("-p", "--prefixFlag",
                      help="Indicate where in command to find a prefix for "
                           "output files. Output extension flags '%%p.ext' "
                           "indicate that the output file will be PREFIX.ext.")
    option_group.add_argument("-t", "--threadsFlag", metavar="FLAG",
                      help="Option to use to tell command how many threads to use")
    option_group.add_argument("-H", "--relative_path", metavar='FLAG',
                            action='append', dest='rel_paths',
                            help="Use the same syntax as -o or -i to indicate " 
                            "elements in the command that are relative "
                            "paths. These will need to be made absolute when "
                            "the command is run in a temp folder.")
    option_group.add_argument("--frag_prep", metavar='COMMAND', default=None,
                            help="command to run on fragments. Can be used "
                            "to fragment databases instead of queries. "
                            " Include place holder '{}' for file name. "
                            "Can optionnally include output suffix to "
                            "prevent for unneeded work when re-using "
                            "fragments. EG: "
                            "'.h3i: hmmpress {}' for fragmenting an HMM db")

    # ways to customize batch behavior
    addFragmentingOptions(parser)

    # way to customize execution
    option_group = parser.add_argument_group('Fragment Execution Options')

    option_group.add_argument("-X", "--queue", default=AUTO, 
                      choices=[SGE,SLURM,LOCAL,AUTO],
                      help="What type of scheduler to use: 'sge', 'slurm', or 'local' (just use threads and command-line). By default, it will check for the qsub (from sge)  and sbatch (from slurm) binaries in the PATH and go with the first one found.")
    option_group.add_argument("-R", "--throttle", default=0, type=int,
            help="Limit number of simultaneously executing fragemnts to given number. Default: 0 => unlimited")
    option_group.add_argument("-w", "--wait", default=False, action='store_true',
                      help="Wait for jobs to finish before exiting script (only when using SGE)")
    option_group.add_argument("-c","--cwd",dest='cwd',
                      default=False,action='store_true',
                      help='Run fragment in current directory, otherwise it will be executed in the tmp location (as configured in Python, usu a local disk like /tmp) of the node (which is good for programs like softnerry that create lots of temporary files)')
    option_group.add_argument("-r","--retries",default=-1, type=int,
                      help='number of times to resubmit failed tasks. Less than zero (the default) means continue as long as some new results come through each time')

    option_group = parser.add_argument_group('Advanced Fragment Execution Options')
    option_group.add_argument("-j", "--jobName", metavar="STRING",
                      help="String for naming queued tasks")
    option_group.add_argument("-d", "--tmp_dir", dest='tmpDir',
                      help="Temporary directory for files")
    option_group.add_argument("-S", "--submitOptions", default=None, dest='sgeOptions',
                      help="Option string to add to SGE or SLURM command")
    option_group.add_argument("-n", "--priority", default=0, type=int,
                      help="Adjust priority of sub-tasks, only applies to SGE")

    # primarily used when calling self
    option_group = parser.add_argument_group('Internal Use Only')
    option_group.add_argument("-f", "--frag_base", dest="fragBase", default=None,
                     help=("naming base for input file fragments"
                           "('file.fragment')"))
    option_group.add_argument("--frag_dir", dest="frag_dir", default=None,
                     help=("folder with input fragments"))
    option_group.add_argument("--frag_suffix", dest="fragSuff", default=None,
                     help=("naming suffix for input file fragments"))
    option_group.add_argument("-m", "--mode", default='launch', metavar='MODE',
                      choices=['launch','run','cleanup'],
                      help="Only used internally to lauch tasks in SGE or SLURM")
    option_group.add_argument("-Z","--taskType",default=None, choices=list(taskTypePatterns.keys()),
                      help="only for task running: what type of task is this? Will be ignored in inital call")

    parser.add_argument("command", nargs="*")

    options = parser.parse_args()
    cmdargs = options.command

    if options.version:
        print (version)
        exit(0)

    if options.about:
        print (description)
        exit(0)

    # if running launcher or cleanup, standard logging
    if options.verbose==0:
        loglevel=logging.ERROR
    elif options.verbose==1:
        loglevel=logging.WARN
    elif options.verbose==2:
        loglevel=logging.INFO
    elif options.verbose>=3:
        loglevel=logging.DEBUG

    if options.logToFile and options.mode != 'run':
        logStream=io.StringIO()
    else:
        logStream=sys.stderr
    logging.basicConfig(stream=logStream, level=loglevel)

    # If running a sub-job, do so now
    if options.mode=='run':
        exitcode=runFragment(options, cmdargs)
        sys.exit(exitcode)

    # cleanup only?
    if options.mode=='cleanup':
        logging.warning("Starting cleanup at: %s" % (datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S UTC %Y")))
        exitcode=cleanup(options, cmdargs)
        sys.exit(exitcode)

    # get SGE options from user
    sgeOptions=[]
    if options.sgeOptions is not None:
        sgeOptions.extend(shlex.split(options.sgeOptions))

    # figure out what the command is and apply some defaults
    if not options.ignore:
        # Inspect command, get any special options for task runner and update sgeOptions
        options.taskType = checkCommand(cmdargs,sgeOptions,options)
        logging.debug('Task type is: %s' % (str(options.taskType)))

        # if any flag not supplied by user, check for defaults based on taskType
        # first, there must be a task type:
        #if options.taskType is None:
        #    if (options.inputFlag is None or not(options.outputFlags)):
        #        parser.error("I'm sorry the command is not recognized. Please use the -i and -o options to specify the input and output flags in your command.")

    else:
        logging.warning('User chose to run command without modification')

    # cerate job name if not supplied by user
    if options.jobName is None:
        options.jobName = getJobName(cmdargs[0])

    # use user input and taskType to get input file and setup command
    (infile,outfile) = prepareCommandForBatch(cmdargs, options)

    # now that we have an output file, redirect logs if asked
    if options.logToFile:
        if outfile is None:
            logFile="%s.launch.err" % options.jobName
        else:
            logFile="%s.launch.err" % outfile
        newLogStream=open(logFile,'w')
        sys.stderr.write("Logging to file: %s\n" % logFile)
        newLogStream.writelines(logStream.getvalue())
        # get old handlers
        oldHandlers = list(logging.getLogger().handlers)
        # add new handler
        logging.getLogger().addHandler(logging.StreamHandler(newLogStream))
        # remove old handlers
        for handler in oldHandlers:
            logging.getLogger().removeHandler(handler)
        logStream.close()
        logStream=newLogStream

    # allow input from STDIN only if chunk given
    if infile is None and options.chunk is None:
        parser.error("We can use STDIN only if a chunk size is given!")

    if options.queue != LOCAL:
        if options.queue == AUTO:
            # auto-detect
            options.queue = lookForQueueingCommands()

        options.sgeOptions = sgeOptions
        logging.debug("Adding the following submit options to task array command: '%s'" % sgeOptions)

        # must use 'wait' to print to stdout
        if not(options.outputFlags) or '%stdout' in options.outputFlags:
            if not(options.wait):
                parser.error("Piping output to stdout only works if we wait for the jobs to finish. Please use wait (-w) or fix your output flags (-o)")

    # get temporary dir
    options.tmpDir = checkTmpDir(options.tmpDir,options.jobName)

    # create fragmented version of input in temp folder
    options.splits=fragmentInput(infile, options, options.fragBase)

    # launch task array
    launchJobs(options,cmdargs,errStream=logStream)

    # launch cleanup
    if options.wait or options.queue == LOCAL:
        # if wait, then we already waited for jobs to finish, just run cleanup
        if options.verbose<=3:
            exitcode=cleanup(options, cmdargs)
        else:
            logging.warning("Debugging level is high, so skipping cleanup and leaving temporary files in place")
            sys.exit(0)
        logging.debug("Done!")
        sys.exit(exitcode)
    else:
        # otherwise, launch cleanup task to wait for processing to finish
        launchCleanup(options, cmdargs, errStream=logStream)

def lookForQueueingCommands():
    """
    Look for qsub and sbatch in the execution path and return the first found queueing system
    """
    for queue, binary in queueBinaryMap.items():
        if checkForBinary(binary):
            return queue
    else:
        raise Exception("Cannot locate a queueing system. None of these executables were found in your PATH: %s" % (queueBinaryMap.values(),))

def checkForBinary(binary):
    """
    Look for binary with `which`.
    """
    try:
        fullPath = subprocess.check_output(['which',binary])
        return True
    except subprocess.CalledProcessError as e:
        return False

def getJobName(programPath=None):
    """
    generate a job name for task array
    """
    if programPath is None:
        return "Bat_%d" % os.getpid()
    else:
        return "Bat_%s" % (os.path.split(programPath)[1])

def prepareCommandForBatch(cmdargs,options):
    """
    If (in/out/thread) flags not defined in options, set from taskType defaults
    Scan command for input flags to get input file name
    return input file name (or None for stdin)

    If input/output files are set with positional arguments by default, translate
    into absolute position in command string.

    if no prefix flag used in command, remove any %p.ext flags from outputFlags.
    """

    infile=None
    prefix=None
    outfile=None
    # remember user specified flags
    user_flags = get_all_flags(options)

    # look up options and flags
    taskType=options.taskType
    if taskType is not None:
        relativePositions={}
        if options.inputFlag is None:
            options.inputFlag = programOptionMap[taskType]['in']
            if isinstance(options.inputFlag,int):
                relativePositions['in']=True
        if not options.outputFlags:
            options.outputFlags = programOptionMap[taskType]['out']
            for flag in options.outputFlags:
                if isinstance(flag,int):
                    relativePositions['out']=True
        if options.prefixFlag is None:
            options.prefixFlag = programOptionMap[taskType].get('prefix',None)
            if isinstance(options.prefixFlag,int):
                relativePositions['prefix']=True
        if options.threadsFlag is None:
            options.threadsFlag = programOptionMap[taskType].get('threads',None)
            if isinstance(options.threadsFlag,int):
                relativePositions['threads']=True
        if not options.rel_paths:
            options.rel_paths = programOptionMap[taskType].get('rel', [])
            for flag in options.rel_paths:
                if isinstance(flag, int):
                    relativePositions['rel'] = True

        if len(relativePositions):
            translatePositionalArgs(options,cmdargs,relativePositions)

    (positionalArgs,flaggedArgs)=getOptionHashes(options)
    logging.debug("positional args: %s from %s" % (str(positionalArgs),options.inputFlag))

    unsupportedFlags=unsupportedOptions.get(taskType,{}).keys()

    # scan command
    nextArgument = None
    usedFlags=[]
    for i in range(len(cmdargs)):

        if nextArgument is None:
            # is this a flag we can't handle?
            if cmdargs[i] in unsupportedFlags:
                raise Exception("Unsupported Flag: %s\n%s" % (cmdargs[i],unsupportedOptions[taskType][cmdargs[i]]))

            # check for positional arguments
            nextArgument=positionalArgs.get(i,None)
            if nextArgument is not None:
                usedFlags.append(i)

        # is this something we're looking for
        if nextArgument is None:
            # no. Well,is it a flag for something?
            arg=cmdargs[i]
            nextArgument=flaggedArgs.get(arg,None)
            if nextArgument is not None:
                usedFlags.append(arg)
            continue

        if nextArgument=='in':
            infile=cmdargs[i]
            logging.debug("Found infile (%s) in arg %d" % (str(infile),i))
        elif nextArgument=='prefix':
            prefix=cmdargs[i]
            logging.debug("Found prefix (%s) in arg %d" % (prefix,i))
        elif nextArgument=='out':
            outfile=cmdargs[i]
            logging.debug("Found output (%s) in arg %d" % (outfile,i))
        elif nextArgument=='threads':
            pass
        elif nextArgument=='rel':
            pass
        else:
            # something went wrong
            raise Exception("Unrecognized nextArgument: %s" % (nextArgument))

        # reset
        nextArgument=None

    # Quit if the user specified flags or any pos arguments are not found
    logging.debug("Used flags: \n%s\npos"
                  "args:\n%s\nuser_flags:\n%s\nflaggedArgs:\n%s",
                  usedFlags, positionalArgs, user_flags, flaggedArgs)
    unused_flags = set(positionalArgs).difference(usedFlags)\
                     .union(set(user_flags).difference(usedFlags))
    logging.debug("unused flags: \n %s", unused_flags)
    if len(unused_flags)>0:
        raise Exception("Unable to find the following arguments in the " + \
                         "command: '%s'" % "', '" \
                          .join("pos "+str(f) if isinstance(f,int) else f \
                                 for f in unused_flags))

    # remove any prefix based output flags if no pref
    if prefix is None and options.outputFlags is not None:
        i=0
        while i<len(options.outputFlags):
            if len(options.outputFlags[i])>1 and options.outputFlags[i][:2]=='%p':
                options.outputFlags.pop(i)
            else:
                i+=1

    return (infile,outfile)


def get_all_flags(options):
    """ return list of flags explicitly chose by user """
    flags = []
    if options.inputFlag:
        flags.append(try_to_int(options.inputFlag))
    if options.outputFlags:
        for flag in options.outputFlags:
            flags.append(try_to_int(flag))
    return flags

def try_to_int(flag):
    try:
        return int(flag)
    except ValueError:
        return flag


def translatePositionalArgs(options,cmdargs,checkMap):
    """
    using list of flags that take arguments for given taskType:
        translate argument positions from relative to absolute

    EG. in 'grep -v -H -f FILE1 FILE2', FILE2 is the first positional argument since FILE1 is the argument to the -f flag. So 1 (the positional argument number) would get translated to 5 (its index in the cmdargs array)
    """
    logging.info('translating positional args: %s' % checkMap)
    logging.debug('Command: %s' % cmdargs)

    # build map from relative positions to absolute positions
    if cmdargs[0] in scriptingLanguages:
        start=2
    else:
        start=1
    positionMap=[]
    skipNext=False
    for i in range(start,len(cmdargs)):
        if skipNext:
            skipNext=False
            continue
        arg=cmdargs[i]
        if len(arg)>0 and arg[0]=='-':
            if arg in flagsWithArguments[options.taskType]:
                skipNext=True
        else:
            positionMap.append(i)

    logging.debug("positions: %s" % positionMap)
    logging.debug(repr(options.inputFlag))

    # scan options and translate
    if 'in' in checkMap:
        if isinstance(options.inputFlag,int):
            options.inputFlag=positionMap[options.inputFlag-1]
    if 'out' in checkMap:
        for i in range(len(options.outputFlags)):
            if isinstance(options.outputFlags[i],int):
                options.outputFlags[i]=positionMap[options.outputFlags[i]-1]
    if 'rel' in checkMap:
        for i in range(len(options.rel_paths)):
            if isinstance(options.rel_paths[i], int):
                options.rel_paths[i] = positionMap[options.rel_paths[i]-1]
    if 'prefix' in checkMap:
        if isinstance(options.prefixFlag,int):
            options.prefixFlag=positionMap[options.prefixFlag-1]
    if 'threads' in checkMap:
        if isinstance(options.threadsFlag,int):
            options.threadsFlag=positionMap[options.threadsFlag-1]

def checkCommand(command,sgeOptions,batchOptions):
    """
    parse command string and do some basic checks:
        if command is recognized, check for required options

    returns tast type
    """
    if len(command)==0:
        raise Exception("No command given. Add ' -- ' followed by the command you'd like to run!")

    logging.debug("Command: %r" % (repr(command)))

    # get program name
    program=command[0]
    if program in scriptingLanguages:
        # if first word in command is perl or python (etc), take second word
        if len(command)==1:
            raise Exception("%s won't do much without a script name" % (program))
        program = command[1]

    # do we know this program?
    for (taskType,taskProgRE) in taskTypePatterns.items():
        m=taskProgRE.search(program)
        if m:
            logging.info("Program recognized as type: %s" % (taskType))
            # run any task specific processing on command
            #  use aplyDefaults() if nothing specific for this task
            if taskType in inspectCommandForTaskType:
                inspectCommandForTaskType[taskType](command,taskType,sgeOptions,m.group(1),batchOptions)
            else:
                applyDefaultsToCommand(command,taskType)
            return taskType

    # if not recognized, just leave as is
    return None

def launchLocalJobs(options, cmdargs, errStream):
    logging.debug("Launching command array: %r" % ({'tmpDir':options.tmpDir,'splits':options.splits,'fragName':options.fragBase,'cmd':cmdargs,'loglevel':options.verbose, 'type':options.taskType, 'throttle':options.throttle}))
    
    numTasks = options.splits
    numThreads = min(numTasks,options.throttle) if options.throttle>0 else numTasks

    # Create thread safe object to report errors to
    errorQueue = queue.Queue()

    # create tasks
    taskQueue = queue.Queue()
    for i in range(1,numTasks+1):
        taskQueue.put(FragmentTask(options, cmdargs, i, errorQueue))

    # create and start threads
    for i in range(numThreads):
        thread = FragmentThread(taskQueue)
        logging.debug("Starting thread %d" % (i+1))
        thread.start()

    logging.debug("Started %d threads" % (numThreads))
    # Wait for all tasks to be done
    taskQueue.join()

    # Check for any errors
    try:
        exitcode=errorQueue.get(False)
    except queue.Empty:
        # no errors
        return

    # There was at least one error
    # add this one to the remaining to get count
    errorCount = errorQueue.qsize()+1
    logging.warning("%d of %d tasks reported non zero exit codes!" % (errorCount, numTasks))
    sys.exit(exitcode)

def launchJobs(options, cmdargs, errStream=sys.stdin):
    """
    launch the SGE task array. Return the job ID

    splits can be a single interger or a pair of integers
    this will launch
         the indicated task number
         the range of tasks numbers indicated by the pair (inclusive)
    """

    if options.queue == LOCAL:
        launchLocalJobs(options,cmdargs,errStream)
        return

    logging.debug("Launching task array: %r" % ({'tmpDir':options.tmpDir,'splits':options.splits,'fragName':options.fragBase,'cmd':cmdargs,'sgeOpts':options.sgeOptions,'job':options.jobName,'priority':options.priority,'loglevel':options.verbose,'wait':options.wait, 'type':options.taskType}))
    
    # SGE or SLURM submission prefix
    command = getSubmissionCommandPrefix(options)

    # batch_runner command
    command.append(BATCHLAUNCHER)
    command+=["--mode","run","--tmp_dir",options.tmpDir,"--frag_base",
              options.fragBase, "--frag_dir", options.frag_dir, "--frag_suffix", options.fragSuff, "--loglevel", str(options.verbose), "--queue", options.queue]
    if options.inputFlag is not None:
        command+=['-i',str(options.inputFlag)]
    if options.prefixFlag is not None:
        command+=['-p',str(options.prefixFlag)]
    if options.threadsFlag is not None:
        command+=['-t',str(options.threadsFlag)]
    if options.outputFlags is not None:
        for flag in options.outputFlags:
            command+=['-o',str(flag)]
    if options.taskType is not None:
        command+=['--taskType',options.taskType]
    if options.cwd:
        command.append('--cwd')
    command.append('--')
    command+=cmdargs

    # redirect qsub output to std, silence if vebose is 0
    #if options.verbose==0:
    #    qsubOuts=open(os.devnull,'w')
    #else:
    #    qsubOuts=errStream
    
    # run command
    logging.debug('Launching task array: %s' % (formatCommand(command)))
    try:
        submissionOutput=subprocess.check_output(command)
        if options.verbose>0:
            errStream.write("Submission Output: " + submissionOutput.decode())
    except subprocess.CalledProcessError as error:
        if options.wait and options.queue != SLURM:
            # when using -sync y, the exit code may come from a task
            #  (which cleanup will handle)
            logging.warning("qsub returned an error code of: %d" 
                    % error.returncode)
        else:
            raise error

    # get job id
    try:
        jobid = re.search(r'(\d+)\s*$',submissionOutput).group(1)
        options.jobid = jobid
    except:
        if options.queue==SLURM:
            logging.error("Cannot parse SLURM job id from '%s'" % (submissionOutput))
            raise

    # SLURM doesn't allow waiting for completion on array jobs, so we hack:
    #  use srun to start a dummy job that will wait for our job array
    if options.wait and options.queue==SLURM:
        waitForSlurmArray(options, errStream)

def waitForSlurmArray(options, errStream):

    # redirect qsub output to std, silence if vebose is 0
    if options.verbose==0:
        qsubOuts=open(os.devnull,'w')
    else:
        qsubOuts=errStream

    command=["srun",]
    command+=["--dependency=afterany:%s" % (options.jobid),]
    #if isinstance(options.sgeOptions,str):
    #    command+=shlex.split(options.sgeOptions)
    #else:
    #    command+=options.sgeOptions

    command += ["true"]

    logging.debug('Launching command to wait: %s' % (formatCommand(command)))
    exitcode=subprocess.call(command, stdout=qsubOuts)

    if exitcode!=0:
        raise subprocess.CalledProcessError(exitcode, command)


def getSubmissionCommandPrefix(options, cleanupFile=None):
    """
    Generate the SGE or SLURM submission prefix (as a list of strings). 
    """
    submitData={}
    #priority
    priority=options.priority
    if cleanupFile is not None:
        priority = min(0,priority)
    else:
        priority = min(-1,priority-10)
    submitData['priority']=str(priority)
    submitData['interpreter']=PYTHON

    # job name and output files
    if cleanupFile is not None:
        submitData['jobName']="%s.cleanup" % (options.jobName)
        submitData['output']="%s.cleanup.out" % (cleanupFile)
        submitData['error']="%s.cleanup.err" % (cleanupFile)
        submitData['waitfor']=options.jobName if options.queue != SLURM else options.jobid
    else:
        submitData['jobName']=options.jobName
        submitData['output']="%s%stasks.out" %(options.tmpDir, os.sep)
        submitData['error']="%s%stasks.err" % (options.tmpDir, os.sep)
        if isinstance(options.splits,int):
            submitData['tasks'] = "1-%d" % options.splits
        else:
            if options.splits[0]==options.splits[1]:
                submitData['tasks']=str(options.splits[0])
            else:
                submitData['tasks']="%d-%d" % tuple(options.splits)
        if options.throttle>0:
            submitData['throttle']=str(options.throttle)
        if options.wait:
            submitData['sync']=True
    if isinstance(options.sgeOptions,str):
        submitData['options']=shlex.split(options.sgeOptions)
    else:
        submitData['options']=options.sgeOptions

    if options.queue==SGE:
        return getSGECommandPrefix(submitData)
    else:
        return getSLURMCommandPrefix(submitData)

def getSLURMCommandPrefix(submitData):
    command=["sbatch",]
    #, "-hard", "-cwd", "-V"]
    if 'waitfor' in submitData:
        command+=["--dependency=afterany:%s" % (submitData['waitfor']),]
    #command+=["-p",submitData['priority']]
    command+=["--job-name=%s" % (submitData['jobName']),
              "--output=%s" % (submitData['output']),
              "--error=%s" % (submitData['error'])]
    if 'tasks' in submitData:
        taskSpec = submitData['tasks']
        if 'throttle' in submitData:
            taskSpec+="%" + submitData['throttle']
        command+=["--array=%s" % (taskSpec),]
    #if 'sync' in submitData:
    #    command+=['-sync','y']
    if 'options' in submitData:
        command+=submitData['options']
    return command

def getSGECommandPrefix(submitData):
    command=["qsub", "-hard", "-cwd", "-V"]
    command+=["-S",submitData['interpreter']]
    command+=["-p",submitData['priority']]
    if 'waitfor' in submitData:
        command+=["-hold_jid", submitData['waitfor']]
    command+=["-N", submitData['jobName'],
              "-o", submitData['output'],
              "-e", submitData['error']]
    if 'tasks' in submitData:
        command+=["-t",submitData['tasks']]
    if 'throttle' in submitData:
        command+=["-tc", submitData['throttle']]
    if 'sync' in submitData:
        command+=['-sync','y']
    if 'options' in submitData:
        command+=submitData['options']
    return command

def launchCleanup(options, cmdargs, errStream=sys.stderr):
    """
    launch the cleanup script, which is this script with the --cleanup flag
    """
    logging.debug("Launching cleanup: %s" % ({'tmpDir':options.tmpDir,'splits':options.splits,'fragBase':options.fragBase,'out':options.outputFlags,'job':options.jobName}))

    # name outputfiles
    outFileByFlag = getOutputFiles(options, cmdargs)
    logging.debug("Outfiles: %s" % (outFileByFlag))
    fileNameBase = getFileNameBase(options.outputFlags,outFileByFlag,options.jobName)
    logging.debug("File Name Base: %s" % (fileNameBase))

    # build command
    # sge
    command = getSubmissionCommandPrefix(options, cleanupFile=fileNameBase)
    
    # cleanup
    command+=[BATCHLAUNCHER,
              '--tmp_dir', options.tmpDir,
              '--frag_base', options.fragBase,
              '--frag_dir', options.frag_dir,
              '--frag_suffix', options.fragSuff,
              '--mode', 'cleanup',
              '--numChunks', str(options.splits), 
              '--loglevel', str(options.verbose),
              '--retries', str(options.retries),
              '--jobName', options.jobName,
              '--chunk', str(options.chunk),
              '--queue', options.queue]
    if options.splitOnSize:
        command.append('--splitOnSize')
    if options.inputFlag is not None:
        command+=['-i',str(options.inputFlag)]
    if options.prefixFlag is not None:
        command+=['-p',str(options.prefixFlag)]
    if options.threadsFlag is not None:
        command+=['-t',str(options.threadsFlag)]
    if options.outputFlags is not None:
        for flag in options.outputFlags:
            command+=['-o',str(flag)]
    if options.taskType is not None:
        command+=['--taskType',options.taskType]
    if options.pattern is not None:
        command+=['--pattern', options.pattern]
    if options.numLines is not None:
        command+=['--recordLines', options.numLines]
    elif options.infileType is not None:
        command+=['--infileType',options.infileType]
    if options.cwd:
        command.append('--cwd')
    if options.sgeOptions:
        command.extend(['--submitOptions'," ".join(options.sgeOptions)])
    command.append('--')
    command+=cmdargs

    # redirect qsub output to std, silence if vebose is 0
    if options.verbose==0:
        qsubOuts=open(os.devnull,'w')
    else:
        qsubOuts=errStream

    # run command
    logging.debug('Launching cleanup: %s' % (formatCommand(command)))
    subprocess.check_call(command, stdout=qsubOuts)

def cleanup(options, cmdargs, errStream=sys.stdin):
    """
    Concatenate all files in tmpDir of the form fragBaseName.###.out into outFile where
### is a number between 1 and splits, inclusive.

    If there is an output file missing, save the output.
    Finally, delete tmpDir and its contents
    """

    logging.debug("Cleanup: retries=%d" % options.retries)
    exitcode=0

    # get list of output flags
    outFileByFlag = getOutputFiles(options, cmdargs)
    logging.debug("Outfiles: %s" % (outFileByFlag))

    # name outputfiles
    fileNameBase = getFileNameBase(options.outputFlags, 
                                   outFileByFlag, 
                                   options.jobName)

    # remove old output files
    errStreamFile="%s.stderr" % fileNameBase
    failureStreamFile="%s.failures" % fileNameBase
    for file in errStreamFile, failureStreamFile:
        if os.path.exists(file):
            logging.debug('Removing previous file: %s' % file)
            os.remove(file)
    for file in outFileByFlag.values():
        if file is not None:
            if os.path.exists(file):
                logging.debug('Removing previous file: %s' % file)
                os.remove(file)

    # set up copy method (some task types might do some filtering)
    copyFilesToStream=taskSpecificCopy.get(options.taskType,addFilesToStream)

    # loop until everything is node or we give up
    taskIds=list(range(1,options.splits+1))
    errStream = None
    failureStream = None
    while True:
        logging.debug("starting to scan for fraqgments: %r (retries: %d)" % (taskIds,options.retries))
        # if no output file specified, add STDOUT
        if len(outFileByFlag)==0:
            outFileByFlag['%stdout']=None
        # Change filenames to (filename,None) tuples that can be populated with streams
        for flag in outFileByFlag.keys():
            if flag == '%stdout':
                # special case for STDOUT
                outFileByFlag[flag]=sys.stdout
            elif isinstance(outFileByFlag[flag],list):
                # reset tuples leftover from previous loop
                outFileByFlag[flag][1]=None
            else:
                # Change filenames to (filename,None) tuples that can be populated with streams
                outFileByFlag[flag]=[outFileByFlag[flag],None]

        # keep track of things to resubmit
        failedTasks=[]
        anySuccess=False
        missingRecords={}
        # look for files
        for i in taskIds:
            # look for output
            fragName = getFragmentName(options.fragBase, i, options.fragSuff)
            prefix = getFragmentPrefix(options.fragBase,i)
            frag = "%s%s%s" % (options.tmpDir, os.sep, fragName)
            fragerr = "%s.exitcode" % (frag)
            outfrag = "%s.stdout" % (frag)
            errfrag = "%s.stderr" % (frag)
            logfrag = "%s.log" % (frag)
            outfragmap={}

            # For each configured output file, map fragment to final
            for (flag, flagOutFile) in outFileByFlag.items():
                if flag=='%stdout':
                    outfragmap[outfrag]=flagOutFile
                else:
                    (tmpDir,otheroutfrag) = getOutputFromFlag(flag,fragName,prefix,options.tmpDir,options.tmpDir)
                    outfragmap["%s%s%s" % (tmpDir,os.sep,otheroutfrag)]=flagOutFile

            if not(os.path.exists(fragerr)):

                # copy results
                try :
                    anyFile=copyFilesToStream(outfragmap,i,frag)

                    # save log,stdout, and stderr if loglevel is high
                    if os.path.exists(logfrag):
                        anyFile=True
                        if options.verbose>=2:
                            if errStream is None:
                                create_parent_dir(errStreamFile)
                                errStream = open(errStreamFile, 'w')
                            addFileToStream(logfrag,errStream,header="## LOGGING from fragment %d:" % (i))
                            if outfrag not in outfragmap:
                                addFileToStream(outfrag,errStream,header="## STDOUT from fragment %d:" % (i))
                            addFileToStream(errfrag,errStream,header="## STDERR from fragment %d:" % (i))

                    # delete files (input, error, outputs)
                    for f in [frag, outfrag, errfrag, logfrag] \
                            + list(outfragmap.keys()):
                        if os.path.exists(f):
                            anyFile=True
                            os.remove(f)

                    if anyFile:
                        anySuccess=True
                        continue

                except FailedFragmentException as ffe:
                    if len(ffe.records) < options.chunk:
                        anySuccess=True
                    logging.info("Task %d has missing records" % i)
                    missingRecords[i]=ffe

            else:
                # there was an error
                logging.info("Task %d failed" % i)
                failedTasks.append(i)

            ## If we got here, there was an error!

            # make sure error streams are open
            if errStream is None:
                errStream = open(errStreamFile, 'w')
            if failureStream is None:
                failureStream = open(failureStreamFile, 'w')

            # append to error streams
            if os.path.exists(logfrag):
                addFileToStream(logfrag,errStream,header="## LOGGING from fragment %d:" % (i))
            else:
                errStream.write("## LOGGING not found for fragment %d!\n" % (i))
            if os.path.exists(errfrag):
                addFileToStream(errfrag,errStream,header="## STDERR from fragment %d:" % (i))
            else:
                errStream.write("## STDERR not found for fragment %d!\n" % (i))
            if outfrag not in outfragmap:
                if os.path.exists(outfrag):
                    addFileToStream(outfrag,errStream,header="## STDOUT from fragment %d:" % (i))
                else:
                    errStream.write("## STDOUT not found for fragment %d!\n" % (i))

            # save failed records to file
            for failfrag in outfragmap:
                if os.path.exists(failfrag):
                    if os.path.isdir(failfrag):
                        failureStream.write("## FAILURES: fragment %d failed." % (i))
                        # TODO: do something with the failed output
                    else:
                        addFileToStream(failfrag,failureStream,header="## FAILURES: %s from fragment %d:" % (failfrag,i))
                        os.remove(failfrag)
                else:
                    failureStream.write("## FAILURES: %s not found for fragment %d!\n" % (failfrag,i))

            # delete files (exitcode, error, outputs) (save input for re-queueing)
            for f in [fragerr, outfrag, errfrag,]:
                if os.path.exists(f):
                    os.remove(f)

        # Finished scanning fragments
        logging.info("Cleanup is done scanning output files: rtrs: %d, aS: %s, fT: %d, mR: %d" % (options.retries, anySuccess, len(failedTasks), len(missingRecords)))

        # close output streams
        for outstream in outfragmap.values():
            if outstream is sys.stdout:
                continue
            if isinstance(outstream,list):
                if outstream[1] is not None:
                    outstream[1].close()

        # If conditions are right, resubmit any failures:
        if anySuccess and options.retries!=0:
            options.retries-=1
            logging.info("Cleanup is checking for anything that needs to be restarted")

            # get the next available task number  (i will still be set from loop)
            nextTaskNum=0

            # build new input fragment from afiled and missed fragments in
            # subdirectory of tmpDir

            # first check tasks that failed completely
            # rename tasks to make them consecutive
            if len(failedTasks)>0:
                nextTaskNum+=reFragmentMissedTasks(failedTasks, options)

            # then, if we were able to identify missing records
            if len(missingRecords)>0:
                # build new fragments out of the missed records
                nextTaskNum=buildMissedRecordFragments(missingRecords, options.tmpDir, options.fragBase, nextTaskNum, options.chunk)

            # rerun any missed records and failed fragments
            if nextTaskNum>0:
                # finish setting up tmp dir
                options.splits=nextTaskNum-1
                moveNewFragmentsToTmpDir(options,nextTaskNum)

                # re-process tmp dir
                options.wait=True
                logging.info("Cleanup will restart tasks: %s" % (options.splits))
                launchJobs(options, cmdargs, errStream=errStream)

                # set up list of fragments to check on next cleanup pass
                taskIds = list(range(1,nextTaskNum))

            else:
                # everything is complete
                logging.debug("All tasks were successful")

                # TODO:
                # Remove failures file if it exists

                break
        else:
            # either everything failed or we give up: exit loop
            logging.debug("Cleanup will not re-start any tasks.")
            exitcode=1
            break


    logging.info("Final cleanup")
    # check contesnts of tasks.err and tasks.out in options.tmpDir
    logging.debug("collecting stderr and stdout from fragments")
    commonerr="%s%stasks.err"%(options.tmpDir,os.sep)
    commonout="%s%stasks.out"%(options.tmpDir,os.sep)
    # if not empty, add to errStream (make sure it's open)
    if os.path.exists(commonerr):
        if os.path.getsize(commonerr)>0:
            if errStream is None:
                errStream = open(errStreamFile, 'w')
            addFileToStream(commonerr,errStream,header="## Uncaptured error output from all tasks:")
        os.remove(commonerr)
    if os.path.exists(commonout):
        if os.path.getsize(commonout)>0:
            if errStream is None:
                errStream = open(errStreamFile, 'w')
            addFileToStream(commonout,errStream,header="## Uncaptured standard output from all tasks:")
        os.remove(commonout)

    # warn if any files left
    logging.debug("Checking for leftover files")
    leftoverFiles=os.listdir(options.tmpDir)
    if len(leftoverFiles)>0:
        if errStream is None:
            errStream = open(errStreamFile, 'w')
        errStream.write("Files left in %s: %r" % (options.tmpDir, leftoverFiles))
        for f in leftoverFiles:
            leftoverFilePath=os.sep.join([options.tmpDir,f])
            if os.path.isdir(leftoverFilePath):
                errStream.write("Cannot delete directory: %s" % (f))
            else:
                os.remove(leftoverFilePath)

    if errStream is not None:
        errStream.close()
    if failureStream is not None:
        failureStream.close()
    # delete directory
    logging.debug("removing tmp folder %s", options.tmpDir)
    os.rmdir(options.tmpDir)
    logging.debug("cleanup is complete")
    return exitcode

def reFragmentMissedTasks(missedTasks, options):
    """
    combine all missed reads into one and re-fragment into the same number of pieces
    """
    options.chunk=1+(options.chunk*len(missedTasks)/options.splits)
    temporaryLocation="%s%stmp"%(options.tmpDir,os.sep)
    os.makedirs(temporaryLocation)

    fileType = getFileType(options, None)

    # create a fileHandle-like object that will read all missed fragments
    inputsToReFragment=[getFragmentPath(options.tmpDir, options.fragBase, i) for i in missedTasks]
    logging.info("Restarting fragments: %s" % missedTasks)
    logging.debug("Restarting fragments: %s" % inputsToReFragment)
    failedRecordStream = fileinput.input(inputsToReFragment)

    # create new fragments in temporaryLocation
    newFragNum=fragmentInputStreamBySize(failedRecordStream, temporaryLocation,
                                         options.chunk, fileType,
                                         options.fragBase,
                                         splitOnSize=options.splitOnSize,
                                         suffix=options.fragSuff)

    # remove old fragments
    for i in missedTasks:
        frag = getFragmentPath(options.tmpDir, options.fragBase, i)
        os.remove(frag)

    return newFragNum+1

def moveNewFragmentsToTmpDir(options,nextTaskNum):
    """
    Move new fragments from tmpDir/tmp/ into tmpDir/
    """
    for i in range(1,nextTaskNum):
        frag = getFragmentPath(options.tmpDir, options.fragBase, i)
        newfrag = getFragmentPath("%s%stmp" % (options.tmpDir, os.sep), options.fragBase, i)
        os.rename(newfrag,frag)
    os.rmdir("%s%stmp" % (options.tmpDir, os.sep))

def getFileNameBase(outFlags,outFileByFlag,default):
    " what do we name all the log files? "
    if outFlags:
        logging.debug("Output Flags: \n%r", outFlags)
        logging.debug("Output Files \n%r", outFileByFlag)
        # base output file on first output file selected that looks good
        # This is a hack to fix fallout from the secondary file code
        best_file = None
        for flag in outFlags:
            out_file = outFileByFlag.get(flag, None)
            if out_file is None:
                continue
            if (not out_file.startswith('-')) and (not
                                                   out_file.startswith('./-')):
                # this is a good file, use it and stop looking
                best_file = out_file
                logging.debug("bestest file: %s", best_file)
                break
            if best_file is None:
                # if we haven't seen anything else yet, save this
                best_file = out_file
                logging.debug("best file: %s", best_file)
        if best_file is not None:
            return best_file
        else:
            # has been translated to integer, take #1 (first in command, not in flags)
            return outFileByFlag[1]

    return default

def getMissedRecordIterator(taskType,records,infile):
    """
    return an iterator that will iterate over records in a given file.
    each record is returned as a list of lines
    """
    # get generator method by task type
    generator=missedRecordGenerator.get(taskType,reBasedRecordGenerator)
    # call generator to get iterator
    generator(taskType,records,infile)

def reBasedRecordGenerator(taskType,records,infile):
    """
    simple generator that assumes first line of each record has the record ID
    There must be a regularExpression stored in recordIDRE for the given task type
    """
    recordIDRE=recordIDRE[taskType]
    currentRecord=None
    currentData=[]
    with open(infile,'rt') as f:
        for line in f:
            m=recordIDRE.match(line)
            if m:
                if currentRecord is not None:
                    yield currentData
                record=m.group(1)
                if record in records:
                    currentRecord=record
                    del currentData[:]
                    currentData.append(line)
                else:
                    currentRecord=None
            elif currentRecord is not None:
                currentData.append(line)
    if currentRecord is not None:
        yield currentData

def buildMissedRecordFragments(missingRecords, tmpDir, fragBaseName, nextFragNum, chunk):
    """
    Using:
        the missing record map which maps frgment numbers to missed record IDs
        the pattern for finding and IDing records
        tmpDir and fragBaseName to find input fragments by number
    TODO:
        scan fragment files in tmpdir to find the next fragment number
        create new fragment file(s) by pulling missed records from other fragment files
        keep new fragment files at least as small as fragment 1
    RETURN:
        tuple with range of new fragment numbers
    """

    # get missed records as list of strings
    recordCount=chunk
    outStream=None
    for ffe in missingRecords.values():
        # iterator will return a list of lines at each step
        missedRecordIterator=getMissedRecordIterator(ffe.taskType,ffe.records,ffe.inputFragment)
        for record in missedRecordIterator:
            if recordCount%chunk==0:
                if outStream is not None:
                    outStream.close()
                recordCount=0
                outStream=open(getInFragmentName(tmpDir,fragBaseName,nextFragNum),'w')
                nextFragNum+=1
            for line in record:
                outStream.write(line)
            recordCount+=1

    outStream.close()
    return nextFragNum

class FailedFragmentException(Exception):
    def __init__(taskType,inputFragment,recordList):
        self.records=recordList
        self.taskType=taskType
        self.inputFragment=inputFragment

def addBlastFilesToStream(fileMap,fragIndex,inputFragment):
    """
    copy fragment output to given srteam, but first, make sure the results are all there
    """
    # TODO: return False if output files are missing completely

    # get list of reads from input fragment
    records={}
    fastaRecordRE=fileTypeMap['fasta'].sepRE
    with open(inputFragment,'rt') as f:
        for line in f:
            m=fastaRecordRE.match(line)
            if m:
                records[m.group(1)]=False

    # assume only one output file for blast
    outFragment=next(iter(fileMap.keys()))
    outStream=getStream(fileMap[outFragment])

    # scan file and write good records to output
    currentRecord=None
    currentRecordOutput=[]
    state=START
    with open(outFragment,'rt') as f:
        line=f.next()
        if blastHeaderRE.match(line):
            state=HEADER
            currentRecordOutput.append(line)
        else:
            # rps blast is weird
            blastHeaderRE=blastRecordRE
        for line in f:
            if state==START:
                if blastHeaderRE.search(line):
                    state=HEADER
                    writeArrayToStream(currentRecordOutput,outStream)
                    del currentRecordOutput[:]
                    if blastHeaderRE is not blastRecordRE:
                        # in most cases, just start saving next record
                        currentRecordOutput.append(line)
                        continue
                    # in the rps case, stay on this line and do header processing
            if state==HEADER:
                m=blastRecordRE.match(line)
                if m:
                    if currentRecord is not None:
                        del records[currentRecord]
                    currentRecord=m.group(1)
                    state=HITS
                else:
                    m=blastHeaderRE.search(line)
                    if m:
                        currentRecord=None
                        del currentRecordOutput[:]
                        line=line[m.start():]
                        state=HEADER
            elif state==HITS:
                if foundHitsRE.match(line):
                    state=START
                elif noHitsRE.match(line):
                    state=START
                else:
                    m=blastHeaderRE.search(line)
                    if m:
                        records[currentRecord]=False
                        currentRecord=None
                        del currentRecordOutput[:]
                        line=line[m.start():]
                        state=HEADER
            currentRecordOutput.append(line)

    if state!=START:
        # ended part way through record
        records[currentRecord]=False
    else:
        del records[currentRecord]

    if len(records)!=0:
        logging.debug("Missing %d(%d bad) records from fragment %s" % (len(records),len(records)-sum(records.values()),fragIndex))
        raise FailedFragmentException(BLAST,inputFragment,records.keys())

def addFilesToStream(fileMap,fragIndex,inputFragment):
    """
    given map of file fragment to file object it needs to be appended to:
        copy data from each fragment to appropriate stream
    """
    logging.debug("addFilesToStream(%r,...)" % (fileMap))
    anyFileExists=False
    for fragfile in fileMap:
        if os.path.exists(fragfile):
            anyFileExists=True
            # copy output fragment contents to appropriate stream or directory
            addFileToStream(fragfile, fileMap[fragfile], 
                            outputIsDir=os.path.isdir(fragfile))

    return anyFileExists

def getStream(outstream):
    """
    Returns an open file handle for writing.

    The 'outstream' parameter can be a an open file handle for writing (which will simply be returned as is) or a [fileName,fielHandle] pair. The handle in the pair may also be set to None, if the file is not yet open for writing. 
    """

    if isinstance(outstream,list):
        # return second element (should be an open handle)
        if outstream[1] is None:
            # open handle if needed
            create_parent_dir(outstream[0])
            outstream[1]=open(outstream[0],'a')
        return outstream[1]
    else:
        return outstream

def addFileToStream(filename, outstream, header=None, outputIsDir=False):
    """
    Copy contents of file into an output stream.
    outstream should be either a fileobject or a two element list: [filename, srtream]
     the strem can be None, it will be opeed here
    Optionally print header before.
    """
    if outputIsDir:
        # Special case if output is a directory
        copyFragmentOutputDir(filename, outstream)
        return

    # if outstream is file,stream pair, get stream
    outstream=getStream(outstream)

    if header is not None:
        outstream.write(header)
        outstream.write('\n')

    with open(filename, 'rt') as f:
        for line in f:
            outstream.write(line)

def copyFragmentOutputDir(fragmentOutputDir, finalOutputData):
    """
    Copy the data in fragment output dir to the location specified by finalOutputData
    """
    # First make sure the final output location exists as a directory
    if not isinstance(finalOutputData,str):
        # if it was a list/tuple, grab the fileName from the first position
        finalOutputData=finalOutputData[0]

    if os.path.exists(finalOutputData):
        # It already exists, make sure it's a directory
        if not os.path.isdir(finalOutputData):
            raise Exception("The output of this fragment is a directory (%s), but the final output already exists as a file (%s)." % (fragmentOutputDir, finalOutputData))
    else:
        # Create output directory 
        os.makedirs(finalOutputData)

    # put fragment directory in final output directory as is, let the user deal with it (TODO: smart merging of files)
    shutil.move(fragmentOutputDir,
                finalOutputData)

def runFragment(options, cmdargs, taskNum=None):
    """
    Execute the command on a fragment of the input file. In the SGE model, this should
    only be run on a compute node.
    The fragment number will be pulled from the SGE_TASK_ID env variable

    When fragmenting locally, the task number will be set explicitly in the
    call.
    """
    t0=time.time()
    exitcode=-3

    # create local tmp dir
    if options.queue == LOCAL:
        localDir = tempfile.mkdtemp(suffix='batchRunner', dir=options.tmpDir)
    else:
        localDir = tempfile.mkdtemp(suffix='batchRunner')
    logging.debug("Local dir (%s): %s" % (os.path.exists(localDir), localDir))

    #if options.queue != LOCAL:
    if taskNum is None:
        # get the task number from env
        taskNum=getTaskId(options)

    hostname=getHostName()
    logFile="%s%stask%05d.%s.log" % (localDir,os.sep,taskNum,hostname)
    errStream=open(logFile,'w')
    # if called from SGE, we will have to:
    if options.queue != LOCAL:
        # set up logging
        if options.verbose==0:
            loglevel=logging.ERROR
        elif options.verbose==1:
            loglevel=logging.WARN
        elif options.verbose==2:
            loglevel=logging.INFO
        elif options.verbose>=3:
            loglevel=logging.DEBUG
        logging.basicConfig(stream=errStream, level=loglevel)

    # set up file names
    infragmentName = getFragmentName(options.fragBase, taskNum,
                                     options.fragSuff)
    fragment_dir = options.frag_dir
    prefix = getFragmentPrefix(options.fragBase, taskNum)
    infragment = "%s%s%s" % (fragment_dir, os.sep, infragmentName)
    stdoutFragment=os.path.join(options.tmpDir, infragmentName) + ".stdout"
    stderrFragment=os.path.join(options.tmpDir, infragmentName) + ".stderr"
    if options.queue == LOCAL:
        stdoutLocal=stdoutFragment
        stderrLocal=stderrFragment
    else:
        stdoutLocal = "%s%stask%s.%s.stdout" % ( localDir, os.sep, taskNum, hostname )
        stderrLocal = "%s%stask%s.%s.stderr" % ( localDir, os.sep, taskNum, hostname )

    try:
        logging.debug("Begin runFragment: type: %s, host: %s" % (options.taskType, hostname))
        # check arguments
        if options.tmpDir is None:
            logging.error("batch_runner needs --tmp_dir parameter!")
            parser.error("MUST supply --temp_dir")

        # if cwd is False, use local tmp directory, otherwise use current directory
        if options.cwd:
            subprocwd=None
        else:
            localDir = os.path.abspath(localDir)
            subprocwd=localDir
            infragment=os.path.abspath(infragment)
            # translate any relative paths in command to absolute
            makeRelativePathsAbsolute(cmdargs)

        # any command specific changes
        if options.taskType in finalCommandProcessingForTask:
            finalCommandProcessingForTask[options.taskType](cmdargs,localDir)

        # modify command
        outLocal = "%s%s%s.output" % (localDir, os.sep, infragmentName)
        (foundI, outputFlags) = prepareCommandForFragment(options, infragment, prefix, outLocal, cmdargs, hostname, errStream)
        logging.debug("ready ro RUN!")

        # fragmented HMMdbs need to be compiled
        prep_input_fragment(infragment,
                            fragment_prep_for_task.get(options.taskType,
                                                       options.frag_prep))

        # setup to run command
        # I/O
        if foundI:
            spin=None
        else:
            spin=open(infragment,'rt')
        spout=open(stdoutLocal,'w')
        sperr=open(stderrLocal,'w')
        # environment
        path=os.environ['PATH']
        #path=os.pathsep.join([BINDIR,os.environ['PATH']])

        # run it
        try:
            logging.debug("Command:\n%s\nwd=%s\npath=%s" % (formatCommand(cmdargs),subprocwd,path))
            t1=time.time()
            # Run the command
            exitcode=subprocess.call(cmdargs, stdin=spin, stdout=spout,
                                     stderr=sperr, cwd=subprocwd)
            #exitcode=subprocess.call(cmdargs, stdin=spin, stdout=spout, stderr=sperr, cwd=subprocwd, env={"PATH":path})
            t2=time.time()
            logging.info("Command took %.2f seconds" % (t2-t1))
        except:
            errStream.write("Exception executing command: %s\n" % (cmdargs))
            errStream.write('-'*60+'\n')
            traceback.print_exc(file=errStream)
            errStream.write('-'*60+'\n')
            exitcode=-1

        if not foundI:
            spin.close()
        spout.close()
        sperr.close()

        # Report any error code
        if exitcode != 0:
            logging.error("Error code %d from command:\n%s" % (exitcode,cmdargs))

        if options.queue != LOCAL:
            # copy stderr and stdout
            logging.info("Copying %s to %s" % (stdoutLocal,stdoutFragment))
            shutil.copy(stdoutLocal,stdoutFragment)
            logging.info("Copying %s to %s" % (stderrLocal,stderrFragment))
            shutil.copy(stderrLocal,stderrFragment)

        # move output files
        for flag in outputFlags:
            # check for other files
            logging.debug("Looking for file from output flag: %s" % (flag))
            # skip stdout
            if flag=='%stdout':
                # we've already grabbed it
                continue

            (outputdir,output) = getOutputFromFlag(flag,infragmentName,prefix,localDir,options.tmpDir)
            logging.debug("file should be: %s in %s" % (output,outputdir))
            if outputdir is localDir:
                localOutput = "%s%s%s" % (outputdir,os.sep,output)
                tmpOutput = "%s%s%s" % (options.tmpDir, os.sep, output)
                if os.path.exists(localOutput):
                    logging.info("Copying %s to %s" % (localOutput, tmpOutput))
                    create_parent_dir(tmpOutput)
                    shutil.move(localOutput,tmpOutput)
                else:
                    logging.warning("Could not find output: %s" % output)
                    logging.debug("Flag: %s" % (flag))

    except:
        exitcode=2
        errStream.write("Exception running fragment %d:\n" % (taskNum))
        errStream.write('-'*60+'\n')
        traceback.print_exc(file=errStream)
        errStream.write('-'*60+'\n')

    # Do some final cleanup:
    if exitcode!=0:
        errCodeFile=os.path.join(options.tmpDir, 
                                 "%s.exitcode" % (infragmentName))
        ecStream=open(errCodeFile,'w')
        ecStream.write(str(exitcode))
        ecStream.close()

    t2=time.time()
    logging.info("Fragment took %.2f seconds" % (t2-t0))

    # copy err to shared dir
    if options.queue != LOCAL:
        logging.shutdown()
    errStream.close()
    if options.queue == LOCAL:
        shutil.move(logFile, os.path.join(options.tmpDir,
                                          "%s.log" % (infragmentName)))
    else:
        shutil.copy(logFile, os.path.join(options.tmpDir,
                                          "%s.log" % (infragmentName)))

    # remove local temorary dir (if not debugging)
    if logging.getLogger().level > logging.DEBUG:
        logging.debug("I AM removing %s", localDir)
        shutil.rmtree(localDir)
    else:
        logging.debug("NOT removing %s", localDir)

    return exitcode


def prep_input_fragment(infragment, frag_prep):
    """ fragmented HMMdbs need to be compiled"""
    if frag_prep:
        m = re.search(r'^([\.a-zA-Z0-9_]+):\s+(\S.+)$', frag_prep)
        if m:
            ext, command_template = m.groups()
            # don't do antthing if output exists and is newer
            if os.path.exists(infragment + ext) and  \
                (os.path.getmtime(infragment) < \
                 os.path.getmtime(infragment + ext)):
                return
        else:
            command_template = frag_prep
        logging.info("Running %s on %s", command_template, infragment)
        command = shlex.split(command_template.format(infragment))
        log_file = infragment + ".prep.log"
        logging.debug("Logging command: %r to %s", command, log_file)
        with open(log_file, 'w') as err_stream:
            subprocess.check_call(command,
                                  stdout=err_stream,
                                  stderr=err_stream,
                                 )


def getOptionHashes(options):
    """
    based on flags set by user and the taskType,
    build hashes that map option flags to what type of
    information follows them, eg {'-i': 'in', '-o': 'out', ...}
    """
    positionalArgs={}
    flaggedArgs={}
    #if options.inputFlag is None and options.taskType is not None:
    #    options.inputFlag=programOptionMap[options.taskType].get('in',None)
    if options.inputFlag is not None:
        try:
            positionalArgs[int(options.inputFlag)]='in'
        except ValueError:
            flaggedArgs[options.inputFlag]='in'
        except TypeError:
            for flag in options.inputFlag:
                flaggedArgs[flag]='in'
    #if not(options.outputFlags) and options.taskType is not None:
    #    options.outputFlags=programOptionMap[options.taskType].get('out',[])
    if options.outputFlags is not None:
        for outputFlag in options.outputFlags:
            try:
                positionalArgs[int(outputFlag)]='out'
            except ValueError:
                flaggedArgs[outputFlag]='out'
            except TypeError:
                for flag in outputFlag:
                    flaggedArgs[flag]='out'
    #if not(options.threadsFlag) and options.taskType is not None:
    #    options.threadsFlag=programOptionMap[options.taskType].get('threads',None)
    if options.threadsFlag is not None:
        try:
            positionalArgs[int(options.threadsFlag)]='threads'
        except ValueError:
            flaggedArgs[options.threadsFlag]='threads'
        except TypeError:
            for flag in options.threadsFlag:
                flaggedArgs[flag]='threads'
    if options.prefixFlag is not None:
        try:
            positionalArgs[int(options.prefixFlag)]='prefix'
        except ValueError:
            flaggedArgs[options.prefixFlag]='prefix'
        except TypeError:
            for flag in options.prefixFlag:
                flaggedArgs[flag]='prefix'
    if options.rel_paths is not None:
        for rel_path_flag in options.rel_paths:
            try:
                positionalArgs[int(rel_path_flag)]='rel'
            except ValueError:
                flaggedArgs[rel_path_flag]='rel'
 
    return (positionalArgs,flaggedArgs)

def getOutputFiles(options, cmdargs):
    """
    return map of (final) output files from this command
    ...keyed on the flag that indicated
    """
    outFileByFlag={}

    logging.debug("Command:\n%s\ninfile: %s\noutfile: %s" % (cmdargs, options.inputFlag, options.outputFlags))

    # setup hashes to look for options
    # returns hash mapping position or flag to 'in' or 'out'
    (positionalArgs,flaggedArgs)=getOptionHashes(options)
    logging.debug("positionalArgs: %s\nflaggedArgs: %s" % (positionalArgs,flaggedArgs))

    # scan command
    namedOutCount=0
    nextArgument = None
    usedFlags=[]
    infile=None
    prefix=None
    for i in range(len(cmdargs)):
        # check for positional arguments if there wasn't a flag in the previous spot
        if nextArgument is None:
            nextArgument=positionalArgs.get(i,None)
            usedFlags.append(i)

        # is this something we're looking for
        if nextArgument is None:
            # no. Well,is it a flag for something?
            arg=cmdargs[i]
            nextArgument=flaggedArgs.get(arg,None)
            if nextArgument is not None:
                usedFlags.append(arg)
            continue

        if nextArgument=='in':
            infile=cmdargs[i]
            logging.debug("Found infile (%s) in arg %d" % (infile,i))
        elif nextArgument=='prefix':
            prefix=cmdargs[i]
            logging.debug("Found prefix (%s) in arg %d" % (prefix,i))
        elif nextArgument=='out':
            if cmdargs[i] != '/dev/null':
                namedOutCount+=1
                outFileByFlag[namedOutCount]=cmdargs[i]
                logging.debug("Found outfile number %d in arg %d" % (namedOutCount,i))
        elif nextArgument=='threads':
            pass
        elif nextArgument=='rel':
            pass
        else:
            # something went wrong
            raise Exception("Unrecognized nextArgument: %s" % (nextArgument))

        # reset
        nextArgument=None

    # get path for infile
    if infile is not None:
        (inpath,infilename) = os.path.split(infile)
    else:
        (inpath,infilename) = ('.','stdin')

    # were there any outputFlags not in the command? These indicate secondary files
    if options.outputFlags is not None:
        for flag in options.outputFlags:
            # is it a positional argument (IE is it an integer)?
            try:
                position=int(flag)
                # Yes, it is. Skip it
                continue
            except ValueError:
                # No, it is not. Move on
                pass

            # is it an option flag (starts with dash)?
            if flag.startswith('-'):
                # I guess this could be a file, but let's ignore anyway
                continue

            if flag not in usedFlags:
                if flag=='%stdout':
                    outFileByFlag[flag] = None
                else:
                    (fOutDir,fOutFile) = getOutputFromFlag(flag,infilename,prefix,os.path.curdir,inpath)
                    outFileByFlag[flag] = os.sep.join((fOutDir,fOutFile))
    else:
        outFileByFlag['%stdout']=None

    return outFileByFlag

def prepareCommandForFragment(options, infragment, prefix, outLocal, cmdargs, hostname, errStream):
    """
    Replace input and output file names (and thread count) with fragment versions

    Translate any other (non i/o) paths to absolute paths.

    return:
        boolean indicating if input file found or not
        list of outputFileFlags. Each is one of:
            a file name (with no path)
            a modification of the input file (starts with %, e.g. %e.ext)
            an integer indicating order in cmdargs
    """

    foundI=False
    logging.debug("Command:\n%s\ninfile: %s\noutfile: %s\ninflag: %r"
                  "\noutflags: %r\nprefixFlags: %r\nrel_paths: %r",
                  cmdargs,
                  infragment,
                  outLocal,
                  options.inputFlag,
                  options.outputFlags,
                  options.prefixFlag,
                  options.rel_paths)

    # setup hashes to look for options
    (positionalArgs,flaggedArgs)=getOptionHashes(options)
    logging.debug("Looking for args: \n%r\n%r" % (positionalArgs,flaggedArgs))

    # scan command
    usedFlags=[]
    namedOutCount=0
    nextArgument = None
    #logging.debug("%d args" % (len(cmdargs)))
    for i in range(len(cmdargs)):
        #logging.debug("%d:%s)" % (i,cmdargs[i]))
        # check for positional arguments
        if nextArgument is None:
            nextArgument=positionalArgs.get(i,None)
            usedFlags.append(i)

        # is this something we're looking for
        if nextArgument is None:
            # no. Well,is it a flag for something?
            arg=cmdargs[i]
            nextArgument=flaggedArgs.get(arg,None)
            if nextArgument is not None:
                usedFlags.append(arg)
            continue

        if nextArgument=='in':
            foundI=True
            cmdargs[i]=infragment
            #logging.debug("Replaced infile at arg %d" % (i))
        elif nextArgument=='prefix':
            cmdargs[i]=prefix
            #logging.debug("Replaced prefix at arg %d" % (i))
        elif nextArgument=='out':
            if cmdargs[i] != '/dev/null':
                namedOutCount+=1
                cmdargs[i]="%s.%s" % (outLocal,namedOutCount)
                #logging.debug("Replaced outfile at arg %d" % (i))
        elif nextArgument=='rel':
            old_path = cmdargs[i]
            new_path = os.path.abspath(old_path)
            logging.debug("changing %s to %s", 
                          old_path, new_path)
            cmdargs[i] = new_path
        elif nextArgument=='threads':
            if options.queue != LOCAL:
                threads = getThreadCountForNode(hostname,errStream)
                logging.debug("%s can take %d threads" % (hostname,threads))
                cmdargs[i]=str(threads)
        else:
            # something went wrong
            logging.warning("Unrecognized nextArgument: %s" % (nextArgument))
            raise Exception("Unrecognized nextArgument: %s" % (nextArgument))

        # reset
        nextArgument=None

    # were there any outputFlags not in the command? These indicate secondary files
    #logging.debug("Done with positional arguments in prepareCommandForFragment")
    otherOutputs=list(range(1,namedOutCount+1))
    if options.outputFlags is not None:
        for flag in options.outputFlags:
            if isinstance(flag, str) and flag not in usedFlags:
                otherOutputs.append(flag)
    else:
        otherOutputs.append('%stdout')

    return (foundI, otherOutputs)

def getOutputFromFlag(flag, inputName, prefix, localDir, fragDir):
    if flag is None or flag=="":
        raise Exception("cannot have an empty output flag!")
    if isinstance(flag,int):
        # numbered output, return corresponding file name
        return (localDir,"%s.output.%d" % (inputName,flag))
    if flag[0]=='%':
        if flag=='%stdout':
            # code calling this should do something special for STDOUT
            return None

        # anything else starting with % is a translation of the input file name
        (name,ext)=os.path.splitext(inputName)
        if len(flag)<=2:
            raise Exception("illegal output flag: %s" % (flag))
        if flag[1]=='e':
            newext=flag[2:]
            return (fragDir,name+newext)
        if flag[1]=='E':
            newext=flag[2:]
            return (localDir,name+newext)
        if flag[1]=='a':
            newext=flag[2:]
            return (fragDir,inputName+newext)
        if flag[1]=='A':
            newext=flag[2:]
            return (localDir,inputName+newext)
        if flag[1]=='p':
            newext=flag[2:]
            return (localDir,prefix+newext)
        if flag[1]=='s':
            raise Exception("Sorry, not yet implemented!")
    else:
        # it's just a file name
        return (localDir, flag)

def getDBType(database):
    for ext in ('pal','nal','pin','nin','psd','nsd'):
        if os.path.exists("%s.%s" % (database,ext)):
            return '%sal' % (ext[0])
    else:
        return None

def getThreadCountForNode(hostname,errStream,queue=SGE):
    """
    figure out how many threads the job can use

    Calls system specific (SGE or SLURM) version of function.
    """
    if queue==SGE:
        return getThreadCountForSGENode(hostname, errStream)
    elif queue==SLURM:
        return getThreadsCountForSLURMNode(hostname, errStream)
    else:
        logging.warning("Unrecognized queue (%s), using 8 cores" % (queue))
        return 8

def getThreadCountForSLRUMNode(hostname, errStream):
    """
    figure out how many threads the job can use
    """
    qhcmd = ["sinfo","-n", hostname, "-o", '"%15N %10c"']
    process=subprocess.Popen(qhcmd,stdout=subprocess.PIPE)
    sinfoCPUsRE = re.compile(r'^\S+\s+(\d+)')
    qhout=""
    for line in process.stdout:
        qhout+=line
        m=sinfoCPUsRE.search(line)
        if m:
            slots = int(m.group(1))
            logging.debug("Node %s has %d slots" % (hostname, slots))
            break
    else:
        slots=8
        logging.warning("Could not parse sinfo output:\n%s" % (qhout))
    return slots

def getThreadCountForSGENode(hostname, errStream):
    """
    figure out how many threads the job can use
    """
    qhcmd = ["qhost","-h", hostname]
    process=subprocess.Popen(qhcmd,stdout=subprocess.PIPE)
    qhostSlotsRE = re.compile(r'^\S+\s+\S+\s+(\d+)\s+')
    qhout=""
    for line in process.stdout:
        qhout+=line
        m=qhostSlotsRE.search(line)
        if m:
            slots = int(m.group(1))
            logging.debug("Node %s has %d slots" % (hostname, slots))
            break
    else:
        slots=8
        logging.warning("Could not parse qhost output:\n%s" % (qhout))
    return slots

def getTaskId(options):
    if options.queue==SGE:
        return getSGETaskId()
    else:
        return getSLURMTaskId()

def getSLURMTaskId():
    logging.debug("SLURM_ARRAY_TASK_ID: %s" % os.environ.get('SLURM_ARRAY_TASK_ID',None))
    logging.debug("SLURM_ARRAY_JOB_ID: %s" % os.environ.get('SLURM_ARRAY_JOB_ID',None))
    logging.debug("SLURM_JOBID: %s" % os.environ.get('SLURM_JOBID',None))
    return int(os.environ.get('SLURM_ARRAY_TASK_ID',1))

def getSGETaskId():
    taskid = int(os.environ.get('SGE_TASK_ID',1))
    logging.debug("getting SGE task id (%d)" % (taskid))
    return taskid

def getHostName():
    try:
        return os.uname()[1].split('.')[0]
    except:
        try:
            return os.uname()[1]
        except:
            logging.warning("Could not parse hostname from %s" % (str(os.uname())))
            return 'NODE'

##################
# task specific code
##################
## FR-HIT
def inspectFrHitCommand(command,taskType,sgeOptions,commandBin,batchOptions):
    """
    get db size of input database and add as parameter to commmand (if not already there)
    add sge option to throttle this task since it's all over network disks
    """

    logging.info("Looking for reference db")
    nextWordIs=None
    refDB=None
    refDBSize=None
    defaultValues=defaultsForTask[taskType]
    for word in command:
        logging.debug("Word is %s" % word)
        if nextWordIs is None:
            if word=='-d':
                nextWordIs='db'
            if word=='-R':
                nextWordIs='dbsize'
            elif word in defaultValues:
                defaultValues.pop(word)
        else:
            if nextWordIs=='db':
                refDB=word
            elif nextWordIs=='dbsize':
                refDBSize=word
            nextWordIs=None
        logging.debug("next word is: %s" % nextWordIs)

    # apply anydefaults not already in command
    for kvPair in defaultValues.items():
        command.extend(kvPair)

    # get total bases in reference db
    if refDB is None:
        raise Exception("You must supply a database to run fr-hit")

    if refDBSize is not None:
        logging.warning("You supplied ref DB size of %s. If you omit the -R option batch_launcher will calculate the db size for you." % (refDBSize))
    else:
        dbInfo = countBasesInFasta(refDB)
        logging.info("Reference db (%s) has %s bases in %s records" % (refDB,dbInfo['bases'],dbInfo['records']))
        command.extend(('-R',str(dbInfo['records']),'-B',str(dbInfo['bases'])))

        # while we know the db size, lets calculate chunk size
        if batchOptions.chunk is None:
            # if the user hasn't set the chunk size, always size chunks by bases
            batchOptions.splitOnSize=True
            dbsize=dbInfo['bases']
            if batchOptions.splits is None:
                # set chunk to max for node RAM (and calculate splits)
                batchOptions.splits = ceil(float(dbsize)/DEFAULT_FRHIT_CHUNK)
                # next, re-adjust chunk so that fragments are similar sizes
            batchOptions.chunk = calculateChunkSize(dbInfo['bases'],batchOptions.splits)
        else:
            if not batchOptions.splitOnSize:
                logging.warning("Are you sure you want to split on number of records? It usually is a good idea to split on number of bases (-s)")

def inspectLastCommand(command,taskType,sgeOptions,commandBin,batchOptions):
    # apply the defaults
    applyDefaultsToCommand(command,taskType,prepend=True)

## Blast
def inspectDCMBCommand(command,taskType,sgeOptions,blastBin,batchOptions):
    # replace first element of command string ('dcmegabalst') with 'blastn -task ...'
    command[0:1]=['blastn','-task','dc-megablast']
    inspectBlastCommand(command,taskType,sgeOptions,blastBin,batchOptions)

def inspectBlastCommand(command,taskType,sgeOptions,blastBin,batchOptions):
    """
    if using blastall, make sure prog/alg name is ok (one of: blastn, blastx, etc)
    """

    logging.info("Looking for blast options")
    nextWordIs=None
    dbs=[]
    program=None
    defaultValues=defaultsForTask[taskType]
    if blastBin in ['rpsblast','rpstblastn']:
        # this is at best ignored and sometimes causes errors, remove it
        defaultValues.pop('-num_threads')

    for word in command:
        logging.debug("Word is %s" % word)
        if nextWordIs is None:
            if word=='-d' or word=='-db':
                nextWordIs='db'
            elif word=='-p':
                nextWordIs='program'
            elif word in defaultValues:
                defaultValues.pop(word)
        else:
            if nextWordIs=='db':
                dbs.append(word)
            elif nextWordIs=='program':
                program=word
            nextWordIs=None
        logging.debug("next word is: %s" % nextWordIs)

    # apply anydefaults not already in command
    for kvPair in defaultValues.items():
        command.extend(kvPair)

    if dbs:
        pass
    else:
        raise Exception("You must supply a database to run blastall")

    if blastBin=='blastall':
        if program in ('rpsblast','rpsblastx','rpstblastn'):
            raise Exception("Use the rpsblast binary, not blastall!")
        if program not in ('blastx','blastn','blastp','tblastn','tblastx'):
            raise Exception("Blast algorithm %s is not recognized" % (program))
    elif blastBin in ['rpstblastn','rpsblast']:
        pass

def makeRelativePathsAbsolute(cmdargs):
    """
    Translate any relative paths (starting with ./ or ../) to absolute
    """
    for i in range(len(cmdargs)):
        if relativePathRE.match(cmdargs[i]):
            cmdargs[i]=os.path.abspath(cmdargs[i])

def mergeDBs(cmdargs, localDir):
    """
    Finds database arguments in a blast(blastall, rpsblast, blast+) command and:
        * creates an alias file to search against if more than one DB given
            NB: this could fail if there is a nuc and pro db with same base name
            alias file is created in localtmp dir and will get deleted with it
    Also enforce defauls for b,v, and a
    """
    logging.debug("mergeDBs.")
    # extract databses from command
    dbs=[]
    nextIsDB=False
    i=0
    dbFlag=None
    foundExtensions=[]
    while i < len(cmdargs):
        if not nextIsDB:
            if cmdargs[i] in ('-d','-db'):
                nextIsDB=True
                dbFlag=cmdargs.pop(i)
                i-=1
        else:
            nextIsDB=False
            # remove db from args
            db=cmdargs.pop(i)
            dbs.append(db)
            dbtype=getDBType(db)
            if dbtype is not None:
                foundExtensions.append(dbtype)
            i-=1
        i+=1

    if len(dbs)>1:
        # determine db type (nuc/prot) based on found extensions of first DB
        if foundExtensions:
            ext = foundExtensions[0]
        else:
            ext = getDBType(dbs[0])

        # create temporary alias
        tmpAlias = "%s%sTemporaryBlastAlias" % (localDir,os.sep)
        tmpAliasFile = "%s.%s" % (tmpAlias,ext)
        tmpHandle=open(tmpAliasFile,'w')
        tmpHandle.write("""#
# Alias file created %s
#
#
TITLE TemporaryBlastAlias
#
DBLIST %s
""" % (datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S UTC %Y"),' '.join(dbs)))

        tmpHandle.close()
        dbs=[tmpAlias,]
        logging.info("Created database alias: %s" % (tmpAlias))

    # restore modified db arguments to command
    cmdargs.extend((dbFlag,dbs[0]))

def applyDefaultsToCommand(command,taskType,prepend=False):
    logging.info("applying defaults to %s" % taskType)
    defaultValues=defaultsForTask.get(taskType,None)
    if defaultValues is None:
        return

    # scan elements of command
    for word in command:
        logging.debug("Word is %s" % word)
        if word in defaultValues:
            defaultValues.pop(word)

    # apply anydefaults not already in command
    if prepend:
        index=1
        for kvPair in defaultValues.items():
            command.insert(index,kvPair[1])
            command.insert(index,kvPair[0])
            index+=2
    else:
        for kvPair in defaultValues.items():
            command.extend(kvPair)

##################
# Constants
##################
BATCHLAUNCHER=os.path.abspath(__file__)
PYTHON=sys.executable
#BATCHLAUNCHER=os.sep.join([BINDIR,'batch_launcher.py'])
scriptingLanguages=['perl','python','ruby','sh', 'bash']
# task types
BLAST='blast'
DCMB='dcmegablast'
BLASTPLUS='blast+'
METARNA='meta_rna'
FGENESB='fgenesb'
GLIMMERMG='glimmermg'
FRHIT='frhit'
LAST='lastal'
HMMER='hmmer'
VIRSORTER='virsorter'
programOptionMap={BLAST:{'in':'-i',
                         'out':['-o'],
                         'rel':['-d'],
                         'threads':'-a',
                        },
                  DCMB:{'in':'-query',
                        'out':['-out'],
                        'threads':'-num_threads'
                       },
                  BLASTPLUS:{'in':'-query',
                         'out':['-out'],
                             'rel':['-db'],
                         'threads':'-num_threads'
                        },
                  VIRSORTER:{'in': '-f',
                             'prefix': '--wdir',
                             'out': ['%p/VIRSorter_global-phage-signal.csv',
                                     '%p/Predicted_viral_sequences/VIRSorter_cat-1.fasta',
                                     '%p/Predicted_viral_sequences/VIRSorter_cat-2.fasta',
                                     '%p/Predicted_viral_sequences/VIRSorter_cat-3.fasta',
                                     '%p/Predicted_viral_sequences/VIRSorter_prophages_cat-4.fasta',
                                     '%p/Predicted_viral_sequences/VIRSorter_prophages_cat-5.fasta',
                                     '%p/Predicted_viral_sequences/VIRSorter_prophages_cat-6.fasta',
                                     '%p/Metric_files/VIRSorter_phage_signal.tab',
                                     '%p/logs/err',
                                     '%p/logs/out',],
                             'threads': '--ncpu',
                            },
                  METARNA:{'in':'-i',
                           'out':['%p.coord','%p.mask','%p.seq','%e.coord','%e.mask','%e.seq'],
                           'prefix':'-o',
                           'threads':'-p',
                          },
                  LAST:{'in':2,
                        'rel':[1,],
                        'out':['%stdout']},
                  GLIMMERMG:{'in':1,
                             'out':['%p.predict','%E.predict'],
                             'prefix':'-o',
                             'threads':'-p'},
                  FGENESB:{'in':'-i',
                           'out':['-o'],
                          },
                  FRHIT:{'in':'-d',
                         'out':['-o'],
                         'threads':'-T'},
                  HMMER:{'in':1,
                         'rel':[2,],
                         'out':['-o','--domtblout','--tblout'],
                         'threads':'--cpu'},
                 }
blastPlusProgRE=re.compile(
        r'(blastn|blastp|blastx|tblastx|tblastn|rpsblast|rpstblastn)$')
blastProgRE=re.compile(r'(blastall|megablast)$')
dcmbProgRE=re.compile(r'^(dcmegablast)$')
mrnaProgRE=re.compile(r'(rna_hmm3\.py)$')
fgbProgRE=re.compile(r'(softberry_wrapper\.pl|fgenesb)$')
glimmermgRE=re.compile(r'(glimmer-mg\.py)$')
frHitRE=re.compile(r'(fr-hit)$')
lastRE=re.compile(r'(lastal)$')
hmmerRE=re.compile(r'(hmmsearch|hmmscan)$')
virsorterProgRE = \
    re.compile(r'^(wrapper_phage_contigs_sorter_iPlant.pl|virsorter)$')
taskTypePatterns={BLAST:blastProgRE,
                  BLASTPLUS:blastPlusProgRE,
                  VIRSORTER:virsorterProgRE,
                  DCMB:dcmbProgRE,
                  METARNA:mrnaProgRE,
                  FGENESB:fgbProgRE,
                  GLIMMERMG:glimmermgRE,
                  FRHIT:frHitRE,
                  LAST:lastRE,
                  HMMER:hmmerRE,
                 }
fragment_prep_for_task={HMMER: ".h3i: hmmpress {}",}
inspectCommandForTaskType={BLAST:inspectBlastCommand,
                           BLASTPLUS:inspectBlastCommand,
                           DCMB:inspectDCMBCommand,
                           FRHIT:inspectFrHitCommand,
                           LAST:inspectLastCommand,
                          }
finalCommandProcessingForTask={BLAST:mergeDBs,
                               BLASTPLUS:mergeDBs,
                               DCMB:mergeDBs,
                              }
defaultsForTask={BLAST:{'-b':'10','-v':'10','-a':'8'},
                 DCMB:{'-max_target_seqs':'10','-num_threads':'8','-word_size':'12','-template_length':'21','-template_type':'optimal'},
                 BLASTPLUS:{'-max_target_seqs':'10','-num_threads':'8'},
                 METARNA:{'-H':'/common/bin/hmmer-3.0',
                          '-L':'/common/bin/meta_rna/rRNA_hmm_fs_wst_v0/HMM3'},
                 FRHIT:{'-r':'25','-T':'0'},
                 LAST:{'-f':'BlastTab'},
                }
unsupportedOptions={}
flagsWithArguments={GLIMMERMG:['--iter','-p','-t','-i','-q','-r','-s','-u','--fudge','--taxlevel','--minbp_pct'],
                    LAST:['-P','-a','-b','-c','-d','-e','-F','-p','-q','-r','-x','-y','-z','-f','-k','-l','-m','-n','-s','-i','-u','-t','-j','-Q','-g','-G','-o'],
                    HMMER:["-o", "-A", "--tblout", "--domtblout",
                           "--pfamtblout", "--textw", "-E", "-T", "--domE",
                           "--domT", "--incE", "--incT", "--incdomE",
                           "--incdomT", "--F1", "--F2", "--F3", "-Z", "--domZ",
                           "--seed", "--tformat", "--cpu"],
                    VIRSORTER:["-d", "--cp", "--db", "--ncpu", "--data_dir",
                               "-f"],
                   }

taskSpecificCopy={}
# code is untested. Turning off here.
#BLAST:addBlastFilesToStream,
#                  BLASTPLUS:addBlastFilesToStream}
missedRecordGenerator={}
recordIDRE={BLAST:fileTypeMap['fasta'].sepRE,
            BLASTPLUS:fileTypeMap['fasta'].sepRE,
            DCMB:fileTypeMap['fasta'].sepRE,
           }

DEFAULT_SPLITS=400
DEFAULT_FRHIT_CHUNK=250000000

# constants for parsing blast
START='start'
HEADER='header'
HITS='hits'
blastHeaderRE=re.compile(r'^[T]?BLAST[NXP]\s+\d+\.\d+\.\d')
blastRecordRE=re.compile(r'^\s*Query=\s+(\S+)')
foundHitsRE=re.compile(r'^(?:Sequences\sproducing|>\S+)')
noHitsRE=re.compile(r'^\s+\*+\s+No\shits')

##############
# Classes
##############
class FragmentThread( threading.Thread ):
    """
    For local multi threading

    This might be better in multiprocessing, but since we're just spawning
    processes, the global lock isn't too much of an issue
    """
    def __init__(self, queue):
        self.queue = queue
        threading.Thread.__init__ ( self )

    def run(self):
        while True:
            # get next fragment from queue
            try:
                task = self.queue.get(True,10)
                task.run()
                self.queue.task_done()
            except queue.Empty:
                # quit if the queue is empty
                break
            

class FragmentTask():
    def __init__(self, options, cmd, taskNum, errorQueue):
        self.options=options
        self.cmd=list(cmd)
        self.taskNum=taskNum
        self.errorQueue=errorQueue

    def run(self):
        logging.debug("Launching job: %s" % (self.cmd))
        exitcode=runFragment(self.options,self.cmd,taskNum=self.taskNum)
        if exitcode!=0:
            logging.warning("Task %d exited with %d" % (self.taskNum,
                exitcode))
            self.errorQueue.put(exitcode)

##################
# Expressions
################
relativePathRE=re.compile(r'^\.\.?\/')

if __name__ == '__main__':
    main()

