#!/usr/bin/python
#$ -S /usr/bin/python
"""
Run (almost) any command accross an SGE cluster as long as the input file can be broken into pieces (using regex). The resulting output fragments will simply be concatenated.

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
#SGEBIN='/common/sge/bin/darwin'
SGEBIN='/usr/bin'
QSUB=SGEBIN+'/qsub'
#BINDIR='/Users/jmeppley/work/delong/projects/scripts'
BINDIR='/slipstream/home/jmeppley/work/delong/projects/scripts'
#BINDIR='/slipstream/opt/scripts'
#BINDIR='/common/bin/scripts'

from optparse import OptionParser
import sys, re, logging, os, tempfile, subprocess, shutil, traceback, datetime, shlex, time, fileinput, io, threading, Queue
sys.path.append(BINDIR)
from numpy import ceil

def main():
    ## set up CLI
    usage = "usage: %prog [options] -- [command to run on cluster]"
    description = """
Given an input file with multiple records, an output file, and a command: process the records in parallel across an SGE cluster.
    """

    parser = OptionParser(usage, description=description)
    # options for processing command
    parser.add_option("-G", "--ignore", default=False, action='store_true',
    help="Do not perform checks based on the program name (in the command). Batch_launcher will skip need you to indicate input/output/threads with -i,-o, and -t, and it will skip any program specific processing (e.g. db checks for blast)")
    parser.add_option("-i", "--inputFlag", metavar="FLAG",
                      help="Indicate where in command to find the input file to fragment. The value can be an option flag (eg '-i' or '-fasta') or a number indicating a positional argument. Only needed if program in command is not recognized.")
    parser.add_option("-o", "--outputFlags", action='append',metavar='FLAG',
                      help="""Indicate where in command (or not) to find the output file.
The value can be an option flag or positional argument (as with -i, above), a file
name (if the command always produces the same output file), '%e.EXT' if the output file
is the input file with a new extension (e.g. %e.faa), '%a.ext' if the output file is the
input fie with an additional extension, or '%s/foo/bar/' if the output file is the input
file with a regexp substitution. Multiple values are permitted to account for multiple
outputs.""")
    parser.add_option("-p", "--prefixFlag",
                      help="Indicate where in command to find a prefix for output files. Output extension flags '%p.ext' indicate that the output file will be PREFIX.ext.")
    parser.add_option("-t", "--threadsFlag", metavar="FLAG",
                      help="Option to use to tell command how many threads to use")

    # ways to customize batch behavior
    addFragmentingOptions(parser)
    parser.add_option("-X", "--local", default=False, action="store_true",
                      help="""Don't use SGE, just fragment input and use
                      subprocess to spawn jobs""")
    parser.add_option("-R", "--throttle", default=0, type='int',
            help="Limit number of simultaneously executing fragemnts to given number. Default: 0 => unlimited")
    parser.add_option("-B", "--bathybius", default=False, action="store_true",
			help="Use bathybius specific SGE flags for throttling jobs and reserving resources (rps, greedy, bigd, etc)")
    parser.add_option("-w", "--wait", default=False, action='store_true',
                      help="Wait for jobs to finish before exiting script")
    parser.add_option("-c","--cwd",dest='cwd',
                      default=False,action='store_true',
                      help='Run fragment in current directory, otherwise it will be executed in the tmp location (as configured in Python, usu a local disk like /tmp) of the node (which is good for programs like softnerry that create lots of temporary files)')
    parser.add_option("-r","--retries",default=-1,type='int',
                      help='number of times to resubmit failed tasks. Less than zero (the default) means continue as long as some new results come through each time')

    parser.add_option("-j", "--jobName", metavar="STRING",
                      help="String for naming SGE tasks")
    parser.add_option("-d", "--tmp_dir", dest='tmpDir',
                      help="Temporary directory for files")
    parser.add_option("-g", "--greedy", type='int', default=None,
                      help="How many resources does each fragment need (1-24, default=24). A good rule of thumb is 1 per 250MB of RAM needed. Also using the maximum claims the whole node which is good for multithreaded programs that can use all the processors.")
    parser.add_option("-S", "--sgeOptions", default=None, dest='sgeOptions',
                      help="Option string to add to SGE command")
    parser.add_option("-n", "--priority", default=0, type='int',
                      help="Adjust priority of sub-tasks")

    # primarily used when calling self via SGE
    parser.add_option("-f", "--frag_base", dest="fragBase", default="file.fragment",
                     help="naming base for input file fragments")
    parser.add_option("-m", "--mode", default='launch', metavar='MODE',
                      choices=['launch','run','cleanup'],
                      help="Only used internally to lauch tasks in SGE")
    parser.add_option("-Z","--taskType",default=None, choices=taskTypePatterns.keys(),
                      help="only for task running: what type of task is this? Will be ignored in inital call")

    # other
    parser.add_option("-l", "--logToFile", default=False, action="store_true",
                      help="print logging messages to file")
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbose", default=1,
                      help="Print log messages. Use twice for debugging")
    parser.add_option("-q", '--quiet', dest='verbose',
                      action="store_const", const=0,
                      help="Suppress warnings. Only print fatal messages")
    parser.add_option("-V","--loglevel", type='int', dest="verbose",
                     help="Shortcut to set verbosity directly")
    parser.add_option("-A", "--about",
                  action="store_true", dest="about", default=False,
                  help="Print description")


    (options, cmdargs) = parser.parse_args()

    if options.about:
        print description
        exit(0)

    # If running a sub-job, do so before logging is set up
    if options.mode=='run':
        exitcode=runFragment(options, cmdargs)
        sys.exit(exitcode)

    # if running launcher or cleanup, standard logging
    if options.verbose==0:
        loglevel=logging.ERROR
    elif options.verbose==1:
        loglevel=logging.WARN
    elif options.verbose==2:
        loglevel=logging.INFO
    elif options.verbose>=3:
        loglevel=logging.DEBUG

    if options.logToFile:
        logStream=io.BytesIO()
    else:
        logStream=sys.stderr
    logging.basicConfig(stream=logStream, level=loglevel)

    # cleanup only?
    if options.mode=='cleanup':
        logging.warn("Starting cleanup at: %s" % (datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S UTC %Y")))
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
        logging.debug('Task type is: %s' % (options.taskType))

        # if any flag not supplied by user, check for defaults based on taskType
        # first, there must be a task type:
        #if options.taskType is None:
        #    if (options.inputFlag is None or not(options.outputFlags)):
        #        parser.error("I'm sorry the command is not recognized. Please use the -i and -o options to specify the input and output flags in your command.")

    else:
        logging.warn('User chose to run command without modification')

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

    if not options.local:
        if options.bathybius:
            # if one of the task specific methods hasn't set greedy:
            for sgeOpt in sgeOptions:
                if len(sgeOpt)>7 and sgeOpt[:7]=='greedy=':
                    logging.debug("greedy already set: %s" % sgeOpt)
                    options.greedy=int(sgeOpt[7:])
                    break
            else:
                # set greedy
                logging.debug('greedy not set in %s' % (sgeOptions))
                if options.greedy is None:
                    # if task type has a greedy value, use it, otherwise use default
                    options.greedy=greedyByTask.get(options.taskType,GREEDY_DEFAULT)
                sgeOptions.extend(['-l', 'greedy=%d'%options.greedy])

            # limit total number of jobs if more than one job can run on a single node
            if options.greedy <= GREEDY_DEFAULT/2:
                # but not if rps is already set
                for sgeOpt in sgeOptions:
                    if (len(sgeOpt)>4 and sgeOpt[:4]=='rps=') or (len(sgeOpt)>8 and sgeOpt[:8]=='rpsjobs='):
                        break
                else:
                    sgeOptions.extend(['-l','rps=%d'%(options.greedy)])

        options.sgeOptions = sgeOptions
        logging.debug("Adding the following to SGE task array command: '%s'" % sgeOptions)

        # must use 'wait' to print to stdout
        if not(options.outputFlags) or '%stdout' in options.outputFlags:
            if not(options.wait):
                parser.error("Piping output to stdout only works if we wait for the jobs to finish. Please use wait (-w) or fix your output flags (-o)")


    fileType = getFileType(options,infile)
    # if we had to get from file name, save:
    if options.infileType is None and options.pattern is None and options.numLines is None:
        options.infileType = fileType.name

    # get temporary dir
    options.tmpDir = checkTmpDir(options.tmpDir,options.jobName)

    # get fragment size if we haven't been told
    if options.chunk is None:
        if options.splits is None:
            options.splits=DEFAULT_SPLITS
        options.chunk = getSizePerChunk(infile,options.splits,fileType,splitOnSize=options.splitOnSize)

    # fragment input file
    options.splits=fragmentInputBySize(infile, options.tmpDir, options.chunk, fileType, options.fragBase, splitOnSize=options.splitOnSize)

    # launch task array
    launchJobs(options,cmdargs,errStream=logStream)

    # launch cleanup
    if options.wait or options.local:
        # if wait, then we already waited for jobs to finish, just run cleanup
        exitcode=cleanup(options, cmdargs)
        sys.exit(exitcode)
    else:
        # otherwise, launch cleanup task to wait for processing to finish
        launchCleanup(options, cmdargs, errStream=logStream)

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

        if len(relativePositions):
            translatePositionalArgs(options,cmdargs,relativePositions)

    (positionalArgs,flaggedArgs)=getOptionHashes(options)
    logging.debug("positional args: %s form %s" % (positionalArgs,options.inputFlag))

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
            logging.debug("Found infile (%s) in arg %d" % (infile,i))
        elif nextArgument=='prefix':
            prefix=cmdargs[i]
            logging.debug("Found prefix (%s) in arg %d" % (prefix,i))
        elif nextArgument=='out':
            outfile=cmdargs[i]
            logging.debug("Found output (%s) in arg %d" % (outfile,i))
        elif nextArgument=='threads':
            pass
        else:
            # something went wrong
            raise Exception("Unrecognized nextArgument: %s" % (nextArgument))

        # reset
        nextArgument=None

    # remove any prefix based output flags if no pref
    if prefix is None and options.outputFlags is not None:
        i=0
        while i<len(options.outputFlags):
            if len(options.outputFlags[i])>1 and options.outputFlags[i][:2]=='%p':
                options.outputFlags.pop(i)
            else:
                i+=1

    return (infile,outfile)

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
    for i in xrange(start,len(cmdargs)):
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
        if command is blastall and dbs NOT in /common/data, return big disk flag

    returns list of SGE options
    """
    if len(command)==0:
        raise Exception("No command given. Add ' -- ' followed by the command you'd like to run!")

    logging.debug("Command: %r" % (command))

    # get program name
    program=command[0]
    if program in scriptingLanguages:
        # if first word in command is perl or python (etc), take second word
        if len(command)==1:
            raise Exception("%s won't do much without a script name" % (program))
        program = command[1]

    # do we know this program?
    for (taskType,taskProgRE) in taskTypePatterns.iteritems():
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

def areDBsLocal(dbList):
    """
    return True only if all DBs are in /common/data
    """
    for db in dbList:
        if not commonDataRE.match(db):
            return False
    else:
        return True

def launchLocalJobs(options, cmdargs, errStream):
    logging.debug("Launching command array: %r" % ({'tmpDir':options.tmpDir,'splits':options.splits,'fragName':options.fragBase,'cmd':cmdargs,'loglevel':options.verbose, 'type':options.taskType, 'throttle':options.throttle}))
    
    numTasks = options.splits
    numThreads = min(numTasks,options.throttle) if options.throttle>0 else numTasks

    # Create thread safe object to report errors to
    errorQueue = Queue.Queue()

    # create tasks
    taskQueue = Queue.Queue()
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
    except Queue.Empty:
        # no errors
        return

    # There was at least one error
    # add this one to the remaining to get count
    errorCount = errorQueue.qsize()+1
    logging.warn("%d of %d tasks reported non zero exit codes!" % (errorCount, numTasks))
    sys.exit(exitcode)

def launchJobs(options, cmdargs, errStream=sys.stdin):
    """
    launch the SGE task array. Return the job ID

    splits can be a single interger or a pair of integers
    this will launch
         the indicated task number
         the range of tasks numbers indicated by the pair (inclusive)
    """

    if options.local:
        launchLocalJobs(options,cmdargs,errStream)
        return

    priority=options.priority

    logging.debug("Launching task array: %r" % ({'tmpDir':options.tmpDir,'splits':options.splits,'fragName':options.fragBase,'cmd':cmdargs,'sgeOpts':options.sgeOptions,'job':options.jobName,'priority':priority,'loglevel':options.verbose,'wait':options.wait, 'type':options.taskType}))

    # sge command
    command=[QSUB, "-hard", "-cwd", "-V"]
    command+=["-t", "This placeholder will be replaced later"]
    taskIndex=5
    if options.throttle>0:
        command+=["-tc",str(options.throttle)]
    command+=["-N",options.jobName]
    outfiles=["-o","%s%stasks.out" %(options.tmpDir, os.sep),"-e","%s%stasks.err" % (options.tmpDir, os.sep)]
    command+=outfiles
    priority-=10
    if priority>=0:
        priority=-1
    command+=["-p",str(priority)]
    if options.wait:
        command+=['-sync','y']
    if isinstance(options.sgeOptions,str):
        command+=shlex.split(options.sgeOptions)
    else:
        command+=options.sgeOptions

    # batch_runner command
    command.append(BATCHLAUNCHER)
    command+=["--mode","run","--tmp_dir",options.tmpDir,"--frag_base", options.fragBase, "--loglevel", str(options.verbose)]
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
    if options.verbose==0:
        qsubOuts=open(os.devnull,'w')
    else:
        qsubOuts=errStream

    if isinstance(options.splits,int):
        command[taskIndex]="1-%d" % options.splits
    else:
        if options.splits[0]==options.splits[1]:
            command[taskIndex]=str(options.splits[0])
        else:
            command[taskIndex]="%d-%d" % tuple(options.splits)

    # run command
    logging.debug('Launching task array: %s' % (formatCommand(command)))
    exitcode=subprocess.call(command, stdout=qsubOuts)
    if exitcode!=0:
        if options.wait:
            # when using -sync y, the exit code may come from a task
            #  (which cleanup will handle)
            logging.warn("qsub returned an error code of: %d" % exitcode)
        else:
            raise subprocess.CalledProcessError(exitcode, command)

def launchCleanup(options, cmdargs, errStream=sys.stderr):
    """
    launch the cleanup script, which is this script with the --cleanup flag
    """
    logging.debug("Launching cleanup: %s" % ({'tmpDir':options.tmpDir,'splits':options.splits,'fragBase':options.fragBase,'out':options.outputFlags,'job':options.jobName}))

    # name outputfiles
    outFileByFlag = getOutputFiles(options,cmdargs)
    logging.debug("Outfiles: %s" % (outFileByFlag))
    fileNameBase = getFileNameBase(options.outputFlags,outFileByFlag,options.jobName)
    logging.debug("File Name Base: %s" % (fileNameBase))

    # build command
    # sge
    if options.priority>0:
        options.priority=0
    command=[QSUB, "-cwd", "-N", "%s.cleanup" % options.jobName, "-hold_jid", options.jobName,"-o","%s.cleanup.out" %(fileNameBase),"-e","%s.cleanup.err" % (fileNameBase), '-p', str(options.priority)]
    # cleanup
    command+=[BATCHLAUNCHER,'--tmp_dir',options.tmpDir,'--frag_base',options.fragBase,'--mode','cleanup','--numChunks',str(options.splits), '--loglevel', str(options.verbose), '--retries', str(options.retries),'--jobName',options.jobName, '--chunk', str(options.chunk)]
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
        command.extend(['--sgeOptions'," ".join(options.sgeOptions)])
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
    outFileByFlag = getOutputFiles(options,cmdargs)

    # name outputfiles
    fileNameBase = getFileNameBase(options.outputFlags,outFileByFlag,options.jobName)

    # remove old output files
    errStreamFile="%s.stderr" % fileNameBase
    failureStreamFile="%s.failures" % fileNameBase
    for file in errStreamFile, failureStreamFile:
        if os.path.exists(file):
            logging.debug('Removing previous file: %s' % file)
            os.remove(file)
    for file in outFileByFlag.itervalues():
        if file is not None:
            if os.path.exists(file):
                logging.debug('Removing previous file: %s' % file)
                os.remove(file)

    # set up copy method (some task types might do some filtering)
    copyFilesToStream=taskSpecificCopy.get(options.taskType,addFilesToStream)

    # loop until everything is node or we give up
    taskIds=xrange(1,options.splits+1)
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
            fragName = getFragmentName(options.fragBase,i)
            prefix = getFragmentPrefix(options.fragBase,i)
            frag = "%s%s%s" % (options.tmpDir, os.sep, fragName)
            fragerr = "%s.exitcode" % (frag)
            outfrag = "%s.stdout" % (frag)
            errfrag = "%s.stderr" % (frag)
            logfrag = "%s.log" % (frag)
            outfragmap={}

            # For each configured output file, map fragment to final
            for (flag, flagOutFile) in outFileByFlag.iteritems():
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
                                errStream = open(errStreamFile, 'w')
                            addFileToStream(logfrag,errStream,header="## LOGGING from fragment %d:" % (i))
                            if outfrag not in outfragmap:
                                addFileToStream(outfrag,errStream,header="## STDOUT from fragment %d:" % (i))
                            addFileToStream(errfrag,errStream,header="## STDERR from fragment %d:" % (i))

                    # delete files (input, error, outputs)
                    for f in [frag, outfrag, errfrag, logfrag] + outfragmap.keys():
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
        for outstream in outfragmap.itervalues():
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
                nextTaskNum+=reFragmentMissedTasks(failedTasks,options)

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
                taskIds = xrange(1,nextTaskNum)

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
    os.rmdir(options.tmpDir)
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
    newFragNum=fragmentInputStreamBySize(failedRecordStream, temporaryLocation, options.chunk, fileType, options.fragBase, splitOnSize=options.splitOnSize)

    # remove old fragments
    for i in missedTasks:
        frag = getFragmentPath(options.tmpDir, options.fragBase, i)
        os.remove(frag)

    return newFragNum+1

def moveNewFragmentsToTmpDir(options,nextTaskNum):
    """
    Move new fragments from tmpDir/tmp/ into tmpDir/
    """
    for i in xrange(1,nextTaskNum):
        frag = getFragmentPath(options.tmpDir, options.fragBase, i)
        newfrag = getFragmentPath("%s%stmp" % (options.tmpDir, os.sep), options.fragBase, i)
        os.rename(newfrag,frag)
    os.rmdir("%s%stmp" % (options.tmpDir, os.sep))

def getFileNameBase(outFlags,outFileByFlag,default):
    # what do we name all the log files?
    if outFlags:
        # base output file on first output file selected
        firstFlag=outFlags[0]
        if firstFlag in outFileByFlag:
            if isinstance(outFileByFlag[firstFlag],str):
                return outFileByFlag[firstFlag]
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
    with open(infile,'rU') as f:
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
    with open(inputFragment,'rU') as f:
        for line in f:
            m=fastaRecordRE.match(line)
            if m:
                records[m.group(1)]=False

    # assume only one output file for blast
    outFragment=fileMap.keys()[0]
    outStream=getStream(fileMap[outFragment])

    # scan file and write good records to output
    currentRecord=None
    currentRecordOutput=[]
    state=START
    with open(outFragment,'rU') as f:
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

    with open(filename, 'rU') as f:
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
    if options.local:
        localDir = tempfile.mkdtemp(suffix='batchRunner', dir=options.tmpDir)
    else:
        localDir = tempfile.mkdtemp(suffix='batchRunner')
    logging.debug("Local dir (%s): %s" % (os.path.exists(localDir), localDir))

    #if not options.local:
    if taskNum is None:
        # get the task number from env
        taskNum=getSGETaskId()

    hostname=getHostName()
    logFile="%s%stask%05d.%s.log" % (localDir,os.sep,taskNum,hostname)
    errStream=open(logFile,'w')
    # if called from SGE, we will have to:
    if not options.local:
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
    infragmentName = getFragmentName(options.fragBase, taskNum)
    prefix = getFragmentPrefix(options.fragBase, taskNum)
    infragment = "%s%s%s" % (options.tmpDir, os.sep, infragmentName)
    stdoutFragment="%s.stdout"%infragment
    stderrFragment="%s.stderr"%infragment
    if options.local:
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
        else:
            makeCommonFilesLocal(cmdargs,localDir)

        # modify command
        outLocal = "%s%s%s.output" % (localDir, os.sep, infragmentName)
        (foundI, outputFlags) = prepareCommandForFragment(options, infragment, prefix, outLocal, cmdargs, hostname, errStream)
        logging.debug("ready ro RUN!")

        # setup to run command
        # I/O
        if foundI:
            spin=None
        else:
            spin=open(infragment,'rU')
        spout=open(stdoutLocal,'w')
        sperr=open(stderrLocal,'w')
        # environment
        if options.bathybius:
            path=os.pathsep.join(['/common/bin',BINDIR,os.environ['PATH']])
        else:
            path=os.pathsep.join([BINDIR,os.environ['PATH']])

        # run it
        try:
            logging.debug("Command:\n%s\nwd=%s\npath=%s" % (formatCommand(cmdargs),subprocwd,path))
            t1=time.time()
            # Run the command
            exitcode=subprocess.call(cmdargs, stdin=spin, stdout=spout, stderr=sperr, cwd=subprocwd, env={"PATH":path})
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

        if not options.local:
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
                    shutil.move(localOutput,tmpOutput)
                else:
                    logging.warn("Could not find output: %s" % output)
                    logging.debug("Flag: %s" % (flag))

    except:
        exitcode=2
        errStream.write("Exception running fragment %d:\n" % (taskNum))
        errStream.write('-'*60+'\n')
        traceback.print_exc(file=errStream)
        errStream.write('-'*60+'\n')

    # Do some final cleanup:
    if exitcode!=0:
        errCodeFile="%s.exitcode" % (infragment)
        ecStream=open(errCodeFile,'w')
        ecStream.write(str(exitcode))
        ecStream.close()

    t2=time.time()
    logging.info("Fragment took %.2f seconds" % (t2-t0))
    # copy err to shared dir
    if not options.local:
        logging.shutdown()
    errStream.close()
    if options.local:
        shutil.move(logFile, "%s.log" % (infragment))
    else:
        shutil.copy(logFile, "%s.log" % (infragment))

    # remove local temorary dir
    for f in os.listdir(localDir):
        logging.debug("File: "+f)
    shutil.rmtree(localDir)

    return exitcode

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

    return (positionalArgs,flaggedArgs)

def getOutputFiles(options,cmdargs):
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
            namedOutCount+=1
            outFileByFlag[namedOutCount]=cmdargs[i]
            logging.debug("Found outfile number %d in arg %d" % (namedOutCount,i))
        elif nextArgument=='threads':
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
    return:
        boolean indicating if input file found or not
        list of outputFileFlags. Each is one of:
            a file name (with no path)
            a modification of the input file (starts with %, e.g. %e.ext)
            an integer indicating order in cmdargs
    """

    foundI=False
    logging.debug("Command:\n%s\ninfile: %s\noutfile: %s" % (cmdargs, infragment, outLocal))

    # setup hashes to look for options
    logging.debug("inflag: %r\noutflags: %r\nprefixFlags: %r" % (options.inputFlag, options.outputFlags, options.prefixFlag))
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
            namedOutCount+=1
            cmdargs[i]="%s.%s" % (outLocal,namedOutCount)
            #logging.debug("Replaced outfile at arg %d" % (i))
        elif nextArgument=='threads':
            if not options.local:
                threads = getThreadCountForNode(hostname,errStream)
                logging.debug("%s can take %d threads" % (hostname,threads))
                cmdargs[i]=str(threads)
        else:
            # something went wrong
            logging.warn("Unrecognized nextArgument: %s" % (nextArgument))
            raise Exception("Unrecognized nextArgument: %s" % (nextArgument))

        # reset
        nextArgument=None

    # were there any outputFlags not in the command? These indicate secondary files
    #logging.debug("Done with positional arguments in prepareCommandForFragment")
    otherOutputs=range(1,namedOutCount+1)
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

def getThreadCountForNode(hostname,errStream):
    """
    figure out how many threads the job can use
    """
    qhcmd = ["%s/qhost"%SGEBIN,"-h", hostname]
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
        logging.warn("Could not parse qhost output:\n%s" % (qhout))
    return slots

def getSGETaskId():
    return int(os.environ.get('SGE_TASK_ID',1))

def getHostName():
    try:
        return os.uname()[1].split('.')[0]
    except:
        try:
            return os.uname()[1]
        except:
            logging.warn("Could not parse hostname from %s" % (str(os.uname())))
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
    for kvPair in defaultValues.iteritems():
        command.extend(kvPair)

    # assume all files are on networked disks
    #sgeOptions.extend(['-l', 'bigd=1'])

    # get total bases in reference db
    if refDB is None:
        raise Exception("You must supply a database to run fr-hit")

    if refDBSize is not None:
        logging.warn("You supplied ref DB size of %s. If you omit the -R option batch_launcher will calculate the db size for you." % (refDBSize))
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
                logging.warn("Are you sure you want to split on number of records? It usually is a good idea to split on number of bases (-s)")

def inspectLastCommand(command,taskType,sgeOptions,commandBin,batchOptions):
    # set rps to 1 if not set:
    # but not if rps is already set
    #for sgeOpt in sgeOptions:
    #    if (len(sgeOpt)>4 and sgeOpt[:4]=='rps=') or (len(sgeOpt)>8 and sgeOpt[:8]=='rpsjobs='):
    #        break
    #else:
    #    sgeOptions.extend(['-l','rps=1'])

    # apply the defaults
    applyDefaultsToCommand(command,taskType,prepend=True)

## Blast
def inspectDCMBCommand(command,taskType,sgeOptions,blastBin,batchOptions):
    # replace first element of command string ('dcmegabalst') with 'blastn -task ...'
    command[0:1]=['blastn','-task','dc-megablast']
    inspectBlastCommand(command,taskType,sgeOptions,blastBin,batchOptions)

def inspectBlastCommand(command,taskType,sgeOptions,blastBin,batchOptions):
    """
    find any blast databseses in the command and check if they are in common data:
        if any are not, add sge option to throttle this task
    also, if using blastall, make sure prog/alg name is ok (one of: blastn, blastx, etc)
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
    for kvPair in defaultValues.iteritems():
        command.extend(kvPair)

    if dbs:
		#if not areDBsLocal(dbs):
        #    # if any DB will not be found in /locadata, slow down this job
        #    sgeOptions.extend(['-l', 'bigd=1'])
        pass
    else:
        raise Exception("You must supply a database to run blastall")

    if blastBin=='blastall':
        if program in ('rpsblast','rpsblastx','rpstblastn'):
            raise Exception("Use the rpsblast binary, not blastall!")
        if program not in ('blastx','blastn','blastp','tblastn','tblastx'):
            raise Exception("Blast algorithm %s is not recognized" % (program))
    elif blastBin in ['rpstblastn','rpsblast']:
        # override greedy
        #logging.debug("overriding greedy for rpsblast")
        #sgeOptions.extend(['-l','greedy=3'])
		pass

def makeRelativePathsAbsolute(cmdargs):
    """
    Translate any relative paths (starting with ./ or ../) to absolute
    """
    for i in xrange(len(cmdargs)):
        if relativePathRE.match(cmdargs[i]):
            cmdargs[i]=os.path.abspath(cmdargs[i])

def makeCommonFilesLocal(cmdargs,localdir):
    """
    Translate anything starting with /common/data to /localdata if it
    exists in both places
    """
    for i in xrange(len(cmdargs)):
        m = commonDataRE.match(cmdargs[i])
        if m:
            if os.path.exists(cmdargs[i]):
                newPath=commonDataRE.sub('/localdata/',cmdargs[i])
                if os.path.exists(newPath):
                    cmdargs[i]=newPath


def makeDBsLocal(cmdargs, localDir):
    """
    Finds database arguments in a blast(blastall, rpsblast, blast+) command and:
        * changes /common/data to /localdata if a localdata version is found
        * creates an alias file to search against if more than one DB given
            NB: this could fail if there is a nuc and pro db with same base name
            alias file is created in localtmp dir and will get deleted with it
    Also enforce defauls for b,v, and a
    """
    logging.debug("makeDBsLocal.")
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
            # get local version (db will be unchanged if not in /common/data)
            localdb=commonDataRE.sub('/localdata/',db)
            dbtype=getDBType(localdb)
            if dbtype is not None:
                # use local version if it exsists
                logging.debug("changed %s to %s" % (db,localdb))
                dbs.append(localdb)
                foundExtensions.append(dbtype)
            else:
                # use common version if not
                logging.debug("Using unchanged: %s" % (db))
                dbs.append(db)
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
        for kvPair in defaultValues.iteritems():
            command.insert(index,kvPair[1])
            command.insert(index,kvPair[0])
            index+=2
    else:
        for kvPair in defaultValues.iteritems():
            command.extend(kvPair)

##################
# Constants
##################
BATCHLAUNCHER=os.sep.join([BINDIR,'batch_launcher.py'])
scriptingLanguages=['perl','python','ruby','sh']
# task types
BLAST='blast'
DCMB='dcmegablast'
BLASTPLUS='blast+'
METARNA='meta_rna'
FGENESB='fgenesb'
GLIMMERMG='glimmermg'
FRHIT='frhit'
LAST='lastal'
programOptionMap={BLAST:{'in':'-i',
                         'out':['-o'],
                         'threads':'-a',
                        },
                  DCMB:{'in':'-query',
                        'out':['-out'],
                        'threads':'-num_threads'
                       },
                  BLASTPLUS:{'in':'-query',
                         'out':['-out'],
                         'threads':'-num_threads'
                        },
                  METARNA:{'in':'-i',
                           'out':['%p.coord','%p.mask','%p.seq','%e.coord','%e.mask','%e.seq'],
                           'prefix':'-o',
                           'threads':'-p',
                          },
                  LAST:{'in':2,
                        'out':['-o']},
                  GLIMMERMG:{'in':1,
                             'out':['%p.predict','%E.predict'],
                             'prefix':'-o',
                             'threads':'-p'},
                  FGENESB:{'in':'-i',
                           'out':['-o'],
                          },
                  FRHIT:{'in':'-d',
                         'out':['-o'],
                         'threads':'-T'}
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
taskTypePatterns={BLAST:blastProgRE,
                  BLASTPLUS:blastPlusProgRE,
                  DCMB:dcmbProgRE,
                  METARNA:mrnaProgRE,
                  FGENESB:fgbProgRE,
                  GLIMMERMG:glimmermgRE,
                  FRHIT:frHitRE,
                  LAST:lastRE,
                 }
inspectCommandForTaskType={BLAST:inspectBlastCommand,
                           BLASTPLUS:inspectBlastCommand,
                           DCMB:inspectDCMBCommand,
                           FRHIT:inspectFrHitCommand,
                           LAST:inspectLastCommand,
                          }
finalCommandProcessingForTask={BLAST:makeDBsLocal,
                               BLASTPLUS:makeDBsLocal,
                               DCMB:makeDBsLocal,
                              }
defaultsForTask={BLAST:{'-b':'10','-v':'10','-a':'8'},
                 DCMB:{'-max_target_seqs':'10','-num_threads':'8','-word_size':'12','-template_length':'21','-template_type':'optimal'},
                 BLASTPLUS:{'-max_target_seqs':'10','-num_threads':'8'},
                 METARNA:{'-H':'/common/bin/hmmer-3.0',
                          '-L':'/common/bin/meta_rna/rRNA_hmm_fs_wst_v0/HMM3'},
                 FRHIT:{'-r':'25','-T':'0'},
                 LAST:{'-b':'1','-f':'0',"-n":'10'},
                }
GREEDY_DEFAULT=24
greedyByTask={BLAST:GREEDY_DEFAULT,
              BLASTPLUS:GREEDY_DEFAULT,
              DCMB:GREEDY_DEFAULT,
              METARNA:GREEDY_DEFAULT,
              FGENESB:1,
              LAST:4,
              GLIMMERMG:GREEDY_DEFAULT,
              FRHIT:GREEDY_DEFAULT,
             }
unsupportedOptions={}
flagsWithArguments={GLIMMERMG:['--iter','-p','-t','-i','-q','-r','-s','-u','--fudge','--taxlevel','--minbp_pct'],
                    LAST:['-a','-b','-c','-d','-e','-F','-p','-q','-r','-x','-y','-z','-f','-k','-l','-m','-n','-s','-i','-u','-t','-j','-Q','-g','-G','-o'],
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
            except Queue.Empty:
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
            logging.warn("Task %d exited with %d" % (self.taskNum,
                exitcode))
            self.errorQueue.put(exitcode)

##################
# Expressions
################
commonDataRE=re.compile(r'/common/data/')
relativePathRE=re.compile(r'^\.\.?\/')

if __name__ == '__main__':
    main()

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

def fragmentInput(infile, options, tmpdir, fragmentBase, suffix='.in'):
    """
    Wraps the following methods into one:
        getFileType(options, infile)
        getSizePerChunk(infile, options.splits, fileType, splitOnSize=options.splitOnSize)
        fragmentInputBySize(infile, tmpdir, options.chunk, fileType, fragmentBase,   splitOnSize=options.splitOnSize, suffix=)
    """
    fileType=getFileType(options, infile)
    if options.chunk is None:
        if options.splits is None:
            sys.exit("Please tell me how many chunks (-N) or how big (-C)!")
        options.chunk=getSizePerChunk(infile, options.splits, fileType, splitOnSize=options.splitOnSize)
    return fragmentInputBySize(infile, tmpdir, options.chunk, fileType, fragmentBase, splitOnSize=options.splitOnSize, suffix=suffix)

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
        logging.warn("""Cannot guess type of %s!

The input file does not have a recognized extension:
%s

You must manually set the file type with '-T TYPE' where TYPE is in:
%s
OR set the record separator pattern with '-p PATTERN' where PATTERN is a regular expression. E.G. '^>" for fasta, or '^LOCUS' for gbk.
""" % (filename, fileExtensionMap.keys(), fileExtensionMap.keys()))
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
    for i in xrange(1,len(recordLines)):
        size+=len(recordLines[i].strip())
    return size

def fastqRecordSizer(recordLines):
    """
    Returns the number of charaters in the lines between the sequence header (@) and the quality header (@) excluding whitespace at the start and end of lines
    """
    size=0
    for i in xrange(1,len(recordLines)):
        line=recordLines[i]
        if len(line)>0 and line[0]=='+':
            return size
        size+=len(line.strip())

    logging.warn("Did not find quality line in record: %s" % recordLines)
    return size

def gbRecordSizer(recordLines):
    """
    uses biopython to parse record
    """
    from Bio import SeqIO
    record = SeqIO.read(recordLines,'genbank')
    return len(record)

def fragmentInputBySize(infile, tmpdir, chunk, fileType, fragmentBase, splitOnSize=True, suffix='.in'):
    """
    Break up input into files of size chunk in tmpdir. Return number of fragments.
    """
    logging.debug("Fragmenting input: %r" % ({'infile':infile,'tmpDir':tmpdir,'chunk':chunk,'base':fragmentBase}))
    inhandle = openInputFile(infile)
    num = fragmentInputStreamBySize(inhandle, tmpdir, chunk, fileType, fragmentBase, splitOnSize=splitOnSize, suffix=suffix)
    if infile is not None:
        inhandle.close()
    return num

def fragmentInputStreamBySize(inhandle, tmpdir, chunk, fileType, fragmentBase, splitOnSize=True, suffix='.in'):
    if splitOnSize:
        # get a custom function that returns the size of this type of record
        recordSizer=fileType.sizer
    else:
        # just return 1 for each record
        recordSizer=lambda x: 1

    count=0
    num=1
    tmpFileName=getFragmentPath(tmpdir,fragmentBase,num,suffix)
    #logging.debug('Writing fragment (%d,%d,%d): %s' % (chunk,count,num,tmpFileName))
    tmpFile = open(tmpFileName, 'w')
    for record in fileType.recordStreamer(inhandle):
        recordSize=recordSizer(record)
        count+=recordSize

        if count>chunk:
            # close previous chunk and open new one
            tmpFile.close
            num+=1
            tmpFileName=getFragmentPath(tmpdir,fragmentBase,num,suffix)
            #logging.debug('Writing fragment (%d,%d,%d): %s' % (chunk,count,num,tmpFileName))
            tmpFile = open(tmpFileName, 'w')
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

def addFragmentingOptions(parser,defaults={"splits":400}):
    parser.add_option("-L", "--recordLines", metavar="NUMLINES", dest='numLines', default=None, type="int",
                       help="Number of lines per record")
    parser.add_option("-P", "--pattern", metavar="PATTERN", dest='pattern', default=None,
                       help="Regular expression to split records")
    parser.add_option("-T","--infileType", dest='infileType', default=None,
                      choices=fileTypeMap.keys(),
                      help='Type of input file. Otherwise, choose by extension. Known types are: %choices')
    parser.add_option("-C", "--chunkSize", type="int", dest='chunk',  metavar="FRAG_SIZE",
                      help="The number of records per fragment. Overrides NUM_FRAGS")
    default=defaults.get("splits",None)
    parser.add_option("-N", "--numChunks", dest='splits', type='int', metavar="NUM_FRAGS", default=default,
                      help="The number of fragments to create (defaults to %default)")
    parser.add_option("-s", "--splitOnSize",  default=False, action='store_true',
                      help="create chunks based on record size, not number of records. For known sequence types (fasta, fastq, gb), sequence length is used, otherwize the full size of the record text is counted")

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
# Constants
##################
FASTA=FragmentableFileType('fasta',sizer=fastaRecordSizer,sepRE=re.compile(r'^>(\S+)'))
FASTQ=FragmentableFileType('fastq',sizer=fastqRecordSizer,numLines=4)
GENBANK=FragmentableFileType('gb',sizer=gbRecordSizer,sepRE=re.compile(r'^LOCUS'))
TABLE=FragmentableFileType('table',sepRE=re.compile(r'^'))
fileTypeMap={FASTA.name:FASTA,
             FASTQ.name:FASTQ,
             GENBANK.name:GENBANK,
             TABLE.name:TABLE,
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
                 }

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
            return open(infile,'rU')
    else:
        return infile
