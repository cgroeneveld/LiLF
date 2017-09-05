import os, sys, re, pickle, random, shutil, logging, glob

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
#import lsmtool

#logger = logging.getLogger("PiLL")


def printLineBold(line):
    """
    Print a line of text 'line' using boldface.
    """

    boldStart = "\033[1m"
    boldEnd   = "\033[0;0m"
    print (boldStart + line + boldEnd)


def distanceOnSphere(RAs1, Decs1, RAs2, Decs2):
    """
    Returns the distances on the sphere from the set of points '(RAs1, Decs1)' to the
    set of points '(RAs2, Decs2)' using the spherical law of cosines.

    It assumes that all inputs are given in degrees, and gives the output in degrees, too.

    Using 'numpy.clip(..., -1, 1)' is necessary to counteract the effect of numerical errors, that can sometimes
    incorrectly cause '...' to be slightly larger than 1 or slightly smaller than -1. This leads to NaNs in the arccosine.
    """

    return np.degrees(np.arccos(np.clip(
           np.sin(np.radians(Decs1)) * np.sin(np.radians(Decs2)) +
           np.cos(np.radians(Decs1)) * np.cos(np.radians(Decs2)) *
           np.cos(np.radians(RAs1 - RAs2)), -1, 1)))


def check_rm(regexp):
    """
    Check if file exists and remove it
    Handle reg exp of glob and spaces
    """
    filenames = regexp.split(' ')
    for filename in filenames:
        # glob is used to check if file exists
        for f in glob.glob(filename):
            os.system("rm -r " + f)


def run_losoto(s, c, h5s, parsets):
    """
    s : scheduler
    c : cycle name, e.g. "final"
    h5s : lists of H5parm files
    parsets : lists of parsets to execute
    """

    logging.info("Running LoSoTo...")

    # concat
    check_rm("cal-" + c + ".h5")
    s.add('H5parm_append.py -v -c freq -s sol000 -o cal-'+c+'.h5 '+' '.join(h5s), log='losoto-'+c+'.log', cmd_type='python', processors='max')
    s.run(check = True)

    check_rm('plots')
    os.makedirs('plots')
    
    for parset in parsets:
        logging.debug('-- executing '+parset+'...')
        s.add('losoto -v cal-'+c+'.h5 '+parset, log='losoto-'+c+'.log', log_append=True, cmd_type='python', processors='max')
        s.run(check=True)

    check_rm('plots-' + c)
    os.system('mv plots plots-' + c)
    

class Scheduler():
    def __init__(self, qsub = None, max_threads = None, max_processors = None, log_dir = 'logs', dry = False):
        """
        qsub:           if true call a shell script which call qsub and then wait 
                        for the process to finish before returning
        max_threads:    max number of parallel processes
        dry:            don't schedule job
        max_processors: max number of processors in a node (ignored if qsub=False)
        """
        self.cluster = self.get_cluster()
        self.qsub    = qsub
        # if qsub/max_thread/max_processors not set, guess from the cluster
        # if they are set, double check number are reasonable
        if self.qsub == None:
            if self.cluster == 'Hamburg': self.qsub = True
            else: self.qsub = False
        else:
            if (self.qsub == False and self.cluster == 'Hamburg') or \
               (self.qsub == True and (self.cluster == 'Leiden' or self.cluster == 'CEP3')):
                logging.critical('Qsub set to %s and cluster is %s.' % (str(qsub), self.cluster))
                sys.exit(1)

        if max_threads == None:
            if self.cluster == 'Hamburg': self.max_threads = 32
            elif self.cluster == 'Leiden': self.max_threads = 64
            elif self.cluster == 'CEP3': self.max_threads = 40
            else: self.max_threads = 12
        else:
            self.max_threads = max_threads

        if max_processors == None:
            if self.cluster == 'Hamburg': self.max_processors = 6
            elif self.cluster == 'Leiden': self.max_processors = 64
            elif self.cluster == 'CEP3': self.max_processors = 40
            else: self.max_processors = 12
        else:
            self.max_processors = max_processors

        self.dry = dry
        logging.info("Scheduler initialized for cluster " + self.cluster + " (Nproc: " + str(self.max_threads) + ", multinode: " +
                     str(self.qsub) + ", max_processors: " + str(self.max_processors) + ").")

        self.action_list = []
        self.log_list    = [] # list of 2-length tuple of the type: (log filename, type of action)

        if (not os.path.isdir(log_dir)):
            logging.info("Creating log dir '" + log_dir + "'.")
            os.makedirs(log_dir)
        self.log_dir = log_dir


    def get_cluster(self):
        """
        Find in which computing cluster the pipeline is running       
        """
        import socket
        hostname = socket.gethostname()
        if (hostname == 'lgc1' or hostname == 'lgc2'):
            return 'Hamburg'
        elif ('leidenuniv' in hostname):
            return 'Leiden'
        elif (hostname[0 : 3] == 'lof'):
            return 'CEP3'
        else: 
            logging.error('Hostname %s unknown.' % hostname)
            return 'Unknown'


    def add(self, cmd = '', log = '', log_append = False, cmd_type = '', processors = None):
        """
        Add a command to the scheduler list
        cmd:        the command to run
        log:        log file name that can be checked at the end
        log_append: if True append, otherwise replace
        cmd_type:   can be a list of known command types as "BBS", "DPPP", ...
        processors: number of processors to use, can be "max" to automatically use max number of processors per node
        """
        if (log != ''):
            log = self.log_dir + '/' + log
        if (log != '' and not log_append):
            cmd += " > " + log + " 2>&1"
        if (log != '' and log_append):
            cmd += " >> " + log + " 2>&1"

        if processors != None and processors == 'max': processors = self.max_processors

        if self.qsub:
            # if number of processors not specified, try to find automatically
            if processors == None:
                processors = 1 # default use single CPU
                if "calibrate-stand-alone" == cmd[:21]: processors = 1
                if "NDPPP" == cmd[:5]: processors = 1
                if "wsclean" == cmd[:7]: processors = self.max_processors
                if "awimager" == cmd[:8]: processors = self.max_processors
            if processors > self.max_processors: processors = self.max_processors
            self.action_list.append([str(processors),'\''+cmd+'\''])
        else:
            self.action_list.append(cmd)

        if log != '':
            self.log_list.append((log,cmd_type))


    def add_casa(self, cmd='', params={}, wkd=None, log='', log_append=False, processors=None):
        """
        Run a casa command pickling the parameters passed in params
        NOTE: running casa commands in parallel is a problem for the log file, better avoid
        alternatively all used MS and CASA must be in a separate working dir

        wkd = working dir (logs and pickle are in the pipeline dir)
        """

        if processors != None and processors == 'max': processors = self.max_processors
        if processors == None: processors=self.max_processors # default use entire node

        # since CASA can run in another dir, be sure log and pickle are in the pipeline working dir
        if log != '': log = os.getcwd()+'/'+self.log_dir+'/'+log
        pfile = os.getcwd()+'/casaparams_'+str(random.randint(0, 1e9))+'.pickle'
        pickle.dump( params, open( pfile, "wb" ) )

        # exec in the script dir?
        if wkd == None: casacmd = 'casa --nogui --log2term --nologger -c '+cmd+' '+pfile
        elif os.path.isdir(wkd):
            casacmd = 'cd '+wkd+'; casa --nogui --log2term --nologger -c '+cmd+' '+pfile
        else:
            logging.error('Cannot find CASA working dir: '+wkd)
            sys.exit(1)

        if self.qsub:
            if log != '' and not log_append: casacmd = str(processors)+' \''+casacmd+' > '+log+' 2>&1'
            elif log != '' and log_append: casacmd = str(processors)+' \''+casacmd+' >> '+log+' 2>&1'
            else: casacmd = str(processors)+' \''+casacmd

            # clean up casa remnants in Hamburg cluster
            if self.cluster == 'Hamburg':
                self.action_list.append(casacmd+'; killall -9 -r dbus-daemon Xvfb python casa\*\'')
                if processors != self.max_processors:
                    logging.error('To clean annoying CASA remnants no more than 1 CASA per node is allowed.')
                    sys.exit(1)
            else:
                self.action_list.append(casacmd+'\'')
        else:
            if log != '' and not log_append: self.action_list.append(casacmd+' > '+log+' 2>&1')
            elif log != '' and log_append: self.action_list.append(casacmd+' >> '+log+' 2>&1')
            else: self.action_list.append(casacmd)

        if log != '':
            self.log_list.append((log,'CASA'))


    def run(self, check=False, max_threads=None):
        """
        If check=True then a check is done on every log in the log_list
        if max_thread != None, then it overrides the global values, useful for special commands that need a lower number of threads
        """
        from threading import Thread
        from Queue import Queue
        import subprocess

        def worker(queue):
            for cmd in iter(queue.get, None):
                if self.qsub and self.cluster == 'Hamburg':
                    # run in priority nodes
                    #cmd = 'salloc --job-name LBApipe --reservation=important_science --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                    #        ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                    # run on all cluster
                    cmd = 'salloc --job-name LBApipe --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                            ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                subprocess.call(cmd, shell = True)
    
        # limit threads only when qsub doesn't do it
        if max_threads != None: max_threads_run = min(max_threads, self.max_threads)
        else: max_threads_run = self.max_threads

        q       = Queue()
        threads = [Thread(target=worker, args=(q,)) for _ in range(max_threads_run)]
    
        for i, t in enumerate(threads): # start workers
            t.daemon = True
            t.start()
    
        for action in self.action_list:
            if self.dry: continue # don't schedule if dry run
            q.put_nowait(action)
        for _ in threads: q.put(None) # signal no more commands
        for t in threads: t.join()

        # check outcomes on logs
        if check:
            for log, cmd_type in self.log_list:
                self.check_run(log, cmd_type)

        # reset list of commands
        self.action_list = []
        self.log_list    = []


    def check_run(self, log = '', cmd_type = ''):
        """
        Produce a warning if a command didn't close the log properly i.e. it crashed
        NOTE: grep, -L inverse match, -l return only filename
        """
        import subprocess

        if not os.path.exists(log):
            logging.warning('No log file found to check results: ' + log)
            return 1

        if cmd_type == 'BBS':
            out = subprocess.check_output('grep -L success '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('BBS run problem on:\n'+out.split("\n")[0])
                return 1

        elif cmd_type == 'NDPPP':
            out = subprocess.check_output('grep -L "Finishing processing" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "Exception" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "**** uncaught exception ****" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('NDPPP run problem on:\n'+out.split("\n")[0])
                return 1

        elif cmd_type == 'CASA':
            out = subprocess.check_output('grep -l "[a-z]Error" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "An error occurred running" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "\*\*\* Error \*\*\*" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('CASA run problem on:\n'+out.split("\n")[0])
                return 1

        elif cmd_type == 'wsclean':
            out = subprocess.check_output('grep -l "exception occurred" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -L "Cleaning up temporary files..." '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('WSClean run problem on:\n'+out.split("\n")[0])
                return 1

        elif cmd_type == 'python':
            out = subprocess.check_output('grep -l "Traceback (most recent call last):" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -i -l "Error" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -i -l "Critical" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('Python run problem on:\n'+out.split("\n")[0])
                return 1

        elif cmd_type == 'general':
            out = subprocess.check_output('grep -l -i "error" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('Run problem on:\n'+out.split("\n")[0])
                return 1

        else:
            logging.warning('Unknown command type for log checking: "'+cmd_type+'"')
            return 1

        return 0