import os, sys, logging, time

import warnings
from astropy.wcs import FITSFixedWarning
# Suppress FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)

class _ColorStreamHandler(logging.StreamHandler):

    DEFAULT = '\x1b[0m'
    RED     = '\x1b[31m'
    GREEN   = '\x1b[32m'
    YELLOW  = '\x1b[33m'
    CYAN    = '\x1b[36m'

    CRITICAL = RED
    ERROR    = RED
    WARNING  = YELLOW
    INFO     = GREEN
    DEBUG    = CYAN

    @classmethod
    def _get_color(cls, level):
        if level >= logging.CRITICAL:  return cls.CRITICAL
        elif level >= logging.ERROR:   return cls.ERROR
        elif level >= logging.WARNING: return cls.WARNING
        elif level >= logging.INFO:    return cls.INFO
        elif level >= logging.DEBUG:   return cls.DEBUG
        else:                          return cls.DEFAULT

    def __init__(self, stream=None):
        logging.StreamHandler.__init__(self, stream)

    def format(self, record):
        color = self._get_color(record.levelno)
        record.msg = color + record.msg + self.DEFAULT
        return logging.StreamHandler.format(self, record)

class Logger():

    def __init__(self, pipename):# logfile = "pipeline.logging", log_dir = "logs"):

        # hopefully kill other loggers
        logger = logging.getLogger()
        logger.propagate = False
        logger.handlers = []
 
        timestamp = time.strftime('%Y-%m-%d_%H:%M:%S', time.localtime())
        self.logfile = pipename+'_'+timestamp+'.logger'
        self.log_dir = 'logs_'+pipename+'_'+timestamp
        os.makedirs(self.log_dir)
        #self.backup(logfile, log_dir)
        self.set_logger(self.logfile)
        # Here we create symlinks for easier debugging and checking of the newest logs
        # we create temp-symlinks and rename them so that we can overwrite existing symlinks.
        os.symlink( self.log_dir, f'logs_{pipename}.tmp')
        os.rename(f'logs_{pipename}.tmp',f'logs_{pipename}')
        os.symlink(self.logfile, pipename+'.logger.tmp')
        os.rename(pipename+'.logger.tmp', pipename+'.logger')


#    def backup(self, logfile, log_dir):
#
#        # bkp old log dir
#        if os.path.isdir(log_dir):
#            current_time = time.localtime()
#            log_dir_old = time.strftime(log_dir+'_bkp_%Y-%m-%d_%H:%M', current_time)
#            os.system('mv %s %s' % ( log_dir, log_dir_old ))
#        os.makedirs(log_dir)
#
#        # bkp old log file
#        if os.path.exists(logfile):
#            current_time = time.localtime()
#            logfile_old = time.strftime(logfile+'_bkp_%Y-%m-%d_%H:%M', current_time)
#            os.system('mv %s %s' % ( logfile, logfile_old ))
            

    def set_logger(self, logfile):
      
        logger = logging.getLogger("LiLF")
        logger.setLevel(logging.DEBUG)

        # create file handler which logs even debug messages
        handlerFile = logging.FileHandler(logfile)
        handlerFile.setLevel(logging.DEBUG)
        
        # create console handler with a higher log level
        handlerConsole = _ColorStreamHandler(stream=sys.stdout)
        handlerConsole.setLevel(logging.INFO)
        
        # create formatter and add it to the handlers
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
        handlerFile.setFormatter(formatter)
        handlerConsole.setFormatter(formatter)
        
        # add the handlers to the logger
        logger.addHandler(handlerFile)
        logger.addHandler(handlerConsole)

        logger.info('Logging initialised in %s (file: %s)' % (os.getcwd(), logfile))


# this is used by all libraries for logging
logger = logging.getLogger("LiLF")
