
import sys
import datetime
import logging
import os

# global hack
sysStdOut = sys.stdout
sysStdErr = sys.stderr
timeStart = None

#----------------------------------------------------------------------
# stream handler object that redirects to Python logging facility
#----------------------------------------------------------------------
class RedirectToLogger(object):
   def __init__(self,is_stderr,logfilename):
      self.is_stderr=is_stderr
      #if is_stderr:
      #  self.log_level=logging.ERROR
      #else:
      self.log_level=logging.DEBUG
      self.logger = logging.getLogger(logfilename)
      #self.logger = logging.getLogger()
      handler=logging.FileHandler(logfilename,mode="w")
      self.logger.addHandler(handler)
 
   def write(self, buf):
      for line in buf.rstrip().splitlines():
        self.logger.log(self.log_level,line.rstrip())
         
   def flush(self):
      pass
 
#----------------------------------------------------------------------
# initialize the logging
#----------------------------------------------------------------------
def init(logFilePrefix,cfg):

    # format a log file name
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y.%m.%d_%H.%M.%S")
    logFileName = os.path.join(cfg.outputPath, logFilePrefix + ".run_log." + timestamp + ".txt")
    # time the entire pipeline
    global timeStart
    timeStart = datetime.datetime.now()
    print("We will write the log file in the path",logFileName)
    print("started_log0 at " + str(timeStart))

    # set up logging to a log file
    logDateFormat = "%Y-%m-%d %H:%M:%S"
    logFormat = "%(asctime)s.%(msecs)03d %(message)s"
    logging.basicConfig(level=logging.DEBUG,format=logFormat,datefmt=logDateFormat)
                     #filename=logFileName,filemode="w")
  
    # create a stream redirection object for stdout, so print() goes to logger
    sys.stdout = RedirectToLogger(False,logFileName)
    sys.stderr = RedirectToLogger(True,logFileName)
    print("sys.argv:",sys.argv)
    print("pid:",os.getpid())
    print("started_log at " + str(timeStart))
    cfg.print_data(logFilePrefix)





#----------------------------------------------------------------------
# close log file and restore stdout and stderr
#----------------------------------------------------------------------
def close():
    # log run completion
    timeEnd = datetime.datetime.now()
    print("total run time: " + str(timeEnd - timeStart))
    print("Execution ended at",timeEnd)
    hdlrtonow=sys.stdout
    logger=hdlrtonow.logger
    # restore stdout and stderr
    sys.stdout = sysStdOut
    sys.stderr = sysStdErr
   
    # close log file, clean up
    #logger = logging.getLogger()
    for (ihandler,handler) in enumerate(logger.handlers):
        handler = logger.handlers[0]
        handler.stream.close()
        logger.removeHandler(handler)

    print("Logger closed at",timeEnd)
