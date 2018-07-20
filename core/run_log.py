from __future__ import print_function
import sys
import datetime
import logging

# global hack
sysStdOut = sys.stdout
sysStdErr = sys.stderr
timeStart = None

#----------------------------------------------------------------------
# stream handler object that redirects to Python logging facility
#----------------------------------------------------------------------
class RedirectToLogger(object):
   def __init__(self,is_stderr):
      self.is_stderr=is_stderr
      if is_stderr:
        self.log_level=logging.ERROR
      else:
        self.log_level=logging.DEBUG
      self.logger = logging.getLogger()
 
   def write(self, buf):
      for line in buf.rstrip().splitlines():
        self.logger.log(self.log_level,line.rstrip())
         
   def flush(self):
      pass
 
#----------------------------------------------------------------------
# initialize the logging
#----------------------------------------------------------------------
def init(logFilePrefix):

   # format a log file name
   now = datetime.datetime.now()
   timestamp = now.strftime("%Y.%m.%d_%H.%M.%S")
   logFileName = logFilePrefix + ".run_log." + timestamp + ".txt"
   print("We will write the log file in the path",logFileName)

   # set up logging to a log file
   logDateFormat = "%Y-%m-%d %H:%M:%S"
   logFormat = "%(asctime)s.%(msecs)03d %(message)s"
   logging.basicConfig(level=logging.DEBUG,
                       format=logFormat,
                       datefmt=logDateFormat,
                       filename=logFileName,
                       filemode="w")

   # create a stream redirection object for stdout, so print() goes to logger
   sys.stdout = RedirectToLogger(False)
   sys.stderr = RedirectToLogger(True)

   # time the entire pipeline
   global timeStart
   timeStart = datetime.datetime.now()
   print("started at " + str(timeStart))

#----------------------------------------------------------------------
# close log file and restore stdout and stderr
#----------------------------------------------------------------------
def close():
   # log run completion
   timeEnd = datetime.datetime.now()
   print("total run time: " + str(timeEnd-timeStart))
   
   # restore stdout and stderr
   sys.stdout = sysStdOut
   sys.stdout = sysStdErr
   
   # close log file, clean up
   logger = logging.getLogger()
   handler = logger.handlers[0]
   handler.stream.close()
   logger.removeHandler(handler)
   
