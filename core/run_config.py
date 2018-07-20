from __future__ import print_function
import ConfigParser, os
from multiprocessing.dummy import cpu_count as cpu_count

class ConfigObj(object):
   """Class intended for its instances to include the configuration values as atrributes-values, and based on the values to do some things like calculating output paths.
   This was created after the object layout that was present in the run function:
   cfg = lambda:0
   cfg.__dict__["readSet"] = readSet
   for section in ("general", readSet):
      for (paramName, paramVal) in parser.items(section):
         if paramName in cfg.__dict__:
            raise Exception("Config file contains duplicate specification of parameter: " + paramName)
         cfg.__dict__[paramName] = paramVal
         print(paramName, paramVal)
   So that, for example when the object is called it returns 0. I don't know if it is of utility, or was an easy way to create an object with free attributes
   """
   def __call__(self):
      return 0
   def file_path(self,filename):
      outputPath=getattr(self,"outputPath",None)
      if outputPath is None:
         #There is no outputpath specified
         return filename
      else:
         if not os.path.exists(outputPath):
            os.makedirs(outputPath)
         return os.path.join(outputPath,filename)

#--------------------------------------------------------------------------------------
def run(readSet,paramFile):

   # read parameter file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)

   # copy all options to a config object - both general params, and params for this readSet
   #cfg = lambda:0
   cfg=ConfigObj()
   for section in ("general", readSet):
      for (paramName, paramVal) in parser.items(section):
         if paramName in cfg.__dict__:
            raise Exception("Config file contains duplicate specification of parameter: " + paramName)
         cfg.__dict__[paramName] = paramVal
         print(paramName, paramVal)

   cfg.__dict__["readSet"] = cfg.file_path(readSet)
   if str(cfg.numCores) == '0':      
      # use all cores if numCores = 0
      cfg.numCores = str(cpu_count())
   
   # convert some params to boolean
   cfg.deleteLocalFiles = cfg.deleteLocalFiles.lower() == "true"

 
   # return config object
   return cfg
