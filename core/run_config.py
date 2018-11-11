from __future__ import print_function
import ConfigParser, os, argparse, sys, traceback
from multiprocessing.dummy import cpu_count as cpu_count

class ArgumentFileParserAction(argparse.Action):
     def __init__(self, option_strings, dest, nargs=None, **kwargs):
         #print("calledact",option_strings, dest, nargs, kwargs)
         if nargs is not None:
             raise ValueError("nargs not allowed")
         super(ArgumentFileParserAction, self).__init__(option_strings, dest, **kwargs)
     def __call__(self, parser, namespace, values, option_string=None):
         print ('action called: %r %r %r' % (namespace, values, option_string))
         #setattr(namespace, self.dest, values)
         # read parameter file
         parser = ConfigParser.SafeConfigParser()
         parser.optionxform = str
         parser.read(paramFile)

         # copy all options to a config object - both general params, and params for this readSet
         #cfgobj = lambda:0
         #cfgobj=ConfigObj()
         #(readsetpath,readsetbasename)=os.path.split(readSet)
         #for section in ("general", readSet):
         for section in ("general", readsetbasename):
            for (paramName, paramVal) in parser.items(section):
               if paramName in namespace.__dict__:
                  raise Exception("Config file contains duplicate specification of parameter {0}, first time defining it to {1} and second one defining to {2}".format(
                     paramName,getattr(namespace,paramName),paramVal))
               setattr(namespace,paramName,paramVal)
               print("Parameter_actioncalled",paramName, "==",paramVal)

         setattr(namespace,"readSet",cfgobj.file_path(readSet))
         if str(cfgobj.numCores) == '0':      
            # use all cores if numCores = 0
            cfgobj.numCores = str(cpu_count())
         
         # convert some params to boolean
         cfgobj.deleteLocalFiles = cfgobj.deleteLocalFiles.lower() == "true"
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
            try:
               assert "{0}" not in outputPath
            except AssertionError as asserr:
                  print("Occurred the exception: {0}, tb: {1}".format(asserr,traceback.format_exc()))

            os.makedirs(outputPath)
         return os.path.join(outputPath,filename)
   def print_data(self,readsetbasename):
      for (paramName,paramVal) in self.__dict__.items():
         print("Parameter {0}=={1}".format(paramName,paramVal))


def parse_command_line_arguments():
   """This function has the arguments definitions needed to process the command line arguments using the argparse standard library."""
   parser = argparse.ArgumentParser(prog="python "+sys.argv[0])
   parser.add_argument("paramFile",metavar="<param_file>",help="Parameter file. This file may contain several sections with names of read sets. Only the parameters the sections 'general', 'smconunter' and with the name of the read set will take effect.")
   parser.add_argument("vc",help="version of smcounter",choices=("v1","v2"))
   parser.add_argument("analysis",help="Type of analysis",choices=("single","tumor-normal"))
   parser.add_argument("outputPath",help="Path where the output files will be created. If it does not exist, program will create it.")
   parser.add_argument("readSet",help="Read set(s)",nargs="+")
   args = parser.parse_args()
   #args=run(" ".join(args.readSet),args.paramFile,args)
 
   return args
#--------------------------------------------------------------------------------------
def run(readSet,paramFile,outputPath):

   # read parameter file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)

   # copy all options to a config object - both general params, and params for this readSet
   #cfgobj = lambda:0
   cfgobj=ConfigObj()
   cfgobj.outputPath=outputPath
   (readsetpath,readsetbasename)=os.path.split(readSet)
   for section in ("general", readsetbasename):
      for (paramName, paramVal) in parser.items(section):
         if paramName in cfgobj.__dict__:
            raise Exception("Config file contains duplicate specification of parameter: " + paramName)
         cfgobj.__dict__[paramName] = paramVal
         print("Parameter",paramName, "==",paramVal)

   cfgobj.__dict__["readSet"] = cfgobj.file_path(readSet)
   if str(cfgobj.numCores) == '0':      
      # use all cores if numCores = 0
      cfgobj.numCores = str(cpu_count())
   
   # convert some params to boolean
   cfgobj.deleteLocalFiles = cfgobj.deleteLocalFiles.lower() == "true"

 
   # return config object
   return cfgobj
