import ConfigParser
import os
import sys

# our modules
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
sm_counter_v1 = __import__("qiaseq-smcounter-v1.sm_counter")
sm_counter_v2 = __import__("qiaseq-smcounter-v2.sm_counter_v2")

#---------------------------------------------------------------------------------
def makeLowPIFile(readSet,smCounterThreshold):
   ''' Make a file that contains variants with PI below-threshold but above 12
   Function used only for smCounterv1
   :param str readSet: The sample name
   :param int smCounterThreshold: Minimum prediction index for a variant to be called.
   '''
   fileout = open(readSet + ".smCounter.GT12PI.txt","w")
   firstLine = True
   idxPI = None
   for line in open(readSet + ".smCounter.all.txt", "r"):
      vals = line.strip().split("\t")
      if firstLine:
         firstLine = False
         idxPI = vals.index("PI")
         fileout.write("read set\t")
         fileout.write(line)
         continue
      predictionIndex = vals[idxPI]
      predictionIndex = int(float(predictionIndex)) if len(predictionIndex) > 0 else 0
      if 12 <= predictionIndex < smCounterThreshold:
         fileout.write(readSet)
         fileout.write("\t")
         fileout.write(line)
   fileout.close()

def run(cfg, paramFile, vc):
   # get read set name
   (path,readSet) = os.path.split(cfg.readSet)

   # get standard smCounter parameters from the main run-params.txt file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)
   cfgSmCounter = {}
   for (paramName, paramVal) in parser.items("smCounter"): 
      cfgSmCounter[paramName] = paramVal
      
   # set up config dictionary to pass to smCounter
   cfgSmCounter["outPrefix"] = readSet
   cfgSmCounter["bamFile"  ] = readSet + ".bam"
   cfgSmCounter["bedTarget"] = cfg.roiBedFile
   cfgSmCounter["rpb"      ] = cfg.readsPerUmi  # this comes from metrics.umi_frags module
   cfgSmCounter["nCPU"     ] = cfg.numCores
   cfgSmCounter["refGenome"] = cfg.genomeFile

   if vc == 'v1':
      cfgSmCounter["mtDepth"] = cfg.umiDepthMean # this comes from metrics.umi_depths module   
      cfgSmCounter["runPath"] = os.path.normpath(path)
      # run smCounter variant caller
      smCounterThreshold = sm_counter_v1.sm_counter.main(cfgSmCounter)      
   else:
      cfgSmCounter["runPath"] = os.path.normpath(path)
      sm_counter_v2.sm_counter_v2.main(cfgSmCounter)
      smCounterThreshold = 6

   # create low PI file for v1
   if vc == 'v1':
      makeLowPIFile(cfg.readSet,smCounterThreshold)
   
   # write smCounter threshold to disk file, for main summary table
   fileout = open(cfg.readSet + ".smCounter.summary.txt", "w")
   fileout.write("{}\tsmCounter variant calling threshold\n".format(smCounterThreshold))
   fileout.close()
   # return number of primitive variants called
   numVariants = -1
   for line in open(cfg.readSet + ".smCounter.cut.txt","r"):
      numVariants += 1
   return numVariants
   
