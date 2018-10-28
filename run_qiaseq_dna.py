from __future__ import print_function
import ConfigParser
import sys, os, argparse
import multiprocessing
# our modules
import core.run_log
import core.run_config
import core.prep
import core.prep_trim_duplex
import core.align
import core.umi_filter
import core.umi_mark
import core.umi_merge
import core.primer_clip
import core.samtools
import core.tumor_normal
import core.sm_counter_wrapper
import metrics.sum_specificity
import metrics.sum_uniformity_primer
import metrics.sum_primer_umis
import metrics.sum_all
import metrics.umi_frags
import metrics.umi_depths
import misc.process_ion
import misc.tvc
import annotate.vcf_complex
import annotate.vcf_annotate

#--------------------------------------------------------------------------------------
# call input molecules, build consenus reads, align to genome, trim primer region
#--------------------------------------------------------------------------------------
def run(args):
   readSet, paramFile, vc, outputPath = args

   # read run configuration file to memory
   cfg = core.run_config.run(readSet,paramFile,outputPath)
   fullReadSetPath=cfg.readSet

   # initialize logger
   core.run_log.init(fullReadSetPath,cfg)

   if cfg.platform.lower() == "illumina":

      if cfg.duplex.lower() == "true": ## Duplex sequencing run
         core.prep_trim_duplex.trimDuplex(cfg)
      else:
         # trim 3' ends of both reads, and extract UMI sequence
         core.prep.run(cfg)
   
      # align trimmed reads to genome using BWA MEM
      readFileIn1 = fullReadSetPath + ".prep.R1.fastq"
      readFileIn2 = fullReadSetPath + ".prep.R2.fastq"
      bamFileOut  = fullReadSetPath + ".align.bam"
      core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
   else:
      misc.process_ion.trimIon(cfg)
      misc.process_ion.alignToGenomeIon(cfg)

   # call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
   bamFileIn  = fullReadSetPath+ ".align.bam"
   core.umi_filter.run(cfg, bamFileIn)
   core.umi_mark.run(cfg)   
   metrics.umi_frags.run(cfg)   
   metrics.umi_depths.run(cfg)   
   core.umi_merge.run(cfg, bamFileIn)
   
   # soft clip primer regions from read alignments
   bamFileIn  = fullReadSetPath + ".umi_merge.bam"
   bamFileOut = fullReadSetPath + ".primer_clip.bam"
   core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)

   # additional metrics to generate
   metrics.sum_primer_umis.run(cfg) # primer-level umi and read metrics
   metrics.sum_specificity.run(cfg) # priming specificity
   metrics.sum_uniformity_primer.run(cfg) # primer-level uniformity

   # sort the final BAM file, to prepare for downstream variant calling
   bamFileIn  = fullReadSetPath + ".primer_clip.bam"
   bamFileOut = fullReadSetPath + ".bam"
   core.samtools.sort(cfg,bamFileIn,bamFileOut)   
  
   if cfg.duplex.lower() == "false": # do not run smCounter for duplex reads
      if cfg.platform.lower() != "illumina": # ion reads
         misc.tvc.run(cfg)

      # run smCounter variant calling
      numVariants = core.sm_counter_wrapper.run(cfg, paramFile, vc)
      
      if cfg.platform.lower() != "illumina":
         numVariants = misc.tvc.smCounterFilter(cfg,vc)

      # create complex variants, and annotate using snpEff
      if numVariants > 0:
         # convert nearby primitive variants to complex variants
         bamFileIn  = fullReadSetPath + ".bam"
         vcfFileIn  = fullReadSetPath + ".smCounter.cut.vcf"
         vcfFileOut = fullReadSetPath + ".smCounter.cplx.vcf"
         annotate.vcf_complex.run(cfg, bamFileIn, vcfFileIn, vcfFileOut, vc)
         
         # annotate variants in the VCF file
         vcfFileIn  = fullReadSetPath + ".smCounter.cplx.vcf"
         vcfFileOut = fullReadSetPath + ".smCounter.anno.vcf"
         annotate.vcf_annotate.run(cfg, vcfFileIn, vcfFileOut,vc)        
   
   # aggregate all summary metrics
   metrics.sum_all.run(cfg)
   
   # close log file
   core.run_log.close()


def run_tumor_normal(readSet,paramFile,vc,outputPath):
   ''' Wrapper around run() for tumor-normal analysis
   '''
   # 2 read set names which are space delimited
   readSets = filter(None,readSet.split(" "))
   assert len(readSets) == 2, "Tumor-Normal Analysis requires exactly 2 read sets !"

   # read parameter file
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(paramFile)

   tumor = None
   normal = None
   for section in parser.sections():
      if section not in ['general','smCounter']:
         for (paramName, paramVal) in parser.items(section):
            if paramName == 'sampleType' and paramVal.lower() == 'normal':
               normal = section
            elif paramName == 'sampleType' and paramVal.lower() == 'tumor':
               tumor = section


   assert tumor!=None and normal!=None, "Could not sync read set names supplied with config file !"
  
   run((normal,paramFile,vc,outputPath))
   run((tumor,paramFile,vc,outputPath))
   ## Additional analysis steps
   cfg = core.run_config.run(tumor,paramFile,outputPath)
   core.tumor_normal.removeNormalVariants(cfg)
   core.tumor_normal.runCopyNumberEstimates(cfg)

#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   cfg=core.run_config.parse_command_line_arguments()#This function handles the command line parsing and validtion, and, in its case, the prinitng of the USAGE screen

   miscfileparts=os.path.split(misc.__file__)
   miscparentparts=os.path.split(miscfileparts[0])#Gets the current path of this file
   os.environ["PATH"]=os.environ["PATH"]+":"+miscparentparts[0]#Adds the pah of the current file to the environment
   #Anterior command line: "\nRun as : python run_qiaseq_dna.py <param_file> <v1/v2> <single/tumor-normal> <readSet(s)>\n"
   print("pid:",os.getpid())
   #referenceUmiFiles = getattr(cfg,'refUmiFiles',cfg.readSet+'.sum.primer.umis.txt').split(",")
   #assert 0,referenceUmiFiles

   if cfg.analysis.lower() == "tumor-normal":
      for (iread,read) in enumerate(cfg.readSet):
         if "{1}" in cfg.outputPath:
            outputpath=cfg.outputPath.format(read,iread)
         elif "{0}" in cfg.outputPath:
            outputpath=cfg.outputPath.format(read)
         else:
            outputpath=cfg.outputPath
         #print("paruntn",read,cfg.paramFile,cfg.vc,outputpath)
         run_tumor_normal(read,cfg.paramFile,cfg.vc,outputpath)
   else: # Single sample, might still need to run quandico
      for (iread,read) in enumerate(cfg.readSet):
         if "{1}" in cfg.outputPath:
            outputpath=cfg.outputPath.format(read,iread)
         elif "{0}" in cfg.outputPath:
            outputpath=cfg.outputPath.format(read)
         else:
            outputpath=cfg.outputPath
         #print("parun",read,cfg.paramFile,cfg.vc,outputpath)
         run((read,cfg.paramFile,cfg.vc,outputpath))
         runcfg = core.run_config.run(read,cfg.paramFile,cfg.outputPath)
         core.tumor_normal.runCopyNumberEstimates(runcfg)
