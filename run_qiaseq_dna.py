#!/usr/bin/env python
from __future__ import print_function
import ConfigParser
import sys, os, argparse, traceback, socket, warnings
import multiprocessing
import shutil
# our modules
import core.run_log
import core.run_config
import core.prep
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
import metrics.duplex_summary
import metrics.sum_primer_duplex
import metrics.fraglen_by_rpu
import misc.process_ion
import misc.tvc
import annotate.vcf_complex
import annotate.vcf_annotate

#--------------------------------------------------------------------------------------
# trim primer and adapters, align to genome, mark UMIs, call variants using smCounter
#--------------------------------------------------------------------------------------
def run_old(args):
   readSet, paramFile, vc, outputPath = args
   
   # read run configuration file to memory
   cfg = core.run_config.run(readSet,paramFile,outputPath)
   fullReadSetPath=cfg.readSet

   # initialize logger
   core.run_log.init(fullReadSetPath,cfg)

   try:
      lastfilegenerated = cfg.readSet + ".vcf_complex.summary.txt"
      if os.path.exists(lastfilegenerated):
         print("Last file listed to be generated was found: {0}, with stats {1}".format(lastfilegenerated,os.stat(lastfilegenerated)))
         core.tumor_normal.runCopyNumberEstimates(cfg)
         # close log file
         core.run_log.close()
         return
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
      core.tumor_normal.runCopyNumberEstimates(cfg)

   except Exception as exc:
      traceback.print_exc()
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
def run(args,tumorNormal):
    readSet, paramFile, vc = args
    # initialize logger
    if not tumorNormal:
        core.run_log.init(readSet)
 
    # read run configuration file to memory
    cfg = core.run_config.run(readSet,paramFile)

    # trim adapters , umi and primers (this module spawns multiple processes)
    core.prep.run(cfg)

    readFileIn1 = readSet + ".prep.R1.fastq"
    readFileIn2 = readSet + ".prep.R2.fastq"
    bamFileOut  = readSet + ".align.bam"    
    if cfg.platform.lower() == "illumina":
        # align trimmed reads to genome using BWA MEM
        core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
    else: # use tmap for ion torrent reads        
        misc.process_ion.alignToGenomeIon(cfg, readFileIn1, bamFileOut)
  
    # call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side    
    bamFileIn  = readSet + ".align.bam"     
    core.umi_filter.run(cfg, bamFileIn)
    core.umi_mark.run(cfg)   
    metrics.umi_frags.run(cfg)
       
    metrics.umi_depths.run(cfg,vc)   
    core.umi_merge.run(cfg, bamFileIn)
    
    # soft clip primer regions from read alignments
    bamFileIn  = readSet + ".umi_merge.bam"
    bamFileOut = readSet + ".primer_clip.bam"
    core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)
 
    # additional metrics to generate
    metrics.sum_primer_umis.run(cfg) # primer-level umi and read metrics
    metrics.sum_specificity.run(cfg) # priming specificity
    metrics.sum_uniformity_primer.run(cfg) # primer-level uniformity

    if cfg.duplex: # additional metrics for Duplex reads
        metrics.duplex_summary.run(cfg)
        metrics.sum_primer_duplex.run(cfg)
        metrics.fraglen_by_rpu.run(cfg)

    # sort the final BAM file, to prepare for downstream variant calling
    bamFileIn  = readSet + ".primer_clip.bam"
    bamFileOut = readSet + ".bam"
    core.samtools.sort(cfg,bamFileIn,bamFileOut)
    
    if cfg.platform.lower() != "illumina": # ion reads
        misc.tvc.run(cfg)

    # run smCounter variant calling
    numVariants = core.sm_counter_wrapper.run(cfg, paramFile, vc)

    if cfg.platform.lower() != "illumina":
        numVariants = misc.tvc.smCounterFilter(cfg,vc)

    # create complex variants, and annotate using snpEff
    if not tumorNormal:
        post_smcounter_work(numVariants, readSet, cfg, tumorNormal=False)
        # close log file
        core.run_log.close()
        
def post_smcounter_work(numVariants, readSet, cfg, tumorNormal):
    ''' Additional Steps after smCounter
    :param int numVariants
    :param str readSet
    :param lambda obj cfg
    :param bool tumorNormal 
    '''
    if numVariants > 0:
        # convert nearby primitive variants to complex variants
        bamFileIn  = readSet + ".bam"
        vcfFileIn  = readSet + ".smCounter.cut.vcf"
        vcfFileOut = readSet + ".smCounter.cplx.vcf"
        annotate.vcf_complex.run(cfg, bamFileIn, vcfFileIn, vcfFileOut, vc)
            
        # annotate variants in the VCF file
        vcfFileIn  = readSet + ".smCounter.cplx.vcf"
        vcfFileOut = readSet + ".smCounter.anno.vcf"
        annotate.vcf_annotate.run(cfg, vcfFileIn, vcfFileOut, vc, tumorNormal)

    elif numVariants == 0: # create a header only anno.vcf from cut.vcf
        vcfFileIn  = readSet + ".smCounter.cut.vcf"
        vcfFileOut = readSet + ".smCounter.anno.vcf"
        shutil.copyfile(vcfFileIn,vcfFileOut)
    
    #else: numVariants == -1 , duplex, emptyBam
        
    # aggregate all metrics
    metrics.sum_all.run(cfg)

def run_tumor_normal(readSet,paramFile,vc):
    ''' Wrapper for tumor-normal analysis
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
    core.run_log.init(tumor)
    run((tumor,paramFile,vc),tumorNormal=True)
    print("--"*20)
    run((normal,paramFile,vc),tumorNormal=True)

    ## Compare Tumor Normal variants and update filter
    if vc == 'v2': # use new TN filter
        print("--"*20)
        cfg = core.run_config.run(tumor,paramFile)
        core.tumor_normal.tumorNormalVarFilter(cfg, normal, tumor)
        print("--"*20)

    ## Create cplx,anno.txt/vcf and sum.all files
    numVariants = core.tumor_normal.getNumVariants(normal)
    cfg = core.run_config.run(normal,paramFile)
    post_smcounter_work(numVariants, normal, cfg, tumorNormal=True)
    print("--"*20)
    cfg = core.run_config.run(tumor,paramFile)
    numVariants = core.tumor_normal.getNumVariants(tumor)
    post_smcounter_work(numVariants, tumor, cfg, tumorNormal=True)
    print("--"*20)

    ## Run old variant substraction code if using v1
    if vc == 'v1':
        print("--"*20)
        print('Warning: Doing naive substraction of normal variants from tumor. Please use smCounter-v2  for newer Tumor-Normal variant Filter')
        cfg = core.run_config.run(tumor,paramFile)
        core.tumor_normal.removeNormalVariants(cfg, normal, tumor)

    ## Run Quandico for CNV ; Note - cfg corresponds to tumor here
    core.tumor_normal.runCopyNumberEstimates(cfg)

    # close log file
    core.run_log.close()

#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   cfg=core.run_config.parse_command_line_arguments()#This function handles the command line parsing and validtion, and, in its case, the prinitng of the USAGE screen
   parser = ConfigParser.SafeConfigParser()
   parser.optionxform = str
   parser.read(cfg.paramFile)
   miscfileparts=os.path.split(misc.__file__)
   miscparentparts=os.path.split(miscfileparts[0])#Gets the current path of this file
   os.environ["PATH"]=os.environ["PATH"]+":"+miscparentparts[0]#Adds the pah of the current file to the environment
   #Anterior command line: "\nRun as : python run_qiaseq_dna.py <param_file> <v1/v2> <single/tumor-normal> <readSet(s)>\n"
   print("pid:",os.getpid())
   if "{0}" in cfg.outputPath or "{1}" in cfg.outputPath:
      cfg.outputPathTemplate=cfg.outputPath
   else:
      if len(cfg.readSet)>1:
         warnings.warn("The directory was not provided in the form of a name template and there is more than one sample. The sample name will be "+
            "append to the end of the name provided (separated by an underscore).")
         cfg.outputPathTemplate=cfg.outputPath+"_{0}"
      #else:
      #   outputpath=cfg.outputPath

   outptemplate=getattr(cfg,"outputPathTemplate",None)
   if cfg.analysis.lower() == "tumor-normal":
      for (iread,read) in enumerate(cfg.readSet):
         if outptemplate is None:
            pass#The outputpath attribute is already with its final value
         elif "{0}" in outptemplate or "{1}" in outptemplate:
            cfg.outputPath=outptemplate.format(read,iread+1)
         run_tumor_normal(read,cfg.paramFile,cfg.vc,cfg.outputPath)
   else: # Single sample, might still need to run quandico
      if len(cfg.readSet)==0:
         readsets=[sec for sec in parser.sections() if sec not in ("general","smCounter")]
      else:
         readsets=cfg.readSet
      for (iread,read) in enumerate(readsets):
         if outptemplate is None:
            pass#The outputpath attribute is already with its final value
         elif "{0}" in outptemplate or "{1}" in outptemplate:
            cfg.outputPath=outptemplate.format(read,iread+1)
         run((read,cfg.paramFile,cfg.vc,cfg.outputPath))
         #runcfg = core.run_config.run(read,cfg.paramFile,cfg.outputPath)

    #if len(sys.argv) > 6 :
    #    print "\nRun as : python run_qiaseq_dna.py <param_file> <v1/v2> <single/tumor-normal> <readSet(s)>\n"
    #    sys.exit(-1)
  
    #paramFile = sys.argv[1]
    #vc = sys.argv[2]
    #analysis = sys.argv[3]
    #readSet   = " ".join(sys.argv[4:]) # 2 readSets in case of tumor-normal
 
    #if analysis.lower() == "tumor-normal":      
    #    run_tumor_normal(readSet,paramFile,vc)
    #else: # Single sample, might still need to run quandico
    #    run((readSet,paramFile,vc),tumorNormal = False)
    #    cfg = core.run_config.run(readSet,paramFile)
    #    core.tumor_normal.runCopyNumberEstimates(cfg)

