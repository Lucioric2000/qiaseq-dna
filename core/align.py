
import subprocess
import os
import os.path

#-------------------------------------------------------------------------
# align reads to refernce genome using BWA MEM
#-------------------------------------------------------------------------
def run(cfg, readFile1, readFile2, bamFileOut):
   #run_through_tempfile(cfg, readFile1, readFile2, bamFileOut)
   #return
   # report start
   print("align: starting read alignment")

   # get some parameters from config
   readSet          = cfg.readSet
   genomeFile       = cfg.genomeFile
   bwaDir           = cfg.bwaDir
   samtoolsDir      = cfg.samtoolsDir
   deleteLocalFiles = cfg.deleteLocalFiles
   numCores         = cfg.numCores
   
   # make file names for local BWA and samtools log files
   #logFileBase     = os.path.basename(bamFileOut)
   tmpOutFile      = bamFileOut.replace(".bam",".out.tmp")
   logFileBwa      = bamFileOut.replace(".bam",".bwa.log")
   logFileSamtools = bamFileOut.replace(".bam",".samtools.log")

   # align reads to reference genome using BWA-MEM, and convert to BAM format
   #cmd = bwaDir + "bwa mem -t " +  numCores \
   #cmd = bwaDir + "bwa mem -C -t " +  numCores \
   #+ " " + genomeFile   \
   #+ " " + readFile1    \
   #+ " " + readFile2    \
   #+ " 2>" + logFileBwa \
   #+ " | " + samtoolsDir + "samtools view -1 -@ " + numCores \
   #+ " - " \
   #+ " 1> " + bamFileOut \
   #+ " 2> " + logFileSamtools
   cmd = bwaDir + "bwa mem -C -t " +  numCores \
   + " " + genomeFile   \
   + " " + readFile1    \
   + " " + readFile2    \
   + " 2>" + logFileBwa \
   + " | " + samtoolsDir + "samtools view -S -1 -@ " + numCores \
   + " - " \
   + " 1> " + bamFileOut \
   + " 2> " + logFileSamtools
   print("Executing bwa and samtools commands:",cmd,"in", cfg.outputPath, "with configuration", cfg)
   #Hasta el momento se han quitado -Cen bwa y -1 en samtools
   subprocess.check_call(cmd, shell=True, cwd=cfg.outputPath)

   # delete local fastq files if not needed anymore, to conserve local disk space
   if deleteLocalFiles:
      for readFile in (readFile1, readFile2):
         if len(os.path.dirname(readFile)) == 0:
            os.remove(readFile)
   #os.remove(tmpOutFile)
   # report completion
   print("align: done aligning reads to genome, bam file: " + bamFileOut)
def run_through_tempfile(cfg, readFile1, readFile2, bamFileOut):
    # report start
    print("align: starting read alignment")

    # get some parameters from config
    readSet          = cfg.readSet
    genomeFile       = cfg.genomeFile
    bwaDir           = cfg.bwaDir
    samtoolsDir      = cfg.samtoolsDir
    deleteLocalFiles = cfg.deleteLocalFiles
    numCores         = cfg.numCores
   
    # make file names for local BWA and samtools log files
    #logFileBase     = os.path.basename(bamFileOut)
    tmpOutFile      = bamFileOut.replace(".bam",".out.tmp")
    logFileBwa      = bamFileOut.replace(".bam",".bwa.log")
    logFileSamtools = bamFileOut.replace(".bam",".samtools.log")

    ## make file names for local BWA and samtools log files
    #logFileBase     = os.path.basename(bamFileOut)
    #logFileBwa      = logFileBase.replace(".bam",".bwa.log")
    #logFileSamtools = logFileBase.replace(".bam",".samtools.log")


    # align reads to reference genome using BWA-MEM, and convert to BAM format
    #cmd1 = bwaDir + "bwa mem -t " +  numCores \
    cmd1 = bwaDir + "bwa mem -C -t " +  numCores \
    + " " + genomeFile   \
    + " " + readFile1    \
    + " " + readFile2    \
    + " 1>" + tmpOutFile \
    + " 2>" + logFileBwa
    #-S is new
    cmd2= samtoolsDir + "samtools view -S -1 -@ " + numCores \
    + " "+ tmpOutFile \
    + " 1> " + bamFileOut \
    + " 2> " + logFileSamtools
    print("Executing bwa command:",cmd1,"and samtools command:",cmd2)
    #Hasta el momento se han quitado -Cen bwa y -1 en samtools
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)

    # delete local fastq files if not needed anymore, to conserve local disk space
    if deleteLocalFiles:
        for readFile in (readFile1, readFile2):
            if len(os.path.dirname(readFile)) == 0:
                os.remove(readFile)
    #os.remove(tmpOutFile)
    # report completion
    print("align: done aligning reads to genome, bam file: " + bamFileOut)
