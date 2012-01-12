#!/usr/bin/env python2.7

# Francesco Favero, CBS DTU Jan 2012

# This might be more improved by you! Feel free to contribute...

import argparse, os, sys, getpass, multiprocessing, subprocess


def check_dir(directory):
   '''
   Check the given directory and return the whole path
   or set the default to working directory.
   '''
   if directory == None:
      directory = os.getcwd()
   if os.path.isdir(directory):
      directory = os.path.abspath(directory)
   else:
      sys.exit("Not existing directory, please create the directory first")
   return directory

def walklevel(directory, level=1):
    directory = directory.rstrip(os.path.sep)
    assert os.path.isdir(directory)
    num_sep = directory.count(os.path.sep)
    for root, dirs, files in os.walk(directory):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]

def one_tophat(tophat,reference,directory,innerdist,ncpu):
   '''
   Assemble all commands to perform tophat on a single 
   directory. This will be used as a core module in a 
   more generic function to handle directories in the
   input.
   '''
   # Check the fastq.gz files contained in the directory
   fastqlist  = []
   for dirname, dirnames, filenames in walklevel(directory,level=1):
      for filename in filenames:
         if filename.endswith(".fastq.gz"):
            fastqlist.append(os.path.join(dirname, filename))
   # Check single/paired end.. eg: if there are *_1.fastq.gz
   # *_2.fastq we have paired end.. otherwise is single end.
   first_end  = []
   second_end = []
   for fastqfile in fastqlist:
      if fastqfile.endswith("_1.fastq.gz"):
         first_end.append(fastqfile)
      if fastqfile.endswith("_2.fastq.gz"):
         second_end.append(fastqfile)
   if first_end == second_end == []:
      paired_end = False
   else:
      if first_end == []:
         sys.exit('Error. In a paired end you need both end. You missing the *_1.fastq files!')
      elif second_end == []:
         sys.exit('Error. In a paired end you need both end. You missing the *_2.fastq files!')
      else:
         paired_end = True
   # Decompress fastq.gz to fastq task.. tophat have
   # to look into le file so can't do it with pipe 
   fastq_gunzip = ['gunzip']
#   if not paired_end:
#      fastq_gunzip.append(' '.join(fastqlist).strip())
   fastq_gunzip.append(' '.join(fastqlist).strip())
#   else:
#      fastq_gunzip.append(' '.join(first_end).strip())
#      fastq_gunzip.append('&')
#      fastq_gunzip.append('gunzip')
#      fastq_gunzip.append(' '.join(second_end).strip())
   # Creating the pipeline for the tophat command switching
   # from single or pired end if needed.
   cmd = [tophat]
   # Add number of CPU
   cmd.append('--num-threads')
   cmd.append(ncpu)
   # Add otput arguments
   cmd.append('--output-dir')
   cmd.append(os.path.join(directory,'tophat_out'))
   if paired_end:
      if innerdist:
         cmd.append('--mate-inner-dist')
         cmd.append(innerdist)
      else:
         sys.exit('Error. In paired end you MUST specify the --mate-inner-dist value. See tophat manual: http://tophat.cbcb.umd.edu/manual.html')
   cmd.append(reference)
   fastq_regzip = ['gzip']
   if paired_end:
      first_end_unzip  = []
      second_end_unzip = []
      for fastqfile in first_end:
         first_end_unzip.append(fastqfile[:-3])
         fastq_regzip.append(fastqfile[:-3])
      for fastqfile in second_end:
         second_end_unzip.append(fastqfile[:-3])
         fastq_regzip.append(fastqfile[:-3])
      cmd.append(','.join(first_end_unzip).strip())
      cmd.append(','.join(second_end_unzip).strip())
   else:
      fastqlist_unzip = []
      for fastqfile in fastqlist:
         fastqlist_unzip.append(fastqfile[:-3])
         fastq_regzip.append(fastqfile[:-3])
      cmd.append(','.join(fastqlist_unzip).strip())
   user      = getpass.getuser()
   last_dir  = os.path.split(directory)[1] 
   toph_proc = {'gunzip' : {user + '_gunzip_' + last_dir : fastq_gunzip},'tophat' : {user + '_tophat_' + last_dir : cmd}, 'regzip' : {user + '_regzip_' + last_dir : fastq_regzip}} 
   return toph_proc

def one_htseq(htseq,samtools,gtf,directory,idattr,mode):
   '''
   Perform the htseq script for the generated BAM
   file from TopHat.. just for one directory. It 
   will be used as a core module of others  
   functions
   '''
   bam_file = os.path.join(directory,'tophat_out/accepted_hits.bam')
   cmd = [samtools]
   cmd.append('view')
   cmd.append(bam_file)
   cmd.append('|')
   cmd.append(htseq)
   if mode:
      cmd.append('--mode')
      cmd.append(mode)
   if idattr:
      cmd.append('--idattr')
      cmd.append(idattr)
   cmd.append('--quiet')
   cmd.append('-')
   cmd.append(gtf)
   return cmd

def pyXmsub(cmd,name,output_err='/dev/null',output_out='/dev/null',dir_exec=None,depends=None,queue=None,partition=None,walltime=None,ncpu=None,mem=None):
   '''
   Helps with the xmsub queue system... nothing special
   '''
   if not name:
      sys.exit('Any xmsub processes might have a name... something went wrong')
   lparms = []
   if mem:
      lparms.append('mem='+str(mem))
   if walltime:
      lparms.append('walltime='+str(walltime))
   if depends:
      lparms.append('depend='+str(depends))
   if ncpu:
      lparms.append('procs='+str(ncpu))
   if partition:
      lparms.append('partition='+str(partition))
      if partition == 'uv':
         lparms.append('flags=sharedmem')
   xmsub = ['/usr/local/bin/xmsub']
   xmsub.append('-N')
   xmsub.append(str(name))
   if queue:
      xmsub.append('-q')
      xmsub.append(str(queue))
   if dir_exec:
      xmsub.append('-d')
      xmsub.append(str(dir_exec))
   xmsub.append('-o') # no STDIN output file in working dir
   xmsub.append('/dev/null') # no STDIN output file in working dir
   xmsub.append('-e') # no STDIN error file in working dir
   xmsub.append('/dev/null') # no STDIN error file in working dir
   xmsub.append('-ro')
   xmsub.append(output_out)
   xmsub.append('-re')
   xmsub.append(output_err)
   if lparms != []:
      xmsub.append('-l')
      xmsub.append(','.join(lparms))  
   xmsub.append('-r')
   xmsub.append('y')
   xmsub.append('-de')
   xmsub.append('"' + cmd + '"')
   return xmsub

def dir_wrapper(directory,prefix):
   '''
   Check if the directory contains subdirectory, and it return a list of path
   (or a single directory name if it doesn't contain any)
   '''
   dirlist = []
   for root, dirnames, filenames in walklevel(directory,level=0):
      for name in dirnames:
         if name.startswith(str(prefix)):
            dirlist.append(os.path.join(root, name))
   if not dirlist:
      dirlist = [directory]
   return dirlist

def main():
   '''
   Execute main function taking arguments
   '''
   parser = argparse.ArgumentParser(description='Cooking a dataset of fastq.gz with Tophat trough the msub system at CBS')
   parser.add_argument('data_dir', metavar='DIR_IN',
                   help='Directory holding the samples (in GEO for example GSE is the dataset containing the GSMs samples)')
   tophat = parser.add_argument_group(title='TopHat',description='Argument for TopHat')
   tophat.add_argument('--reference',dest='reference',
                   help='Specify the reference genome to use. Note, it must be in bowtie format')
   tophat.add_argument('-r', '--mate-inner-dist',dest='innerdist',
                   help=' from TopHat manual: This is the expected (mean) inner distance between mate pairs. For, example, for paired end runs with fragments selected at 300bp, where each end is 50bp, you should set -r to be 200. There is no default, and this parameter is required for paired end runs.')
   htseq  = parser.add_argument_group(title='HTSeq',description='Argument for the htseq-count script')
   htseq.add_argument('--gtf', dest='gtf',
                   help='Give the GTF file tom map the reference genome. The genome reference and the GTF have to be related to make it sense.')
   htseq.add_argument('--idattr', dest='idattr',
                   help='Tells HTSeq wich ID to take from the GTF; eg: gene_id, ens_genese... ens_transcript (but depends on the contents of the gtf file)')
   htseq.add_argument('--mode', dest='mode',
                   help='Is the HTSeq-count mode. possible option are: union, intersection-strict and intersection-nonempty, default is union.')
   output = parser.add_argument_group(title='Output',description='Argument that involve the output destination')
   output.add_argument('-d','--dir', dest='dir_out',
                   help='Set root directory to save files (default working directory)')
   msub = parser.add_argument_group(title='Msub',description='Optional argument that control the Msub queue system')
   msub.add_argument('--queue', dest='queue',
                   default='cbs',
                   help='Queue to use, default is cbs')
   msub.add_argument('--ncpu', dest='ncpu',
                   default=None,
                   help='Number of CPU to use')
   msub.add_argument('--mem', dest='mem',
                   default=None,
                   help='Amount of RAM to reserve for the process')
   msub.add_argument('--partition', dest='partition',
                   default='uv',
                   help='Cluster to use, default is uv')

   tweaks = parser.add_argument_group(title='Tweaks',description='Optional argument that might effect the processes')
   tweaks.add_argument('--samples_prefix', dest='samples_prefix',
                   default='GSM',
                   help='Set the prefix of the directory with the sample files, default is GSM')
   tweaks.add_argument('--tophat', dest='toph',
                   default='/tools/bin/tophat',
                   help='Path of the tophat executable')
   tweaks.add_argument('--htseq', dest='htseq',
                   default='/home/projects3/sbge_cancer_seq/sequence/RNAseq/src/myPy/bin/htseq-count',
                   help='Path of the tophat executable')
   tweaks.add_argument('--samtools', dest='samtools',
                   default='/usr/cbs/bio/bin/samtools',
                   help='Path of the samtools executable')
   tweaks.add_argument('--noqueue', dest='noqueue',action='store_true',
                   help='Do not use the msub sytem, just run it on the shell (use screen is suggested)')
   args      = parser.parse_args()
   data_dir  = check_dir(args.data_dir)
   last_dir  = check_dir(args.dir_out)
   reference = args.reference
   gtf       = args.gtf
   # Chech the complete path for the reference
   if reference:
      reference_path = os.path.split(reference)
      reference_base = check_dir(reference_path[0])
      reference      = os.path.join(reference_base,reference_path[1])
   else:
      sys.exit('Error. You MUST specify a reference genome. see --help or tophat manual.')
   # Chech the complete path for the GTF file
   if gtf:
      gtf_path = os.path.split(gtf)
      gtf_base = check_dir(gtf_path[0])
      gtf      = os.path.join(gtf_base,gtf_path[1])
   else:
      sys.exit('For HTSeq you have to specify the .gtf file to index the genome')
   #print one_tophat(args.toph,reference,data_dir,args.innerdist)
   #hts = one_htseq(args.htseq,args.samtools,args.gtf,data_dir,args.idattr,args.mode)
   #print " ".join(pyXmsub(' '.join(hts),'HTSeq1','tmp/err.txt','tmp/out.txt',dir_exec='./',queue='cbs',partition='uv',ncpu=124))
   samples = dir_wrapper(data_dir,args.samples_prefix)
   for sample in samples:
      tophat     = one_tophat(args.toph,reference,sample,args.innerdist,args.ncpu)
      htsq       = one_htseq(args.htseq,args.samtools,gtf,sample,args.idattr,args.mode)
      gunzip_id  = tophat['gunzip'].keys()[0]
      gunzip_out = os.path.join(sample,'logs',gunzip_id + '_o.txt')
      gunzip_err = os.path.join(sample,'logs',gunzip_id + '_e.txt')
      gunzip_cmd = ' '.join(tophat['gunzip'][gunzip_id])
      gunzip_xms = " ".join(pyXmsub(str(gunzip_cmd),str(gunzip_id),output_err=gunzip_err,output_out=gunzip_out,dir_exec=sample,queue=args.queue,partition=args.partition,ncpu=args.ncpu))
      tophat_id  = tophat['tophat'].keys()[0]
      tophat_out = os.path.join(sample,'logs',tophat_id + '_o.txt')
      tophat_err = os.path.join(sample,'logs',tophat_id + '_e.txt')
      tophat_cmd = ' '.join(tophat['tophat'][tophat_id])
      tophat_xms = " ".join(pyXmsub(str(tophat_cmd),str(tophat_id),output_err=tophat_err,output_out=tophat_out,dir_exec=sample,queue=args.queue,partition=args.partition,ncpu=args.ncpu,mem=args.mem,depends=gunzip_id))
      regzip_id  = tophat['regzip'].keys()[0]
      regzip_out = os.path.join(sample,'logs',regzip_id + '_o.txt')
      regzip_err = os.path.join(sample,'logs',regzip_id + '_e.txt')
      regzip_cmd = ' '.join(tophat['regzip'][regzip_id])
      regzip_xms =  " ".join(pyXmsub(str(regzip_cmd),str(regzip_id),output_err=regzip_err,output_out=regzip_out,dir_exec=sample,queue=args.queue,partition=args.partition,ncpu=args.ncpu,depends=str(tophat_id)))
      last_dir   = os.path.split(sample)[1] 
      htseq_id   = getpass.getuser() + '_htsqCounts_' + last_dir
      htseq_out  = os.path.join(sample,last_dir + '_htseq_o.txt')
      htseq_err  = os.path.join(sample,last_dir + '_htseq_e.txt')
      htseq_cmd  = ' '.join(htsq)
      htseq_xms  = " ".join(pyXmsub(str(htseq_cmd),str(htseq_id),output_err=htseq_err,output_out=htseq_out,dir_exec=sample,queue=args.queue,partition=args.partition,ncpu=args.ncpu,mem=args.mem,depends=str(regzip_id)))
      if args.noqueue:
         print gunzip_cmd
         work0      = subprocess.Popen(gunzip_cmd,stdout=subprocess.PIPE,shell=True)
         out0       = work0.communicate()[0]
         print tophat_cmd
         work1      = subprocess.Popen(tophat_cmd,stdout=subprocess.PIPE,shell=True)
         out1       = work1.communicate()[0]
         print regzip_cmd
         work2      = subprocess.Popen(regzip_cmd,stdout=subprocess.PIPE,shell=True)
         out2       = work2.communicate()[0]
         print htseq_cmd
         with open(htseq_out,"wb") as out:
            with open(htseq_err,"wb") as err:
               work3      = subprocess.Popen(htseq_cmd,stdout=out,stderr=err,shell=True)
               out3       = work3.communicate()[0]
      else:
         print gunzip_xms
         work0      = subprocess.Popen(gunzip_xms,stdout=subprocess.PIPE,shell=True)
         out0       = work0.communicate()[0]
         print tophat_xms
         work1      = subprocess.Popen(tophat_xms,stdout=subprocess.PIPE,shell=True)
         out1       = work1.communicate()[0]
         print regzip_xms
         work2      = subprocess.Popen(regzip_xms,stdout=subprocess.PIPE,shell=True)
         out2       = work2.communicate()[0]
         print htseq_xms
         work3      = subprocess.Popen(htseq_xms,stdout=subprocess.PIPE,shell=True)
         out3       = work3.communicate()[0]
         print 'unzip '+ last_dir +' job in queue id: ' + out0.strip()
         print 'Tophat for ' + last_dir + ' job in queue id: ' + out1.strip()
         print 'gzip files ' + last_dir + ' job in queue id: ' + out2.strip()
         print 'HTSeq count '+ last_dir +' job in queue id: ' + out3.strip()





if __name__ == "__main__":
    main()
