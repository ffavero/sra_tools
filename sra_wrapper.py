#!/usr/bin/env python2.7

# Francesco Favero, CBS DTU Nov 2011

# Wrapper to retrive information about NGS public dataset on NCBI
# GEO and download the whole dataset in fastq.gz format (for paired
# and single-end).

# This might be more improved by you! Feel free to contribute...

from xml.etree import ElementTree as ET
from urllib2 import urlopen, Request
from urllib import urlencode
from urlparse import urlparse
import ftplib, contextlib, argparse, os, sys, socket, getpass, multiprocessing, subprocess

def sra_tree(host,user,directory):
   '''
   Get all the .sra file given a directory.
   wrapping the sra directory from a GSM in 
   GEO, will give all the sra files for the 
   requested sample.
   '''
   # The server can be reeeally slow just
   # increment the timeout to 2 minutes
   socket.setdefaulttimeout(120)
   ftp = ftplib.FTP(host)
   ftp.login(user)
   dirs  = []
   files = []
   try:
      dirs = ftp.nlst( directory )
   except ftplib.error_perm, resp:
      if str(resp) == "550 No files found":
         print "no files in this directory"
      else:
         raise
   for d in dirs:
      tmp = []
      try:
         tmp = ftp.nlst(d)
      except ftplib.error_perm, resp:
         if str(resp) == "550 No files found":
            print "no files in this directory"
         else:
            raise
      ftp.sendcmd("TYPE i")
      for f in tmp:
         size = ftp.size(f)/1048576
         if size != None:
            if f.endswith('.sra'):
               files.append({'name':f,'size':size})
   ftp.quit()
   return files

def get_sra_gsm(gsm_id):
   '''
   Get the sra url referring to the gsm_id
   '''
   if gsm_id.startswith('GSM'):
      url   = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'
      NS = '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}'
      path  = './/' + NS + 'Supplementary-Data'
      values = {'acc' : gsm_id,
                'form' : 'xml',
                'view' : 'brief',
                'targ' : 'gsm' }
      data = urlencode(values)
      req = Request(url, data)
      try:
         with contextlib.closing(urlopen(req)) as xml:
            geo = ET.parse(xml)
         elts = geo.findall(path)
         sra_url = []
         for elem in elts:
            content = elem.get('type')
            if content == "SRA Experiment":
               sra_url.append(elem.text)
         if not sra_url:
            sys.exit('Sorry the dataset does not contains SRA information.')
         else:
            return sra_url
      except ET.ParseError:
         sys.exit('Error parsing the XML, probably the file does not exists')
   else:
      sys.exit(str(gsm_id)+' is not a GSM id.')

def get_gsm_list(gse_id):
   '''
   Get all the GSM of a GSE serie
   '''
   if gse_id.startswith('GSE'):
      url   = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'
      NS = '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}'
      path  = './/' + NS + 'Sample'
      values = {'acc' : gse_id,
                'form' : 'xml',
                'view' : 'brief',
                'targ' : 'gse' }
      data = urlencode(values)
      req = Request(url, data)
      try:
         with contextlib.closing(urlopen(req)) as xml:
            geo = ET.parse(xml)
         elts = geo.findall(path)
         gsm_list = []
         for elem in elts:
            content = elem.get('iid')
            if content.startswith('GSM'):
               gsm_list.append(content)
         return gsm_list
      except ET.ParseError:
         sys.exit('Error parsing the XML, probably the file does not exists')
   else:
      sys.exit(str(gse_id)+' is not a GSE id.')

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

def gse_parse(gse_id,tmp_dir,logs,sra,paired):
   '''
   just the case in which a GSE is specified
   '''
   gsm_list = get_gsm_list(gse_id)
   geo = GEOobj(gse_id,gsm_list,tmp_dir,logs,sra,paired)
   return geo
 
def gsm_parse(gsm_id,tmp_dir,logs,sra,paired):
   '''
   just the case in which a GSM is specified
   '''
   gsm_list = [gsm_id]
   geo = GEOobj(gsm_id,gsm_list,tmp_dir,logs,sra,paired)
   return geo

class GEOobj:
   '''
   Just create an object with the information related to
   the given argument(s). If a GSE is specified, the    
   object wiull contain all the relative GSMs, if a GSM 
   is specified the object will contains the given GSM. 
   '''
   def __init__(self,data,gsm_list,tmp_dir,logs,sra,paired):
      self.parent   = data
      self.gsm_list = gsm_list
      self.tmp_dir  = tmp_dir
      self.logs     = logs 
      self.sra      = sra
      self.paired   = paired
      gsm_dict = {}
      for gsm in gsm_list:
         gsm_dict[gsm]          = {}
         gsm_dict[gsm]['files'] = []
         urls      = get_sra_gsm(gsm)
         tmpurl    = urls[0]
         tmp_parts = urlparse(tmpurl.strip())
         for url in urls:
            url_parts = urlparse(url.strip())
            if tmp_parts.netloc == url_parts.netloc:
               tmp = sra_tree(url_parts.netloc,'anonymous',url_parts.path)
               for f in tmp:
                  gsm_dict[gsm]['files'].append(f)
            else:
               sys.exit('Cannot handle file hosted in differents repositories')
         gsm_dict[gsm]['host'] = tmp_parts.netloc
      self.gsm_dict = gsm_dict
   def summary(self):
      '''
      Summary of the files, parts and size 
      of the dataset
      '''
      data        = self.gsm_dict
      samples     = {}
      tot_size    = 0
      tot_files   = 0
      for sample in data:
         samples[sample] = {}
         size = 0
         for d in data[sample]['files']:
            size += d['size']
         n = len(data[sample]['files'])
         tot_size       += size
         tot_files      += n
         samples[sample] = {'size': size,'n': n}
      tot_samples = len(samples)
      summary     = '\n\t####\t\n'
      summary    += '\t####\t' + self.parent + ' is composed by ' + str(tot_samples) + ' samples \n'
      summary    += '\t####\ttotal number of SRA files to download: ' + str(tot_files) + '\n'   
      summary    += '\t####\toverall download size: ' + str(float(tot_size/1000)) + 'GB \n'
      summary    += '\t####\t\n'
      print summary
   def tree(self,root_directory,fastq):
      '''
      Generate the directory structure: one root directory
      (the GEO_ID passed to the args) and the relative
      subforlders (the GSMs). Each of these subfolder contains
      a metadata, with the various .sra files information on it
      '''
      data      = self.gsm_dict
      root_path = os.path.join(root_directory,self.parent)
      if os.path.isdir(root_path):
         print root_path + ' Is already there....'
      else:
         os.mkdir(root_path)
      gsm_deps = None
      for gsm in self.gsm_list:
         gsm_path = os.path.join(root_path,gsm)
         if os.path.isdir(gsm_path):
            print gsm_path + ' Is already there....'
         else:
            os.mkdir(gsm_path)
         with open(os.path.join(gsm_path,'metainfo.txt'),'wb') as metainfo:
            host      = data[gsm]['host']
            file_list = data[gsm]['files']
            for file_line in file_list:
               metainfo.write('ftp://'+host+file_line['name']+'\n')
         if fastq:
            sra_attrs = data[gsm]['files']
            gsm_deps = sra_to_fastq(gsm_path,self.tmp_dir,sra_attrs,'anonymous',self.sra,self.paired,self.logs,gsm_deps)

def sra_to_fastq(directory,tmp_dir,sra_attrs,user,sra_sdk,paired,log,gsm_deps):
   '''
   Reads the metainfo.txt file in the directory
   Download the file in the tmpdir and write the 
   output in fastq.gz in the directory.
   It sends the request one by one in the queue 
   system
   '''
   last = None
   with open(os.path.join(directory,'metainfo.txt'),'rb') as metainfo:
      url_deps = None
      for metaline in metainfo:
         ftp_url    = metaline.strip()
         ftp_parts  = urlparse(ftp_url)
         ftp_path   = ftp_parts.path
         sra_size   = 0
         sra_status = False
         for sra_attr in sra_attrs:
            if sra_attr['name'] == ftp_path:
               sra_size = sra_attr['size']
               sra_status = True
         if not sra_status:
            sys.exit('Weird... the file '+ ftp_path +' has not being found in the parent object.')
         try:
            downtime   = int(sra_size/0.1) # assuming low download speed at 0.1 M/s
         except:
            sys.exit('Something went wrong on retriving the size of '+ ftp_path)
         if downtime == 0:
            sys.exit('Something went wrong on retriving the size of '+ ftp_path)
         filename   = os.path.split(ftp_path)[1]
         file_out   = ''
         if '.sra' in filename[-4:]:
            file_out = filename[:-4] + '.fastq.gz'
         else:
            sys.exit('Strange the file '+ filename +' has not the sra extension')
         user         = getpass.getuser()
         proc_dw_sra  = user + '_' + 'download_'+ filename
         proc_sra     = user + '_' + filename
         proc_rm_sra  = user + '_' + 'rm_' + filename
         last         = proc_rm_sra
         pipe_sra     = sra_fastq_pipe(sra_sdk,os.path.join(tmp_dir,filename),paired)
         ncpu         = multiprocessing.cpu_count()
         resource0    = '-l walltime='+str(downtime)
         if gsm_deps:
            resource0    = '-l walltime='+str(downtime)+',depend='+str(gsm_deps)
         if url_deps:
            resource0    = '-l walltime='+str(downtime)+',depend='+str(url_deps)
         url_deps     = proc_dw_sra
         resource1    = '-l procs='+str(ncpu)+',depend='+str(proc_dw_sra)
         proc0        = []
         proc0.append('/usr/local/bin/xmsub')
         proc0.append('-N')
         proc0.append(proc_dw_sra)
         proc0.append('-q')
         proc0.append('cbs')            
         proc0.append('-d')
         proc0.append(directory)
         proc0.append('-o') # no STDIN output file in working dir
         proc0.append('/dev/null') # no STDIN output file in working dir
         proc0.append('-e') # no STDIN error file in working dir
         proc0.append('/dev/null') # no STDIN error file in working dir
         if log:      
            msg_dir = os.path.join(directory,'logs')
            if os.path.isdir(msg_dir):
               pass
            else:
               os.mkdir(msg_dir)
            proc0.append('-re')
            proc0.append(os.path.join(msg_dir,proc_dw_sra+'.err'))
         proc0.append('-de')
         proc0.append('-ro')
         proc0.append((os.path.join(tmp_dir,filename)))
         proc0.append(resource0)
         proc0.append('/usr/bin/curl')
         proc0.append('--connect-timeout')
         proc0.append('120')
         proc0.append(ftp_url)
         line0 = ''
         for proc in proc0:
            line0 += str(proc) + ' '
         line0 = line0.strip()
         proc1        = []
         proc1.append('/usr/local/bin/xmsub')
         proc1.append('-N')
         proc1.append(proc_sra)
         proc1.append('-q')
         proc1.append('cbs')
         proc1.append('-d')
         proc1.append(directory)
         proc1.append('-o') # no STDIN output file in working dir
         proc1.append('/dev/null') # no STDIN output file in working dir
         proc1.append('-e') # no STDIN error file in working dir
         proc1.append('/dev/null') # no STDIN error file in working dir
         if log:      
            msg_dir = os.path.join(directory,'logs')
            if os.path.isdir(msg_dir):
               pass
            else:
               os.mkdir(msg_dir)
            proc1.append('-re')
            proc1.append(os.path.join(msg_dir,proc_sra+'.err'))
            if paired:
               proc1.append('-ro')
               proc1.append(os.path.join(msg_dir,proc_sra+'.out'))  
         proc1.append('-de')
         if not paired:
            proc1.append('-ro')
            proc1.append(os.path.join(directory,file_out))         
         proc1.append(resource1)
         proc1.append('"'+pipe_sra+'"')
         line1 = ''
         for proc in proc1:
            line1 += str(proc) + ' '
         line1 = line1.strip()
         proc2        = []
         proc2.append('/usr/local/bin/xmsub')
         proc2.append('-N')
         proc2.append(proc_rm_sra)
         proc2.append('-q')
         proc2.append('cbs')
         proc2.append('-d')
         proc2.append(directory)
         proc2.append('-o') # no STDIN output file in working dir
         proc2.append('/dev/null') # no STDIN output file in working dir
         proc2.append('-e') # no STDIN error file in working dir
         proc2.append('/dev/null') # no STDIN error file in working dir
         if log:      
            msg_dir = os.path.join(directory,'logs')
            if os.path.isdir(msg_dir):
               pass
            else:
               os.mkdir(msg_dir)
            proc2.append('-re')
            proc2.append(os.path.join(msg_dir,proc_rm_sra+'.err'))
         proc2.append('-de')
         proc2.append('-ro')
         proc2.append('/dev/null')
         proc2.append('-l')
         proc2.append('depend='+proc_sra)
         proc2.append('/bin/rm')
         proc2.append('-f')
         proc2.append(os.path.join(tmp_dir,filename))
         line2 = ''
         for proc in proc2:
            line2 += str(proc) + ' '
         line2 = line2.strip()
#         print line0+'\n'
#         print line1+'\n'
#         print line2+'\n'
         work0        = subprocess.Popen(line0,stdout=subprocess.PIPE,shell=True)
         out0         = work0.communicate()[0]
         work1        = subprocess.Popen(line1,stdout=subprocess.PIPE,shell=True)
         out1         = work1.communicate()[0]
         work2        = subprocess.Popen(line2,stdout=subprocess.PIPE,shell=True)
         out2         = work2.communicate()[0]
         print 'Downlod '+filename+' job in queue id: ' + out0.strip()
         print 'SRA to '+file_out+' job in queue id: ' + out1.strip()
         print 'Remove '+filename+' job in queue id: ' + out2.strip()
   return last

def sra_fastq_pipe(sra_sdk,sra_file,paired):
   '''
   Create the commands pipe to convert .sra to
   fastg.gz
   '''
   fastq_dump  = os.path.join(sra_sdk,'fastq-dump')
   sdk_help    = subprocess.Popen([fastq_dump,'--help'],stdout=subprocess.PIPE)
   sdk_help    = sdk_help.communicate()[0]
   sdk_help    = sdk_help.split('\n')
   contains_Z  = False
   for help_line in sdk_help:
      if '--stdout' in help_line:
         contains_Z = True
   if not contains_Z:
      sys.exit('The sra_sdk does not have the --stdout (-Z) option, try to use another version (eg.: 2.1.7)')
   pipe        = ''
   cmd         = []
   cmd.append(fastq_dump)
   if not paired:
      cmd.append('--stdout')
      cmd.append(sra_file)
      cmd.append('|')
      cmd.append('gzip')
   else:
      cmd.append('--split-files')
      cmd.append('--gzip')
      cmd.append(sra_file)
   for ar in cmd:
      pipe +=ar + ' '
   return pipe.strip()
   
def main():
   '''
   Execute main function taking arguments
   '''
   parser = argparse.ArgumentParser(description='Wrapper to process the SRA of a GEO experiment. It is allowed to pass either a GSE experiment id or a single GSM sample id (--gsm). The script will provide info about the content of the SRA included in the specified item (--print) or can generate a tree directory structure containing the urls of the various .sra files (--just-tree). If nothing is specified, the various .sra file will be downloaded and transformed in fastq.gz using the xmsub queue system. Each job_id submitted to the queue will be also printed in STDOUT.')
   parser.add_argument('geo_id', metavar='GEO_ID',
                   help='A GEO id, a series id (GSE) or a sample id (GSM)')
   switches = parser.add_argument_group(title='Actions',description='Action to perform, if nothing is specified it will downloads all the files.')
   switches.add_argument('-p','--print', dest='summary', action='store_true',
                   help='Print a summary of the dataset, numbers of .sra files and total size')
   switches.add_argument('-t','--just-tree', dest='justtree', action='store_true',
                   help='Build ONLY the directory structure and .sra urls info.')
   output = parser.add_argument_group(title='Input/Output',description='Argument that involve the output destination')
   output.add_argument('--paired', dest='paired',action='store_true', 
                   help='Tell sra tolkit that it is a paired-end experiment')
   output.add_argument('--gsm', dest='geo',action='store_const', 
                   const=gsm_parse, default=gse_parse,
                   help='Process a GSM and not a GSE')
   output.add_argument('-d','--dir', dest='dir',
                   help='Set root directory to save files (default working directory)')
   output.add_argument('--quiet', dest='no_logs', action='store_true',
                   help='Do not write the xmsub logs')
   tweaks = parser.add_argument_group(title='Tweaks',description='Optional argument that might effect the processes')
   tweaks.add_argument('--srasoft', dest='sra',
                   default='/home/people/favero/src/sratoolkit.2.1.7/',
                   help='Path of the SRA_SDK directory (>=2.1.7)')
   tweaks.add_argument('--tmp_dir', dest='tmp_dir', default='/scratch/',
                   help='Temp directory where downloads the .sra files')
   args      = parser.parse_args()
   directory = check_dir(args.dir)
   tmp_dir   = check_dir(args.tmp_dir)
   logs      = True
   if args.no_logs:
      logs = False
   if args.summary:
      args.geo(args.geo_id,tmp_dir,logs,args.sra,args.paired).summary()
   else:
      if args.justtree:
         args.geo(args.geo_id,tmp_dir,logs,args.sra,args.paired).tree(directory,False)
         #sys.exit('Done\n')
      else:
         args.geo(args.geo_id,tmp_dir,logs,args.sra,args.paired).tree(directory,True)

if __name__ == "__main__":
    main()
