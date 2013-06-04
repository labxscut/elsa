#!/usr/bin/env python
# submit jobs in $1 while total cores<300, total mem<300G
# current limit is 299
# current mem limit is 300G

import os, shutil, subprocess, argparse, fnmatch, time
core_max=63
mem_max=300
uname="lixia"
qname="main"

def mem_size(mem):
  if mem[-2:]=='mb':
    return int(mem[:-2])/1000
  elif mem[-2:]=='gb':
    return int(mem[:-2])
  else:
    return -1

def peek_current(uname):
  tmp=subprocess.Popen("qstat -u %s | awk '{print $8}'" % uname, shell=True, stdout=subprocess.PIPE).communicate()
  mems=tmp[0].split('\n')
  #print mems
  tmp=subprocess.Popen("qstat -u %s | awk '{print $6}'" % uname, shell=True, stdout=subprocess.PIPE).communicate()
  cores=tmp[0].split('\n')
  #print cores
  tmp=subprocess.Popen("qstat -u %s | awk '{print $5}'" % uname, shell=True, stdout=subprocess.PIPE).communicate()
  sessions=tmp[0].split('\n')
  #print sessions
  hskip=5
  tskip=1
  total_core=0
  total_mem=0
  #print queue
  #assert len(mems) == len(cores)

  for i in range(hskip+1,len(sessions)-tskip):
    try:
      session=int(sessions[i])
    except ValueError:
      return True #if seesion is not started, just wait

  for i in range(hskip+1,len(cores)-tskip):
    total_core+=int(cores[i])
  
  for i in range(hskip+1,len(mems)-tskip):
    if mem_size(mems[i])>0:
      total_mem+=mem_size(mems[i])
    else: #Let's just wait
      total_core=core_max
      total_mem=mem_max
      break

  full_status=(total_core>=core_max) or (total_mem>=mem_max)
  return full_status

def submit( pbsFile ):
  tmp=subprocess.Popen("qsub -q %s %s" % (qname, pbsFile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
  submitted=True
  if tmp[1] == '':
    print tmp[0]
    print pbsFile, "submitted"
  else:
    print tmp[1]
    print pbsFile, "error"
    submitted=False
  return submitted

def main():

  parser = argparse.ArgumentParser(description="MCB Queue Checking and Submission Tool")
  parser.add_argument("pbsFile", metavar="pbsFile", help="single pbs file to be submitted")
  arg_namespace = parser.parse_args()
  pbsFile = vars(arg_namespace)['pbsFile']
  pbsFiles = set([pbsFile])

  #print pbsFiles

  while( len(pbsFiles) != 0 ):
    if peek_current(uname):  #full
      time.sleep(1)
    else:
      pbsFile=pbsFiles.pop()
      #print len(pbsFiles)
      if not submit(pbsFile):
        pbsFiles.add(pbsFile)

if __name__=="__main__":
  main()
