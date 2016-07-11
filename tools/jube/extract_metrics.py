###################################################################
#                extract_metrics.py                               #
#  Created : March 2016                                           #
#                                                                 #
#  Author:                                                        #
#     Matthieu Haefele                                            #
#     matthieu.haefele@maisondelasimulation.fr                    #
#     Maison de la Simulation USR3441                             #
#  Extract metrics from a trace file                              #
###################################################################

import os
import subprocess
import re
import math
import sys
import json
import argparse

try:
  import numpy
except ImportError:
  pass

def get_output(cmd, stdout=True):
  lines=[]
  if stdout:
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    for m in proc.stdout:
      lines.append(m[:-1])
  else:
    proc= subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    for m in proc.stderr:
      lines.append(m[:-1])
                                              
  return lines

def extract(cmd, res_type='f'):
  answer = get_output(cmd)
  for s in answer:
    if re.search(r"\(id=0\)",s):
        if res_type == 'f':
          values = map(float,[i for i in re.split('\s+',s)[1:] if i != ""])
          return sum(values)/len(values)
        elif res_type == 'i':
          values = map(int,[i for i in re.split('\s+',s)[1:] if i != ""])
          return int(sum(values)/len(values))
        else:
          print 'Error extract_timing: res_type %s not handled'
          return None

  
def extract_metrics_scalasca(args):
  """Returns a list of numpy arrays that contains the data in the file file_to_open. Parameter names and value for this particular experiment are passed as parameters if it can help the file parsing """

  #Check if collective and point 2 point timings are well extracted
  res={}
  metrics_incl = ['time', 'mpi', 'mpi_nxn_completion', 'delay_mpi', 'delay_omp', 'bytes', 'critical_path_imbalance', 'critical_path', 'omp_time']
  for i in metrics_incl: 
    cmd='cube_dump -m %s -x incl -c 0 -z incl %s'%(i, args.trace_file)
    res[i] = extract(cmd)

  metrics_incl_int = ['comms_p2p', 'comms_coll']
  for i in metrics_incl_int: 
    cmd='cube_dump -m %s -x incl -c 0 -z incl %s'%(i, args.trace_file)
    res[i] = extract(cmd, 'i')

  metrics_excl = ['mpi_point2point', 'mpi_collective']
  for i in metrics_excl: 
    cmd='cube_dump -m %s -x excl -c 0 -z incl %s'%(i, args.trace_file)
    res[i] = extract(cmd)

  dump_json(res, args.o)
  
def extract_right_line(filename):
  f = open(filename, 'r')
  l = f.readlines()
  f.close()
  pattern = re.compile('IdrisMemMPI')
  for i in reversed(l):
    res = pattern.search(i)
    if res != None:
      return i
  return None


def extract_time(args):
   f = open(args.stderr_file, 'r')
   l = f.readlines()
   f.close()
   pattern = re.compile('real\s+?(.+?)m(.+?)[.,](.+?)s')
   match = None
   for i in reversed(l):
     match = pattern.search(i)
     if match is not None:
        break
   key = "time_" + args.tag
   res={key: 'Fail'}
   if match is not None:
     res[key] = int(match.group(1))*60 + int(match.group(2)) + int(match.group(3)) / 1000.
   else:
     print 'Warning: No timing information'
   dump_json(res, args.o)

def extract_time2(args):
  l = extract_right_line(args.stderr_file)  
  key = "time_" + args.tag
  res={key: 'Fail'}

  if l != None:
    res[key] = float(l.split(',')[0].split(':')[1].split('=')[1].split(' ')[1])
  else:
    print 'Warning: No timing information'
  dump_json(res, args.o)
  
def extract_mem(args):
  l = extract_right_line(args.stderr_file)  
  key = "mem"
  res={key: 'Fail'}

  if l != None:
    to_return = l.split(',')[2].split('=')[1].split(' ') 
    to_return = to_return[1] + ' ' + to_return[2]
    res[key] = to_return
  else:
    print 'Warning: No memory footprint information'
  dump_json(res, args.o)

def extract_metrics_darshan(args):
  f = open(args.darshan_metric_file, 'r')
  l = f.readlines()
  f.close()
  pattern = {"total_mb_read" : (re.compile('Total MB read: (.+?)$'),float),
             "total_mb_written" : (re.compile('Total MB written: (.+?)$'),float),
             "read_calls" : (re.compile('posix calls: (.+?) '),int),
             "write_calls" : (re.compile('posix calls: (?:.+? ){1}(.+?) '),int),
             "open_calls" : (re.compile('posix calls: (?:.+? ){2}(.+?) '),int),
             "stat_calls" : (re.compile('posix calls: (?:.+? ){3}(.+?) '),int),
             "seek_calls" : (re.compile('posix calls: (?:.+? ){4}(.+?) '),int),
             "mmap_calls" : (re.compile('posix calls: (?:.+? ){5}(.+?) '),int),
             "fsync_calls" : (re.compile('posix calls: (?:.+? ){6}(.+?)$'),int),
             "ind_read_time" : (re.compile('cumul_avg_io_time: (.+?) '),float),
             "ind_write_time" : (re.compile('cumul_avg_io_time: (?:.+? ){1}(.+?) '),float),
             "ind_meta_time" : (re.compile('cumul_avg_io_time: (?:.+? ){2}(.+?) '),float),
             "sh_read_time" : (re.compile('cumul_avg_io_time: (?:.+? ){3}(.+?) '),float),
             "sh_write_time" : (re.compile('cumul_avg_io_time: (?:.+? ){4}(.+?) '),float),
             "sh_meta_time" : (re.compile('cumul_avg_io_time: (?:.+? ){5}(.+?) '),float),
             "io_time" : (re.compile('io_time: (.+?)$'),float),
             "avg_io_ac_size" : (re.compile('avg_io_ac_size: (.+?)$'),float)
             }

  res = dict()
  for i in l:
    for pat in pattern:
        match = pattern[pat][0].search(i)
        if match is not None and pat not in res:
            res[pat] = pattern[pat][1](match.group(1))
  if len(res) == 0:
    print 'Warning: No darshan information'
  dump_json(res, args.o)

def dump_json(dic, file_name):
  to_write = json.dumps(dic)
  f = open(file_name, "w")
  f.write(to_write)
  f.close()

def build_parser():
  # create the top-level parser
  parser = argparse.ArgumentParser(prog='extract_metrics')
  parser.add_argument('-o', type=str, default="metrics.json", help='name of the output file (default metrics.json)')
  subparsers = parser.add_subparsers(help='sub-command help')
 
  # create the parser for the "scalasca" command
  parser_sca = subparsers.add_parser('scalasca', help='Extracts metrics from the trace file from scoreP')
  parser_sca.add_argument('trace_file', type=str, help='ScoreP trace file')
  parser_sca.set_defaults(func=extract_metrics_scalasca)
 
  # create the parser for the "scalasca" command 
  parser_sca = subparsers.add_parser('darshan', help='Extracts metrics from the darshan metric file')
  parser_sca.add_argument('darshan_metric_file', type=str, help='Darshan metric file')
  parser_sca.set_defaults(func=extract_metrics_darshan)

  # create the parser for the "time" command
  parser_time = subparsers.add_parser('time', help='Extracts time from the std.err file')
  parser_time.add_argument('stderr_file', type=str, help='stderr file')
  parser_time.add_argument('-tag', type=str, help='Tag added to a metric name to differentiate the context it has been obtained')
  parser_time.set_defaults(func=extract_time)
  
  # create the parser for the "time" information given by libmem  command
  parser_time = subparsers.add_parser('time2', help='Extracts time from the std.err file')
  parser_time.add_argument('stderr_file', type=str, help='stderr file')
  parser_time.add_argument('-tag', type=str, help='Tag added to a metric name to differentiate the context it has been obtained')
  parser_time.set_defaults(func=extract_time2)
  
  # create the parser for the "mem" command
  parser_mem = subparsers.add_parser('mem', help='Extracts memory ffotprint from the std.err file')
  parser_mem.add_argument('stderr_file', type=str, help='stderr file')
  parser_mem.set_defaults(func=extract_mem)

  return parser


if __name__ == '__main__':
  p =  build_parser()
  args = p.parse_args(sys.argv[1:])
  args.func(args)

