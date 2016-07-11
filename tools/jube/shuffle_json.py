###################################################################
#                Shuffle_json.py                                  #
#  Created : March 2016                                           #
#                                                                 #
#  Author:                                                        #
#     Matthieu Haefele                                            #
#     matthieu.haefele@maisondelasimulation.fr                    #
#     Maison de la Simulation USR3441                             #
#                                                                 #
#  Manipulates dictionaries stored in json format that            #
#  contains performance metrics and genrates a Tex tabular        #
###################################################################

import sys
import argparse
import json

def load_dic(file_name):
  f = open(file_name, "r")
  d = json.loads(f.read())
  f.close()
  return d

def agg(args):
  metrics={}
  for f_name in args.json_file:
    metrics.update(load_dic(f_name))
  f = open(args.o, "w")
  f.write(json.dumps(metrics))
  f.close()

def write_group(group_name, group, experiment, to_write):
  to_write = ''
  nb_item = len(group)
  first = True
  for m in group:
    if first:
      to_write += "\\multirow{%d}{*}{\\rotatebox{90}{%s}} & %s "%(nb_item, group_name, m)
      first = False
    else:
      to_write += " & %s "%(m)

    for e in experiment:
      if m in e:
        if type(e[m]) == float:
          to_write += '& %.1f '%e[m]
        elif type(e[m]) == int:
          to_write += '& %d '%e[m]
        else:
          to_write += '& %s '%str(e[m])
      else:
        to_write += '& N.A. '
    to_write += '\\\\\n\\cline{2-%d}\n'%(len(experiment)+2)
  return to_write

def keys_here(keys_to_check, dic):
  for k in keys_to_check:
    if k not in dic:
      return False

  return True  

def build_final_metrics(res):
  res2={}
  keys = ['time_ref']
  if keys_here(keys, res):
    res2['Total Time (s)'] = res[keys[0]]

  keys = ['mpi']
  if keys_here(keys, res):
   res2['Time MPI (s)'] = res[keys[0]]

  keys = ['comms_p2p']
  if keys_here(keys, res):
    res2['P2P Calls (nb)'] = res[keys[0]]

  keys = ['mpi_point2point']
  if keys_here(keys, res):
    res2['P2P Calls (s)'] = res[keys[0]] #?

  keys = ['comms_coll']
  if keys_here(keys, res):
    res2['Collective Calls (nb)'] = res[keys[0]]

  keys = ['mpi_collective','mpi_nxn_completion']
  if keys_here(keys, res):
    res2['Collective Calls (s)'] = res[keys[0]] + res[keys[1]] #?

  keys = ['delay_mpi']
  if keys_here(keys, res):
    res2['Synchro / Wait MPI (s)'] = res[keys[0]]

  keys = ['delay_mpi','mpi']
  if keys_here(keys, res):
    if res[keys[1]] > 0.0:
      res2['Ratio Synchro / Wait MPI'] = res[keys[0]] / res[keys[1]] * 100
    else:
      res2['Ratio Synchro / Wait MPI'] = 0.0

  keys = ['bytes','comms_coll','comms_p2p']
  if keys_here(keys, res):
    if res[keys[1]] > 0:
      res2['Message Size (kB)'] = res[keys[0]] / (res[keys[1]]+res[keys[2]]) / 1000.

  keys = ['critical_path_imbalance','critical_path']
  if keys_here(keys, res):
    if res[keys[1]] > 0:
      res2['Load Imbalance MPI'] = res[keys[0]] / (res[keys[1]]) * 100

  keys = ['omp_time','time','mpi','Time IO (s)'] #TODO handle properly IOs
  if keys_here(keys, res):
    if res[keys[1]] > 0:
      res2['Ratio OpenMP'] = res[keys[0]] / (res[keys[1]] - res[keys[2]] - res2[keys[3]]) * 100

  keys = ['delay_mpi','omp_time']
  if keys_here(keys, res):
    if res[keys[1]] > 0.0:
      res2['Ratio Synchro / Wait OpenMP'] = res[keys[0]] / res[keys[1]] * 100
    else:
      res2['Ratio Synchro / Wait OpenMP'] = 0.0

  keys = ["mem"]
  if keys_here(keys, res):
    res2["Memory Footprint (B)"] = res[keys[0]]

  keys = ['total_mb_written', 'total_mb_read']
  if keys_here(keys, res):
    res2["IO Volume (MB)"] = res[keys[0]] + res[keys[1]]

  keys = ["io_time"]
  if keys_here(keys, res):
    res2['Time IO (s)'] = res[keys[0]]

  keys = ['open_calls', 'seek_calls', 'mmap_calls', 'write_calls', 'fsync_calls', 'stat_calls', 'read_calls']
  if keys_here(keys, res):
    res2["Calls (nb)"] = res[keys[0]] + res[keys[1]] + res[keys[2]] + res[keys[3]] + res[keys[4]] + res[keys[5]] + res[keys[6]]

  keys = ["io_time", 'total_mb_written', 'total_mb_read']
  if keys_here(keys, res):
    if res2['Time IO (s)'] > 1e-9:
        res2["Throughput (MB/s)"] = res2["IO Volume (MB)"] / res2['Time IO (s)']
    else:
        res2["Throughput (MB/s)"] = 0.0

  keys = ["avg_io_ac_size"]
  if keys_here(keys, res):
    res2["Individual IO Access (kB)"] = res[keys[0]] / 1024.

  keys = ["time_no-vec"]
  if keys_here(keys, res):
    res2["Runtime without vectorisation (s)"] = res[keys[0]]

  keys = ["time_no-fma"]
  if keys_here(keys, res):
    res2["Runtime without FMA (s)"] = res[keys[0]]

  keys = ["time_no-fma","time_ref"]
  if keys_here(keys, res):
    res2["FMA efficiency"] = res[keys[0]] / res[keys[1]]

  keys = ["time_no-vec","time_ref"]
  if keys_here(keys, res):
    res2["Vectorisation efficiency"] = res[keys[0]] / res[keys[1]]

  keys = ["time_compact","time_scatter"]
  if keys_here(keys, res):
    res2["Memory vs Compute Bound"] = res[keys[0]] / res[keys[1]]
  return res2

def tex(args):
  global_metrics=["Golbal", "Total Time (s)", "Time IO (s)", "Time MPI (s)", "Memory vs Compute Bound"]

  io_metrics=["IO", "IO Volume (MB)", "Calls (nb)", "Throughput (MB/s)", "Individual IO Access (kB)"]
  mpi_mterics=["MPI", "P2P Calls (nb)", "P2P Calls (s)", "Collective Calls (nb)", "Collective Calls (s)", "Synchro / Wait MPI (s)", "Ratio Synchro / Wait MPI", "Message Size (kB)", "Load Imbalance MPI"]

  node_metrics=["Node", "Ratio OpenMP", "Load Imbalance OpenMP", "Ratio Synchro / Wait OpenMP"]

  memory_metrics=["Mem", "Memory Footprint (B)", "Cache Usage Intensity", "RAM Avg Throughput (GB/s)"]
 
  cpu_metrics=["Core", "IPC", "Runtime without vectorisation (s)", "Vectorisation efficiency", "Runtime without FMA (s)", "FMA efficiency"]
  
  to_write = '\\begin{tabular}{|l|c|'
  for i in range(len(args.json_file)):
    to_write += 'c|'
  to_write += '}\n'
  to_write += "\\hline\n & Metric name "
  experiment = []
  for f_name in args.json_file:
    to_write += '& ' + f_name + ' '
    experiment.append(build_final_metrics(load_dic(f_name)))
  to_write += '\\\\\n\\hline\n'

  #to_write += write_group(global_metrics, experiment, to_write)
  for g in [global_metrics, io_metrics, mpi_mterics, node_metrics, memory_metrics, cpu_metrics]:
    to_write += write_group(g[0], g[1:], experiment, to_write)
    to_write += '\n\\hline\n'

  to_write += '\\end{tabular}\n'
  f = open(args.o, "w")
  f.write(to_write)
  f.close()

def build_parser():
  # create the top-level parser
  parser = argparse.ArgumentParser(prog='shuffle_json')
  #parser.add_argument('--foo', action='store_true', help='foo help')
  subparsers = parser.add_subparsers(help='sub-command help')
 
  # create the parser for the "aggregate" command
  parser_agg = subparsers.add_parser('agg', help='aggregate the json files passed in parameter')
  parser_agg.add_argument('-o', type=str, default="out.json", help='name of the output file (default out.json)')
  parser_agg.add_argument('json_file', type=str, nargs='+', help='json files to aggregate')
  parser_agg.set_defaults(func=agg)
 
  # create the parser for the "tex" command
  parser_tex = subparsers.add_parser('tex', help='generates a TeX array, each json file passed in parameter is a column of the array')
  parser_tex.add_argument('-o', type=str, default="out.tex", help='name of the output file (default out.tex)')
  parser_tex.add_argument('json_file', type=str, nargs='+', help='json files to put in the table')
  parser_tex.set_defaults(func=tex)
  return parser


if __name__ == '__main__':
  p =  build_parser()
  args = p.parse_args(sys.argv[1:])
  args.func(args)

