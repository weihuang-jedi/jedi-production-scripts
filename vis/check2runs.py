import getopt
import os, sys
import subprocess

import matplotlib.pyplot as plt
import numpy as np

#import time
#import datetime
from datetime import *
from datetime import timedelta

#=========================================================================
def hour2date(hour, sdate=None):
  if(sdate is None):
    st = datetime(2020, 1, 1, 0, 0, 0)
  else:
    st = sdate
  dt = timedelta(hours=hour)
  ct = st + dt

  date = ct.strftime("%Y%m%d%H")

 #print('st = ', st)
 #print('ct = ', ct)
 #print('dt = ', dt)
 #print('date = ', date)

  return date

#########################################################################
""" TimeSeries """
class TimeSeries:
  """ Constructor """
  def __init__(self, debug=0, casename=None, output=0,
               firstdir=None, seconddir=None, interval=6, linear=0):

    """ Initialize class attributes """
    self.debug = debug
    self.firstdir = firstdir
    self.seconddir = seconddir
    self.casename = casename
    self.interval = interval
    self.output = output
    self.linear = linear

    if(firstdir is None):
      print('firstdir not defined. Exit.')
      sys.exit(-1)
    if(seconddir is None):
      print('seconddir not defined. Exit.')
      sys.exit(-1)

    self.colorlist1 = ['red', 'green', 'blue', 'magenta', 'cyan', 'skyblue']
    self.colorlist2 = ['firebrick', 'lime', 'navy', 'purple', 'forestgreen', 'steelblue']

    self.linestyle1 = ['solid', 'solid', 'solid', 'solid', 'solid', 'solid']
    self.linestyle2 = ['dashed', 'dashed', 'dashed', 'dashed', 'dashed', 'dashed']

    self.funclist = ['util::Timers::Total',
                     'oops::GETKFSolver::computeHofX',
                     'oops::GETKFSolver::computeWeights',
                     'oops::GETKFSolver::measurementUpdate',
                     'oops::State::State',
                     'oops::ObsFilter']

   #self.funclist = ['util::Timers::Total',
   #                 'oops::GETKFSolver::computeHofX',
   #                 'oops::GETKFSolver::computeWeights',
   #                 'oops::GETKFSolver::measurementUpdate',
   #                 'oops::State::State',
   #                 'oops::GetValues',
   #                 'oops::ObsError',
   #                 'oops::ObsFilter',
   #                 'oops::ObsSpace',
   #                 'oops::ObsVector',
   #                 'oops::ObsOperator',
   #                 'oops::VariableChange']

    self.stats1 = []
    self.times1 = []

    self.stats2 = []
    self.times2 = []
      
  def set_linear(self, linear=1):
    self.linear = linear

  def set_output(self, output=1):
    self.output = output

  def process(self):
    self.stats1, self.times1 = self.process_dir(self.firstdir)
    self.stats2, self.times2 = self.process_dir(self.seconddir)

    self.ntime = len(self.times1)
    if(self.ntime > len(self.times2)):
      self.ntime = len(self.times2)

  def process_dir(self, workdir):
    hour = 12
    search_more = True
    nf = 0
    timelist = []
    filelist = []
    while(search_more):
      rdate = hour2date(hour)
      flnm = '%s/%s/stdoutNerr/stdout.00000000' %(workdir, rdate)

      if(os.path.exists(flnm)):
        if(self.debug):
          print('Processing No %d file: %s' %(nf, flnm))

        timelist.append(rdate)
        filelist.append(flnm)

        hour += self.interval
        nf += 1
      else:
        print('file: %s does not exist.' %(flnm))
        search_more = False
        break

    nf -= 1
    slist = []
    tlist = []
    for n in range(nf):
      flnm = filelist[n]
      available, stats = self.process_file(flnm)
      if(available):
        slist.append(stats)
        tlist.append(timelist[n])

    return slist, tlist
  
  def process_file(self, flnm):
    available = False
    stats = {}
    with open(flnm) as fp:
      lines = fp.readlines()
     #line = fp.readline()
      num_lines = len(lines)
     #print('Total number of lines: ', num_lines)

      nl = 0
      while(nl < num_lines):
        if(lines[nl].find('Parallel Timing Statistics') > 0):
         #if(self.debug):
         #  print('lines[%d]: %s' %(nl, lines[nl]))
          nl, stats = self.parallel_time_stats(lines, nl)
          nl, wcm = self.wall_clock_time(lines, nl)
          stats['Runtime'] = wcm
          available = True
          break
        nl += 1

    return available, stats
  
  def wall_clock_time(self, lines, nl):
    going = 1
    wcm = {}
    ns = nl
    while(going):
      line = lines[ns].strip()
     #print('Line ' + str(ns) + ': ' + line)
      ns += 1
      if(line.find('OOPS_STATS Run end') >= 0):
        idx_Runtime = line.find('Runtime:')
        idx_total = line.find('Memory: total:')
        idx_min = line.find('min =')
        idx_max = line.find('max =')
        maxmemstr = line[idx_max+5:].strip()
        minmemstr = line[idx_min+5:idx_max].strip()
        ttlmemstr = line[idx_total+14:idx_min].strip()
        rtstr = line[idx_Runtime+8:idx_total].strip()

       #print('maxmemstr:', maxmemstr)
        maxitem = maxmemstr.split(' ')
        wcm['max'] = float(maxitem[0])
        wcm['max unit'] = maxitem[1][0:2]

       #print('minmemstr:', minmemstr)
        minitem = minmemstr.split(' ')
        wcm['min'] = float(minitem[0])
        wcm['min unit'] = minitem[1][0:2]

       #print('ttlmemstr:', ttlmemstr)
        totalitem = ttlmemstr.split(' ')
        wcm['total'] = float(totalitem[0])
        wcm['total unit'] = totalitem[1][0:2]

       #print('rtstr:', rtstr)
        runtimeitem = rtstr.split(' ')
        wcm['runtime'] = float(runtimeitem[0])
        wcm['runtime unit'] = runtimeitem[1][0:3]
        
        going = 0
        break

   #print('wcm: ', wcm)

    return ns, wcm

  def parallel_time_stats(self, lines, nl):
    going = 1
    ns = nl + 3
    stats = {}
    for name in self.funclist:
      stats[name] = {}
      stats[name]['min'] = 0.0
      stats[name]['max'] = 0.0
      stats[name]['avg'] = 0.0

    while(going):
      line = lines[ns].strip()
     #if(self.debug):
     #  print('lines[%d]: %s' %(ns, line))

      ns += 1
      if(line.find('Parallel Timing Statistics') > 0):
        going = 0
        break

      item = line.split(': ')
     #print(item)
      nlist = item[0].strip().split(' ')
      name = nlist[1]

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      nlist = tstr.split(' ')
      tmin = float(nlist[0])
      tmax = float(nlist[1])
      tavg = float(nlist[2])

      for i in range(len(self.funclist)):
        if(name.find(self.funclist[i]) >= 0):
          stats[self.funclist[i]]['min'] += tmin
          stats[self.funclist[i]]['max'] += tmax
          stats[self.funclist[i]]['avg'] += tavg
      
   #print('stats: ', stats)

    return ns, stats

  def plot(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = '%s Timing' %(self.casename)

    n = self.ntime
    x = np.zeros((n), dtype=float)
    y = np.zeros((n), dtype=float)
    z = np.zeros((n), dtype=float)
    xlabels = []
    for k in range(n):
      x[k] = float(k*self.interval)/24.0
    k = 0
    xp = []
   #print('self.times = ', self.times1)
    while(k < n):
      lbl = self.times1[int(k)][0:8]
      xlabels.append(lbl)
      xp.append(float(k*self.interval)/24.0)
      k += 24/self.interval

    pmin = 1.0/16.0
    pmax = 256.0
    ylabels = []
    yp = []
    pcur = pmin/2.0
    while(pcur < pmax):
      pcur *= 4.0
      lbl = '%6.2f' %(pcur)
      ylabels.append(lbl)
      yp.append(pcur)
     #print('yp = ', yp)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.yscale('linear')
    else:
     #plt.yscale('log', base=10)
      plt.yscale('log', base=2)
     #plt.yscale('log', basey=2)
      plt.yticks(yp, ylabels)

    plt.xscale('linear')
   #plt.xticks(xp, xlabels)
   #plt.xticks(xp, xlabels, rotation ='vertical')
    plt.xticks(xp, xlabels, rotation = 45)

    vmin = 1.0e20
    vmax = 0.0

    for i in range(len(self.funclist)):
      name = self.funclist[i]
      for k in range(n):
        y[k] = 0.001*self.stats1[k][name]['max']/60.0
        if(vmin > y[k]):
          vmin = y[k]
        if(vmax < y[k]):
          vmax = y[k]
        if(y[k] < pmin):
          y[k] = pmin
       #ymin = 0.001*self.stats1[k][name]['min']/60.0
       #yavg = 0.001*self.stats1[k][name]['avg']/60.0
      ax.plot(x, y, color=self.colorlist1[i], linewidth=1,
              linestyle=self.linestyle1[i], alpha=0.9)

      for k in range(n):
        y[k] = 0.001*self.stats2[k][name]['max']/60.0
        if(vmin > y[k]):
          vmin = y[k]
        if(vmax < y[k]):
          vmax = y[k]
        if(y[k] < pmin):
          y[k] = pmin
       #ymin = 0.001*self.stats2[k][name]['min']/60.0
       #yavg = 0.001*self.stats2[k][name]['avg']/60.0
      ax.plot(x, y, color=self.colorlist2[i], linewidth=2,
              linestyle=self.linestyle2[i], alpha=0.9)

    plt.grid()

    print('vmin: %f, vmax: %f' %(vmin, vmax))

    pmin = 1.0/16.0
    pmax = 256.0
   
    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
   #title = '%s Timing (in minutes), min: %8.2f, max: %8.2f' %(self.casename, vmin, vmax)
    title = '%s Timing (in minutes)' %(self.casename)
   #plt.suptitle(title, fontsize=13, fontweight=0, color='black', style='italic', y=1.02)
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
   #plt.subplots_adjust(bottom=0.2, right=0.70, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.675, top=0.8)
    plt.subplots_adjust(bottom=0.2, right=0.65, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.5, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

   #bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (minutes)', labelpad=20)

    lblist = self.funclist
   #for func in self.funclist:
   #  lblist.append(func)

   #Create the legend
    fig.legend(ax, labels=lblist, 
           loc="center right",   # Position of legend
           fontsize=6,
           borderpad=1.2,
           handlelength=1.5)

   #Adjust the scaling factor to fit your legend text completely outside the plot
   #(smaller value results in more space being made for the legend)

    if(self.linear):
      imgname = 'lin_timing_%s.png' %(self.casename)
    else:
      imgname = 'log_timing_%s.png' %(self.casename)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  casename = 'GSI-QC-Filter'
  firstdir = '/work2/noaa/da/weihuang/cycling/jedi_C96_lgetkf_sondesonly'
  seconddir = '/work2/noaa/gsienkf/weihuang/gsi/jedi_C96_lgetkf_sondesonly'
  output = 0
  linear = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'firstdir=',
                             'seconddir=', 'output=', 'casename=', 'linear='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--firstdir'):
      firstdir = a
    elif o in ('--seconddir'):
      seconddir = a
    elif o in ('--casename'):
      casename = a
    elif o in ('--output'):
      output = int(a)
    elif o in ('--linear'):
      linear = int(a)
    else:
      assert False, 'unhandled option'

  pr = TimeSeries(debug=debug, output=output, firstdir=firstdir,
                  seconddir=seconddir, casename=casename, linear=linear)
  pr.process()
  pr.plot()

