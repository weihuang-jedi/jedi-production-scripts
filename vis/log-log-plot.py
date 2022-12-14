#########################################################################
#$Id: bld.py 28 2021-01-21 15:10:31Z whuang $
#$Revision: 28 $
#$HeadURL: file:///Users/whuang/.wei_svn_repository/trunk/jedi-build-tools/bld.py $
#$Date: 2021-01-21 08:10:31 -0700 (Thu, 21 Jan 2021) $
#$Author: whuang $
#########################################################################

import getopt
import os, sys
import subprocess
import time
import datetime

import matplotlib.pyplot as plt
import numpy as np

#--------------------------------------------------------------------------------
class Profiler:
  """ Constructor """
  def __init__(self, debug=0, output=0, linear=0):

    """ Initialize class attributes """
    self.debug = debug
    self.output = output
    self.linear = linear

    self.nodelist = [1, 2, 4, 8]
    self.corelist = 36 * self.nodelist

    self.shortest_time = 0.1

    self.colorlist = ['red', 'green', 'cyan', 'blue', 'magenta',
                      'firebrick', 'lime']
    self.linestyle = ['solid', 'solid', 'solid', 'solid', 'solid',
                      'dashed', 'dashed']

    self.halo_develop = [10038.5, 9454.27, 8080.81, 6975.5]
    self.halo_anna1 = [8560.52, 7829.97, 6141.86, 6653.33]
    self.halo_anna2 = [8102.79, 7121.7, 6283.36, 6531.6]

    self.RR_develop = [7277.24, 11000.0, 4937.51, 4623.36]
    self.RR_anna1 = [11000.0, 11000.0, 5462.4, 7104.29]
    self.RR_anna2 = [11000.0, 11000.0, 5111.65, 5538.25]

    self.RR_develop_observer = [2350.22, 2069.85, 1758.18, 1555.16]
    self.RR_anna1_observer = [1882.66, 1927.42, 2394.42, 4210.54]
    self.RR_anna2_observer = [1592.92, 1389.98, 1245.70, 1180.93]

    self.Halo_develop_solver = [4927.02, 11000.0, 3179.33, 3068.20]
    self.Halo_anna1_solver = [11000.0, 11000.0, 3067.98, 2893.75]
    self.Halo_anna2_solver = [11000.0, 11000.0, 3865.95, 4357.32]

    self.namelist = ['halo_develop', 'halo_anna1', 'halo_anna2',
                     'RR_develop', 'RR_anna1', 'RR_anna2']

#--------------------------------------------------------------------------------
  def set_linear(self, linear=1):
    self.linear = linear

  def set_output(self, output=1):
    self.output = output

  def set_minmax(self, min=0.0, max=1.0):
    self.min = min
    self.max = max

#--------------------------------------------------------------------------------
  def plot(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = 'Halo Distr vs RR Observer + Halo Solver Timing'

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

    if(self.linear == 0):
      pmin = 64
      pmax = 256

      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur *= 2
    else:
      pmin = 60
      pmax = 180
      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur += 10
   #print('yp = ', yp)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.xscale('log', base=2)
      plt.yscale('linear')
    else:
     #plt.xscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)
     #plt.xscale('log', basex=2)
     #plt.yscale('log', basey=2)
      plt.xticks(x, xlabels)
      plt.yticks(yp, ylabels)

    ylist = [self.halo_develop, self.halo_anna1, self.halo_anna2,
             self.RR_develop, self.RR_anna1, self.RR_anna2]

    for i in range(len(ylist)):
      n = 0
      for v in ylist[i]:
        y[n] = v/60.0
        n += 1
      print('i = ', i)
      print('y = ', y)
      ib = 0
      if(i > 2):
        ib = 2
      ax.plot(x[ib:], y[ib:], color=self.colorlist[i], linewidth=2,
              linestyle=self.linestyle[i], alpha=0.9)

    plt.grid()

    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.2, right=0.65, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.5, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (minutes)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=self.namelist,
           loc="center right",   # Position of legend
           fontsize=6,
           borderpad=1.2,
           handlelength=1.5)

    if(self.linear):
      imgname = 'lin_timing_all.png'
    else:
      imgname = 'log_timing_all.png'

    plt.savefig(imgname)
    plt.show()

#--------------------------------------------------------------------------------
  def plot_halo_all(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = 'Halo All Timing'

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

    if(self.linear == 0):
      pmin = 64
      pmax = 256

      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur *= 2
    else:
      pmin = 60
      pmax = 180
      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur += 10
   #print('yp = ', yp)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.xscale('log', base=2)
      plt.yscale('linear')
    else:
     #plt.xscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)
     #plt.xscale('log', basex=2)
     #plt.yscale('log', basey=2)
      plt.xticks(x, xlabels)
      plt.yticks(yp, ylabels)

    ylist = [self.halo_develop, self.halo_anna1, self.halo_anna2]

    ib = 0
    for i in range(len(ylist)):
      n = 0
      for v in ylist[i]:
        y[n] = v/60.0
        n += 1
      print('i = ', i)
      print('y = ', y)
      ax.plot(x[ib:], y[ib:], color=self.colorlist[i], linewidth=2,
              linestyle=self.linestyle[i], alpha=0.9)

    plt.grid()

    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.2, right=0.65, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.5, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (minutes)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=['Develop', 'Anna 1', 'Anna 2'],
           loc="center right",   # Position of legend
           fontsize=6,
           borderpad=1.2,
           handlelength=1.5)

    if(self.linear):
      imgname = 'lin_timing_JEDI_all.png'
    else:
      imgname = 'log_timing_JEDI_all.png'

    plt.savefig(imgname)
    plt.show()

#--------------------------------------------------------------------------------
  def plot_observer(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = 'Observer Timing'

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

    if(self.linear == 0):
      pmin = 8
      pmax = 128

      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur *= 2
    else:
      pmin = 10
      pmax = 120
      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur += 10
   #print('yp = ', yp)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.xscale('log', base=2)
      plt.yscale('linear')
    else:
     #plt.xscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)
     #plt.xscale('log', basex=2)
     #plt.yscale('log', basey=2)
      plt.xticks(x, xlabels)
      plt.yticks(yp, ylabels)

    ylist = [self.RR_develop_observer, self.RR_anna1_observer, self.RR_anna2_observer]

    ib = 0
    for i in range(len(ylist)):
      n = 0
      for v in ylist[i]:
        y[n] = v/60.0
        n += 1
      print('i = ', i)
      print('y = ', y)
      ax.plot(x[ib:], y[ib:], color=self.colorlist[i], linewidth=2,
              linestyle=self.linestyle[i], alpha=0.9)

    plt.grid()

    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.2, right=0.65, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.5, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (minutes)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=['Develop', 'Anna 1', 'Anna 2'],
           loc="center right",   # Position of legend
           fontsize=6,
           borderpad=1.2,
           handlelength=1.5)

    if(self.linear):
      imgname = 'lin_timing_Observer.png'
    else:
      imgname = 'log_timing_Observer.png'

    plt.savefig(imgname)
    plt.show()

#--------------------------------------------------------------------------------
  def plot_solver(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = 'Solver Timing'

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

    if(self.linear == 0):
      pmin = 8
      pmax = 128

      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur *= 2
    else:
      pmin = 10
      pmax = 120
      ylabels = []
      yp = []
      pcur = pmin
      while(pcur <= pmax):
        lbl = '%d' %(pcur)
        ylabels.append(lbl)
        yp.append(pcur)
        pcur += 10
   #print('yp = ', yp)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.xscale('log', base=2)
      plt.yscale('linear')
    else:
     #plt.xscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)
     #plt.xscale('log', basex=2)
     #plt.yscale('log', basey=2)
      plt.xticks(x, xlabels)
      plt.yticks(yp, ylabels)

    ylist = [self.Halo_develop_solver, self.Halo_anna1_solver, self.Halo_anna2_solver]

    ib = 2
    for i in range(len(ylist)):
      n = 0
      for v in ylist[i]:
        y[n] = v/60.0
        n += 1
      ax.plot(x[ib:], y[ib:], color=self.colorlist[i], linewidth=2,
              linestyle=self.linestyle[i], alpha=0.9)

    plt.grid()

    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.2, right=0.65, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.5, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (minutes)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=['Develop', 'Anna 1', 'Anna 2'],
           loc="center right",   # Position of legend
           fontsize=6,
           borderpad=1.2,
           handlelength=1.5)

    if(self.linear):
      imgname = 'lin_timing_Solver.png'
    else:
      imgname = 'log_timing_Solver.png'

    plt.savefig(imgname)
    plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  casename = 'sondes'
  nodelist = [1, 2, 4, 8]
  output = 0
  linear = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'linear='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--linear'):
      linear = int(a)
    else:
      assert False, 'unhandled option'

  pr = Profiler(debug=debug, output=output, linear=linear)
 #pr.set_output(output=output)
  pr.plot()
  pr.plot_halo_all()
  pr.plot_observer()
  pr.plot_solver()

