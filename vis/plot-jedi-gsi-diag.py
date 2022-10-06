import getopt
import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

#from dateutil.rrule import *
#from dateutil.parser import *
from datetime import *
from datetime import timedelta

import netCDF4 as nc4

#=========================================================================
def hour2date(hour):
  st = datetime(1970, 1, 1, 0, 0, 0)
  dt = timedelta(hours=hour)
  ct = st + dt

  date = ct.strftime("%Y-%m-%d:%H")

 #print('st = ', st)
 #print('ct = ', ct)
 #print('dt = ', dt)
 #print('date = ', date)

  return date

#=========================================================================
def read_stats(filename):
  fin = open(filename,'r')
  lines = fin.readlines()
  fin.close()

  stats = {}

 #fin.write('# %s %s %s %s-%s\n' % (datapath,self.runid,self.hem,self.sdate,self.edate))
  item = lines[0].strip().split(' ')
 #print('item = ', item)
  stats['datapath'] = item[1]
  stats['runid'] = item[2]
  stats['hem'] = item[3]
  dateitem = item[4].split('-')
  stats['sdate'] = dateitem[0]
  stats['edate'] = dateitem[1]

 #fin.write('# press wind_count wind_rms temp_count temp_rms temp_bias humid_rmsx1000 humid_biasx1000\n')
  item = lines[1].strip().split(' ')
  stats['stats_name'] = item[1:]
 #print("stats['stats_name'] = ", stats['stats_name'])

 #fin.write('# 1000-0 %10i %7.4f %10i %7.4f %7.4f %10i %7.4f %7.4f\n' % (count_wind.sum(), rms_wind.mean(), count_temp.sum(), rms_temp.mean(), bias_temp.mean(), count_humid.sum(), 1000*rms_humid.mean(), 1000*bias_humid.mean()))
  sumline = lines[2]
  while(sumline.find('  ') >= 0):
    sumline = sumline.replace('  ', ' ')
  item = sumline.strip().split(' ')
 #print('item = ', item)
  stats['count_wind_sum'] = float(item[2])
  stats['rms_wind_mean'] = float(item[3])
  stats['count_temp_sum'] = float(item[4])
  stats['rms_temp_mean'] = float(item[5])
  stats['bias_temp_mean'] = float(item[6])
  stats['count_humid_sum'] = float(item[7])
  stats['rms_humid_mean'] = float(item[8])
  stats['bias_humid_mean'] = float(item[9])

  stats['p'] = []
  stats['count_wind'] = []
  stats['rms_wind'] = []
  stats['count_temp'] = []
  stats['rms_temp'] = []
  stats['bias_temp'] = []
  stats['count_humid'] = []
  stats['rms_humid'] = []
  stats['bias_humid'] = []
  for n in range(3, len(lines)):
    line = lines[n].strip()
    while(line.find('  ') >= 0):
      line = line.replace('  ', ' ')
    item = line.strip().split(' ')
    if(item[0] != '#'):
      stats['p'].append(float(item[0]))
      stats['count_wind'].append(float(item[1]))
      stats['rms_wind'].append(float(item[2]))
      stats['count_temp'].append(float(item[3]))
      stats['rms_temp'].append(float(item[4]))
      stats['bias_temp'].append(float(item[5]))
      stats['count_humid'].append(float(item[6]))
      stats['rms_humid'].append(float(item[7]))
      stats['bias_humid'].append(float(item[8]))
    else:
      stats['rms_wind_mean'] = float(item[1])
      stats['rms_temp_mean'] = float(item[2])
      stats['rms_humid_mean'] = float(item[3])

 #print('stats = ', stats)
  return stats

#=========================================================================
def plot_lines(plevs, gsirms, jedirms, header='temp', output=0):
  try:
    plt.close('all')
    plt.clf()
    plt.cla()
  except Exception:
    pass

  title = header + ' RMS'

 #print('plevs = ', plevs)
 #print('gsirms = ', gsirms)
 #print('jedirms = ', jedirms)

  print('\n')
  print('Variable name: ', header)
  print('pressure   gsirms    jedirms')
  k = len(plevs)
  for n in range(len(plevs)):
    k -= 1
    print('%4.0f %7.4f %7.4f' %(plevs[k], gsirms[k], jedirms[k]))

  pmin = 0.0
 #pmin = np.min(jedirms)
 #gmin = np.min(gsirms)
 #if(pmin > gmin):
 #  pmin = gmin

  pmax = np.max(jedirms)
  gmax = np.max(gsirms)
  if(pmax < gmax):
    pmax = gmax
  vmax = 0.0
 
  x = [0]
  xlabels = ['0']
  k = 0
  while(vmax < pmax):
    k += 1
    vmax += 1.0
    lbl = '%d' %(k)
    xlabels.append(lbl)
    x.append(k)

  y = np.linspace(0,1000,11)
  ylabels = []
  for v in y:
    lbl = '%d' %(1000-v)
    ylabels.append(lbl)

  fig = plt.figure()
  ax = plt.subplot()

  k = 0
  yp = []
  xd = []
  for pres in plevs:
    yp.append(1000-pres)
    jmg = jedirms[k] - gsirms[k]
    xd.append(jmg)
    k += 1

  ax.plot(gsirms, yp, color='blue', linewidth=2, alpha=0.9)
  ax.plot(jedirms, yp, color='red', linewidth=2, alpha=0.9)

  plt.xscale('linear')
 #plt.xscale('log', base=2)
 #plt.yscale('log', base=2)
 #plt.yscale('log', base=10)
  plt.xticks(x, xlabels)
 #plt.xticks(x, xlabels, rotation ='vertical')
  plt.yticks(y, ylabels)

 #ax.plot(xd, yp, color='cyan', linewidth=2, alpha=0.9)

  plt.grid()

 #Same limits for everybody!
 #print('pmin: %f, pmax: %f' %(pmin, pmax))
  plt.xlim(0, vmax)
  plt.ylim(0, 1000)
 
 #general title
  title = 'GSI and JEDI %s rms' %(header)
  plt.suptitle(title, fontsize=16, fontweight=1, color='black')

 #Create a big subplot
  bs = fig.add_subplot(111, frameon=False)
  plt.subplots_adjust(bottom=0.2, right=0.70, top=0.8)

 #hide tick and tick label of the big axes
  plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

  bs.set_xlabel('(m/s)', labelpad=10)
  bs.set_ylabel('Pressure (hPa)', labelpad=20)

 #labels = ['GSI', 'JEDI', 'JEDI - GSI']
  labels = ['GSI', 'JEDI']

 #Create the legend
  fig.legend(ax, labels=labels,
         loc="center right",   # Position of legend
         fontsize=8,
         borderpad=1.2,
         labelspacing=1.2,
         handlelength=1.5
         )

 #Adjust the scaling factor to fit your legend text completely outside the plot
 #(smaller value results in more space being made for the legend)

  imgname = '%s_rms.png' %(header)

  if(output):
    plt.savefig(imgname)
  else:
    plt.show()

#=========================================================================
class Plot_JEDI_GSI_Diag():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

  def plot_contour(self, times, plevs, data=[]):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    nrows = 1
    ncols = len(data)

    y = []
    for pres in plevs:
      y.append(1000-pres)

    ylp = np.linspace(0,1000,11)
    yt = []
    ylabels = []
    for pv in ylp:
      lbl = '%d' %(1000-pv)
      yt.append(pv)
      ylabels.append(lbl)

    x = (times - times[0])/24.0
    xv = 0
    xt = [0]
    xlabels = [0]
    while(xv < x[-1]):
      xv += 1
      lbl = '%d' %(xv)
      xlabels.append(lbl)
      xt.append(xv)

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            figsize=(11,8.5))
 
   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
     #axs[i].set_global()

      pvar = data[i]

     #print('Plot No. ', i)
     #print('\tpvar.shape = ', pvar.shape)

      vmin = np.min(pvar)
      vmax = np.max(pvar)

      if(i < 2):
       #self.cmapname = 'coolwarm'
       #self.cmapname = 'rainbow'
        self.cmapname = 'jet'
        self.clevs = np.arange(0.0, 10.1, 0.1)
        self.cblevs = np.arange(0.0, 11.0, 1.0)
      else:
        self.cmapname = 'bwr'
       #self.clevs = np.arange(-2.0, 2.1, 0.1)
       #self.cblevs = np.arange(-2.0, 3.0, 1.0)
        self.clevs = np.arange(-0.5, 0.52, 0.02)
        self.cblevs = np.arange(-0.5, 0.6, 0.2)

     #if((vmax - vmin) > 1.0e-5):
     #  self.clevs, self.cblevs = get_plot_levels(pvar)

      xm, ym = np.meshgrid(x, y)
      cs=axs[i].contourf(xm, ym, pvar,
                         levels=self.clevs, extend=self.extend,
                         alpha=self.alpha, cmap=self.cmapname)

      title = '%s min: %5.2f, max: %5.2f' %(self.runname[i], vmin, vmax)
      axs[i].set_title(title)

     #axs[i].set_xlabel('Hours', fontsize=10)
      axs[i].set_xlabel('Days', fontsize=10)
      if(i < 1):
        axs[i].set_ylabel('hPa', fontsize=10)

      axs[i].set_xticks(xt, xlabels)
      axs[i].set_yticks(yt, ylabels)

     #minor_ticks_top=np.linspace(0,60,13)
     #minor_ticks_top=np.linspace(-60,0,13)
     #axs[i].set_yticks(minor_ticks_top,minor=True)

     #Draw the colorbar
      cbar=plt.colorbar(cs, ax=axs[i], pad=self.pad,
                        ticks=self.cblevs,
                        orientation='horizontal')

     #cbar.set_label(self.label, rotation=90)
      cbar.set_label(self.label)

   #Add a big title at the top
    hour = float(times[0])
    date = hour2date(hour)
    title = '%s start from %s' %(self.title, date)
    plt.suptitle(title)

    fig.canvas.draw()
   #plt.tight_layout()

    if(self.output):
      if(self.imagename is None):
        imagename = 'sample.png'
      else:
        imagename = self.imagename
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

  def set_default(self):
    self.imagename = 'sample.png'

    self.runname = ['GSI', 'JEDI', 'JEDI - GSI']

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
   #self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
    self.cmapname = 'jet'

    self.clevs = np.arange(-0.2, 0.21, 0.01)
    self.cblevs = np.arange(-0.2, 0.3, 0.1)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Unit (C)'
    self.title = 'Temperature Increment'

  def set_label(self, label='Unit (C)'):
    self.label = label

  def set_title(self, title='Temperature Increment'):
    self.title = title

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_imagename(self, imagename):
    self.imagename = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

#--------------------------------------------------------------------------------
def get_plot_levels(var):
  vmin = np.min(var)
  vmax = np.max(var)

 #print('\tget_plot_lelves: vmin = %f, vmax = %f' %(vmin, vmax))

  vcen = 0.5*(vmax + vmin)

  if(vcen >= 1.0):
    vcen = 10.0*int(vcen/10.0)
    pmax = 10.0*int(0.5 + (vmax-vcen)/10.0) + vcen
    pmin = vcen - 10.0*int(0.5 + (vmax-vcen)/10.0)
  else:
    pcen = vcen
    fact = 1.0
    while(pcen < 1.0):
      fact *= 10.0
      pcen *= 10.0
    vcen = 10.0*int(pcen/10.0)/fact
    pmax = 10.0*int(0.5 + fact*(vmax-vcen)/10.0)/fact + vcen
    pmin = vcen - 10.0*int(0.5 + fact*(vmax-vcen)/10.0)/fact

 #print('\tget_plot_lelves: pmin = %f, pcen = %f, pmax = %f' %(pmin, vcen, pmax))

  delt = (pmax - pmin)/200.0
  clevs = np.arange(pmin, pmax + delt, delt)

  delt = (pmax - pmin)/20.0
  cblevs = np.arange(pmin, pmax + delt, delt)

 #print('\tclevs: ', clevs[::10])
 #print('\tcblevs: ', cblevs)

  return clevs, cblevs

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  type = 'bfg'
  gsifile = 'gsi_stats.nc'
  jedifile = 'jedi_stats.nc'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'type=',
                                                'jedifile=', 'gsifile='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--jedifile'):
      jedifile = a
    elif o in ('--type'):
      type = a
    elif o in ('--gsifile'):
      gsifile = a
    else:
      assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
  pjgd = Plot_JEDI_GSI_Diag(debug=debug, output=output)

  gsi_stats = read_stats('gsi_stats')
  jedi_stats = read_stats('jedi_stats')

 #stats['p'].append(float(item[0]))
 #stats['count_wind'].append(float(item[1]))
 #stats['rms_wind'].append(float(item[2]))
 #stats['count_temp'].append(float(item[3]))
 #stats['rms_temp'].append(float(item[4]))
 #stats['bias_temp'].append(float(item[5]))
 #stats['count_humid'].append(float(item[6]))
 #stats['rms_humid'].append(float(item[7]))
 #stats['bias_humid'].append(float(item[8]))
  p = gsi_stats['p']
  gsirms = gsi_stats['rms_temp']
  jedirms = jedi_stats['rms_temp']

 #print('len(p) = ', len(p))
 #print('p = ', p)
 #print('len(gsirms) = ', len(gsirms))
 #print('gsirms = ', gsirms)

  plot_lines(p, gsirms, jedirms, header='temp', output=output)

  gsirms = gsi_stats['rms_wind']
  jedirms = jedi_stats['rms_wind']

  plot_lines(p, gsirms, jedirms, header='wind', output=output)

  gsirms = gsi_stats['rms_humid']
  jedirms = jedi_stats['rms_humid']

  for n in range(len(jedirms)):
    if(jedirms[n] < 0.0):
      jedirms[n] = 0.0
    if(gsirms[n] < 0.0):
      gsirms[n] = 0.0

  plot_lines(p, gsirms, jedirms, header='humidity', output=output)

#-----------------------------------------------------------------------------------------
  ncgsi = nc4.Dataset(gsifile, 'r')
  ncjedi = nc4.Dataset(jedifile, 'r')

  times = ncgsi.variables['times']
  plevs = ncgsi.variables['plevs']

#-----------------------------------------------------------------------------------------
  clevs = np.arange(-1.0, 1.01, 0.01)
  cblevs = np.arange(-1.0, 1.1, 0.1)

  pjgd.set_clevs(clevs=clevs)
  pjgd.set_cblevs(cblevs=cblevs)

#-----------------------------------------------------------------------------------------
  varlist = ['omf_rmswind', 'omf_rmstemp', 'omf_rmshumid']
  unitlist = ['Unit (m/s)', 'Unit (C)', 'Unit (g/kg)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist)):
    jedivar = ncjedi.variables[varlist[n]][:, :]
    gsivar = ncgsi.variables[varlist[n]][:, :]

    if('omf_rmshumid' == varlist[n]):
      jedivar = np.where(jedivar < 0.0, 0.0, 1.0e7*jedivar)
      gsivar = np.where(gsivar < 0.0, 0.0, 1.0e7*gsivar)

    ntime, nlev = jedivar.shape
   #print('jedivar.shape = ', jedivar.shape)
   #print('gsivar.shape = ', gsivar.shape)

    pjgd.set_label(unitlist[n])

    v0 = gsivar.transpose()
    v1 = jedivar.transpose()
   #v0 = gsivar
   #v1 = jedivar
    v2 = v1 - v0

    data = [v0, v1, v2]

    title = varlist[n]
    pjgd.set_title(title)

    print('Plotting ', title)
   #print('\tv0.shape = ', v0.shape)
   #print('\tv1.shape = ', v1.shape)
   #print('\tv2.shape = ', v2.shape)
 
   #print('\tv0.max: %f, v0.min: %f' %(np.max(v0), np.min(v0)))
   #print('\tv1.max: %f, v1.min: %f' %(np.max(v1), np.min(v1)))
   #print('\tv2.max: %f, v2.min: %f' %(np.max(v2), np.min(v2)))

    imagename = 'diag_%s.png' %(varlist[n])
    pjgd.set_imagename(imagename)

    pjgd.plot_contour(times, plevs, data=data)

#-----------------------------------------------------------------------------------------
  ncjedi.close()
  ncgsi.close()

