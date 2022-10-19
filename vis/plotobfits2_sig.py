import getopt
#import dateutils, sys
import sys
from netCDF4 import Dataset
from scipy import stats
from datetime import *
from datetime import timedelta

import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

#=========================================================================
def ttest(data1, data2, inflate=False):
    # calculate means
    mean1 = data1.mean(axis=0); mean2 = data2.mean(axis=0)
    # number of paired samples
    n = data1.shape[0]
    # sum squared difference between observations
    d1 = ((data1-data2)**2).sum(axis=0)
    # sum difference between observations
    d2 = (data1-data2).sum(axis=0)
    # standard deviation of the difference between means
    sd = np.sqrt((d1 - (d2**2 / n)) / (n - 1))
    # standard error of the difference between the means
    inflation = 1.0
    if inflate:
        # inflation to represent autocorrelation (see Geer 2016 appendix, Wilks 2006)
        x = data1-data2
        r1 = np.empty(data1.shape[1])
        r2 = np.empty(data1.shape[1])
        for i in range(data1.shape[1]):
            r1[i] = np.corrcoef(x[:-1,i], x[1:,i],rowvar=False)[0,1]
            r2[i] = np.corrcoef(x[:-2,i], x[2:,i],rowvar=False)[0,1]
        #r2 = r1 # AR(1)
        phi1 = r1*(1.-r2)/(1.-r1**2)
        phi2 = (r2-r1**2)/(1.-r1**2)
        rho1  = phi1/(1.-phi2)
        rho2 = phi2 + phi1**2/(1.-phi2)
        inflation = np.sqrt((1.-rho1*phi1-rho2*phi2)/(1.-phi1-phi2)**2)
        inflation = np.where(inflation < 1.0, 1.0, inflation)
    sed = inflation*sd / np.sqrt(n)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # return the p-values
    return 1.-(1.-stats.t.cdf(abs(t_stat), n-1)) * 2.0 # two sided

#=========================================================================
def ttest2(data1,data2):
    t, p = stats.ttest_rel(data1,data2)
    return 1.-p

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
class Plot_JEDI_GSI_Diag():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()
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

    self.levtop = 125
    self.levbot = 925 # plot limits
    self.sigthresh = 0.99       # significance threshold (p-value)

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

  def get_stats(self, expt=None):
    nc_expt = Dataset(expt+'.nc','r')
    levels = nc_expt['plevs'][:]
    levels_up = nc_expt['plevs_up'][:]
    levels_down = nc_expt['plevs_down'][:]
    times = nc_expt['times'][:]
    print(times)

   #dates = [dateutils.hrstodate(time) for time in times] 
    dates = [hour2date(time) for time in times] 
    wind_counts = nc_expt['wind_obcounts'][:]
    temp_counts = nc_expt['temp_obcounts'][:]
    temp_fits = nc_expt['omf_rmstemp'][:]
    temp_bias = nc_expt['omf_biastemp'][:]
    wind_fits  = nc_expt['omf_rmswind'][:]
   #print(wind_fits.mean(axis=0))
   #raise SystemExit

    nc_expt.close()

    stats = {}
    stats['levels'] = levels
    stats['levels_up'] = levels_up
    stats['levels_down'] = levels_down
    stats['times'] = times
    stats['dates'] = dates
    stats['wind_counts'] = wind_counts
    stats['temp_counts'] = temp_counts
    stats['temp_fits'] = temp_fits
    stats['temp_bias'] = temp_bias
    stats['wind_fits'] = wind_fits

    return stats

  def plot(self, stats1=None, stats2=None, lbl1=None, lbl2=None, title=None, region=None):
    color1 = 'r'
    color2 = 'b'
    linewidth1 = 1.0
    linewidth2 = 1.0

    obtypes = 'All Insitu'

    rcParams['figure.subplot.left'] = 0.1 
    rcParams['figure.subplot.top'] = 0.85 
    rcParams['legend.fontsize']=12

    dates = stats1['dates']
    dates_txt = '%s-%s' % (dates[0],dates[-1])
    levels = stats1['levels']
    levels_up = stats1['levels_up']
    levels_down = stats1['levels_down']
    wind_fits1 = stats1['wind_fits']
    wind_fits2 = stats2['wind_fits']

    fig = plt.figure(figsize=(8,6))
    plt.subplot(1,2,1)
    wind_fits1mean = wind_fits1.mean(axis=0)
    wind_fits2mean = wind_fits2.mean(axis=0)
    wind_fits_pval = ttest(wind_fits1,wind_fits2,inflate=True)
    sig = wind_fits_pval >= self.sigthresh
    plt.plot(wind_fits1mean,levels,color=color1,linewidth=linewidth1,marker='o',label=lbl1)
    plt.plot(wind_fits2mean,levels,color=color2,linewidth=linewidth2,marker='o',label=lbl2)
    for n in range(len(levels)):
      if sig[n]:
        plt.axhspan(levels_up[n], levels_down[n], facecolor='lightgoldenrodyellow')
    plt.ylabel('pressure (hPa)')
    plt.title('%s: %s' % (obtypes+' V',region))
    plt.xlabel('RMS (mps)')
    plt.axis('tight')
    if int(dates[1][4:6]) > 10 or int(dates[2][4:6]) < 4:
      plt.xlim(2.5,4.0)
    else:
      plt.xlim(2.2,3.4)
      #plt.xlim(2.2,3.4)
      #plt.xlim(3.0,5.5)
    plt.ylim(self.levbot,self.levtop)
    plt.grid(True)

    temp_fits1 = stats1['temp_fits']
    temp_fits2 = stats2['temp_fits']

    plt.subplot(1,2,2)
    temp_fits1mean = np.sqrt(temp_fits1).mean(axis=0)
    temp_fits2mean = np.sqrt(temp_fits2).mean(axis=0) 
    temp_fits_pval = ttest(temp_fits1,temp_fits2,inflate=True)
    sig = temp_fits_pval >= self.sigthresh
    plt.plot(temp_fits1mean,levels,color=color1,linewidth=linewidth1,marker='o',label=lbl1)
    plt.plot(temp_fits2mean,levels,color=color2,linewidth=linewidth2,marker='o',label=lbl2)
    for n in range(len(levels)):
      if sig[n]:
        plt.axhspan(levels_up[n], levels_down[n], facecolor='lightgoldenrodyellow')
    plt.xlabel('RMS (K)')
    plt.title('%s: %s' % (obtypes+' T',region))
    plt.axis('tight')
    plt.xlim(0.5,1.4)
    plt.xlim(0.75,2.25)       
    plt.ylim(self.levbot,self.levtop)
    plt.grid(True)
    plt.legend(loc=0)

    plt.figtext(0.5,0.93,'%s RMS O-F (%s)' % (title,dates_txt),horizontalalignment='center',fontsize=14)
   #plt.savefig('obfits_%s_%s_%s.png' % (lbl1,lbl2,region))
    plt.savefig('obfits_%s_%s.png' % (title,region))

    plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  title = 'JEDI and GSI'
  region = 'All Insitu'
  fexp = 'gsi_stats'
  sexp = 'jedi_stats'
  lbl1 = 'GSI'
  lbl2 = 'JEDI'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'title=',
                                                'region=', 'fexp=', 'sexp=',
                                                'lbl1=', 'lbl2='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--fexp'):
      fexp = a
    elif o in ('--sexp'):
      sexp = a
    elif o in ('--title'):
      title = a
    elif o in ('--region'):
      region = a
    else:
      assert False, 'unhandled option'

#-----------------------------------------------------------------------------------------
  pjgd = Plot_JEDI_GSI_Diag(debug=debug, output=output)

  stats1 = pjgd.get_stats(expt=fexp)
  stats2 = pjgd.get_stats(expt=sexp)

  pjgd.plot(stats1=stats1, stats2=stats2, lbl1=lbl1, lbl2=lbl2, title=title, region=region)

