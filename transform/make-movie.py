import os
import sys
import getopt
import moviepy.video.io.ImageSequenceClip

#-----------------------------------------------------------------------------------------
debug = 1
dir = '.'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'dir='])

for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--dir'):
    dir = a
  else:
    assert False, 'unhandled option'

fps=1
mni = 1000
varlist = ['u', 'v', 'T', 'delp', 'delz', 'sphum', 'o3mr']

for variable in varlist:
  videoname = 'video_%s.mp4' %(variable)

  image_files = []
  n = 0
  while(n < mni):
    n += 1
    imgname = '%s/%s_lev_%3.3d.png' %(dir, variable, n)
    if(os.path.exists(imgname)):
      image_files.append(imgname)

  if(len(image_files) > 5):
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(videoname)

