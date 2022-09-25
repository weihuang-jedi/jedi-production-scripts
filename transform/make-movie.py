import os
import sys
import getopt
import moviepy.video.io.ImageSequenceClip

#-----------------------------------------------------------------------------------------
debug = 1
dir = '.'
variable = 'T'

#-----------------------------------------------------------------------------------------
opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'dir=', 'variable='])

for o, a in opts:
  if o in ('--debug'):
    debug = int(a)
  elif o in ('--variable'):
    variable = a
  elif o in ('--dir'):
    dir = a
  else:
    assert False, 'unhandled option'

videoname = 'video_%s.mp4' %(variable)

image_files = []
mni = 1000
n = 0
while(n < mni):
  n += 1
  imgname = '%s/%s_lev_%3.3d.png' %(dir, variable, n)
  if(os.path.exists(imgname)):
    image_files.append(imgname)

fps=1
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile(videoname)

