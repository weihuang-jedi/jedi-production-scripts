#!/bin/bash

 set -x

#------------------------------------------------------------------------------
 firstlist=(orig.jedi  orig.gsi gsi  orig.gsi  gsi   orig.gsi)
 secondlist=(jedi      jedi     jedi orig.jedi orig.jedi gsi)
 for j in ${!firstlist[@]}
 do
   first=${firstlist[$j]}
   second=${secondlist[$j]}
   echo "first: ${first}, second: ${second}"
   rm -f ${second}-${first}
 done
 rm -f ~/jg.tar

 exit 0

