'''
fetchAMCFilesFromWeb.py

Quick python script to gather Mocap sensor data from CMU database on the web.
This sensor data comes in two flavors:
   ASF : skeleton files (one per actor) that define joint angles, etc
   AMC : actual motions observed in each sequence

The particular dataset we're interested comes from two subjects
   ID # 13 and 14 on the CMU mocap website
'''
import convertASFtoSkeletonKey as CC
import sys
import os
from subprocess import call

trialInfo = dict()
trialInfo[13] = [29, 30, 31] 
trialInfo[14] = [6, 14, 20]

'''
L=0
for kk in trialInfo.keys():
  print " ==================================== Subj %02d" % (kk )
  print "Acquiring Skeleton file (ASF)"
  urlPath = "http://mocap.cs.cmu.edu:8080/subjects/%02d/%02d.asf" % (kk, kk)
  outPath = "asf/%02d.asf" % (kk)
  return_code = call("wget -nv -O"+outPath+" "+urlPath, shell=True)
  L += len( trialInfo[kk] )
  print "Acquiring motion info for %d sequences (AMC)" % (len(trialInfo[kk]) )
  for trialID in trialInfo[kk]:
    urlPath = "http://mocap.cs.cmu.edu:8080/subjects/%02d/%02d_%02d.amc" % (kk, kk, trialID)
    outPath = "amc/%02d_%02d.amc" % (kk, trialID)
    #    print urlPath
    return_code = call("wget -nv -O"+outPath+" "+urlPath, shell=True)
print '... SUCCESS! %d total sequences acquired and stored in amc/ and asf/' % (L)
'''

# Now, extract human-readable summary of the actual joint angles tracked.
#  only need to do this for one actor sequence,
#  since the skeleton joint names are consistent for ALL trials
kk = 13
asfname = "asf/%02d.asf" % (kk)
amcname = "amc/%02d_%02d.amc" % (kk, trialInfo[kk][0] )
outname = 'amc/SkeletonJoints.key'
CC.convertASFtoSkeletonKey( asfname, amcname, outname )
print '... SUCCESS! Wrote SkeletonJoints.key to amc/'
