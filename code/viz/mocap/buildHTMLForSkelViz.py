#!/bin/py
'''
Create summary visualization for a dataset
'''
import glob
import os.path
import sys
import getopt
import numpy as np
import commands
import urllib2
from BeautifulSoup import BeautifulSoup

vizdir = '/home/mhughes/git/NPBayesHMM/figs/SkelVisBIG/'

MAXbID=33

SubjInfoMap = dict()
TrialInfoMap = dict()

def ExtractInfoFromWeb( subjID, trialID ):
  doSubj = doTrial = True
  if subjID in SubjInfoMap:
    subjMsg = SubjInfoMap[subjID]
    doSubj = False
  if subjID in TrialInfoMap and trialID in TrialInfoMap[subjID]:
    trialMsg = TrialInfoMap[subjID][trialID]
    doTrial = False
  if doSubj or doTrial:
    rawHTML = urllib2.urlopen( 'http://mocap.cs.cmu.edu/search.php?subjectnumber=%d' % (subjID) ).read()
    sloc = rawHTML.index( 'Subject #' )
    sloc = rawHTML.index( '(', sloc )
    eloc = rawHTML.index( ')', sloc )
    subjMsg = rawHTML[ sloc+1:eloc ].strip()
    SubjInfoMap[subjID] = subjMsg
    
    trialMarker =  '>%d</TD><TD>' % (trialID)
    try:
      sloc = rawHTML.index(trialMarker )
      sloc = sloc + len(trialMarker)
      eloc = rawHTML.index( '</TD>', sloc )
      trialMsg = rawHTML[ sloc:eloc ].strip()
    except Exception:
      trialMsg = ''
    if subjID not in TrialInfoMap:
      TrialInfoMap[subjID] = dict()
    TrialInfoMap[subjID][trialID] = trialMsg
  return subjMsg, trialMsg

def GetBasename( filepathstr ):
  baseFName, baseExt = os.path.splitext( filepathstr )
  return os.path.basename( baseFName )

def CreateHTML( bID, nX=3, nY=2 ):
  print 'Building page visualizations for behavior %d' % (bID)
  # Obtain list of all videos in category
  VNames = glob.glob( os.path.join( vizdir, 'Behavior%03d*.png'%(bID) ) )
  VNames = [GetBasename(v)+'.png' for v in VNames]
  
  bIDs = np.arange( MAXbID )
  [iip,iin] = bIDs.take( [bID-1, bID+1], mode='wrap' )
  if iip==0: iip = MAXbID

  origHTMLstr = open( os.path.join( vizdir, 'html', "template.html") ,'r' ).read()
  HTMLstr = origHTMLstr
  HTMLstr = HTMLstr.replace( "CURNAME", "Behavior %02d" % (bID) )
  HTMLstr = HTMLstr.replace( "PREVNAME", "%02d" % (iip) )
  HTMLstr = HTMLstr.replace( "NEXTNAME", "%02d" % (iin) )
  for ii in range(6):
    if ii >= len( VNames): 
      imgpath = ''
      seqtext = ''
    else:
      imgpath = VNames[ii]      
      seqname = VNames[ii][-9:-4]
      Fields = seqname.split('_')
      subjID = int( Fields[0] )
      trialID = int( Fields[1] )
      subjMsg, trialMsg = ExtractInfoFromWeb( subjID, trialID )
      seqtext = 'Subj. %d: %s <br /> Trial %d: %s' % (subjID, subjMsg, trialID, trialMsg)
    HTMLstr = HTMLstr.replace( "IMGPATH%02d"%(ii+1), imgpath )
    HTMLstr = HTMLstr.replace( "SEQNAME%02d"%(ii+1), seqtext )
  outhtmlpath =os.path.join(vizdir, 'html', "%02d.html" % (bID))
  outF = open(outhtmlpath , "w+")
  outF.write( HTMLstr )
  outF.close()

for kk in range( MAXbID ):
  CreateHTML( kk+1 )

