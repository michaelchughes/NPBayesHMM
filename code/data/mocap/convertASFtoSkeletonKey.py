#!/usr/bin/python
'''
Parse Acclaim Skeleton File (.asf) for Motion Capture Data
  and create a simple, easy-to-read summary of the 
  Skeleton joints measured by the capture process.
Author: Mike Hughes (mhughes@cs.brown.edu)
Dependencies: none
Examples:
  Parse "MyASFFile.asf" and save the result in folder "amc/"
    >>> python convertASFtoSkeletonKey.py MyASFFile.asf amc/
'''
import sys
import getopt
import io


# =========================================================
# convertToSkeletonKey( asffilename, outfilename )
#   parses input file (.asf), and writes simple output to 
#    outfilename
# INPUT: string file paths
#    outname = '/data/liv/mhughes/data/MoCapBIG/amc/skeldof'
#    fname   = '/data/liv/mhughes/data/MoCapBIG/asf/01.asf'
#    <or>
#    outname = '/data/liv/mhughes/data/KitchenMocap/amc/skeldof'
#      fname = '/data/liv/mhughes/data/KitchenMocap/asf/S01_Brownie.asf'
# OUTPUT: none
# =========================================================
def convertASFtoSkeletonKey( fname, amcname, outname ):
  fid = open( fname, 'rt' )
  ASFlist = fid.readlines()
  ASFstring = ''.join( ASFlist )
  fid.close()

  name = 'root'
  dof = 'tx ty tz rx ry rz'
  JNames = list()
  JKey = list()
  JNames.append( name )
  JKey.append( "%s %s" % (name, dof)  )

  print "Reading joint measurement info from ASF file"
  sID = 0
  eID = 0
  while sID >= 0:
    sID = ASFstring.find( "begin")
    eID = ASFstring.find( "end" )
    PartSpec = ASFstring[sID:eID]
    Fields = [f.strip() for f in PartSpec.split("\n")]
    dof = ''
    for line in Fields:
      if line.startswith("dof"):
        dof = line[4:]      
      if line.startswith("name"):
        name = line[5:]
    if len( dof ) > 1:
      JNames.append( name )
      JKey.append( "%s %s" % (name, dof)  )
      #outfile.write( "%s %s\n" %(name, dof)  )
    ASFstring = ASFstring[eID+1:]

  print "Aligning with proper order in AMC file"
  # Read in AMC file to ensure field names 
  #   are in correct order
  AMCNames = list()
  amcfid = open( amcname, 'rt' )
  didSeeData = False
  for line in amcfid.readlines()[:100]:
    if line.startswith(":") or line.startswith("#"):
      continue
    line = line.strip()
    Fields = line.split()
    if len( Fields ) == 1 and didSeeData:
      break
    if len( Fields ) > 1:
      didSeeData = True
      AMCNames.append( Fields[0] )
  #print AMCNames
  #print JNames
  sortIDs = list()
  for aname in AMCNames:
    for aa in xrange( 0, len(JNames) ):  
      if aname == JNames[aa]:
        sortIDs.append( aa )
  JKey = [ JKey[ x ] for x in sortIDs ]

  print "Writing to output key: ", outname
  outfile = open( outname, 'w' )
  for line in JKey:
    outfile.write( "%s\n" % (line)  )
  outfile.close()

# =========================================================
# MAIN()
# =========================================================
NUM_ARGS = 3
def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
  except getopt.error, msg:
    print msg
    print "for help use -h OR --help"
    sys.exit(2)
  # ------------------------------------------------ parse opts
  for o, a in opts:
    if o in ("-h", "--help"):
      print __doc__
    sys.exit(0)
  # ------------------------------------------------ parse args
  if len( args ) is not NUM_ARGS:
    print "Usage: convertASFtoSkeletonKey.py [OPTIONS] <asffilename> <amcfilename> <outfilename>"
  else:
    convertASFtoSkeletonKey( args[0], args[1], args[2] )      
    
if __name__ == "__main__":
    main()

