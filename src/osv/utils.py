"""
Jython utilities for optimal surface voting.
Author: Xinming Wu, University of Texas at Austin
Version: 2016.01.28
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../data/"
_pngdir = "../../png/"

#############################################################################
# Setup

def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2
  Example: setupForSubset("pnz")
  """
  global pngDir
  global seismicDir
  global s1,s2,s3
  global n1,n2,n3
  if name=="f3d2d":
    """ 2d fault """
    print "setupForSubset: 2d fault"
    pngDir = _pngdir+"2d/f3d2d/"
    seismicDir = _datdir+"2d/f3d2d/"
    n1,n2 = 222,440 #f3d75s
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="f3d":
    print "setupForSubset: f3d"
    pngDir = _pngdir+"3d/f3d/"
    seismicDir = _datdir+"3d/f3d/"
    n1,n2,n3 = 100,400,420
    d1,d2,d3 = 1.0,1.0,1.0 
    #j1,j2,j3 = 344,0,0
    #d1,d2,d3 = 0.004,0.025,0.024999 # (s,km,km)
    #f1,f2,f3 = 0.472,0.0,0.0 # = 0.000,0.000,0.000
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="petrobras":
    print "setupForSubset: petrobras"
    pngDir = _pngdir+"petrobras/"
    seismicDir = _datdir+"petrobras/"
    n1,n2,n3 = 400,1100,750
    d1,d2,d3 = 1.0,1.0,1.0 
    #j1,j2,j3 = 344,0,0
    #d1,d2,d3 = 0.004,0.025,0.024999 # (s,km,km)
    #f1,f2,f3 = 0.472,0.0,0.0 # = 0.000,0.000,0.000
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="campos":
    print "setupForSubset: campos"
    pngDir = _pngdir+"3d/campos/"
    seismicDir = _datdir+"3d/campos/"
    n1,n2,n3 = 450,1950,1200
    n1,n2,n3 = 300,600,400
    d1,d2,d3 = 1.0,1.0,1.0 
    #j1,j2,j3 = 344,0,0
    #d1,d2,d3 = 0.004,0.025,0.024999 # (s,km,km)
    #f1,f2,f3 = 0.472,0.0,0.0 # = 0.000,0.000,0.000
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="swj":
    print "setupForSubset: swj"
    pngDir = _pngdir+"swj/"
    seismicDir = _datdir+"swj/"
    n1,n2,n3 = 451,371,271
    d1,d2,d3 = 1.0,1.0,1.0 
    #j1,j2,j3 = 344,0,0
    #d1,d2,d3 = 0.004,0.025,0.024999 # (s,km,km)
    #f1,f2,f3 = 0.472,0.0,0.0 # = 0.000,0.000,0.000
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="anadarko":
    print "setupForSubset: anadarko"
    pngDir = _pngdir+"anadarko/"
    seismicDir = _datdir+"anadarko/"
    n1,n2,n3 = 333,480,640
    d1,d2,d3 = 1.0,1.0,1.0 
    #j1,j2,j3 = 344,0,0
    #d1,d2,d3 = 0.004,0.025,0.024999 # (s,km,km)
    #f1,f2,f3 = 0.472,0.0,0.0 # = 0.000,0.000,0.000
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="f3d2d":
    """ 2d fault """
    print "setupForSubset: 2d fault"
    pngDir = _pngdir+"2d/f3d2d/"
    seismicDir = _datdir+"2d/f3d2d/"
    n1,n2 = 222,440 #f3d75s
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  else:
    print "unrecognized subset:",name
    System.exit

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

def getPngDir():
  return pngDir

#############################################################################
# read/write images
def readImageChannels(basename):
  """ 
  Reads three channels of a color image
  """
  fileName = seismicDir+basename+".jpg"
  il = ImageLoader()
  image = il.readThreeChannels(fileName)
  return image
def readColorImage(basename):
  """ 
  Reads three channels of a color image
  """
  fileName = seismicDir+basename+".jpg"
  il = ImageLoader()
  image = il.readColorImage(fileName)
  return image

def readImage(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
def readImageX(m1,m2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(m1,m2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage3DX(m1,m2,m3,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(m1,m2,m3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image


def readImage3D(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage3DL(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def readImageL(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readInts(image)
  ais.close()
  return image

def readImage2L(n1,n2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image


def readImage1D(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage1L(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  #aos.writeBytes(image)
  aos.close()
  return image

def writeImageL(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# read/write fault skins

def skinName(basename,index):
  return basename+("%03i"%(index))
def skinIndex(basename,fileName):
  assert fileName.startswith(basename)
  i = len(basename)
  return int(fileName[i:i+3])

def listAllSkinFiles(basename):
  """ Lists all skins with specified basename, sorted by index. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  return fileNames

def removeAllSkinFiles(basename):
  """ Removes all skins with specified basename. """
  fileNames = listAllSkinFiles(basename)
  for fileName in fileNames:
    File(seismicDir+fileName).delete()

def readSkin(basename,index):
  """ Reads one skin with specified basename and index. """
  return FaultSkin.readFromFile(seismicDir+skinName(basename,index)+".dat")

def readSkins(basename):
  """ Reads all skins with specified basename. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  skins = []
  for iskin,fileName in enumerate(fileNames):
    index = skinIndex(basename,fileName)
    skin = readSkin(basename,index)
    skins.append(skin)
  return skins

def writeSkin(basename,index,skin):
  """ Writes one skin with specified basename and index. """
  FaultSkin.writeToFile(seismicDir+skinName(basename,index)+".dat",skin)

def writeSkins(basename,skins):
  """ Writes all skins with specified basename. """
  for index,skin in enumerate(skins):
    writeSkin(basename,index,skin)

from org.python.util import PythonObjectInputStream
def readObject(name):
  fis = FileInputStream(seismicDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  obj = ois.readObject()
  ois.close()
  return obj
def writeObject(name,obj):
  fos = FileOutputStream(seismicDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(obj)
  oos.close()
