#!/usr/bin/env python

from __future__ import division, print_function

import array
import sys
import struct
import math
import copy

import numpy as np


def read_outchem(outchem):
    lines = list(outchem)

    center = np.loadtxt(lines[9:10], dtype=np.float)
    spacing = float(lines[14].strip())
    gridpoints = np.loadtxt(lines[16:17], dtype=np.float)

    return center, spacing, gridpoints


def load_vdw_grids(outchem, vdwfile):
    with open(outchem) as f:
        center, spacing, gridpoints = read_outchem(outchem)

    scale = 1.0 / spacing
    size = max(*gridpoints)

    vdwA, vdwB = load_vdwgrids(vdwfile, gridpoints=gridpoints, byteswap=True)



    vdw_attractive = vanderwaals(vdwA,
                                 scale=scale,
                                 center=center,
                                 size=size,
                                 name='Attractive VdW')

    vdw_repulsive = vanderwaals(vdwB,
                                 scale=scale,
                                 center=center,
                                 size=size,
                                 name='Repulsive VdW')

    return vdw_attractive, vdw_repulsive


def load_vdwgrids(vdwfile, gridpoints=None, byteswap=True):
    if byteswap:
        vdw = np.fromfile(vdwfile, dtype='>f4')
    else:
        vdw = np.fromfile(vdwfile, dtype='<f4')

    vdwA, vdwB = vdw.reshape((2, -1))[:,1:-1]

    if gridpoints:
        vdwA = vdwA.reshape(*gridpoints, order='f')
        vdwB = vdwB.reshape(*gridpoints, order='f')
    else:
        vdwA = vdwA.reshape(order='f')
        vdwB = vdwB.reshape(order='f')

    return vdwA, vdwB


class vanderwaals(object):

  def __init__(self, array, scale, size, center=(0,0,0), name="VDW Grid"):
    '''reads the phi file from disk'''

    self.gridDimension = size
    self.scale = scale
    self.vdwArray = array
    self.oldmid = center
    self.toplabel = name
    self.title = ""
    self.head = ""
    self.botlabel = ""
    self.__minsmaxs = None
    self.__boundaries = None

  def copyPhi(self):
    '''make a deep copy of the phimap that can be edited without disturbing the
    original.'''
    newVdw = vanderwaals()
    newVdw.oldmid = self.oldmid
    newVdw.scale = self.scale
    newVdw.vdwArray = self.vdwArray
    newVdw.gridDimension = self.gridDimension
    newVdw.__minsmaxs = None
    newVdw.__boundaries = None
    return newVdw

  def write(self, phiFileName=False):
    '''write data to member data structure manually,
    then call this to write to file
    the pad lines reproduce the binary padding of an original
    fortran formatted phi file'''
    if phiFileName:  # do nothing if no filename given
      outArray = copy.deepcopy(self.vdwArray)
      outArray.byteswap()  # switch endianness back, only for writing
      phiFile = open(phiFileName, 'wb')  # b may be unnecessary, have to check
      phiFile.write(struct.pack('4b', 0, 0, 0, 20))  # pad
      phiFile.write(struct.pack('20s', self.toplabel))
      phiFile.write(struct.pack('8b', 0, 0, 0, 20, 0, 0, 0, 70))  # pad
      phiFile.write(struct.pack('10s', self.head))
      phiFile.write(struct.pack('60s', self.title))
      phiFile.write(struct.pack('4b', 0, 0, 0, 70))  # pad, always same
      phiFile.write(struct.pack('>l', len(outArray)*4))  # diff. pad sometimes
      #print "writing this many data points in phimap:", len(outArray)
      outArray.tofile(phiFile)  # array
      phiFile.write(struct.pack('>l', len(outArray)*4))  # diff. pad sometimes
      phiFile.write(struct.pack('4b', 0, 0, 0, 16))  # pad, always same
      phiFile.write(struct.pack('16s', self.botlabel))
      phiFile.write(struct.pack('8b', 0, 0, 0, 16, 0, 0, 0, 16))  # pad
      phiFile.write(struct.pack(
          '>ffff', self.scale, self.oldmid[0], self.oldmid[1], self.oldmid[2]))
      #> on previous line forces big-endian writing
      phiFile.write(struct.pack('4b', 0, 0, 0, 16))  # pad
      phiFile.close()

  def trimVdw(self, newmidIndices, newSize):
    '''for a new center index and a desired cubic grid size, trim the current
    phimap and return the new trimmed phimap'''
    plusMinus = (newSize - 1) / 2  # how many to add or subtract from the center
    newPhi = vanderwaals()
    newPhi.oldmid = self.getXYZlist(newmidIndices)  # only change of these data
    newPhi.toplabel = self.toplabel
    newPhi.head = self.head
    newPhi.title = self.title
    newPhi.botlabel = self.botlabel
    newPhi.scale = self.scale
    #the phiArray does change
    newPhi.phiArray = array.array('f')
    for oldIndexZ in range(
        newmidIndices[2] - plusMinus, newmidIndices[2] + plusMinus + 1):
      for oldIndexY in range(
          newmidIndices[1] - plusMinus, newmidIndices[1] + plusMinus + 1):
        for oldIndexX in range(
            newmidIndices[0] - plusMinus, newmidIndices[0] + plusMinus + 1):
          if oldIndexX >= 0 and oldIndexX < self.gridDimension and \
              oldIndexY >= 0 and oldIndexY < self.gridDimension and \
              oldIndexZ >= 0 and oldIndexZ < self.gridDimension:
            newPhi.phiArray.append(
                self.getValue(oldIndexX, oldIndexY, oldIndexZ))
          else:
            newPhi.phiArray.append(0.0)  # outside the original grid.
    #print "total array size is:", len(newPhi.phiArray)
    return newPhi

  def findVdwCorners(self, newmidIndices, newSize):
    '''for a new center index and a desired cubic grid size, find the new
    corners of the phimap'''
    plusMinus = (newSize - 1) / 2  # how many to add or subtract from the center
    lowerLeft = \
        [newmidIndices[0] - plusMinus, newmidIndices[1] - plusMinus,
         newmidIndices[2] - plusMinus]
    upperRight = \
        [newmidIndices[0] + plusMinus + 1, newmidIndices[1] + plusMinus + 1,
         newmidIndices[2] + plusMinus + 1]
    return lowerLeft, upperRight

  def findNewVdwIndices(self, newmidIndices, newSize):
    '''for a new center index and a desired cubic grid size, return xyz coords
    of each coordinate in the new box'''
    coordList = []
    plusMinus = (newSize - 1) / 2  # how many to add or subtract from the center
    for oldIndexZ in range(
        newmidIndices[2] - plusMinus, newmidIndices[2] + plusMinus + 1):
      for oldIndexY in range(
          newmidIndices[1] - plusMinus, newmidIndices[1] + plusMinus + 1):
        for oldIndexX in range(
            newmidIndices[0] - plusMinus, newmidIndices[0] + plusMinus + 1):
          coordList.append((oldIndexX, oldIndexY, oldIndexZ))
    return coordList

  def getMinsMaxs(self):
    '''finds the positions of the extreme grid corners'''
    if self.__minsmaxs is None:
      mins, maxs = [], []
      for center in self.oldmid:
        mins.append(center - ((self.gridDimension - 1.)/(2. * self.scale)))
        maxs.append(center + ((self.gridDimension - 1.)/(2. * self.scale)))
      self.__minsmaxs = mins, maxs
    return self.__minsmaxs

  def getMinMaxValues(self):
    '''finds the minimum and maximum value'''
    return min(self.vdwArray), max(self.vdwArray)

  def getMeanAbsoluteValues(self):
    '''takes the abs value of each phi value, then the average'''
    sum = 0.0
    for value in self.vdwArray:
      sum += math.fabs(value)
    return sum/float(len(self.vdwArray))

  def getMeanValues(self):
    '''mean of all phi values'''
    sum = 0.0
    for value in self.vdwArray:
      sum += value
    return sum/float(len(self.vdwArray))

  def getMaxValues(self):
    '''just the max'''
    return max(self.vdwArray)

  def countValues(self):
    '''counts the occurence of each value'''
    counts = {}
    for value in self.vdwArray:
      if value in counts:
        counts[value] += 1
      else:
        counts[value] = 1
    return counts

  def histogramValues(self, width=1., useMin=None, useMax=None):
    '''makes a basic histogram'''
    ends = list(self.getMinMaxValues())
    if useMin is not None:
      ends[0] = useMin
    if useMax is not None:
      ends[1] = useMax
    bars = int(math.ceil((ends[1] - ends[0]) / width) + 1)
    counts = [0 for x in range(bars)]
    for value in self.vdwArray:
      if value >= ends[0] and value <= ends[1]:
        counts[int(math.floor((value - ends[0]) / width))] += 1
    return counts

  def getXYZlist(self, xyz):
    '''changes list to x,y,z calls getXYZ'''
    return self.getXYZ(xyz[0], xyz[1], xyz[2])

  def getXYZ(self, xInd, yInd, zInd):
    '''returns the xyz coordinate of the center of the box'''
    mins, maxs = self.getMinsMaxs()
    gap = 1./self.scale
    return mins[0]+(xInd*gap), mins[1]+(yInd*gap), mins[2]+(zInd*gap)

  def getValueList(self, xyz):
    '''changes list into x, y, z then calls getValue'''
    return self.getValue(xyz[0], xyz[1], xyz[2])

  def getValue(self, xInd, yInd, zInd):
    '''for a given set of indices, return the value in the array'''
    index = int(zInd*(self.gridDimension**2.) + yInd*self.gridDimension + xInd)
    return self.vdwArray[index]

  def getValueListCheckBounds(self, xyzList, retValueIfBad=0):
    '''passes to getValueCheckBounds'''
    return self.getValueCheckBounds(
        xyzList[0], xyzList[1], xyzList[2], retValueIfBad)

  def getValueCheckBounds(self, xInd, yInd, zInd, retValueIfBad=0):
    '''does grid bounds checking first, returns retValueIfBad if outside grid,
    otherwise call getvalue'''
    if xInd >= 0 and xInd < self.gridDimension and \
        yInd >= 0 and yInd < self.gridDimension and \
        zInd >= 0 and zInd < self.gridDimension:
      return self.getValue(xInd, yInd, zInd)
    else:
      return retValueIfBad

  def setValueList(self, xyz, value):
    '''calls setValue with expanded xyz into items'''
    self.setValue(xyz[0], xyz[1], xyz[2], value)

  def setValue(self, xInd, yInd, zInd, value):
    '''puts the value into the phi array'''
    index = int(zInd*(self.gridDimension**2.) + yInd*self.gridDimension + xInd)
    self.vdwArray[index] = value

  def transform(self, threshold=6.0, inside=-2.0, outside=-1.0):
    '''for every value in the array, change it to inside or outside,
    destructively overwrites old values'''
    for index in range(len(self.vdwArray)):
      value = self.vdwArray[index]
      if value < threshold:
        where = outside
      else:
        where = inside
      self.vdwArray[index] = where

  def translate(self, newmid):
      oldmin = self.oldmid
      translate = [
          newmid[0] - oldmin[0],
          newmid[1] - oldmin[1],
          newmid[2] - oldmin[2]
      ]
      self.oldmid = newmid
      self.__boundaries = None
      self.__minsmaxs = None
      return translate

  def center(self):
      return self.translate([0, 0, 0])

  def subtract(self, other):
    '''subtract other from self, destructively write over self'''
    self.modify(other, -1)

  def add(self, other):
    '''add other to self, destructively write over self.'''
    self.modify(other, 1)

  def modify(self, other, change):
    '''modify other to self, destructively write over self. allows +-/etc
    presume without checking that grids are compatible (same mid etc)'''
    for index in range(len(self.vdwArray)):
      value = other.vdwArray[index]
      #save = self.phiArray[index]
      self.vdwArray[index] += (value * change)
      #if self.phiArray[index] != 0.0:
      #  print self.phiArray[index], value, save, index

  def findBoundaries(
      self, inside=-2.0, border=2, pointXYZ=None, pointList=None):
    '''finds the extreme x,y,z positions that enclose all inside positions'''
    if self.__boundaries is None:  # need to calculate it
      if pointXYZ is not None:
        self.__boundaries = self.findPointMinsMaxs(pointXYZ, pointList)
      else:
        self.__boundaries = [self.gridDimension, self.gridDimension,
                             self.gridDimension], [0, 0, 0]
      for x in range(self.gridDimension):
        for y in range(self.gridDimension):
          for z in range(self.gridDimension):
            if x < self.__boundaries[0][0] or x > self.__boundaries[1][0] or \
                y < self.__boundaries[0][1] or y > self.__boundaries[1][1] or \
                z < self.__boundaries[0][2] or z > self.__boundaries[1][2]:
              value = self.getValue(x, y, z)
              if value == inside:
                indices = (x, y, z)
                for coord in range(3):
                  self.__boundaries[0][coord] = min(
                      self.__boundaries[0][coord], indices[coord])
                  self.__boundaries[1][coord] = max(
                      self.__boundaries[1][coord], indices[coord])
      for coord in range(3):
        self.__boundaries[0][coord] = max(
            0, self.__boundaries[0][coord] - border)
        self.__boundaries[1][coord] = min(
            self.gridDimension, self.__boundaries[1][coord]+border)
    return self.__boundaries

  def getBoundaryLengths(self, inside=-2.0, border=2):
    '''calls findBoundaries if necessary, returns the lengths (max-min)'''
    if self.__boundaries is None:  # need to calculate it
      self.findBoundaries(inside, border)
    lengths = [self.__boundaries[1][0] - self.__boundaries[0][0],
               self.__boundaries[1][1] - self.__boundaries[0][1],
               self.__boundaries[1][2] - self.__boundaries[0][2]]
    return lengths

  def createFromGrid(
      self, grid, gridSize, defaultValue=0.0, toplabel="",
      head="", title="", botlabel="", lowestGridSize=65):
    '''does grid->phi data structure conversion'''
    self.toplabel = toplabel[:20]  # easy stuff first
    self.head = head[:10]
    self.title = title[:60]
    self.botlabel = botlabel[:16]
    lens = [len(grid), len(grid[0]), len(grid[0][0])]
    #have to expand to valid gridSize
    newGridSize = 0
    for possibleGridSize in self.gridSizes:
      good = True
      if possibleGridSize < lowestGridSize:
        good = False
      for oneLength in lens:
        if oneLength > possibleGridSize:
          good = False
      if good:
        newGridSize = possibleGridSize
    self.gridDimension = newGridSize
    #now take care of the grid
    self.vdwArray = array.array('f')
    for z in range(self.gridDimension):
      for y in range(self.gridDimension):
        for x in range(self.gridDimension):
          if x < lens[0] and y < lens[1] and z < lens[2]:
            self.vdwArray.append(grid[x][y][z][0])
          else:  # outside real grid
            self.vdwArray.append(defaultValue)
    #scale and oldmid are all that is left
    self.scale = 1./gridSize
    for coord in range(3):
      self.oldmid[coord] = grid[0][0][0][coord + 1] \
          - (gridSize / 2.) + (self.gridDimension / self.scale) / 2.
    #data should be ready for writing now

  def findPointMinsMaxs(self, pointXYZ, pointList):
    minsPts = pointXYZ[0][1:]
    maxsPts = pointXYZ[0][1:]
    for point in pointList:
      xyz = pointXYZ[point-1][1:]
      for coord in range(3):
        minsPts[coord] = min(minsPts[coord], xyz[coord])
        maxsPts[coord] = max(maxsPts[coord], xyz[coord])
    newMins = list(self.getIndices(minsPts))
    newMaxs = list(self.getIndices(maxsPts))  # so they initialize to pts
    return newMins, newMaxs

  def getIndices(self, pt):
    '''helper function to find the box a point is in'''
    mins, maxs = self.getMinsMaxs()
    gridSize = 1./self.scale
    xIndex = int(math.floor((pt[0]-mins[0])/gridSize))
    yIndex = int(math.floor((pt[1]-mins[1])/gridSize))
    zIndex = int(math.floor((pt[2]-mins[2])/gridSize))
    #print xIndex, yIndex, zIndex, mins, pt, maxs
    return xIndex, yIndex, zIndex

  def trilinear_interpolation(self, point):
    '''for a given point, find the box it is in, the trilinearly interpolate
    and return the value at that point. this is in kT, as that is what phiMaps
    hold. for usual applications, you want to take this times the charge and
    0.5924 to put it in kcal/mol'''
    ptX, ptY, ptZ = self.getIndices(point)
    values = [0. for count in range(8)]
    values[7] = self.getValue(ptX, ptY, ptZ)
    values[6] = self.getValue(ptX, ptY, ptZ + 1) - values[7]
    values[5] = self.getValue(ptX, ptY + 1, ptZ) - values[7]
    values[4] = self.getValue(ptX + 1, ptY, ptZ) - values[7]
    values[3] = self.getValue(ptX, ptY + 1, ptZ + 1) - values[7] - \
        values[6] - values[5]
    values[2] = self.getValue(ptX + 1, ptY, ptZ + 1) - values[7] - \
        values[6] - values[4]
    values[1] = self.getValue(ptX + 1, ptY + 1, ptZ) - values[7] - \
        values[5] - values[4]
    values[0] = self.getValue(ptX + 1, ptY + 1, ptZ + 1) - values[7] - \
        values[6] - values[5] - values[4] - values[3] - values[2] - values[1]
    gridPoint = self.getXYZ(ptX, ptY, ptZ)
    fraction = [0. for count in range(3)]
    for count in range(3):
      fraction[count] = point[count] - gridPoint[count]
    returnVdwValue = values[0] * fraction[0] * fraction[1] * fraction[2] + \
        values[1] * fraction[0] * fraction[1] + \
        values[2] * fraction[0] * fraction[2] + \
        values[3] * fraction[1] * fraction[2] + values[4] * fraction[0] + \
        values[5] * fraction[1] + values[6] * fraction[2] + values[7]
    #print values, fraction, returnPhiValue
    return returnVdwValue


  def trimToBoxCenterAndSize(self, corners, center, dimensions):
    '''given a box, find the new center and size of a valid phimap based on
    this current phimap'''
    #print corners, center, dimensions
    #print self.scale, self.oldmid
    #find the midpoint and corners
    centerIndices = self.getIndices(center)
    newmid = self.getXYZlist(centerIndices)  # becomes the new oldmid
    onecorner = self.getIndices(corners[0:3])
    twocorner = [coord + 1 for coord in self.getIndices(corners[3:6])]
    #phimap grid can only be cubic
    biggestDimension = 0
    if twocorner[1] - onecorner[1] > twocorner[0] - onecorner[0]:
      biggestDimension = 1
    if (twocorner[2] - onecorner[2] >
        twocorner[biggestDimension] - onecorner[biggestDimension]):
      biggestDimension = 2
    newSize = twocorner[biggestDimension] - onecorner[biggestDimension]
    if 0 == newSize % 2:  # if size is even, that's not allowed, so,
      newSize += 1  # make it odd
    return centerIndices, newSize

  def trimToBox(self, corners, center, dimensions):
    '''given a box (see box.py) trim so that the box is enclosed but not more.
    returns the new trimmed phimap'''
    centerIndices, newSize = self.trimToBoxCenterAndSize(
        corners, center, dimensions)
    return self.trimVdw(centerIndices, newSize), centerIndices, newSize

if __name__ == '__main__':
  pass
