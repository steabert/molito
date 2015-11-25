# chemistry module

import sys
import numpy

class Element():
    def __init__(self, Z, radius, max_bonds, r, g, b):
        self.Z = Z
        self.radius = radius
        self.max_bonds = max_bonds
        self.color = numpy.array([r,g,b], dtype='float32')
    def __str__(self):
        format_string = ' '.join(['{:4d}','{:9f}','{:4d}','{:9f}'*3])
        return format_string.format(self.Z, self.radius, self.max_bonds, *self.color)

class Molecule():
    def __init__(self, size=0):
        self.natoms = size
        # set up a PSE dictionary with element data
        self.PSE = {}
        for line in open('PSE.txt'):
            if line.startswith('#'): continue
            data = line.split()
            self.PSE[data[0]] = Element(*data[1:])

    def readxyz(self, filename, extension='xyz', center=True):
        if extension != 'xyz':
            print('I can only handle xyz format')
            sys.exit(1)
        f = open(filename)
        # read number of atoms and title
        self.natoms = int(next(f))
        self.hdr = next(f)
        # check number of atoms is valid
        if self.natoms < 0:
            print('Invalid number of atoms (negative)')
            sys.exit(1)
        elif self.natoms > 100000:
            print('Molecule too large (> 100 000 atoms)')
            sys.exit(1)
        # read and store atom coordinates
        data = [line.split() for line in f]
        f.close()
        self.atom_types = [d[0] for d in data]
        self.pos = numpy.array([d[1:] for d in data], dtype='float32')
        self.rad = numpy.array([self.PSE[el].radius for el in self.atom_types], dtype='float32')
        self.col = numpy.array([self.PSE[el].color for el in self.atom_types], dtype='float32')
        # get origin and box size
        self.pos = self.pos - self.pos.min(axis=0)
        self.ori = numpy.array([0,0,0],dtype='float32')
        self.box = self.pos.max(axis=0)
        # centrate if requested
        if center:
            R = numpy.array([0.,0.,0.])
            R = sum(self.pos)/self.natoms
            self.pos = self.pos - R
    def genbonds(self):
        # create a distance matrix, rather inefficient since each atoms is checked against
        # all the other atoms (memory + time scaling as natoms^2)
        # TODO: block the distance compution into a loop over 3D blocks of ~(3Å)^3 and
        # only consider distances to atoms in adjacent blocks, bringing the scaling down
        # to n instead of n^2.
        delta = [None]*3
        for i in range(3):
            coord = self.pos[:,i]
            delta[i] = coord - coord.reshape((len(coord),1))
        dist = delta[0]**2 + delta[1]**2 + delta[2]**2 + 10. * numpy.identity(self.natoms)
        rcov = (self.rad + self.rad.reshape((len(self.rad),1)))*2
        # get array of indices where the distances is below the covalent radius
        self.con = numpy.array(numpy.where(dist-rcov<0.)).transpose()
        self.dist = numpy.sqrt(dist[self.con[:,0],self.con[:,1]])
        self.nbonds = len(self.dist)

