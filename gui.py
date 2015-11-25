#!/usr/bin/python

import sys
import numpy
import math

from OpenGL.GL import *
from OpenGL.GLU import gluPerspective, gluNewQuadric, gluSphere, gluCylinder, gluLookAt
from PyQt4.QtCore import Qt
from PyQt4.QtGui import QMainWindow, QWidget
from PyQt4.QtOpenGL import QGLWidget

from chem import Molecule

class MainWindow(QMainWindow):
    def __init__(self, filename, quality, parent=None):
        super(MainWindow, self).__init__(parent)
        self.canvas = OpenGLCanvas(filename, quality)
        self.setCentralWidget(self.canvas)
        self.setWindowTitle("Window")
        self.resize(500,500)
    def keyPressEvent(self, event):
        key = event.key()
        if key == Qt.Key_Plus:
                self.canvas.translateModel(0.05)
        if key == Qt.Key_Minus:
                self.canvas.translateModel(-0.05)
    def mousePressEvent(self, event):
        self.canvas.x, self.canvas.y = self.canvas.lastx, self.canvas.lasty = event.pos().x(), event.pos().y()
    def mouseMoveEvent(self, event):
        self.canvas.lastx, self.canvas.lasty = self.canvas.x, self.canvas.y
        self.canvas.x, self.canvas.y = event.pos().x(), event.pos().y()
        self.canvas.updateGL()
    def wheelEvent(self, event):
        self.canvas.zoomModel(event.delta()*0.001)

class OpenGLCanvas(QGLWidget):
    def __init__(self, filename, quality, parent=None, name=None):
        QGLWidget.__init__(self,parent,name)
        # load molecule and read from xyz file
        self.mol = Molecule()
        self.mol.readxyz(filename)
        self.mol.genbonds()
        # adjust quality to sanity
        # set display quality
        if self.mol.natoms > 10000:
            self.quality = min(2,quality)
        else:
            self.quality = quality
        print(quality)
        self.atom_size = 0.4
        self.bond_size = 0.1
        self.offset = -self.mol.box[0]-self.mol.box[1]
        # set initial coordinates
        self.x,self.y = self.lastx,self.lasty = [0,0]
    def initializeGL(self):
        # set OpenGL global states
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_LIGHT1)
        glEnable(GL_COLOR_MATERIAL)
        glClearColor(0., 0., 0., 1.)
        glClearDepth(1.0)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glCullFace(GL_BACK)
        glEnable(GL_CULL_FACE)
        glEnable(GL_BLEND)
        # load identity matrix and zoom out
        glLoadIdentity()
        glTranslate(0.,0.,self.offset)
        # draw the elements, set up paint callback
        if self.quality > 0:
            self.drawBallAndStick(self.atom_size, self.bond_size)
            self.paintScene = self.paintBallAndStick
        else:
            self.drawWireFrame()
            self.paintScene = self.paintWireFrame
    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        # rotation matrix adopted from molcas gv, needs improvement!!
        w = max(self.w, 1.0)
        h = max(self.h, 1.0)
        xScale = 2.0 / w
        yScale = 2.0 / h
        x = (self.y - self.lasty) * yScale
        y = (self.x - self.lastx) * xScale
        z = (self.x * self.lasty - self.y * self.lastx) * xScale * yScale
        a = numpy.sqrt((x**2+y**2+z**2) * (self.lastx**2+self.lasty**2) / (self.x**2+self.y**2+1)) * 180.
        glRotatef(a, x, y, z)
        # paint the scene
        self.paintScene()
    def resizeGL(self, w, h):
        self.w = w
        self.h = h
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(w)/float(h), 0.1, 1000.0)
        glMatrixMode(GL_MODELVIEW)
    def drawBallAndStick(self, r_ball = 0.4, r_stick = 0.1):
        res = self.quality * 6
        quad = gluNewQuadric()
        glNewList(1,GL_COMPILE)
        for i in range(self.mol.natoms):
            glPushMatrix()
            # set atom color
            r,g,b = self.mol.col[i]
            glColor3f(r,g,b)
            # translate to atom position
            x,y,z = self.mol.pos[i]
            glTranslatef(x,y,z)
            # draw the atom
            gluSphere(quad, r_ball*self.mol.rad[i], res, res)
            glPopMatrix()
        for i,c in enumerate(self.mol.con):
            m = c[0]
            n = c[1]
            glPushMatrix()
            # set bond color
            r,g,b = self.mol.col[m]
            glColor3f(r,g,b)
            # translate to atom position
            x,y,z = self.mol.pos[m]
            glTranslatef(x,y,z)
            # draw the bond
            v = self.mol.pos[n] - self.mol.pos[m]
            angle = math.acos(numpy.dot([0.0,0.0,1.0], v)/self.mol.dist[i])*180./math.pi
            x,y,z = numpy.cross([0.0,0.0,1.0], v)
            glRotatef(angle,x,y,z)
            h = self.mol.dist[i] * self.mol.rad[m] / (self.mol.rad[m] + self.mol.rad[n])
            gluCylinder(quad, r_stick, r_stick, h, res, res)
            glPopMatrix()
        glEndList()
    def paintBallAndStick(self):
        glCallList(1)
    def drawWireFrame(self):
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        glVertexPointer(3, GL_FLOAT, 0, self.mol.pos)
        glColorPointer(3, GL_FLOAT, 0, self.mol.col)
        # experimental VBO usage, currently not working
        #glGenBuffers(2, vboID)
        #glBindBuffer(GL_ARRAY_BUFFER, vboID[0])
        #glBufferData(GL_ARRAY_BUFFER, len(self.mol.pos)*4, self.mol.pos, GL_STATIC_DRAW)
        #glBindBuffer(GL_ARRAY_ELEMENTS_BUFFER, vboID[1])
        #glBufferData(GL_ARRAY_ELEMENTS_BUFFER, len(self.mol.con)*4, self.mol.con, GL_STATIC_DRAW)
    def paintWireFrame(self):
        glDrawElements(GL_POINTS, self.mol.natoms, GL_UNSIGNED_INT, numpy.arange(self.mol.natoms))
        glDrawElements(GL_LINES, 2*self.mol.nbonds, GL_UNSIGNED_INT, self.mol.con)
    def drawSurface(self):
        pass
    def paintSurface(self):
        pass
    def translateModel(self, zoom):
        mat = glGetFloatv(GL_MODELVIEW_MATRIX)
        glLoadIdentity()
        glTranslatef(0.,0.,zoom*abs(self.offset))
        glMultMatrixf(mat)
        mat = glGetFloatv(GL_MODELVIEW_MATRIX)
        self.updateGL()
    def zoomModel(self, zoom):
        self.translateModel(zoom)
