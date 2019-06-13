import math
import trep
from trep import Potential
from trep._trep import _TorsionalSpringPotential
import numpy as np
import numpy.linalg

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    _opengl = True
except:
    _opengl = False

class TorsionalSpring(_TorsionalSpringPotential, Potential):
    """
    Spring implements a fixed-length, fixed-spring-constant
    spring against a configuratoin variable:

    V(q) = 0.5*k*(q - q0)**2
    """
    def __init__(self, system, frame1, frame2, k, q0=0.0, name=None):
        Potential.__init__(self, system, name)
        _TorsionalSpringPotential.__init__(self)
        self._q0 = 0.0
        self._k = 0.0

        if not system.get_frame(frame1):
            raise ValueError("Could not find frame %r" % frame1)
        self._frame1 = system.get_frame(frame1)
        
        if not system.get_frame(frame2):
            raise ValueError("Could not find frame %r" % frame2)
        self._frame2 = system.get_frame(frame2)
        
        self._q0 = q0
        self._k = k

    @property
    def frame1(self):
        return self._frame1

    @property
    def frame2(self):
        return self._frame2

    @property
    def q0(self):
        return self._q0

    @q0.setter
    def q0(self, q0):
        self._q0 = q0

    @property 
    def k(self):
        return self._k

    @k.setter
    def k(self, k):
        self._k = k
            
