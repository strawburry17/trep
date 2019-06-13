:class:`TorsionalSpring` -- Torsional spring between two points
=========================================================

.. currentmodule:: trep.potentials


:class:`TorsionalSpring` creates a spring force between two coordinate frames 
which are coaxial (eg: the x-axis of those two frames are coaxial and point to 
the same direction) in 3D space. (Note if this prerequisite is not satisfied an error
message will be generated indicating that the frames selected are not coaxial):

.. math::

   q = q_1 - q_2

   V(q) = -k(q - q_0)^2


where :math:`q_1` and :math:`q_2` are the angles of two coordinate
frames relative to a fixed coordinate system 'World', math:`q` is the angular difference between two coordinate frames, math:`k` is the spring stiffness and :math:`q_0` is the natural 'length' of the spring. 


   .. table:: **Implemented Calculations**

      ===========   ===========
      Calculation   Implemented
      ===========   ===========
      V                 Y
      V_dq              Y
      V_dqdq            Y
      V_dqdqdq          Y 
      ===========   ===========

.. warning::
   
   The current implementation will fail if :math:`q_1` equals :math:`q_2` and
   :math:`q_0` is nonzero because of a divide by zero problem.

   If the  the angles of two coordinate frames relative to a fixed
   coordinate system 'World' are equal and :math:`q_0` is not zero,
   there should be a force.  But since there is no vector between the
   two frames, the direction of this force is undefined.  When the
   natural 'length' is zero, this problem can be corrected because the
   force also goes to zero.

Examples
--------

We can create a simple 1D harmonic oscillator using
:class:`TorsionalSpring` with a frame that is free to translate:

.. literalinclude:: code_snippets/Torsionalspring-ex1.py

The :class:`TorsionalSpring` works between arbitrary frames, not just
frames connected by a translation.  Here, we create two pendulums and
connect their masses with a spring:

.. literalinclude:: code_snippets/Torsionalspring-ex2.py

Torsional Spring Objects
--------------------

.. class:: TorsionalSpring(system, frame1, frame2, k[, q0=0.0, name=None])

   Create a new torsional spring between *frame1* and *frame2*.  The frames must
   already exist in the system and coaxial at the same time.

.. attribute:: TorsionalSpring.frame1
               
   The coordinate frame at one end of the spring.

   *(read only)*

.. attribute:: TorsionalSpring.frame1
               
   The coordinate frame at the other end of the spring.

   *(read only)*

.. attribute:: TorsionalSpring.q0

   The natural 'length' of the spring.  
                   
.. attribute:: TorsionalSpring.k

   The spring constant of the spring.




