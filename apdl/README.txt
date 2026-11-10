ME5103 LAB 5 - MECHANICAL APDL SCRIPTS
======================================

These command files use N-mm-MPa units, plane stress, and a 1 mm thickness:

  E  = 200000 MPa
  nu = 0.30
  T  = S = 100 MPa
  plate = 300 mm x 300 mm
  hole radius = 10 mm

Files
-----

q1_full_tri.txt
    Question 1, full plate, PLANE182 triangular mesh.

q1_full_quad.txt
    Question 1, full plate, PLANE182 quadrilateral mesh.

q2_quarter_tri.txt
    Question 2, one-quarter model, PLANE182 triangular mesh.

q2_quarter_quad.txt
    Question 2, one-quarter model, PLANE182 quadrilateral mesh.

q1_partitioned_diagonal_tri.txt
    Question 1, full plate with a free triangular mesh. It is the diagonal/
    triangle comparison case, but it is not the recommended final mesh.

q1_partitioned_square_ring_quad.txt
    Question 1, square-ring topology comparison with a mapped quadrilateral
    mesh. This file uses a 20 mm square opening, not the exact circular hole.

q2_partitioned_square_ring_quad.txt
    Question 2, quarter square-ring topology comparison with a mapped
    quadrilateral mesh. This file uses a 10 mm square corner opening, not the
    exact quarter-circle.

How to run
----------

In Mechanical APDL:
  File > Read Input From...
or at the command prompt:
  /INPUT,q1_full_tri,txt

Each file writes a database named lab5_*.db and opens PNG plots for:
  SX, SY, SXY, and principal stress intensity (when available).

The scripts also print the nodal stress extrema and reaction-force totals.
Use /POST1, SET,LAST, PRNSOL,S,COMP and PRRSOL for additional tabulated data.

Important modeling note
-----------------------

The circular hole is represented by straight finite-element edges. Refine the
mesh around the hole and compare at least three meshes. The scripts use
consistent node-force totals on the outer boundaries for the shear case,
because a tangential traction cannot be represented by the ordinary PRES
surface-load label on a 2-D solid edge. For a publication-quality result,
replace this with surface traction elements or a verified SFE/SFL formulation
in the ANSYS version being used.

Partition recommendation
------------------------

For the actual circular-hole answer, use q1_full_quad.txt or q2_quarter_quad.txt
with a sufficiently small ESIZE and local sizing near the hole. The square-ring
files demonstrate the partition topology only; do not compare their numerical
stress concentration directly with the circular-hole analytical solution.
The best production partition is a four-block O-grid whose inner edges are
true circular arcs. That requires creating four arc lines with LARC and is
preferable to a square opening. Use the diagonal case only as a triangular
mesh comparison.
