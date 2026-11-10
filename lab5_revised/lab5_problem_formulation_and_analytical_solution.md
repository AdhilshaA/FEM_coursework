# Lab 5: Cartesian formulation and analytical solutions

This document provides report-ready formulation text for the two plate-with-hole
problems. The FEM model uses N-mm-MPa units:

| Quantity | Value |
|---|---:|
| Plate length | \(L=300\ {\rm mm}\) |
| Plate height | \(H=300\ {\rm mm}\) |
| Hole radius | \(a=10\ {\rm mm}\) |
| Young's modulus | \(E=200000\ {\rm MPa}\) |
| Poisson's ratio | \(\nu=0.30\) |
| Applied traction | \(S=T=100\ {\rm MPa}\) |
| Thickness | \(t=1\ {\rm mm}\) |

The analytical expressions below are the classical **infinite-plate
Kirsch/Airy solutions**. The FEM plate is finite, so agreement is expected to
improve with mesh refinement but is not expected to be exact.

## 1. Common Cartesian formulation

Let \(u(x,y)\) and \(v(x,y)\) be the displacement components in the \(x\) and
\(y\) directions. Under plane stress,

\[
\boldsymbol{\varepsilon} =
\begin{bmatrix}
\varepsilon_x\\
\varepsilon_y\\
\gamma_{xy}
\end{bmatrix}
=
\begin{bmatrix}
\partial u/\partial x\\
\partial v/\partial y\\
\partial u/\partial y+\partial v/\partial x
\end{bmatrix}.
\]

The Cartesian equilibrium equations without body force are

\[
\frac{\partial\sigma_x}{\partial x}
+\frac{\partial\tau_{xy}}{\partial y}=0,
\qquad
\frac{\partial\tau_{xy}}{\partial x}
+\frac{\partial\sigma_y}{\partial y}=0.
\]

For plane stress,

\[
\begin{bmatrix}
\sigma_x\\
\sigma_y\\
\tau_{xy}
\end{bmatrix}
=
\frac{E}{1-\nu^2}
\begin{bmatrix}
1&\nu&0\\
\nu&1&0\\
0&0&(1-\nu)/2
\end{bmatrix}
\begin{bmatrix}
\varepsilon_x\\
\varepsilon_y\\
\gamma_{xy}
\end{bmatrix}.
\]

The hole boundary is traction-free:

\[
\boldsymbol{\sigma}\boldsymbol{n}=0
\quad\text{on }r=a.
\]

The outer boundary is sufficiently far from the hole that the applied
tractions approximate the remote analytical stress state.

The Cartesian-to-polar coordinates are

\[
r=\sqrt{x^2+y^2},\qquad
\theta=\tan^{-1}\left(\frac{y}{x}\right).
\]

For comparison with the FEM Cartesian stresses, use

\[
\sigma_x
=\sigma_{rr}\cos^2\theta
+\sigma_{\theta\theta}\sin^2\theta
-2\tau_{r\theta}\sin\theta\cos\theta,
\]

\[
\sigma_y
=\sigma_{rr}\sin^2\theta
+\sigma_{\theta\theta}\cos^2\theta
+2\tau_{r\theta}\sin\theta\cos\theta,
\]

\[
\tau_{xy}
=(\sigma_{rr}-\sigma_{\theta\theta})\sin\theta\cos\theta
+\tau_{r\theta}(\cos^2\theta-\sin^2\theta).
\]

## 2. Question 1: plate under pure shear

### 2.1 Cartesian problem statement

The remote Cartesian stress tensor is

\[
\boldsymbol{\sigma}^{\infty}
=
\begin{bmatrix}
0&S\\
S&0
\end{bmatrix},
\qquad S=100\ {\rm MPa}.
\]

Therefore, on the four outer edges, the applied traction is

\[
\boldsymbol{t}=\boldsymbol{\sigma}^{\infty}\boldsymbol{n}.
\]

For the edges of the full plate:

| Boundary | Outward normal \(\boldsymbol n\) | Applied traction |
|---|---|---|
| \(x=+150\) | \((1,0)\) | \((0,S)\) |
| \(x=-150\) | \((-1,0)\) | \((0,-S)\) |
| \(y=+150\) | \((0,1)\) | \((S,0)\) |
| \(y=-150\) | \((0,-1)\) | \((-S,0)\) |

The FEM model must remove rigid-body motion without introducing artificial
stress. A convenient constraint set is fixing both displacements at one
corner and fixing one displacement at a second corner.

### 2.2 Analytical infinite-plate solution

For remote pure shear \(S\), the Airy/Kirsch solution in polar coordinates is

\[
\sigma_{rr}
=-S\left(1-\frac{4a^2}{r^2}+\frac{3a^4}{r^4}\right)\sin 2\theta,
\]

\[
\sigma_{\theta\theta}
=S\left(1+\frac{3a^4}{r^4}\right)\sin 2\theta,
\]

\[
\tau_{r\theta}
=S\left(1+\frac{2a^2}{r^2}-\frac{3a^4}{r^4}\right)\cos 2\theta.
\]

At the hole boundary \(r=a\):

\[
\sigma_{rr}=0,\qquad
\tau_{r\theta}=0,\qquad
\sigma_{\theta\theta}=4S\sin 2\theta.
\]

Thus, the maximum circumferential normal stress is

\[
\left|\sigma_{\theta\theta}\right|_{\max}=4S=400\ {\rm MPa}.
\]

The maximum in-plane shear stress is

\[
\tau_{\max}
=\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}.
\]

At the hole boundary, the maximum value is

\[
\tau_{\max}=2S=200\ {\rm MPa},
\qquad
K_t^{(\text{in-plane shear})}
=\frac{200}{100}=2.
\]

The circumferential normal-stress concentration is instead

\[
K_t^{(\text{circumferential normal})}
=\frac{4S}{S}=4.
\]

These are different concentration factors and should not be mixed in the
report.

### 2.3 Section A-A

For the positive \(x\)-axis, \(y=0\), \(\theta=0\), and \(r=x\). The analytical
shear stress becomes

\[
\tau_{xy}(x,0)
=S\left(1+\frac{2a^2}{x^2}-\frac{3a^4}{x^4}\right),
\qquad x\ge a.
\]

The boundary values are

\[
\tau_{xy}(a,0)=0,\qquad
\tau_{xy}(\infty,0)=S.
\]

The FEM section plot should be compared with this curve on the material
segments \(x\le -a\) and \(x\ge a\), not through the hole.

## 3. Question 2: plate under uniaxial tension

### 3.1 Cartesian problem statement

The remote stress tensor is

\[
\boldsymbol{\sigma}^{\infty}
=
\begin{bmatrix}
T&0\\
0&0
\end{bmatrix},
\qquad T=100\ {\rm MPa}.
\]

Only the first quadrant is modeled:

\[
0\le x\le150\ {\rm mm},\qquad
0\le y\le150\ {\rm mm},\qquad
r\ge a.
\]

The symmetry conditions are

\[
u=0\quad\text{on }x=0,
\qquad
v=0\quad\text{on }y=0.
\]

The right boundary \(x=150\) mm carries the tensile traction

\[
\boldsymbol t=(T,0).
\]

The upper boundary and the hole boundary are traction-free, apart from the
symmetry constraints on \(x=0\) and \(y=0\).

### 3.2 Analytical infinite-plate solution

For remote tension \(T\) in the \(x\)-direction, the Kirsch solution is

\[
\sigma_{rr}
=\frac{T}{2}\left(1-\frac{a^2}{r^2}\right)
+\frac{T}{2}\left(1-\frac{4a^2}{r^2}+\frac{3a^4}{r^4}\right)\cos2\theta,
\]

\[
\sigma_{\theta\theta}
=\frac{T}{2}\left(1+\frac{a^2}{r^2}\right)
-\frac{T}{2}\left(1+\frac{3a^4}{r^4}\right)\cos2\theta,
\]

\[
\tau_{r\theta}
=-\frac{T}{2}\left(1+\frac{2a^2}{r^2}-\frac{3a^4}{r^4}\right)\sin2\theta.
\]

At \(r=a\):

\[
\sigma_{rr}=0,\qquad
\tau_{r\theta}=0,
\]

\[
\sigma_{\theta\theta}=T(1-2\cos2\theta).
\]

The maximum occurs at \(\theta=\pi/2\), at the top and bottom of the hole:

\[
\sigma_{\theta\theta,\max}=3T=300\ {\rm MPa}.
\]

Therefore:

\[
K_t=\frac{\sigma_{\max}}{T}=3.
\]

In the first-quadrant model, the corresponding point is the top of the
quarter-hole, \((x,y)=(0,a)\). A mesh that does not resolve this location will
underestimate the peak stress.

### 3.3 Cartesian section comparison

Along the positive \(x\)-axis, \(y=0\), \(\theta=0\), \(r=x\):

\[
\sigma_x(x,0)
=\frac{T}{2}\left[
2-\frac{a^2}{x^2}-\frac{3a^4}{x^4}
\right],
\qquad x\ge a,
\]

\[
\sigma_y(x,0)
=\frac{T}{2}\left[
\frac{a^2}{x^2}+\frac{3a^4}{x^4}
\right],
\]

\[
\tau_{xy}(x,0)=0.
\]

Consequently:

\[
\sigma_x(a,0)=0,\qquad
\sigma_y(a,0)=2T,\qquad
\sigma_x(\infty,0)=T.
\]

The exact shear stress on this symmetry axis is zero. Therefore, if the lab
statement asks for shear stress along \(x=0\) or \(y=0\), a zero or nearly
zero FEM result is expected by symmetry; the normal-stress plot is the more
informative comparison.

## 4. FEM-to-analytical comparison quantities

For each mesh, report:

| Quantity | Question 1 | Question 2 |
|---|---|---|
| Remote traction | \(S=100\) MPa | \(T=100\) MPa |
| Peak measure | maximum in-plane shear | maximum principal/normal stress |
| Analytical peak | \(2S=200\) MPa for in-plane shear | \(3T=300\) MPa |
| Analytical \(K_t\) | \(2\) for in-plane shear | \(3\) |
| Section | \(\tau_{xy}(x,0)\) | \(\sigma_x(x,0)\) |
| Hole boundary | traction-free | traction-free |

Use the same stress definition in Python and APDL. Do not compare Python
\(\tau_{xy}\) directly with APDL maximum principal stress. Also distinguish
element stresses, nodal-averaged stresses, and extrapolated integration-point
stresses, because their peak values can differ even for the same mesh.

## 5. Limitations to state in the report

1. The analytical solution assumes an infinite plate; the FEM plate is finite.
2. A linear CST element has constant stress within each element.
3. A polygonal approximation of the hole changes the exact local geometry.
4. APDL nodal-averaged or extrapolated stresses may be smoother than raw CST
   element stresses.
5. Peak stress should be assessed only after mesh and hole-boundary refinement.
6. A converged FEM value may approach the analytical value without matching it
   exactly because the outer boundary is only 15 hole radii away.
