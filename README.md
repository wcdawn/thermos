# Thermos

Numerical solutions to the one-dimensional heat-conduction equation.

Geometries supported:
- Cartesian (slab)
- Cylindrical (equal-radius and equal-area)
- Spherical (equal-radius and equal-volume)

Numerical methods available:
- Finite Difference Method (FDM)
- Finite Element Method (FEM) (coming soon)

General non-linear heat conduction is supported of the form $k(T)$. A Picard iteration is used to resolve any non-linearity of the problem.

# Usage

To build.

```
git clone <url>
cd thermos
cmake -B build
cmake --build build
```

Then, to run.

```
cd thermos/cases/cyl_lin_klin_area
../../build/thermos.x ./cyl_lin_kklin_area.inp
```
