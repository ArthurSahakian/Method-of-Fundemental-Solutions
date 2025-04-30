# Method of Fundamental Solution for the Biharmonic Steklov Problem

This librairy computes the numerical approximations of the Steklov Biharmonic eigenvalues via the Method of Fundamental Solutions (MFS), used in the thesis Biharmonic eigenvalue problems (click [here](https://www.researchgate.net/publication/379825121_Polyharmonic_eigenvalue_problems_of_the_steklov_type)).

##  Description

This project computes the eigenvalues of two Biharmonic eigenvalue problems of the Steklov kind; the classical Steklov problem and the Variant one. The eigenvalues of these problems vary according to the domain of the object. 

The script **Dimension_two_class** treats domains that are of dimension two. In particular, it treats smooth domains (boundaries of class C2 at least).

## Mathematical Background 

### Support Function — Definition

Let **D** ⊂ ℝ² be a smooth, convex domain. The **support function** of **D** is defined as:

```
h_D(θ) = max_{x ∈ D} <x, u(θ)> = max_{x ∈ D} (x₁ cos(θ) + x₂ sin(θ))
```

where **u(θ) = (cos(θ), sin(θ))** is the unit vector in the direction **θ**.  
This function returns the distance from the origin to the tangent line of **D** in the direction **θ**.

---

### Fourier Approximation of the Support Function

Because **h_D(θ)** is a **2π**-periodic function, it can be approximated by a truncated Fourier series:

```
h_D(θ) ≈ a₀ + Σ (ak cos(nθ) + bₙ sin(nθ)) from n=1 to K
```

This compact representation is useful for numerical shape optimization and inverse problems.

---

### Example: Support Function of a Circle

Consider the unit circle of radius **R** centered at the origin. Its support function is constant:

```
h_D(θ) = R
```

And its Fourier series contains only the **a₀** term and can be written as 

```
(1,0,0,0,...0)
```

### Method of Fundamental Solutions (MFS)

The **Method of Fundamental Solutions (MFS)** is a numerical technique for solving boundary value problems (BVPs) by using source points both on the boundary of the domain and on an artificial boundary slightly outside the domain.

#### Key Steps:
1. **Source Points on Boundary and Artificial Boundary**:  
   We place source points on the boundary of the domain as well as on an artificial boundary that lies just outside the domain.
   
2. **Estimation of Eigenfunctions**:  
   The method estimates the eigenfunctions of the problem by using the source points on the artificial boundary.
   
3. **Computation of Eigenvalues**:  
   The eigenvalues are computed based on the eigenfunctions obtained from the source points on the artificial boundary.

This approach allows for an efficient estimation of the solution, particularly for problems involving eigenvalue computations, as it avoids directly solving the differential equation in the interior of the domain.

## Use Case 

The script **Dimension_two_class** solves the the Steklov and Variant problems for domains in dimension 2.

### Importation of the class and initialization

First we will import initialize the class by feeding it two arguments:

- N: The number of points to be considered on the boundary of the domain **D** and the artificial boundary; default is 80.
- beta: The distance between the boundary of **D** and the artificial boundary associated with **D**, default is 0.4.

```python
from Dimension_two_class import dimension_2
mfs = dimension_2(80, 0.4)
```

### Fitting of the domain

Let us consider the example of the circle and let us set K=3, where K is the order of the Fourrier approximation. Then the circle can be represented by the vector (1,0,0,0,0,0,0).

```python
mfs.fourrier_fit(np.array([1,0,0,0.0,0]))
```

### Computation of the Eigenvalues

Finally, the user chooses one of two options, either "steklov" or "variant", to compute the eigenvalues of the associated problem. 

```python
mfs.real_eigenvalues(problem = "steklov")
```

#### `real_eigenvalue` Function

This function computes the approximate eigenvalues of a problem using the Method of Fundamental Solutions. It can solve either the original or a variant problem and supports different eigenvalue computation methods.

##### Parameters:
- `i` (integer, optional):  
  The rank of the eigenvalue to compute. Default is `0`.

- `method` (string, optional):  
  The method for solving the eigenvalue problem. Options are:
  - `"direct"`: Direct computation.
  - `"inverse"`: Inverse method.
  - `"lu"`: LU decomposition.  
  Default is `"direct"`.

- `problem` (string, optional):  
  Defines the problem to solve. Options are:
  - `"steklov"`: The original Steklov problem.
  - `"variant"`: The variant of the problem.  
  Default is `"steklov"`.

##### Returns:
- `first_eigen` (float):  
  The eigenvalue of the given rank `i`.

- `positive` (list of floats):  
  A sorted list of all the eigenvalues of the generalized problem, only including the positive eigenvalues.

##### Function Overview:
This function generates angles, computes matrices based on the selected problem (`steklov` or `variant`), and calculates the eigenvalues using one of the specified methods (`direct`, `inverse`, `lu`, or QR). It returns the eigenvalue of the given rank `i` and a list of all positive eigenvalues.

The method of calculation (e.g., LU decomposition, inverse, or direct) determines how the eigenvalues are computed from the generalized problem matrices.




