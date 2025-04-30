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

And its Fourier series contains only the **a₀** term and can be written as (1,0,0,0,...0).

```

```

## Use Case 
