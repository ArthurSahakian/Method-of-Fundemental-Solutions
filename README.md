# Method of Fundamental Solution for the Biharmonic Steklov Problem

This librairy computes the numerical approximations of the Steklov Biharmonic eigenvalues via the Method of Fundamental Solutions (MFS), used in the thesis Biharmonic eigenvalue problems (click [here](https://www.researchgate.net/publication/379825121_Polyharmonic_eigenvalue_problems_of_the_steklov_type)).

##  Description

This project computes the eigenvalues of two Biharmonic eigenvalue problems of the Steklov kind; the classical Steklov problem and the Variant one. The eigenvalues of these problems vary according to the domain of the object. 

The script **Dimension_two_class** treats domains that are of dimension two. In particular, it treats smooth domains (boundaries of class C2 at least).

## Mathematical Background 

###  Support Function — Definition

Let \\( D \subset \mathbb{R}^2 \\) be a smooth, convex domain. The **support function** of \\( D \\) is defined as:

\[
h_D(\theta) = \max_{x \in D} \langle x, u(\theta) \rangle = \max_{x \in D} (x_1 \cos\theta + x_2 \sin\theta)
\]

where \\( u(\theta) = (\cos\theta, \sin\theta) \\) is the unit vector in the direction \\( \theta \\).  
This function returns the distance from the origin to the tangent line of **D** in the direction \\( \theta \\).

---

###  Fourier Approximation of the Support Function

Because \\( h_D(\theta) \\) is a \\( 2\pi \\)-periodic function, it can be approximated by a truncated Fourier series:

\[
h_D(\theta) \approx a_0 + \sum_{n=1}^{N} \left( a_n \cos(n\theta) + b_n \sin(n\theta) \right)
\]

This compact representation is useful for numerical shape optimization and inverse problems.

---

### Example: Support Function of a Circle

Consider the unit circle of radius \\( R \\) centered at the origin. Its support function is constant:

\[
h_D(\theta) = R
\]

And its Fourier series contains only the \\( a_0 \\) term:

\[
h_D(\theta) = R + 0\cdot \cos(n\theta) + 0\cdot \sin(n\theta)
\]

---
