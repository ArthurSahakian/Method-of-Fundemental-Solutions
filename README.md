This librairy computes the numerical approximations of the Steklov Biharmonic eigenvalues via the Method of Fundamental Solutions (MFS).
The librairy is the code used in the thesis Biharmonic eigenvalue problems ( click [here](https://www.researchgate.net/publication/379825121_Polyharmonic_eigenvalue_problems_of_the_steklov_type)).

# Project Title

A short description of your project in one sentence.

## 📘 Description

This project solves the following mathematical problem:

> **[Insert Problem Statement]**

For example, we aim to minimize the following loss function:

![Loss Function](https://latex.codecogs.com/png.image?\dpi{120}&space;L(\theta)=\frac{1}{n}\sum_{i=1}^n(y_i-f_\theta(x_i))^2)

The algorithm is implemented in Python and relies on [insert key libraries, e.g., NumPy, PyTorch, etc.].

## 🧠 Mathematical Formulation

Let \\( x \in \mathbb{R}^n \\) be the input features and \\( y \in \mathbb{R} \\) the target. The model seeks to find parameters \\( \theta \\) that minimize:

```math
\min_{\theta} \; \frac{1}{n} \sum_{i=1}^n \left( y_i - f_\theta(x_i) \right)^2
