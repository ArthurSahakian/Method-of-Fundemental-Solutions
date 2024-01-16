# Introduction: #
This module covers the functionalities of the library SelectPoint, used to sample equidistant points on the boundary of a 2D domain or a 3D domain.
The boundaries of closed 2D domains are arcs, and the boundary of 3D domains are surfaces.
## Examples of boundaries: ##
- The circle is the boundary of the disk.
- The 6 faces (squares) form the boundary of the associated cube.

# Main Functionalities of the library: #
## Dimension 2 ##
The domains in dimension two, or the closed arcs, have to be represented by their support function, and ,ore precisely by the Fourier coefficients of the truncated Fourier series associated with their support function. 
###  `point_selector_dimension_two ` Class

A class for point selection in two dimensions.

#### Constructor

```python
obj= =point_selector_dimension_two()
```
### fit_fourier:

Fit a Fourier series to the given data.

#### Parameters

- `N` (int): Number of terms in the Fourier series.
- `fourier_coeff` (list): Fourier coefficients for the series. Should contain an odd number (2M+1) of real numbers of the form (a_0, a_1, ..., a_M, b_1, ..., b_M).
- `domain` (list): Domain of the data, typically [start_value, end_value].

#### Usage Example

```python
obj.fit_fourier(N=10, fourier_coeff=[...], domain=[start_value, end_value])
```

###  `angles_gen_arc ` Method

Generates equidistant angles on an arc, given the parametric equation of its boundary.
By default, the domain of definition of the function is [0, 2pi].

#### Parameters

- `p` (array): Fourier coefficients of the boundary function.
- `domain` (tuple or list): Domain of definition of the function. Default is [0, 2pi].

#### Returns

- `angles` (array): Angles equisitantly repartitioned along the arc.

#### Usage Example

```python
result = obj.angles_gen_arc(p=[...], domain=[start_value, end_value])
```

###  `two_dimensional_equidistant_points ` Method

Computes the boundary points of the domain enclosed by the arc `p`.
This function assumes that the domain is centered around the origin.

#### Parameters

- `include_normals` (boolean, optional): Computes and returns the exterior unit normals at the points `X` if set to True. The default is False.

#### Returns

- `X` (numpy array): Returns a 2D numpy array containing points on the boundary.
- `normals` (numpy array): Returns a 2D numpy array containing the normals vectors. (Returned only if `include_normals` is True)

#### Usage Example

```python
result = obj.two_dimensional_equidistant_points(include_normals=False)
```

###  `two_dimensional_plotte ` Method

Create a 2D plot based on the specified plot type.

#### Parameters

- `plot_type` (str, optional): Type of plot to create ("scatter" or "line"). Default is "scatter".

#### Usage Example

```python
obj.two_dimensional_plotter(plot_type="scatter")
```

## Dimension 3 ## 
Three Dimensional domains, or closed surfaces, can be represented in a lot of ways.
Currently this library supports two types of representations:
- Surfaces reprsented by a function F(theta,phi) = ( u(theta,phi), v(theta,phi), w(theta,phi))
- Surfaces represented by a radius function r(theta,phi) in spherical coordinates.

Remark: It is neccesarry to determine the first order derivatives of the input function as well.

###  `point_selector_dimension_three ` Class

A class for point selection and surface fitting in three dimensions.

#### Constructor

```python
obj = point_selector_dimension_three()
```
### `fit_parametric_funcs` Method

Fit parametric surface functions for a 3D surface.

#### Parameters

- `F` (callable): Parametric function representing the surface.
- `F_theta` (callable): Derivative of `f` with respect to theta.
- `F_phi` (callable): Derivative of `f` with respect to phi.
- `N` (int): The maximal number of points on the largest band of the domain.
- `inter_t` (tuple): Tuple representing the interval for parameter theta (start_t, end_t).
- `inter_s` (tuple): Tuple representing the interval for parameter phi (start_s, end_s).

#### Description

This method fits parametric surface functions for a 3D surface. The parameters `F`, `F_theta`, and `F_phi` are callables representing the parametric surface function, its derivative with respect to theta, and its derivative with respect to phi, respectively. The maximal number of points on the largest band of the domain.

The fitted functions are stored as attributes: `self.F`, `self.F_theta`, `self.F_phi`.

#### Usage Example

```python
obj.fit_parametric_funcs(F, F_theta, F_phi, N, inter_t=(0, 2*np.pi), inter_s=(0, np.pi))
```

### `fit_spherical_funcs` Method

Fit spherical functions for a 3D surface.

#### Parameters

- `r` (callable): The spherical function defining the radial distance.
- `r_theta` (callable): The derivative of the spherical function with respect to theta.
- `r_phi` (callable): The derivative of the spherical function with respect to phi.
- `N` (int): The maximal number of points on the largest band of the domain.
- `inter_t` (tuple): Tuple representing the interval for parameter theta (start_t, end_t).
- `inter_s` (tuple): Tuple representing the interval for parameter phi (start_s, end_s).


#### Description

This method fits spherical functions for a 3D surface based on the provided spherical coordinates. The parameters `r`, `r_theta`, and `r_phi` are callables representing the spherical function, its derivative with respect to theta, and its derivative with respect to phi, respectively. `N` signifies the maximal number of points on each band of the domain.

The fitted functions are stored as attributes: `self.r`, `self.r_theta`, `self.r_phi`.

#### Usage Example

```python
obj.fit_spherical_funcs(r, r_theta, r_phi, N, inter_t=(0, 2*np.pi), inter_s=(0, np.pi))
```

### `main_point_generator` Method

Main method for generating points and normals based on the specified tracker type.

#### Returns

- `X`: Generated points.
- `normals`: Generated normals.

#### Note

Ensure that the required functions (`point_generator`, `f`, `f_theta`, `f_phi`, `r`, `r_theta`, `r_phi`, `creator_r`, `creator_r_theta`, `creator_r_phi`) and interpolations (`inter_t`, `inter_s`) are defined in the class.

#### Raises

- `RuntimeError`: Raised if no function associated with the surface. Call any of the `fit_polar_funcs`, `fit_spherical_funcs`, or `fit_harmonics` methods before `main_initialization`.

#### Usage Example

```python
X, normals = obj.main_point_generator()
```

### `point_generator` Method

Generate 3D points on a parametric surface defined by functions `f`, `f_theta`, and `f_phi`.

#### Parameters

- `f` (callable): Parametric function representing the surface.
- `f_theta` (callable): Derivative of `f` with respect to theta.
- `f_phi` (callable): Derivative of `f` with respect to phi.
- `inter_t` (tuple): Tuple representing the interval for parameter theta (`start_t`, `end_t`).
- `inter_s` (tuple): Tuple representing the interval for parameter phi (`start_s`, `end_s`).

#### Returns

A tuple containing two numpy arrays:

- First array: 3D points representing the generated surface.
- Second array: Normalized normals at each generated point.

#### Description

This method generates equidistant points on a 3D surface defined by parametric functions. The parameters `f`, `f_theta`, and `f_phi` are callables representing the surface function, its derivative with respect to theta, and its derivative with respect to phi, respectively. `inter_t` and `inter_s` are tuples specifying the intervals for theta and phi, respectively.

The method returns a tuple containing two arrays:

- The first array represents the 3D points on the surface.
- The second array contains the normalized normals corresponding to each generated point.

#### Usage Example

```python
complete_pts, normals = obj.point_generator(f, f_theta, f_phi, inter_t, inter_s)
```

### `shape_plotter` Method

3D shape plotter method for visualizing the generated points.

#### Parameters

- `ann` (bool, optional): If `True`, annotate points with indices. Default is `False`.

#### Description

This method creates a 3D plot using the generated points. It uses the `matplotlib` library to visualize the points in a 3D space. If `ann` is set to `True`, each point is annotated with its index.

#### Usage Example

```python
obj.shape_plotter(ann=True)
```
