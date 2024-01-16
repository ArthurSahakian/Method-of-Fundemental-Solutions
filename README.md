# Introduction: #
This module covers the functionalities of the library SelectPoint, used to sample equidistant points on the boundary of a 2D domain or a 3D domain.
The boundaries of closed 2D domains are arcs, and the boundary of 3D domains are surfaces.
## Examples of boundaries: ##
- The circle is the boundary of the disk.
- The 6 faces (squares) form the boundary of the associated cube.

# Main Functionalities of the library: #
## Dimension 2 ##
The domains in dimension two, or the closed arcs, have to be represented by their support function, and ,ore precisely by the Fourier coefficients of the truncated Fourier series associated with their support function. 
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

### <span style="color:blue;">angles_gen_arc</span>:

Generates equidistant angles on an arc, given the parametric equation of its boundary.
By default, the domain of definition of the function is [0, 2pi].

#### Parameters

- `p` (array): Fourier coefficients of the boundary function.
- `domain` (tuple or list): Domain of definition of the function. Default is [0, 2pi].

#### Returns

- `angles` (array): Angles equisitantly repartitioned along the arc.

#### Usage Example

```python
result = obj.<span style="color:blue;">angles_gen_arc</span>(p=[...], domain=[start_value, end_value])
```

### <span style="color:blue;">two_dimensional_equidistant_points</span>:

Computes the boundary points of the domain enclosed by the arc `p`.
This function assumes that the domain is centered around the origin.

#### Parameters

- `include_normals` (boolean, optional): Computes and returns the exterior unit normals at the points `X` if set to True. The default is False.

#### Returns

- `X` (numpy array): Returns a 2D numpy array containing points on the boundary.
- `normals` (numpy array): Returns a 2D numpy array containing the normals vectors. (Returned only if `include_normals` is True)

#### Usage Example

```python
result = obj.<span style="color:blue;">two_dimensional_equidistant_points</span>(include_normals=False)
```

### <span style="color:blue;">two_dimensional_plotter</span>:

Create a 2D plot based on the specified plot type.

#### Parameters

- `plot_type` (str, optional): Type of plot to create ("scatter" or "line"). Default is "scatter".

#### Usage Example

```python
obj.<span style="color:blue;">two_dimensional_plotter</span>(plot_type="scatter")
```

## Dimension 3 ##

