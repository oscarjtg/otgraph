# otgraph
Python matplotlib wrappers

## Instructions

```python
import otgraph.gragh as otg
import otgraph.video as vid
```

## Uses

Use the `CuboidPlot` class to plot data on the faces of a cuboid, 
oriented such that you can see the top, front, and right hand side of the cuboid.

The `CuboidPlot` instance is initialised with arrays representing the data 
to be plotted on each face.

* The top face is a slice through the x-y plane. 
On the top face, the data in the 2d numpy array `arr_top`
will be plotted on a grid delimited by the 1d numpy arrays `x_top` and `y_top`.

* The front face is a slice through the x-z plane.
On the front face, the data in the 2d numpy array `arr_front`
will be plotted on a grid delimited by the 1d numpy arrays `x_front` and `z_front`.

* The side face is a slice through the y-z plane.
On the side face, the data in the 2d numpy array `arr_side`
will be plotted on a grid delimited by the 1d numpy arrays `y_side` and `z_side`.

```python
# Assuming the appropriate arrays are already defined...
oc = otg.CuboidPlot(x_top=x_top, 
                    y_top=y_top, 
                    arr_top=arr_top, 
                    x_front=x_front, 
                    z_front=z_front, 
                    arr_front=arr_front, 
                    y_side=y_side, 
                    z_side=z_side, 
                    arr_side=arr_side)

# Change the colormaps
oc.cmap_top = "RdBu_r"
oc.cmap_front = "RdBu_r"
oc.cmap_side = "RdBu_r"

# Set the colorbar label
oc.cbar_label_top = r"vorticity  [s$^{-1}$]"
# etc...

# Call the draw_plot() method to draw the plot. Returns a matplotlib.pyplot.figure.
oc.draw_plot()

# Save the figure using usual matplotlib syntax.
plt.savefig("plots/cuboidplot.png", dpi=150)
```

## TODO
- CuboidPlots cosmetic tweaks, e.g.
    - better way to show physical dimensions of domain (e.g. arrows)
    - better placement of colorbars and title(s)
