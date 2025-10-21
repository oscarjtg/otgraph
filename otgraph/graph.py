import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import netCDF4 as nc

class CuboidPlot:
    def __init__(self, x_top, y_top, arr_top, x_front, z_front, arr_front, y_side, z_side, arr_side):
        self.x_top = x_top
        self.y_top = y_top
        self.arr_top = arr_top
        self.x_front = x_front
        self.z_front = z_front
        self.arr_front = arr_front
        self.y_side = y_side
        self.z_side = z_side
        self.arr_side = arr_side

    #########
    ## Colourmaps for the three faces.
    #########
    cmap_front = "Purples"
    cmap_top = "BrBG"
    cmap_side = "RdBu_r"

    #########
    ## Factor by which to multiply the true vmax if plotting anomalies centred on zero.
    #########
    vmax_factor_top = 1.0
    vmax_factor_front = 1.0
    vmax_factor_side = 1.0

    #########
    ## Set colorbar limits such that the central colour of the colormap
    ## is centred on zero.
    ## This option is appropriate for plotting data representing residuals or anomalies.
    #########
    centred_front = True
    centred_side = True
    centred_top = True

    #########
    ## The labels for each of the colorbars.
    #########
    cbar_label_top = "cbar_label_top"
    cbar_label_front = "cbar_label_front"
    cbar_label_side = "cbar_label_side"

    ########
    ## Set the title
    ########
    plot_title = "plot title"

    def draw_plot(self, plot_L=12, plot_W=2, plot_H=2, alpha_deg=45):
        """
        Draws plot. Plot is shaped like a cuboid. 
        
        On the top face of the cuboid, we plot the 2d numpy array called `arr_top` 
        on a grid defined by the 1d numpy arrays called `x_top` and `y_top`.

        The arrays plotted on the front and side (right hand side) faces of the cuboid
        follow a similar naming convention, but with `front` or `side` instead of `top`.

        Parameters
        plot_L (float)
            The horizontal length of the cuboid drawn on the figure

        plot_W (float)
            The diagnoal (projection of out-of-plane direction) length 
            of the cuboid drawn on the figure

        plot_H (float)
            The vertical height of the cuboid drawn on the figure

        alpha_deg (float)
            The angle that the diagonal makes with the vertical

        Returns
        -------
        matplotlib.pyplot.figure
            The figure object of the plot that has been drawn.
        
        """
        # Find dimensions of plot from the input arrays.
        Lx = self.x_front[-1] - self.x_front[0]
        Ly = self.y_side[-1] - self.y_side[0]
        Lz = abs(self.z_front[-1] - self.z_front[0])

        def transform_coordinates(x, y, z, X, Y, Z, L=plot_L, W=plot_W, H=plot_H, alpha_deg=alpha_deg):
            """
            Projects 3D coordinates (x, y, z) in a 3D cuboid domain of size (X, Y, Z)
            onto a 2D plane that looks like a cuboid of length L, width W, and height H.

            In other words, this works out how to plot 3D data on a cuboid.

            Parameters
            ---------
            x (float)
                x-coordinate in 3D space

            y (float)
                y-coordinate in 3D space

            z (float)
                z-coordinate in 3D space

            X (float)
                Size of 3D domain

            Y (float)
                Size of 3D domain

            Z (float)
                Size of 3D domain

            L (float)
                Length of cuboid (horizontal direction) drawn on the plane, in inches

            W (float)
                Width of cuboid (diagonal direction) drawn on plane, in inches

            H (float)
                Height of cuboid (vertical direction) drawn on plane, in inches

            alpha_deg (float)
                Angle that the diagonal (out-of-plane) direction makes on the page
                relative to the vertical

            Returns
            -------
            (float, float)
                Tuple of floats giving the x and y coordinate on the plane, in inches.
            
            """
            alpha = alpha_deg * np.pi / 180.0
            sina = np.sin(alpha)
            tana = np.tan(alpha)
            xt = x * (L/X) + (y + Y/2)/Y * (W / sina)
            yt = z * (H/Z) + (y + Y/2)/Y * (W / tana)
            return xt, yt
        
        # Desired pixel dimensions divisible by 16
        width_px = 1504  # divisible by 16
        height_px = 608  # divisible by 16
        dpi = 150

        # Convert to inches
        width_in = width_px / dpi
        height_in = height_px / dpi
        
        # Create figure.
        fig, ax = plt.subplots(figsize=(width_in, height_in))

        # Plot top side.
        X_top, Y_top = np.meshgrid(self.x_top, self.y_top)
        z_top = 0.0
        XX_top, YY_top = transform_coordinates(X_top, Y_top, z_top, Lx, Ly, Lz)
        if self.centred_top:
            vmax = self.vmax_factor_top * np.abs(np.min(self.arr_top))
            pmesh_top = ax.pcolormesh(XX_top, YY_top, self.arr_top, cmap=self.cmap_top, vmin=-vmax, vmax=vmax)
        else:
            pmesh_top = ax.pcolormesh(XX_top, YY_top, self.arr_top, cmap=self.cmap_top)

        #plt.colorbar(pmesh_top, label=self.cbar_label_top, fraction=0.03, pad=0.02)
        #fig.colorbar(pmesh_top, ax=ax, orientation='horizontal', label=self.cbar_label_top)#, #fraction=0.046, pad=0.04)

        # Plot front side.
        X_front, Z_front = np.meshgrid(self.x_front, self.z_front)
        y_front = -Ly/2
        XX_front, YY_front = transform_coordinates(X_front, y_front, Z_front, Lx, Ly, Lz)
        if self.centred_front:
            umax = self.vmax_factor_front * np.max(np.abs(self.arr_front))
            pmesh_front = ax.pcolormesh(XX_front, YY_front, self.arr_front, cmap=self.cmap_front, vmin=-umax, vmax=umax)
        else:
            pmesh_front = ax.pcolormesh(XX_front, YY_front, self.arr_front, cmap=self.cmap_front)
        #plt.colorbar(pmesh_front, label=self.cbar_label_front, fraction=0.03, pad=0.02)
        #fig.colorbar(pmesh_front, ax=ax, orientation='horizontal', label=self.cbar_label_front)#, fraction=0.046, pad=0.04)

        # Plot side.
        Y_side, Z_side = np.meshgrid(self.y_side, self.z_side)
        XX_side, YY_side = transform_coordinates(Lx, Y_side, Z_side, Lx, Ly, Lz)
        if self.centred_side:
            umax = self.vmax_factor_side * np.max(np.abs(self.arr_side))
            pmesh_side = ax.pcolormesh(XX_side, YY_side, self.arr_side, cmap=self.cmap_side, vmin=-umax, vmax=umax)
        else:
            pmesh_side = ax.pcolormesh(XX_side, YY_side, self.arr_side, cmap=self.cmap_side)
        #plt.colorbar(pmesh_side, label=self.cbar_label_side, fraction=0.03, pad=0.02)
        #fig.colorbar(pmesh_side, ax=ax, orientation='horizontal', label=self.cbar_label_side)#, fraction=0.046, pad=0.04)

        ax.set_aspect("equal")
        plt.axis("off")
        ax.plot([0, plot_L], [-plot_H, -plot_H], color="black", lw=1)
        ax.plot([0, plot_L], [0, 0], color="black", lw=1)
        XX, YY = transform_coordinates(np.array([0, Lx]), Ly/2, np.array([0, 0]), Lx, Ly, Lz)
        ax.plot(XX, YY, color="black", lw=1)
        ax.plot([0, 0], [-plot_H, 0], color="black", lw=1)
        ax.plot([plot_L, plot_L], [-plot_H, 0], color="black", lw=1)
        XX, YY = transform_coordinates(np.array([Lx, Lx]), Ly/2, np.array([-Lz, 0]), Lx, Ly, Lz)
        ax.plot(XX, YY, color="black", lw=1)
        XX, YY = transform_coordinates(0, np.array([-Ly/2, Ly/2]), 0, Lx, Ly, Lz)
        ax.plot(XX, YY, color="black", lw=1)
        XX, YY = transform_coordinates(Lx, np.array([-Ly/2, Ly/2]), 0, Lx, Ly, Lz)
        ax.plot(XX, YY, color="black", lw=1)
        XX, YY = transform_coordinates(Lx, np.array([-Ly/2, Ly/2]), -Lz, Lx, Ly, Lz)
        ax.plot(XX, YY, color="black", lw=1)

        # Add text
        downshift = 0.3
        ax.text(0, -plot_H-downshift, "0")
        ax.text(plot_L, -plot_H-downshift, f"{Lx}")

        ax.text(-1.5, -plot_H, f"{Lz}")

        ax.text(1, plot_H, f"{Ly}")

        ax.set_title(self.plot_title, pad=-20)

        
        # Top colorbar
        cax_top = inset_axes(ax, width="50%", height="5%", loc="upper center", borderpad=-2.2)
        cbar_top = fig.colorbar(pmesh_top, cax=cax_top, orientation='horizontal', label=self.cbar_label_top)
        cbar_top.ax.xaxis.set_ticks_position('top')
        cbar_top.ax.xaxis.set_label_position('top')


        # Front colorbar
        cax_front = inset_axes(ax, width="50%", height="5%", loc="lower left", borderpad=-2.2)
        fig.colorbar(pmesh_front, cax=cax_front, orientation='horizontal', label=self.cbar_label_front)

        # Side colorbar
        cax_side = inset_axes(ax, width="25%", height="5%", loc="lower right", borderpad=-2.2)
        fig.colorbar(pmesh_side, cax=cax_side, orientation='horizontal', label=self.cbar_label_side)


        #plt.tight_layout()
        return fig
    


def plot_vertical_profile_over_time_3d(file_path, time_str, z_str, var_str, xlabel_str, ylabel_str, cbar_str, x_str, x_val, subtract_initial, diverging):
    """
    Plots a heatmap depicting how the vertical profile of a variable evolves with time.
    THe vertical profiles are taken at a given x-value, averaged across the y direction, at every time step.

    Parameters
    ----------
    file_path (string)
        A string that gives the absolute or relative path to the NetCDF file containing the data.

    time_str (string)
        A string that gives the name of the "time" variable in the NetCDF file.

    z_str (string)
        A string that gives the name of the "z" variable in the NetCDF file.

    var_str (string)
        A string that gives the name of the variable whose vertical profiles are to be plotted over time,
        as given in the NetCDF file.

    xlabel_str (string)
        A string giving the label for the horizontal axis on the plot.

    ylabel_str (string)
        A string giving the label for the vertical axis on the plot.

    cbar_str (string)
        A string giving the label for the colorbar.

    x_str (string)
        A string giving the name of the "x" variable in the NetCDF file.

    x_val (float)
        A float giving the x-value at which the vertical profiles are to be sampled.

    subtract_initial (bool)
        A boolean value that says whether or not to subtract the profile at t=0 from all profiles.
        If subtract_initial is True, the t=0 profile is subtracted from all the profiles.
        Otherwise, this step is not carried out.

    diverging (bool)
        A boolean value that says whether or not the colorscale should be centred symmetrically on zero 
        and a diverging colormap should be used.
        If divering is True, the vmin and vmax parameters of plt.pcolormesh() are set symmetrically about zero
        as the largest absolute value from all the vertical profiles, and the diverging colormap "RdBu_r" is used.
        If divering is False, the default vmin, vmax, and colorbar are used.
        TODO: add customisation of colorbar.

    Returns
    -------
    (matplotlib.pyplot.figure, matplotlib.pyplot.axis)

    """
    with nc.Dataset(file_path, "r") as ncfile:
        x = np.array(ncfile.variables[x_str][:])
        x_idx = np.argmax(x > x_val)
        time = np.array(ncfile.variables[time_str][:])
        depth = np.array(ncfile.variables[z_str][:])
        var = np.mean(np.array(ncfile.variables[var_str][:, :, :, x_idx]), axis=2)
        if subtract_initial:
            var_init = np.mean(np.array(ncfile.variables[var_str][0, :, :, x_idx]), axis=1)
            var = var - var_init[np.newaxis, :]
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    var = np.transpose(var)
    if diverging:
        vmax = np.max(np.abs(var))
        pmesh = ax.pcolormesh(time, depth, var, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    else:
        pmesh = ax.pcolormesh(time, depth, var)
    cbar = plt.colorbar(pmesh, ax=ax)
    cbar.set_label(cbar_str, fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xlabel(xlabel_str, fontsize=16)
    ax.set_ylabel(ylabel_str, fontsize=16)
    ax.tick_params(labelsize=14)
    return fig, ax

