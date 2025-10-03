import matplotlib.pyplot as plt
import numpy as np

class CuboidPlot:
    def __init__(self, x_top, y_top, arr_top, x_front, y_front, arr_front, x_side, y_side, arr_side):
        self.x_top = x_top
        self.y_top = y_top
        self.arr_top = arr_top
        self.x_front = x_front
        self.y_front = y_front
        self.arr_front = arr_front
        self.x_side = x_side
        self.y_side = y_side
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
        Ly = self.x_side[-1] - self.x_side[0]
        Lz = abs(self.y_front[-1] - self.y_front[0])

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

        plt.colorbar(pmesh_top, label=self.cbar_label_top)

        # Plot front side.
        X_front, Z_front = np.meshgrid(self.x_front, self.y_front)
        y_front = -Ly/2
        XX_front, YY_front = transform_coordinates(X_front, y_front, Z_front, Lx, Ly, Lz)
        if self.centred_front:
            umax = self.vmax_factor_front * np.max(np.abs(self.arr_front))
            pmesh_front = ax.pcolormesh(XX_front, YY_front, self.arr_front, cmap=self.cmap_front, vmin=-umax, vmax=umax)
        else:
            pmesh_front = ax.pcolormesh(XX_front, YY_front, self.arr_front, cmap=self.cmap_front)
        plt.colorbar(pmesh_front, label=self.cbar_label_front)

        # Plot side.
        Y_side, Z_side = np.meshgrid(self.x_side, self.y_side)
        XX_side, YY_side = transform_coordinates(Lx, Y_side, Z_side, Lx, Ly, Lz)
        if self.centred_side:
            umax = self.vmax_factor_side * np.max(np.abs(self.arr_side))
            pmesh_side = ax.pcolormesh(XX_side, YY_side, self.arr_side, cmap=self.cmap_side, vmin=-umax, vmax=umax)
        else:
            pmesh_side = ax.pcolormesh(XX_side, YY_side, self.arr_side, cmap=self.cmap_side)
        plt.colorbar(pmesh_side, label=self.cbar_label_side)

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
        ax.text(0, -plot_H-1, "0")
        ax.text(plot_L, -plot_H-1, f"{Lx}")

        ax.text(-1.5, -plot_H, f"{Lz}")

        ax.text(1, plot_H, f"{Ly}")

        plt.tight_layout()
        return fig
    


        