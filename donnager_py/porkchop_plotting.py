import numpy as np
import matplotlib.pyplot as plt

import donnager

class PorkchopPlotter:
    """
    Class Implementation for Porkchop Plot.
    """

    def __init__(
        self,
        departure_body,
        target_body,
        launch_span,
        arrival_span,
        ax=None,
        is_plot_tof=True,
        is_plot_arrive_v=True,
        max_c3=45.0,
        max_vhp=5
    ):
        """
        
        Parameters
        ----------
        departure_body: 
            Body from which departure is done

        target_body: 
            Body for targetting

        launch_span: 
            Time span for launch

        arrival_span: 
            Time span for arrival

        ax: matplotlib.axes.Axes
            For custom figures

        is_plot_tof: bool
            For plotting time flight contour lines

        vhp: bool
            For plotting arrival velocity contour lines

        max_c3: float
            Sets the maximum C3 value for porkchop

        max_vhp: float
            Sets the maximum arrival velocity for porkchop

        """
        self.departure_body = departure_body
        self.target_body = target_body
        self.launch_span = launch_span
        self.arrival_span = arrival_span
        self.ax = ax
        self.is_plot_tof = is_plot_tof
        self.is_plot_arrive_v = is_plot_arrive_v
        self.max_c3 = max_c3
        self.max_vhp = max_vhp

    def porkchop(self):
        """
        Plots porkchop between two bodies.
        """
        dv_launch, dv_arrival = donnager.interplan.calc_porkchop(
            self.departure_body,
            self.target_body,
            self.launch_span[np.newaxis, :],
            self.arrival_span[:, np.newaxis])

        if self.ax is None:
            fig, self.ax = plt.subplots(figsize=(15, 15))
        else:
            fig = self.ax.figure

        c3_levels = np.linspace(0, self.max_c3, 30)

        c = self.ax.contourf(
            [depart.to_datetime() for depart in self.launch_span],
            [arrive.to_datetime() for arrive in self.arrival_span],
            dv_launch**2,
            c3_levels)

        line = self.ax.contour(
            [depart.to_datetime() for depart in self.launch_span],
            [arrive.to_datetime() for arrive in self.arrival_span],
            dv_launch**2,
            c3_levels,
            colors="black",
            linestyles="solid")

        cbar = fig.colorbar(c)
        cbar.set_label("km2 / s2")
        self.ax.clabel(line, inline=1, fmt="%1.1f", colors="k", fontsize=10)

        if self.is_plot_tof:
            time_levels = np.linspace(100, 500, 5)

            tfl_contour = self.ax.contour(
                [depart.to_datetime() for depart in self.launch_span],
                [arrive.to_datetime() for arrive in self.arrival_span],
                tof,
                time_levels,
                colors="red",
                linestyles="dashed",
                linewidths=3.5)

            self.ax.clabel(
                tfl_contour, 
                inline=1, 
                fmt="%1.1f", 
                colors="r", 
                fontsize=14)

        if self.is_plot_arrive_v:
            vhp_levels = np.linspace(0, self.max_vhp, 5)

            vhp_contour = self.ax.contour(
                [depart.to_datetime() for depart in self.launch_span],
                [arrive.to_datetime() for arrive in self.arrival_span],
                dv_arrival,
                vhp_levels,
                colors="navy",
                linewidths=2.0)

            self.ax.clabel(
                vhp_contour, 
                inline=1, 
                fmt="%1.1f", 
                colors="navy", 
                fontsize=12)

        self.ax.grid()
        fig.autofmt_xdate()

        if not hasattr(self.target_body, "name"):
            self.ax.set_title(
                f"{self.departure_body.name} - Target Body for year {self.launch_span[0].datetime.year}, C3 Launch",
                fontsize=14,
                fontweight="bold")
        else:
            self.ax.set_title(
                f"{self.departure_body.name} - {self.target_body.name} for year {self.launch_span[0].datetime.year}, C3 Launch",
                fontsize=14,
                fontweight="bold")

        self.ax.set_xlabel(
            "Launch date", 
            fontsize=10, 
            fontweight="bold")
        
        self.ax.set_ylabel(
            "Arrival date", 
            fontsize=10, 
            fontweight="bold")
    