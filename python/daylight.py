# Plot contour plot of daylight by lattitude and day of year
import numpy as np
import matplotlib.pyplot as plt

import donnager as dgr

if __name__ == '__main__':
    days = range(365)
    lats = np.linspace(-90, 90, 300)

    mesh_days, mesh_lats = np.meshgrid(days, lats, indexing='xy')
    longs = np.zeros_like(lats)

    hours = dgr.calc_earth_day_length(
        lattitude_deg=mesh_lats,
        longitude_deg=longs,
        days=mesh_days,
    )

    plt.imshow(hours)
    plt.xlabel('Day of year')
    plt.ylabel('Lattitude (deg)')

    plt.colorbar()
    plt.title('Length of day vs lattitude, Earth')
