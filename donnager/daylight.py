# Plot contour plot of daylight by lattitude and day of year
import numpy as np
import matplotlib.pyplot as plt

import donnager as dgr

if __name__ == '__main__':
    days = range(365)
    lats = np.linspace(-90, 90, 300)
    # meshgrid
    hours = dgr.calc_earth_day_length(
        days,
        lats,
        longs
    )
    plt.imshow(

    )
    plt.xlabel('Day of year')
    plt.ylabel('Lattitude (deg)')

    plt.colorbar()
    plt.title('Length of day vs lattitude, Earth')
