from cmath import pi
from datetime import datetime, date, timedelta
from skyfield.api import load, wgs84
from skyfield.api import EarthSatellite
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

start, end = datetime.now(), datetime.now() + timedelta(days=1)
print("Интересуемый диапазон: c {} до {}".format(start.strftime("%d.%m.%Y %H:%M:%S"), end.strftime("%d.%m.%Y %H:%M:%S")))


tle_1 = '1 33591U 09005A   22069.18326583  .00000083  00000-0  69793-4 0  9999'
tle_2 = '2 33591  99.1613 101.1760 0014836 106.3564 253.9239 14.12536020674436'


satellite = EarthSatellite(tle_1, tle_2, "NOAA 19")
print(satellite)

moscow = wgs84.latlon(55.75, 37.62)

ts = load.timescale()
t0 = ts.now() + timedelta(hours = 3)
t1 = t0 + timedelta(days = 1)

print(t0.utc_strftime('%Y %b %d %H:%M:%S'), t1.utc_strftime('%Y %b %d %H:%M:%S'))

ecis = []

for i in range(0, 1440):
    geocentric = satellite.at(t0)
    ecis.append(geocentric.position.km)
    t0 += timedelta(minutes = 1)


a, b = 6378.137, 6356.7523142

ecis = np.array(ecis)
x, y, z = ecis[:, 0], ecis[:, 1], ecis[:, 2]

fig = plt.figure(figsize=(8, 8), dpi = 400)
ax = fig.add_subplot(111, projection='3d')

ax.plot(x, y, z)

rx, ry, rz = a, a, b
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = rx * np.outer(np.cos(u), np.sin(v))
y = ry * np.outer(np.sin(u), np.sin(v))
z = rz * np.outer(np.ones_like(u), np.cos(v))

ax.plot_wireframe(x, y, z, alpha=0.1)


max_radius = max(rx, ry, rz)
for axis in 'xyz':
    getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

plt.savefig("satell.png")


diff = satellite - moscow
t0 = ts.now() + timedelta(hours = 3)
vision = [] # сюда записываются показания видно/не видно и координаты связи
altidudes = [] # сюда только координаты
azimuths = []

for i in range(1440):
    topocentric = diff.at(t0)
    alt, az, distance = topocentric.altaz()
    if alt.degrees > 0:
        print(f"NOAA is VISIBLE: alt({alt.degrees}), az({az.degrees}), the time is {t0.utc_strftime('%Y %b %d %H:%M:%S')}")
        vision.append(f"NOAA is VISIBLE:  {alt.degrees}, {az.degrees}, the time is {t0.utc_strftime('%Y %b %d %H:%M:%S')}")
        altidudes.append(alt.radians)
        azimuths.append(az.radians)
    else:
        #print("NOAA is INVISIBLE, the time is ",t0.utc_strftime('%Y %b %d %H:%M:%S'))
        vision.append(f"NOAA is INVISIBLE, the time is {t0.utc_strftime('%Y %b %d %H:%M:%S')}")
    t0 += timedelta(minutes = 1)


fig = plt.figure(figsize=(8, 6), dpi = 400)
ax1 = fig.add_subplot(111, projection='polar')
plt.plot(azimuths, altidudes, c = 'black')
plt.ylim(0.5*pi, 0)
ax1.grid(True)
plt.savefig("coords.png")