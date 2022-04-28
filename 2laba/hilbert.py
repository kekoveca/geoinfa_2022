from scipy.io.wavfile import read
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import hilbert, chirp
import numpy as np

mpl.rcParams['agg.path.chunksize'] = 10000

input = read("2laba/signal.wav")
audio = input[1]
plt.figure()
plt.plot(audio[4000:4500])
plt.ylabel("Amplitude")
plt.xlabel("Time")
plt.title("Sample Wave")
plt.savefig("2laba/wave.png")

audio = np.array(audio, dtype=float)
n_data = (audio - 128) / 128

analytic_signal = hilbert(n_data)
amplitude_envelope = np.abs(analytic_signal)

n_ndata = amplitude_envelope[0:(len(amplitude_envelope) - 193)]
new_data = np.reshape(n_ndata, (int(len(amplitude_envelope)/5512), 5512))

plt.figure()
# Последний элемент убираем, потому что он (может быть) заполнен не полностью
plt.imshow(new_data[:-1], cmap='gray')
plt.title('Image')
plt.savefig("2laba/image.png")

sinch = new_data[750][1860:1959]
plt.figure()
plt.plot(sinch)
plt.savefig("2laba/sinch.png")

zros = np.zeros(((int(len(amplitude_envelope)/5512), 0)))
new_data_buffered = np.hstack((new_data, zros))

deltas = []

for i in range(np.shape(new_data_buffered)[0]):
    deltas.clear()
    for j in range(1500, 2500):
        delta = (new_data_buffered[i][j:j+len(sinch)] - sinch)**2
        deltas.append(np.sum(delta))
        minimal = deltas.index(min(deltas))
    new_data_buffered[i] = np.roll(new_data_buffered[i], -(1500+minimal))


plt.figure()
# Последний элемент убираем, потому что он (может быть) заполнен не полностью
plt.imshow(new_data_buffered[:-1], cmap='gray')
plt.title('Image_shifted')
plt.savefig("2laba/image_shifted.png")

min = 0
max = 0

for j in range(400, 900):
      for i in range(2650, 2750):
        if new_data_buffered[j][i] > max:
            max = new_data_buffered[j][i]
        if new_data_buffered[j][i] < min:
            min = new_data_buffered[j][i]

convert = np.polyfit([min, max], [0, 1], 1)

for i in range(np.shape(new_data_buffered)[0]):
    for j in range(np.shape(new_data_buffered)[1]):
        new_data_buffered[i][j] = np.polyval(convert, new_data_buffered[i][j])

        if new_data_buffered[i][j] >= 1:
            new_data_buffered[i][j] = 1

        if new_data_buffered[i][j] <= 0:
            new_data_buffered[i][j] = 0

plt.figure()
# Последний элемент убираем, потому что он (может быть) заполнен не полностью
plt.imshow(new_data_buffered[:-1], cmap='gray')
plt.title('Image_shifted_normalized')
plt.savefig("2laba/image_shifted_normalized.png")

resistors = []
vert = 2700
hor = 550
for i in range(30):
    resistors.append(255*new_data_buffered[hor+i*8][vert])

plt.figure()
plt.plot(resistors[13:])
plt.title('Gradient')
plt.savefig("2laba/telemetry.png")

el = 13

c_arr = []

for i in range (4):
    c_arr.append(resistors[el + 9 + i])

## pixel coords
x, y = 600, 1600

c_s = new_data_buffered[600][133]*255
c_e = new_data_buffered[x][y]*255
c_bb = resistors[el+15]

prt_1 = [276.6067, 0.051111, 1.405783*pow(10, -6), 0, 0]
prt_2 = [276.6119, 0.051090, 1.496037*pow(10, -6), 0, 0]
prt_3 = [276.6311, 0.051033, 1.496990*pow(10, -6), 0, 0]
prt_4 = [276.6228, 0.051058, 1.493110*pow(10, -6), 0, 0]

prt = [prt_1, prt_2, prt_3, prt_4]
T = []

for i in range(4):
    temp = 0
    for j in range(5):
        temp = temp + prt[i][j]*pow(c_arr[i], j)
    T.append(temp)

print(f"Temperatures of resistors {T}")

T_bb = np.mean(T)

print(f"Temperature of black body {T_bb}")

A, B = 1.67396, 0.997364

T_bb_shtrih = A + B*T_bb

print(f"Effective temp of BB {T_bb_shtrih}")

u_e, u_c = 2670, 2670

c1, c2 = 1.1910427*pow(10, -5), 1.4387752

N_bb = (c1*u_e**3)/(np.exp((c2*u_c/T_bb_shtrih)) - 1)

print(f"N_bb {N_bb}")

N_e = N_bb * (c_s - c_e) / (c_s - c_bb)

print(f"N_e {N_e}")

T_e_shtrih = (c2*u_c)/np.log(1+((c1*(u_c**3))/(N_e)))

print(f"T*e {T_e_shtrih}")

T_e = (T_e_shtrih - A)/B

print(f"Temperature of pixel at [{x}][{y}] {T_e} K")