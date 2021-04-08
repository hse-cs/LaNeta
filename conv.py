import numpy as np
from scipy.fft import fft, ifft2

def fft_beauty_conv_up(x):

    x = np.array(x)
    if len(x.shape) == 1:
        x = x.reshape(1, -1)

    p = x.shape[1]
    fft_Beauty = np.zeros((2 * p - 1, 2 * p - 1), dtype=complex)
    print('FFT: |', end='')

    # # # # # # # # # # # For Output
    step = len(x)//50
    power = 1
    odd = 0
    if step == 0:
        step = 1
        power = 50//len(x)
        odd = 50%len(x)
    # # # # # # # # # # #

    for i in range(len(x)):
        b = np.append(np.array(x[i]), np.zeros(p - 1))
        B = fft(b)
        for j in range(2*p - 1):
            for k in range(2*p - 1):
                fft_Beauty[j][k] += B[k] * B.conj()[j] * B[j - k]  # при замене k на j - будет (ds, d)
    # # # # # # # # # # # # # # # # # For Output
        if i%(step) == 0:
            print('#'*power, end='')
        print('#'*odd, end='')
    # # # # # # # # # # # # # # # # #

    Beauty = ifft2(fft_Beauty)

    print('| 100%', end='')
    return Beauty.real
