import numpy as np
from scipy.fft import fft, ifft2


import os



def slow_conv(a1, a2):
    conv = []
    out_n = len(a1) + len(a2) - 1
    start_1 = 0  # начальные позиции для обеих последовательностей
    start_2 = 0
    for i in range(0, out_n):
        conv_elem = 0
        for j in range(0, min(len(a1) - start_1, start_2 + 1)):  # Чтобы складывать, нужно знать номера с каких начинать
            conv_elem += a1[start_1 + j] * a2[start_2 - j]
        conv.append(conv_elem)
        if start_2 != len(a2) - 1:  # Когда последовательность 2 стала "ведущей" ее старт не меняется, так как он
            # максимален
            start_2 += 1
        else:
            start_1 += 1
    return conv


def fft_conv(a1, a2):  # нужен NumPy
    A1 = np.append(np.array(a1), np.zeros(len(a2) - 1))
    A2 = np.append(np.array(a2), np.zeros(len(a1) - 1))
    fft_A1 = np.fft.fft(A1)
    fft_A2 = np.fft.fft(A2)
    fft_c = fft_A1 * fft_A2
    conv_arr = np.fft.ifft(fft_c).real
    conv_arr = np.around(conv_arr)  # можно закоментить если что
    conv = list(conv_arr)
    return conv


def fft_tw_eq_circle_conv(a):  # нужен NumPy
    A = np.append(np.array(a), np.zeros(len(a) - 1))
    # A = np.array(a)
    fft_A = np.fft.fft(A)
    fft_c = fft_A * fft_A.conj()
    conv_arr = np.fft.ifft(fft_c).real
    conv_arr = np.around(conv_arr)  # можно закоментить если что
    conv = list(conv_arr)
    return conv


def fft_th_eq_conv(a):
    A = np.append(np.array(a), np.zeros(2 * len(a) - 2))
    fft_A = np.fft.fft(A)
    fft_c = fft_A * fft_A.conj() * fft_A
    conv_arr = np.fft.ifft(fft_c).real
    conv_arr = np.around(conv_arr)  # можно закоментить если что
    conv = list(conv_arr)
    return conv


def fft_th_eq_conv_test(a, j, k):
    n = 10
    A = np.append(np.array(a), np.zeros(n))
    fft_A = np.fft.fft(A)
    fft_c = np.zeros(n, dtype=complex)
    for i in range(n):
        fft_c[i] = fft_A[j] * fft_A.conj()[k] * fft_A[k - j]
    conv_arr = np.fft.ifft(fft_c).real
    conv_arr = np.around(conv_arr)  # можно закоментить если что
    conv = list(conv_arr)
    return conv


def fft_beauty_conv(x):
    p = len(x)
    b = np.append(np.array(x), np.zeros(p - 1))
    B = fft(b)
    fft_Beauty = np.zeros((2 * p - 1, 2 * p - 1), dtype=complex)
    for j in range(2*p - 1):
        for k in range(2*p - 1):
            fft_Beauty[j][k] = B[k] * B.conj()[j] * B[j - k]  # при замене k на j - будет (ds, d)
    Beauty = ifft2(fft_Beauty)
    # Beauty = np.around(Beauty)  # можно закоментить если что
    return Beauty.real

def fft_beauty_conv_up(x):
    x = np.array(x)
    p = x.shape[1]
    fft_Beauty = np.zeros((2 * p - 1, 2 * p - 1), dtype=complex)
    print('FFT: |', end='')
    for i in range(len(x)):



        b = np.append(np.array(x[i]), np.zeros(p - 1))
        B = fft(b)
        for j in range(2*p - 1):
            for k in range(2*p - 1):
                fft_Beauty[j][k] += B[k] * B.conj()[j] * B[j - k]  # при замене k на j - будет (ds, d)
        if i%(len(x)//50):
            print('#', end='')
    Beauty = ifft2(fft_Beauty)

    print('| 100%', end='')
    return Beauty.real
