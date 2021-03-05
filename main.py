import conv
import numpy as np

a1 = [1, 2, 3]
a2 = [4, 5, 6]
a3 = [7, 8, 9]

A1 = [a1, a2, a3]

A = np.zeros((5, 5))
A += conv.fft_beauty_conv(a1)
A += conv.fft_beauty_conv(a2)
A += conv.fft_beauty_conv(a3)

B = conv.fft_beauty_conv_up(A1)

print(A)
print(B)
