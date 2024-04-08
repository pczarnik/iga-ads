import cv2
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

im = cv2.imread("forest_map.png")
imbw = im[:, :, 0]*0.299 + im[:, :, 1]*0.587 + im[:, :, 2]*0.114
imbw = imbw.astype(int)
imbw = imbw.reshape(imbw.shape[0], imbw.shape[1])
# print(type(im))

imbw_normalized = imbw / 255
imbw = imbw_normalized

print(imbw.shape)
print(imbw)

plt.imshow(imbw, cmap='gray')
plt.show()

np.savetxt("image.txt", imbw, fmt='%.5f')

df = pd.DataFrame(imbw)

df.to_csv('image.csv', index=False, header=False)
