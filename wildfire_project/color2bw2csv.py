import cv2
import matplotlib.pyplot as plt

import pandas as pd

im = cv2.imread("forest_map.png")
imbw = im[:, :, 0]*0.0 + im[:, :, 1]*0.9 + im[:, :, 2]*0.0
imbw = imbw.astype(int)
imbw = imbw.reshape(imbw.shape[0], imbw.shape[1])
# print(type(im))

imbw_normalized = imbw / 255
imbw = imbw_normalized

print(imbw.shape)
print(imbw)

plt.imsave('forest_map_gray.png', imbw, cmap='gray')

df = pd.DataFrame(imbw)

df.to_csv('image.csv', index=False, header=False)
