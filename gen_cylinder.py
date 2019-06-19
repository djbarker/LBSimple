"""
Generate a simple domain with a cylinder and two walls. 
"""

#%%
import numpy as np

Nx = 800
Ny = 400

R = 0.05 * Ny
Cx = 0.1 * Nx 

dom = np.ones((Nx, Ny), dtype=np.uint8) * 255

# velocity boundary
dom[0, :] = 3
dom[-1, :] = 3

# wall boundaries
dom[:, 0] = 0
dom[:, -1] = 0

for i in range(Nx):
	for j in range(Ny):
		if ((j-Ny//2)**2 + (i-Cx)**2) <= R**2:
			dom[i, j] = 0

dom = dom.T

with open("cylinder.data", "wb") as fout:
	dom.tofile(fout)

#%%

import matplotlib.pyplot as plt

plt.figure(figsize=(12, 12))
plt.imshow(dom)
plt.show()
