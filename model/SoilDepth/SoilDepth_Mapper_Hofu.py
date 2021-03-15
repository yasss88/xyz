# -*- coding: utf-8 -*-
"""
SoilDept_Mapper_Hohu.py
by Yasushi TANAKA: 2021/3/3
"""

import numpy as np
import matplotlib.pyplot as plt
import time
start = time.time()

X,Y=1536,2358
dem = np.fromfile('katsu01m_SoilDepth.flt', np.float32).reshape([Y,X])
dem[dem==-9999]=np.nan
#plt.imshow(dem, cmap="rainbow")
plt.imshow(dem, cmap="jet", vmin=0, vmax=1)
plt.colorbar()

process_time = time.time() - start
print(process_time)
