import matplotlib.pyplot as plt
import os
import struct
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')

#网格的尺寸
r = 100
c = 100
#文件的时间步间隔
dt = 10
#dat文件目录
datroot = "./result/"
#输出图片文件目录
picroot = "./pic/"

files = os.listdir(datroot)


for i in range(0, len(files)):
    print(i)
    filename = "{}.dat".format(dt*i)
    with open(datroot + filename, 'rb') as fp:
        data = fp.read()
        if len(data) != r*c*8:
            print("Error size")

    arr = struct.unpack('{}d'.format(r*c), data)
    plt.axis('off')
    plt.imshow(np.array(arr).reshape(r, c), origin="lower")

    plt.savefig(picroot+str(i)+'.jpg')
    plt.close('all')
