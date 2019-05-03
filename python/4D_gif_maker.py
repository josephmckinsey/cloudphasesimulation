
from __future__ import division
 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def import_matrix(file_name):
    matrix = {}
    with open(file_name, "r") as f:
        first_line = f.readline()
        first_line = first_line.split(',')
        first_line = [int(x) for x in first_line]
        lines = f.readlines()
        for line in lines:
            values = [int(x) for x in line.split(",")[:-1]]
            values.append(float(line.split(',')[-1]))
            matrix[(values[0], values[1], values[2])] = values[3]
    return (first_line, matrix)

class AnimatedGif:
    def __init__(self, size=(640, 480)):
        self.fig = plt.figure()
        self.fig.set_size_inches(size[0]/20, size[1]/20)
        ax = self.fig.add_axes([0, 0, 1, 1], frameon=False, aspect=1)
        ax.set_xticks([])
        ax.set_yticks([])
        self.images = []
 
    def add(self, image, label=''):
        plt_im = plt.imshow(image, cmap='Blues', vmin=0, vmax=1, animated=True)
        plt_txt = plt.text(10, 310, label, color='red')
        self.images.append([plt_im, plt_txt])
 
    def save(self, filename):
        animation = anim.ArtistAnimation(self.fig, self.images)
        animation.save(filename, writer='ffmpeg', fps=10)

(dims, matrix) = import_matrix("../matlab/u_matrix.csv")
print(dims)
print(matrix[(dims[0]-1,dims[1]-1,dims[2]-1)])
cloud_gif = AnimatedGif(size=(dims[0],dims[1]))
for t in range(dims[2]-1):
    M = np.zeros((dims[0],dims[1]))
    for x in range(dims[0]-1):
        for y in range(dims[1]-1):
            if matrix[(x+1,y+1,t+1)] < 0.5:
                M[x,y] = 1
            else:
                M[x,y] = 0
    cloud_gif.add(M, label=str(t))
cloud_gif.save("cloud.mp4")

