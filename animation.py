import numpy as np
import matplotlib
import os

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import rc
rc('animation', html='jshtml')


gridsize = 5
timesteps = 30
interval = 100


def remove_spaces(data):
    if data ==  '\n':
        return False
    return True

def initial_seed_check(gridsize):
    fig, ax = plt.subplots(figsize=(50,50))
    im = ax.imshow(np.reshape(states[0:gridsize**2], (gridsize,gridsize)), cmap=plt.cm.binary)
    return

f = open("cmake-build-debug/outputs.txt", "r")
spaceddata = np.asarray(list(f.read()))
filterdata = list(filter(remove_spaces,spaceddata))
data = list(float(x) for x in filterdata)
states = np.asarray(data)

def grid_animation(states, gridsize, timesteps, interval):

    fig = plt.figure(figsize = (6,6))
    ax = plt.axes([0, 0, 1, 1], xticks=[], yticks=[], frameon=False)
    im = plt.imshow(np.reshape(states[0:(gridsize**2)], (gridsize,gridsize)), cmap=plt.cm.binary, interpolation = 'nearest')
    plt.close()
    def init():
        im.set_data(np.zeros((gridsize,gridsize)))
        return [im]

    def update_animation(i, gridsize):
        state = states[(gridsize**2)*i:((gridsize**2)*(i+1))]
        state = np.reshape(state,(gridsize,gridsize))
        im.set_array(state)

        return [im]



    anim = animation.FuncAnimation(fig, update_animation,frames = timesteps, fargs = [gridsize],interval=interval)
    anim.save("animation.gif")
    plt.show()
    return anim

grid_animation(states,gridsize,timesteps,interval)