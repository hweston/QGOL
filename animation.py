import numpy as np
import matplotlib 
import os
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import rc

gridsize = 6
timesteps = 30
interval = 200
entanglement = 8

def initial_seed_check(gridsize):
  fig, ax = plt.subplots(figsize=(10,10))
  im = ax.imshow(np.reshape(states[0:gridsize,:], (gridsize,gridsize)), cmap= 'binary')
  return

def read_grid(entanglement):
  f = open("outputs.txt", "r")
  data = pd.read_csv("outputs.txt", header = None)
  states = np.asarray(data)

  if entanglement != 0:
    states = states[:,(gridsize*((2**entanglement))):]

  states = states[:,:-1]
  return states

states = read_grid(entanglement)
initial_seed_check(gridsize)



rc('animation', html='jshtml')
def grid_animation(states, gridsize, timesteps, interval):

  fig = plt.figure(figsize = (6,6))
  ax = plt.axes([0, 0, 1, 1], xticks=[], yticks=[], frameon=False)
  im = plt.imshow(np.reshape(states[0:gridsize,:], (gridsize,gridsize)), cmap='binary', interpolation = 'nearest')
  plt.close()
  def init():
      im.set_data(np.zeros((gridsize,gridsize)))
      return [im]

  def update_animation(i, gridsize):
      state = states[gridsize*i:gridsize*(i+1),:]
      state = np.reshape(state,(gridsize,gridsize))   
      im.set_array(state)

      return [im]



  anim = animation.FuncAnimation(fig, update_animation,frames = timesteps, fargs = [gridsize],interval=interval)
  plt.show()
  #writervideo = animation.FFMpegWriter(fps=20) 
  #anim.save('Name.mp4', writer=writervideo)  
  return anim

grid_animation(states,gridsize,timesteps,interval)
