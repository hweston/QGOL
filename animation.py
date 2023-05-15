import numpy as np
import matplotlib 
import os
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import rc

gridsize = 60
timesteps = 200
interval = 200
entanglement = 0

def initial_seed_check(gridsize):
  fig, ax = plt.subplots(figsize=(10,10))
  im = ax.imshow(np.reshape(states[0:gridsize,:], (gridsize,gridsize)), cmap= 'RdBu', vmax = 1, vmin = -1)
  return

def read_grid(entanglement):
  f = open("outputs.txt", "r")
  data = pd.read_csv("outputs.txt", header = None)
  states = np.asarray(data)
  print(states.shape)



  if entanglement != 0:
    states = states[:,(gridsize*(entanglement)):]


  states = states[:,:-1]
  return states

states = read_grid(entanglement)
initial_seed_check(gridsize)


rc('animation', html='jshtml')
def grid_animation(states, gridsize, timesteps, interval):

  fig = plt.figure(figsize = (6,6))
  ax = plt.axes([0, 0, 1, 1], xticks=[], yticks=[], frameon=False)
  im = plt.imshow(np.reshape(states[0:gridsize,:], (gridsize,gridsize)), cmap='RdBu',vmin = -1, vmax = 1, interpolation = 'nearest')

  plt.axhline(-0.5,color = "black", linewidth = 10)
  plt.axvline(-0.5,color = "black",  linewidth = 10)

  plt.axhline(gridsize - 0.5,color = "black", linewidth = 10)
  plt.axvline(gridsize - 0.5,color = "black", linewidth = 10)

  
  for i in range(0,gridsize):
    plt.axhline(i+0.5,color = "grey", linewidth = 0.2)
    plt.axvline(i+0.5, color = "grey", linewidth = 0.2)
  plt.close()
  s = "0"
  txt = ax.text(-.4,0,s, fontsize =15)
  def init():
      im.set_data(np.zeros((gridsize,gridsize)))
      return [im]

  def update_animation(i, gridsize):
      state = states[gridsize*i:gridsize*(i+1),:]
      state = np.reshape(state,(gridsize,gridsize))   
      im.set_array(state)
      s = "Time = {}".format(i)
      txt.set_text(s)
      return [im]


  anim = animation.FuncAnimation(fig, update_animation,frames = timesteps, fargs = [gridsize],interval=interval)
  plt.show()
  writervideo = animation.FFMpegWriter(fps=5) 
  anim.save("example_grid.mp4", writer=writervideo)  
  return anim

grid_animation(states,gridsize,timesteps,interval)

volume = [2,3,4]
entropy = [0.950, 1.33218,1.60944]
x = [2,4]
y =[0.95,1.9]

fig = plt.figure(figsize = (9,9))
plt.plot(volume,entropy, label = "Observed increase in SQGoL", color = 'blue')
plt.plot(x,y, label = "Volume Law Maximum Increase", linestyle = "--", color = "red")
plt.tick_params(axis='both', which='major', labelsize=15)
plt.title('Comparison of Entropy Increase with Volume Law', fontsize = 20)
plt.xlabel('Subsystem Volume', fontsize = 18)
plt.ylabel('Entropy', fontsize = 18)
plt.xlim([2,4])
plt.ylim([0.9,1.7])
plt.grid()
plt.legend(fontsize = 15)
plt.show()

boundary =[6,8]
entropy = [0.950, 1.60944]
x = [6,8]
y =[0.95,(0.95*4/3)]

fig = plt.figure(figsize = (9,9))
plt.plot(boundary,entropy, label = "Observed increase in SQGoL", color = 'blue')
plt.plot(x,y, label = "Area Law Maximum Increase", linestyle = "--", color = "red")
plt.tick_params(axis='both', which='major', labelsize=15)
plt.xlim([6,8])
plt.ylim([0.9,1.7])
plt.title('Comparison of Entropy Increase with Area Law', fontsize = 20)
plt.xlabel('Boundary Size between Subsystems', fontsize = 18)
plt.ylabel('Entropy', fontsize = 18)
plt.grid()
plt.legend(fontsize = 15)
plt.show()
