import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a new figure for 3D plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plotting the first segment
x1 = [0, 1]
y1 = [0, 1]
z1 = [0, 0]
ax.plot(x1, y1, z1, label=f'Segment 1 ({x1[0],y1[0],z1[0]} to {x1[1],y1[1],z1[1]})', color='b')

# Plotting the second segment
x2 = [-5, -0.5]
y2 = [-5, -0.5]
z2 = [ 1, 1]
ax.plot(x2, y2, z2, label=f'Segment 2 ({x2[0],y2[0],z2[0]} to {x2[1],y2[1],z2[1]})', color='r')

# Adding labels
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.legend()

# Show the plot
plt.show()