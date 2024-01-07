import matplotlib.pyplot as plt

# Define the number of rows and columns in the grid
rows = 8
cols = 40

# Create a figure and axis
fig, ax = plt.subplots(figsize=(8,2))
plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05) 

# Plot the grid of cells
for i in range(rows + 1):
    ax.plot([9,cols], [i, i], color="k")

for i in range(9,cols + 1):
    ax.plot([i,i], [0, rows], color="k")

# plot computation domain boundary
ax.plot([0,9,9,0,0], [0, 0, rows, rows,0], color="k")

# plot boundary condition
ax.plot([40,0,0,40],[0,0,8,8,], color = "r")

# plot forcing condition
for i in range(rows):
    ax.plot([9,10],[i,i+1], color='b')
    ax.plot([10,9],[i,i+1], color='b')

ax.set_aspect('equal')
plt.show()
