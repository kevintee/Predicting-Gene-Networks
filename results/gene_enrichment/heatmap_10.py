# This will do heatmap stuff.
import matplotlib.pyplot as plt
import numpy as np

column_labels = [
    'cell cycle',
    'proliferation',
    'apoptosis',
    'leukocyte',
    'immune',
    'T cell',
    'immune response to tumor cell',
    'Killer cell',
    'DNA repair',
    'DNA replication',
]

row_labels = [
	'sbm',
	'merlin',
]

data_unprocessed = [[4,2], [2,5], [3,1], [4,5], [4,5], [3,4], [1, 0], [1,3], [1,3], [2,3]]

# normalize the data
n_sbm = 10
n_merlin = 20
data_processed = []
for d in data_unprocessed:
	data_processed.append([1.0*d[0]/n_sbm, 1.0*d[1]/n_merlin])
data = np.array(data_processed)
print data
rows = list(column_labels)
columns = list(row_labels)
 
#  Finishing Touches
fig,ax=plt.subplots()
# using the ax subplot object, we use the same
# syntax as above, but it allows us a little
# bit more advanced control
ax.pcolor(data,cmap=plt.cm.Greys_r,edgecolors='k')
ax.set_xticks(np.arange(0,2)+0.5)
ax.set_yticks(np.arange(0,10)+0.5)
 
# Here we put the x-axis tick labels
# on the top of the plot.  The y-axis
# command is redundant, but inocuous.
ax.xaxis.tick_top()
ax.yaxis.tick_left()
# similar syntax as previous examples
ax.set_xticklabels(columns,minor=False,fontsize=20)
ax.set_yticklabels(rows,minor=False,fontsize=20)
 
# Here we use a text command instead of the title
# to avoid collision between the x-axis tick labels
# and the normal title position
plt.text(0.5,1.08,'SBM vs Merlin pathway enrichments',
         fontsize=25,
         horizontalalignment='center',
         transform=ax.transAxes
         )
 
# standard axis elements
plt.ylabel('Y Axis Label',fontsize=20)
plt.xlabel('X Axis Label',fontsize=20)
plt.show()