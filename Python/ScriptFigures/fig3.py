import numpy as np
import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig3'
LEFT = 0.7; TOP = 0.3
fig = fma.plotting.figure(n=0, size=(17.5,10.9))

# example
ax = isau.axes(fig, LEFT, TOP, 8, 3)

# oscillations are more likely in ISA
results = isau.load_results(res_root / 'HpcPfcEntrainment.pkl')
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+0.4, top, 2.6, height)

# ISA entrains



fig.show()
#fma.plotting.saveFigure(fig, output_file, ['png','svg'])