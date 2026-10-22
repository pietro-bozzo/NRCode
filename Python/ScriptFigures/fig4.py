import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig4'
LEFT = 0.8; TOP = 0.3; VSPACE = 0.8
fig = fma.plotting.figure(n=0, size=(17.5,10.9))

# corr
results = isau.load_results(res_root / 'unitLagCorr.pkl')
ax = isau.axes(fig, LEFT, TOP, 3, 3)





fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])