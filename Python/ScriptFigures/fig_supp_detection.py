import numpy as np
import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'FigSuppDetection'
LEFT = 0.9; TOP = 0.3
axs = {}
fig = fma.plotting.figure(n=0, size=(8.5,5.5))
results = isau.load_results(res_root / 'exampleISA.pkl')


# detection method
ax = isau.axes(fig, LEFT, TOP+0.5,7,2)
axs['detection_fr'] = ax
start = results['start_detect']; stop = results['stop_detect']; fr = results['fr_detect'];  spikes = results['spikes_detect']; aval = results['aval_detect']

threshold = np.percentile(fr[:,1], 30)
ax.axhline(threshold, ls='--', color='k', lw=0.5, label='threshold')
fma.plotting.plotIntervals(aval, ax=ax)
fma.plotting.plotXY(fr, color=isau.paperColors('nr'), lw=0.8)
ax.set(xlim=(start,stop), xticks=(), ylabel='population firing\nrate (Hz)', ylim=(0,ax.get_ylim()[1]))

left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, LEFT, top+height, width, height)
axs['detection_raster'] = ax
fma.plotting.plotRaster(spikes, color=isau.paperColors('nr'), compact=True, lw=0.3)
ax.set(ylim=(0,results['n_units_nr']-0.5))
ax.set(xlabel='time', xlim=(start,stop), xticklabels=(), ylabel='NR units',yticks=())
ax.yaxis.labelpad = 7

# timescale
left, top, width, height = isau.pos(ax, fig)
isau.line(fig, [left+width*4/5,left+width], [top+height+0.15,top+height+0.15], color='k', lw=0.5)
isau.text(fig, left+width*4.5/5, top+height+0.2, '50 s', fontsize=4, ha='center', va='top')



fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])
isau.save_axes_figs(axs, fig, output_file)