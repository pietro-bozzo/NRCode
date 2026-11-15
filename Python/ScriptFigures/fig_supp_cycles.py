import numpy as np
import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'FigSuppCycles'
LEFT = 1; TOP = 0.5; HSPACE = 2; VSPACE = 1
axs = {}
fig = fma.plotting.figure(n=0, size=(10,4))

# ON
results = isau.load_results(res_root / 'characterizeISA.pkl')
on = results['on']; off = results['off']
width = 3
ax = isau.axes(fig, LEFT, TOP, width, 2.5)
axs['on_pdf'] = ax
fma.plotting.semPlot(on['nr'].t, on['nr'], color=isau.paperColors('nr'))
ax.set(xlabel='ISA ON duration (s)', ylabel='distribution', ylim=(None,0.2))

ax = isau.axes(fig, width*0.65, 0, 0.55, 0.65, ax, 'wh')
axs['on_mode'] = ax
peak_t_f = lambda a, axis: np.nanmean(on['nr'].t[np.argmax(a, axis=1-axis)])
b, ci = fma.statistics.hierarchicalBootstrap(on['nr'], on['nr'].rat, stat=peak_t_f, return_ci=True)
fma.plotting.boxPlot(b[0], mode='violin', color='gray', facecolor='none')
print(f'ON mode {b[0].mean()}, 95% CI: {ci}')
[ax.axhline(c, color=isau.paperColors('nr'), ls='--', lw=0.6) for c in ci[0]]
ax.set(xticks=(), ylabel='ON mode (s)', yticks=(2,3,4))

# OFF
ax = isau.axes(fig, HSPACE, 0, 1, 1, axs['on_pdf'], 'lwh')
axs['off_pdf'] = ax
fma.plotting.semPlot(off['nr'].t, off['nr'], color=isau.paperColors('nr'))
ax.set(xlabel='ISA OFF duration (s)', ylabel='distribution', ylim=(None,0.2))

ax = isau.axes(fig, width*0.65, 0, 0.55, 0.65, ax, 'wh')
axs['off_mode'] = ax
peak_t_f = lambda a, axis: np.nanmean(off['nr'].t[np.argmax(a, axis=1-axis)])
b, ci = fma.statistics.hierarchicalBootstrap(off['nr'], off['nr'].rat, stat=peak_t_f, return_ci=True)
fma.plotting.boxPlot(b[0], mode='violin', color='gray', facecolor='none')
print(f'OFF mode {b[0].mean()}, 95% CI: {ci}')
[ax.axhline(c, color=isau.paperColors('nr'), ls='--', lw=0.6) for c in ci[0]]
ax.set(xticks=(), ylabel='OFF mode (s)', yticks=(2,3,4))



fig.show()
fma.plotting.saveFigure(fig, output_file, ('png','svg','pdf'))
isau.save_axes_figs(axs, fig, output_file)