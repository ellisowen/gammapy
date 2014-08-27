"""Plot background estimation image.
"""
from source_diffuse_estimation import run_estimation
from aplpy import FITSFigure

_, background = run_estimation()

fig = FITSFigure(background)
fig.show_colorscale(stretch='linear', interpolation='none')
fig.add_colorbar()