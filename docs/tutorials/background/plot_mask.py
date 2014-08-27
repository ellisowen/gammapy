"""Plot source mask.
"""
from source_diffuse_estimation import run_estimation
from aplpy import FITSFigure

mask, _ = run_estimation()

fig = FITSFigure(mask)
fig.show_grayscale(stretch='linear', interpolation='none')
