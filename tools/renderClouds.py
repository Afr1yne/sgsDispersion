#!/usr/bin/env pvpython
"""Render all native final parcels with the same camera and screen marker size."""
from paraview.simple import *
from paraview.simple import _DisableFirstRenderCameraReset
from paraview import servermanager
from pathlib import Path
import os, json, hashlib

root = Path(os.environ['FOAM_RUN'])/'sgsDispersion'
out = root/'results/abl'
names = ['RANS_DRW', 'LES_resolved', 'LES_SGS']
labels = ['(a) RANS + random walk', '(b) LES: resolved tracking', '(c) LES + SGS dispersion']
_DisableFirstRenderCameraReset()
views, evidence = [], []

def label(view, value, position, color=(0, 0, 0), size=34):
    source = Text()
    source.Text = value
    rep = Show(source, view, 'TextSourceRepresentation')
    rep.WindowLocation = position
    rep.FontFamily = 'Times'
    rep.FontSize = size
    rep.Color = color
    rep.Bold = 0

for name, title in zip(names, labels):
    folder = root/'abl'/name/'VTK/lagrangian/kinematicCloud'
    series = json.loads((folder/'kinematicCloud.vtp.series').read_text())
    path = folder/next(x['name'] for x in series['files'] if x['time'] == 440)
    source = XMLPolyDataReader(FileName=[str(path)])
    source.UpdatePipeline()
    assert servermanager.Fetch(source).GetNumberOfPoints() == 6000
    bounds = list(servermanager.Fetch(source).GetBounds())
    vertices = ProgrammableFilter(Input=source)
    vertices.OutputDataSetType = 'vtkPolyData'
    vertices.Script = '''from vtkmodules.vtkFiltersGeneral import vtkVertexGlyphFilter
glyph = vtkVertexGlyphFilter()
glyph.SetInputData(inputs[0].VTKObject)
glyph.Update()
self.GetOutput().ShallowCopy(glyph.GetOutput())
'''
    vertices.UpdatePipeline()
    view = CreateView('RenderView')
    views.append(view)
    view.ViewSize = [1600, 1100]
    view.UseColorPaletteForBackground = 0
    view.Background = [1, 1, 1]
    view.OrientationAxesVisibility = 0
    view.CenterAxesVisibility = 0
    view.CameraParallelProjection = 1
    view.CameraFocalPoint = [164, -12, 48]
    view.CameraPosition = [234, -117, 128]
    view.CameraViewUp = [0, 0, 1]
    view.CameraParallelScale = 25
    for diameter, color in [(20e-6, [0, 0, 1]), (80e-6, [1, 0, 0])]:
        selected = Threshold(Input=vertices)
        selected.Scalars = ['POINTS', 'd']
        selected.ThresholdMethod = 'Between'
        selected.LowerThreshold = diameter*.999
        selected.UpperThreshold = diameter*1.001
        selected.UpdatePipeline()
        assert servermanager.Fetch(selected).GetNumberOfPoints() == 3000
        rep = Show(selected, view, 'GeometryRepresentation')
        rep.Representation = 'Points'
        rep.ColorArrayName = ['POINTS', '']
        rep.DiffuseColor = color
        rep.AmbientColor = color
        rep.PointSize = 3.5
        rep.DisableLighting = 1
        rep.Opacity = 1
    label(view, title, 'Upper Left Corner', size=44)
    label(view, 't = 440 s', 'Upper Right Corner')
    label(view, '20 um', 'Lower Left Corner', (0, 0, 1))
    label(view, '80 um', 'Lower Right Corner', (1, 0, 0))
    Render(view)
    SaveScreenshot(str(out/f'cloud_{name}.png'), view,
                   ImageResolution=[1600, 1100], TransparentBackground=0)
    evidence.append({'case': name, 'input': str(path),
                     'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                     'native_points': 6000, 'subsampling': False,
                     'interpolation': False, 'bounds_m': bounds})

layout = CreateLayout()
layout.SplitVertical(0, 1/3)
layout.SplitVertical(2, .5)
for view, hint in zip(views, [1, 5, 6]):
    AssignViewToLayout(view, layout, hint=hint)
layout.SetSize(1600, 3300)
RenderAllViews()
SaveScreenshot(str(out/'particle_clouds.png'), layout,
               ImageResolution=[1600, 3300], TransparentBackground=0)
SaveState(str(out/'particle_clouds.pvsm'))
(out/'figure_provenance.json').write_text(json.dumps({
    'renderer': 'ParaView 6.1.0 EGL', 'time_s': 440,
    'camera': {'position': [234, -117, 128], 'focal_point': [164, -12, 48],
               'up': [0, 0, 1], 'parallel_scale_m': 25},
    'colors': {'20_um': '#0000FF', '80_um': '#FF0000'},
    'marker_size_px': 3.5, 'marker_size_represents_diameter': False,
    'sources': evidence}, indent=2), encoding='utf-8')
print('PASS native points: 6000 per case; fixed camera; no subsampling or interpolation')
