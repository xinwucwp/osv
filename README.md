## Automatic fault interpretation with optimal surface voting

This repository contains computer programs written and used by 
[Xinming Wu](http://www.jsg.utexas.edu/wu/) 
for 2D and 3D horizon extraction that is discussed in our Geophysics paper 
[Automatic fault interpretation with optimal surface voting](http://www.jsg.utexas.edu/wu/files/wu2018automaticFaultInterpretationWithOptimalSurfaceVotingLow.pdf).

If you find this work helpful in your research, please cite:

    @article{wu2018least,
        author = {Xinming Wu and Sergey Fomel},
        title = {Automatic fault interpretation with optimal surface voting},
        journal = {GEOPHYSICS},
        volume = {83},
        number = {5},
        pages = {O67-O82},
        year = {2018},
        doi = {10.1190/GEO2018-0115.1},
        URL = {https://library.seg.org/doi/abs/10.1190/geo2018-0115.1},
    }

### The source codes will be coming soon. As people are asking me for the test datasets in this paper, we decide to first share the datasets here before the source codes.

---
## Examples

2D and 3D examples published in the [paper](http://www.jsg.utexas.edu/wu/files/wu2018LeastSquaresHorizons.pdf).
### 2D examples

#### 1) Campos in Figure 9 (seismic data was provided by Dr. Michael Hudec)
Dimensions: n1=300, n2=550

Data type:  binary with BIG_ENDIAN

Seismic: ./data/2d/campos/gx373.dat

OSV fault: ./data/2d/campos/fv.dat

Thinned OSV fault: ./data/2d/campos/fvt.dat

<p align="left">
  <img src="png/2d/campos/gx.png" width="445px" height="300px"/>
  <img src="png/2d/campos/el.png" width="445px" height="300px"/>
</p>
<p align="left">
  <img src="png/2d/campos/fl.png" width="445px" height="300px"/>
  <img src="png/2d/campos/fv.png" width="445px" height="300px"/>
</p>

#### 2) Netherlands off-shore F3 (provided by the Dutch Government through TNO and dGB Earth Sciences)
Left: predictive horizons with only local slopes

Center: least-squares horizons with only local slopes

Right: least-squares horizons with both local slopes and multi-grid correlations (proposed)

<p align="left">
  <img src="png/2d/f3d/f3dp.png" width="295px" height="200px"/>
  <img src="png/2d/f3d/f3ds.png" width="295px" height="200px"/>
  <img src="png/2d/f3d/f3dm.png" width="295px" height="200px"/>
</p>

#### 3) Curt (provided by Australian government)
Top row: predictive horizons with only local slopes

Middle row: least-squares horizons with only local slopes

Bottom row: least-squares horizons with both local slopes and multi-grid correlations (proposed)

<p align="left">
  <img src="png/2d/curt/curtp.png" width="885px" height="350px"/>
</p>
<p align="left">
  <img src="png/2d/curt/curts.png" width="885px" height="350px"/>
</p>
<p align="left">
  <img src="png/2d/curt/curtm.png" width="885px" height="350px"/>
</p>

### 3D examples

#### 1) Netherlands off-shore F3 (provided by the Dutch Government through TNO and dGB Earth Sciences)
Top row: least-squares horizons with only local slopes

Bottom row: least-squares horizons with both local slopes and multi-grid correlations (proposed)
<p align="left">
  <img src="png/3d/f3d/surfs1.png" width="445px" height="350px"/>
  <img src="png/3d/f3d/surfs2.png" width="445px" height="350px"/>
</p>
<p align="left">
  <img src="png/3d/f3d/surfm1.png" width="445px" height="350px"/>
  <img src="png/3d/f3d/surfm2.png" width="445px" height="350px"/>
</p>

#### 2) provided by RCRL at BEG (purchased from Australian Government-Geoscience Australia)
A horizon surface extracted using the proposed method with one control point (green point in (b))
<p align="left">
  <img src="png/3d/aust3d/aust.png"/>
</p>

---
Copyright (c) 2018, Xinming Wu. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
