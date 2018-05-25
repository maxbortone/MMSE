# Multi Mode Source Extraction
_Pushing the limits of SExtractor_

This repository contains all of the research I did during my master degree semester project,
which took place at the Harvard-Smithsonian Center of Astrophysics in Cambridge, MA.
My work has been supervised by Sandro Tacchella<sup>[1](#harvard)</sup>, Daniel Eisenstein<sup>[1](#harvard)</sup>, Ben Johnson<sup>[1](#harvard)</sup> and Alexandre Refregier<sup>[2](#eth)</sup>.

## Content
- `data` folder contains:
  - HST extreme deep field fits files for various bands
  - 3D-HST catalog
  - a H band PSF
  - a pickled version of the final catalog
- `*.ipynb` files are Jupyter notebooks that document the research
- `mmse.py` contains the final `MultiModeSourceExtractor` class

## Getting started
Before running any notebooks, download the HST XDF
```
cd ./data/
wget $(cat hlsp_xdf_hst_download_v1.txt)
```

## Jupyter notebooks
- *Overshredding*: analysis on SExtractor's parameters and how to get it to detect sub-structures
- *Kernels*: comparison of the detections obtained using different filtering kernels
- *Asymmetries*: analysis on the difference between peak value and barycentre of detected objects in order to set priors on final positions
- *Multi Mode Source Extractor*: example on how to use MMSE using H and I bands of HST XDF images
- *Stacking*: example on how to increase detection efficiency by stacking multiple H bands together

To see the notebooks, run `jupyter notebook` from the root directory of the project.

## Slides
The *Slides* notebook contains a short presentation about my work. To see the presentation run
```
jupyter nbconvert Slides.ipynb --to slides --post serve
```
from the root directory of the project or run it as a normal notebook and start it with the [RISE](https://github.com/damianavila/RISE) extension.

----
<ol>
  <li><a name="harvard">Harvard-Smithsonian Center of Astrophysics, Cambridge MA</a></li>
  <li><a name="eth">Institute for Particle Physics and Astrophysics, ETH Zürich</a></li>
</ul>
