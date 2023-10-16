# Notes
mac Si: Make sure you're in x86_64/intel mode
huge issues if you try to rename the package NAME in setup.py


# How To prep with pup
``` bash
conda env create
conda activate harp
pip install pup
pip install ./
pup package ./ --icon-path ./harp/gui/HARP_logo.png --nice-name HARP --license-path ./LICENSE
```
