# Aquaponics #

Model and control of an aquaponics system.

## Installation ##

It is *highly* recommended that all installation is performed within a python
virtual environment. The code is designed to run with Python 3 (specifically
Python 3.5.2).

Install core dependencies:

        > pip install -U numpy pandas matplotlib jupyter

Navigate to the project root and install project and any additional
dependencies:

        > cd /path/to/aquaponics
        > python setup.py develop

## Project Structure ##

The core libraries are located at `/path/to/aquaponics/aquaponics`. The sources
from which the models are defined are found at `/path/to/aquaponics/papers/`.

Experiments and results are contained within Jupyter notebooks found at
`/path/to/aquaponics/notebooks`. To access these notebooks:

        > cd /path/to/aquaponics/notebooks
        > jupyter notebook

## Project Writeup ##

Click on the following link to navigate to an editable copy of the final
report for this project:

- [Aquaponics Status Reports (sharelatex.com)](https://www.sharelatex.com/3163581738vwgxshvxxccr)
- [Final Report (sharelatex.com)](https://www.sharelatex.com/7897664638gzxkbzcjxngx)

Milestone reports will be incomplete versions of this document.
