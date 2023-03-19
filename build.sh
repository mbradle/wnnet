# Script to automate build for PyPI.

rm -fr dist
python -m pip install --upgrade build
python -m build
