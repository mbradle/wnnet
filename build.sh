# Script to automate build for PyPI.

rm -fr dist
black --line-length=79 wnnet/*.py
python -m pip install --upgrade build
python -m build
