# Script to automate build for PyPI.

rm -fr dist
cd wnnet
black --line-length=79 *.py
pylint --disable=C0103,W0611 consts.py
pylint __about__.py __init__.py flows.py  graph.py graph_helper.py net.py  nuc.py reac.py zones.py

cd ../.github/workflows/
pytest
cd ../..

python -m pip install --upgrade build
python -m build
python -m pip install --upgrade twine

echo ""
echo "All version numbers must be the same:"
echo ""

grep version wnnet/__about__.py | grep -v ","
grep version CITATION.cff | grep -v "cff-version"
grep Version doc/source/changelog.rst | grep -v Versioning | head -1
grep version pyproject.toml

echo ""
echo "Check the release date:"
echo ""
grep date CITATION.cff
echo ""
