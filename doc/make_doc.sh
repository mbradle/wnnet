rm -f source/wnflows.*.rst
mkdir -p source/_static source/_templates
sphinx-apidoc -M -f -n -o source ../wnflows
make html
