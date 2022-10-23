#!/bin/bash
set -e
GHUSER=gboeing
PACKAGE=osmnx

# ensure script was run from the package root, so cd and rm behave predictably
if [[ ${PWD##*/} != $PACKAGE ]]; then
    echo "You must run this script from the package's root directory."
    exit 1
fi

# update necessary python packaging packages
eval "$(conda shell.bash hook)"
conda activate base
mamba update conda-smithy setuptools twine --yes --no-banner

# get the current package version number
cd ./$PACKAGE
VERSION=$(python -c "from _version import __version__ as v; print(v)")
cd ..

# build the distribution then get its SHA256
rm -rf ./build ./dist
python setup.py sdist bdist_wheel
twine check dist/*
SHA=$(openssl dgst -sha256 -r "./dist/$PACKAGE-$VERSION.tar.gz"  | awk '{print $1}')

# rerender the conda-forge feedstock
cd ../$PACKAGE-feedstock
git pull origin main
git pull upstream main
conda smithy rerender --commit auto

# wait for user to update feedstock recipe
echo ""
echo "Ready to push git tags, upload PyPI distribution, and update conda-forge feedstock."
echo "$PACKAGE version is $VERSION and sha256 is $SHA"
echo "Please update \"$PACKAGE-feedstock/recipe/meta.yaml\" version, sha256, and build number before proceeding!"
echo ""
read -r -n 1 -p "Ready to proceed? (y/n) " INPUT
echo ""
if [[ $INPUT != "Y" && $INPUT != "y" ]]; then
    exit 1
fi

# commit and push feedstock recipe changes
git commit -am "version bump to v$VERSION"
git push origin main

# tag new version in repo, push tags, and upload distribution to pypi
cd ../$PACKAGE
git tag -a "v$VERSION" -m "v$VERSION"
git push origin --tags
twine upload ./dist/*
rm -rf ./build ./dist

# all done
echo ""
echo "Done: git tags pushed, PyPI distribution uploaded, and conda-forge feedstock updated."
echo "Please open a pull request at https://github.com/$GHUSER/$PACKAGE-feedstock/pulls"
echo ""
