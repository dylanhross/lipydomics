#!/bin/bash

python3 -m pydoc -w lipydomics
python3 -m pydoc -w lipydomics.data
python3 -m pydoc -w lipydomics.stats
python3 -m pydoc -w lipydomics.plotting

mv -f lipydomics*.html lipydomics/doc/
