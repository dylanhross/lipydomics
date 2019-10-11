#!/bin/bash

python3 -m pydoc -w lipydomics
python3 -m pydoc -w lipydomics.data
python3 -m pydoc -w lipydomics.stats
python3 -m pydoc -w lipydomics.plotting
python3 -m pydoc -w lipydomics.identification

mv -f lipydomics*.html lipydomics/doc/
