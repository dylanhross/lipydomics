#!/bin/bash

python3 -m pydoc -w lipydomics
python3 -m pydoc -w lipydomics.data
python3 -m pydoc -w lipydomics.stats
python3 -m pydoc -w lipydomics.plotting
python3 -m pydoc -w lipydomics.identification
python3 -m pydoc -w lipydomics.identification.rt_calibration
python3 -m pydoc -w lipydomics.interactive
python3 -m pydoc -w lipydomics.util

mv -f lipydomics*.html lipydomics/doc/
