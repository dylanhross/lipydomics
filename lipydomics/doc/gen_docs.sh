#!/bin/bash

python3 -m pydoc -w lipydomics
python3 -m pydoc -w lipydomics.data
python3 -m pydoc -w lipydomics.stats
python3 -m pydoc -w lipydomics.plotting
python3 -m pydoc -w lipydomics.test
python3 -m pydoc -w lipydomics.test.__main__
python3 -m pydoc -w lipydomics.test.data
python3 -m pydoc -w lipydomics.test.stats
python3 -m pydoc -w lipydomics.test.plotting


mv -f lipydomics*.html lipydomics/doc/
