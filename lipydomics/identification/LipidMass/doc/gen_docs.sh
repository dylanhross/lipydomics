#!/bin/bash

# generate documentation for the module

python3 -m pydoc -w LipidMass
    python3 -m pydoc -w LipidMass.base
    python3 -m pydoc -w LipidMass.monoiso
    python3 -m pydoc -w LipidMass.lipids
        python3 -m pydoc -w LipidMass.lipids.glycerolipids
        python3 -m pydoc -w LipidMass.lipids.glycolipids
        python3 -m pydoc -w LipidMass.lipids.glycerophospholipids
        python3 -m pydoc -w LipidMass.lipids.sphingolipids
        python3 -m pydoc -w LipidMass.lipids.lysoglycerophospholipids
        python3 -m pydoc -w LipidMass.lipids.misc

# move all of the generated documentation into the doc directory
mv -f ./*.html LipidMass/doc/
