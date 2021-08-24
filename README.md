# LiPydomics (version 1.6.8)
_Dylan H. Ross and Jang Ho Cho_
  
  
## Overview
`lipydomics` is a Python package for performing standard lipidomics analysis protocols on data in an efficient and 
reproducible fashion. The primary modules in the package are `data`, `stats`, `plotting`, `identification` and 
`interactive`. The `data` module is responsible for storing lipidomics datasets, along with their associated metadata
and computed statistics. The `stats` module contains functions for computing statistics on a lipydomics `Dataset` 
instance and the `plotting` module has functions for generating various plots using computed statistics or the data
itself. The `identification` module has functions for identifying lipid features using different levels of information. 
The `interactive` module provides a more user-friendly text-based interface for performing lipidomics data analysis
for those who are not as familiar/comfortable with the more flexible Python API. 


## API
Documentation for the `lipydomics` API, organized by module, is available separately in [api.md](api.md). An example of 
a complete scripted analysis is included in [example_standard_analysis.py](example_standard_analysis.py). An additional 
example analysis (executed in a Jupyter notebook) is available [here](notebook_examples/API_Example.ipynb). 



## Interactive
Documentation for the `lipydomics.interactive` module is available separately in [interactive.md](interactive.md). A
complete example demonstrating use of the interactive module (executed in a Jupyter notebook) is also available 
[here](notebook_examples/Interactive_Example.ipynb).


## Support
If you have any questions or suggestions, or notice any bugs please feel free to contact dhross92@uw.edu or 
[file an issue](https://github.com/dylanhross/lipydomics/issues/new).


## License 
This software is made available under the [MIT license](LICENSE). We ask that if you use this package in published
work, please cite [our paper](https://pubs.acs.org/doi/10.1021/acs.analchem.0c02560). 

