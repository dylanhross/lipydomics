## Lipydomics 
__Dylan H. Ross__
  
  
### Overview
`lipydomics` is a Python package for performing standard lipidomics analysis protocols on data in an efficient and modular fashion. 
  
  
### Usage
The library contains two main components: the `Dataset` object and functions from `lipydomics/stats.py`. The `Dataset` object encapsulates a typical lipidomics dataset, with methods for loading data, assigning metadata, 
  
  
## TODOs
`lipydomics/data.py`
- [ ] Implement some sort of serialization or other method for saving a lipidomics dataset (`lipydomics.data.Dataset`) to disk in a way that it can be loaded later

`lipydomics/stats.py`
- [ ] Add a PLS-DA function for analyzing the features that contribute to group separation
  
`lipydomics/test/tests.py`
- [ ] The tests are a mess

