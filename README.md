# Lipydomics 
_Dylan H. Ross and Jang Ho Cho_
  
  
# Overview
`lipydomics` is a Python package for performing standard lipidomics analysis protocols on data in an efficient and 
reproducible fashion. The primary modules in the package are `data`, `stats`, `plotting`, `identification` and 
`interactive`. The `data` module is responsible for storing lipidomics datasets, along with their associated metadata
and computed statistics. The `stats` module contains functions for computing statistics on a lipydomics `Dataset` 
instance and the `plotting` module has functions for generating various plots using computed statistics or the data
itself. The `identification` module has functions for identifying lipid features using different levels of information. 
The `interactive` module provides a more user-friendly text-based interface for performing lipidomics data analysis
for those who are not as familiar/comfortable with the more flexible Python API. Documentation for all of the above
modules is available in HTML format under `lipydomics/doc/lipydomics.html`.
  
  
# API 

## Dataset
The `lipydomics.data` module contains the `Dataset` object, which organizes all of the relevant data from a lipidomics 
analysis together. This includes raw data, normalized data, group assignments, experimental metadata, computed 
statistics, retention time calibration, and compound identifications. 

### Initialization
A `Dataset` instance may either be initialized using a `.csv` file containing feature labels 
(m/z, retention time, CCS) and intensities, or by loading a previously saved `Dataset` instance from serialized binary 
format.

**Loading from `.csv`:**
```python
from lipydomics.data import Dataset

# this data was acquired in positive ESI mode
dset = Dataset('example.csv', esi_mode='pos')
```

Where `example.csv` has the following layout:

| m/z | rt | ccs | sample1 | sample2 | ... | sampleN |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 234.5678 | 3.45 | 123.4 | 45678 | 56789 | ... | 23456 |
| ... | ... | ... | ... | ... | ... | ... |

The first three columns contain the feature identifiers, while the remaining columns contain the feature intensities for
each sample.

**Loading from serialized binary:**
```python
from lipydomics.data import Dataset

# old_dset.pickle is a Dataset that has been saved in a serialized binary format
dset = Dataset.load_bin('old_dset.pickle')
```

### Assigning Groups
Groups can be assigned either explicitly by column index or by providing a list of group names and number of replicates
in each group. Explicit group assignment is more flexible, as it allows you to assign groups multiple ways and not
necessarily in order, but assignment by replicates is easier when the data is organized for it.

**Explicit group assignment:**
```python
# dictionary mapping group identifiers (str) to sample indices (int)
group_idx = {
    'control': [0, 1, 2],
    'knockout': [3, 4, 5],
    'bugbear': [6, 7, 8]
}
dset.assign_groups(group_idx)

# method can be called again to edit existing or add more group assignments
dset.assing_groups({'kobold': [9, 10, 11]})
```

**Group assignment by replicates:**
```python
# has the same effect as the above example
dset.assign_groups_with_replicates(['control', 'knockout', 'bugbear' ,'kobold'], 3)
```

### Data Normalization
When a `Dataset` is first initialized, the intensities are stored in the `Dataset.intensities` instance variable which 
is a `numpy.ndarray` of shape: _(n_features, n_samples)_. A user-defined normalization may be applied to these 
intensities using a vector of weights (one per sample), which is multiplied down the intensities matrix and the result
stored in the`Dataset.normed_intensities` instance variable.

**Applying a normalization:**
```python
from numpy import array

weights = array([0.95, 0.90, 1.10, 1.05, 1.00, 0.90, 1.05, 1.10, 0.95])
dset.normalize(weights)
```
The `weights` vector would be determined prior using _e.g._ some external factor like tissue weight or the signal 
intensities of a known feature (internal standard). Multiple normalizations can be performed by using a weights vector
that is a product of the weights vectors from the individual normalizations.

**Applying multiple normalizations:**
```python
from numpy import array

# based on tissue weight
weights_tissue = array([0.95, 0.90, 1.10, 1.05, 1.00, 0.90, 1.05, 1.10, 0.95])
# based on internal standard
weights_intstd = array([0.88, 0.99, 1.11, 1.00, 1.22, 0.77, 1.33, 1.11, 0.99])

# apply both normalizations
weights = weights_tissue * weights_intstd
dset.normalize(weights)
```

### Saving to File
A `Dataset` instance can be saved in a serialized binary format, retaining group assignments, normalization, computed
statistics, retention time calibration, and lipid identifications. Internally, this is done using Python's `pickle` 
library. 

```python
dset.save_bin('saved_dataset.pickle')
```

## Stats
The `lipydomics.stats` module contains several functions for applying statistical analyses to a `Dataset` instance.

