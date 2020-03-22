# API 
The following is a brief overview of the package API with usage examples, organized by module.

**Complete documentation for all modules is available in HTML format under `lipydomics/doc/lipydomics.html`.**

## Data
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

### Batch Selecting Features
Feature data can be selected and exported to a .csv file using the `Dataset.select_feature_data` method. This method
takes as input the name of a .csv file containing the features to select, the name of the output .csv file to generate,
and optional search tolerances for m/z, retention time, and CCS. All matching features and their corresponding 
intensities (raw and normalized, if available) are written to the output file. This method is useful for selecting out
data for features that you know about ahead of time. If no matching features are found, then the output file is not 
written.

**Example:**
```python
dset.select_feature_data('input.csv', 'output.csv', tolerance=(0.025, 0.25, 2.5))
```

*Structure of `input.csv`*:

| m/z | rt | ccs |
|:---:|:---:|:---:|
| 234.5678 | 3.45 | 123.4 |
| ... | ... | ... |

*Structure of `output.csv`*:

| m/z | retention time | CCS | raw intensities | | ... |  |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 234.5678 | 3.45 | 123.4 | 45678 | 56789 | ... | 23456 |
| ... | ... | ... | ... | ... | ... | ... |


## Stats
The `lipydomics.stats` module contains several functions for applying statistical analyses to a `Dataset` instance. All
stats functions have a similar interface, taking a reference to the `Dataset` instance, a list of group names to use for
the analysis, and an optional boolean indicating whether to use raw or normalized intensities as input. Depending on the
statistic being applied, one or more entries may be added to the `Dataset.stats` instance variable, which is a 
dictionary that maps statistic labels (`str`) to their corresponding data (typically `numpy.ndarray`).


### ANOVA p-value
Adds a single column containing ANOVA p-values calculated using the specified groups.

```python
from lipydomics.stats import add_anova_p

# compute ANOVA between groups A, B, C, and D with normalized data
add_anova_p(dset, ['A', 'B', 'C', 'D'], normed=True)
```
The above example would add an entry to `dset.stats` with the label `ANOVA_A-B-C-D_normed`.


### Pearson Correlation
Computes per-feature Pearson correlation coefficients between two specified groups.

```python
from lipydomics.stats import add_2group_corr

# compute correlation between 'wt' and 'ko' using raw data
add_2group_corr(dset, ['wt', 'ko'])
```
The above example would add an entry to `dset.stats` with the label `2-group-corr_wt-ko_raw`.


### Principle Components Analysis
Computes a 3-component PCA using the specified groups, storing the PCA projections for each sample and the component
loadings for each feature separately.

```python
from lipydomics.stats import add_pca3

# compute the PCA for 4 groups using normalized data
add_pca3(dset, ['sweet', 'sour', 'salty', 'bitter'], normed=True)
```
The above example would add two entries to `dset.stats` with labels `PCA3_sweet-sour-salty-bitter_loadings_normed` and 
`PCA3_sweet-sour-salty-bitter_projections_normed`.


### Partial Least-Squares Discriminant Analysis
Computes a PLS-DA between two specified groups, storing the projections for each sample and the component
loadings for each feature separately.

```python
from lipydomics.stats import add_plsda

# compute PLS-DA between 'wt' and 'ko' using raw data
add_plsda(dset, ['wt', 'ko'])
```

### Partial Least-Squares Regression Analysis
Computes a PLS-RA using specified groups and an external continuous variable, storing the projections for each sample
and the component loadings for each feature separately. 

```python
from lipydomics.stats import add_plsra

# load the external continuous variable
with open('external.txt', 'r') as f:
    external = [float(_.strip()) for _ in f.readlines()]

# compute PLS-RA with multiple groups, using normalized data
add_plsra(dset, ['elf', 'halforc', 'drow', 'aasimar'], external, normed=True)
```

### Two Group Log2(fold-change)
Computes the Log2(fold-change) between the mean intensities of two specified groups for all lipid species. If lipid
identifications have also been made, this data can be used to produce heatmaps of these fold-changes organized by lipid
class.

```python
from lipydomics.stats import add_log2fc

# compute Log2(fold-change) with raw data
add_log2fc(dset, ['Grog', 'Scanlan'])
```

## Plotting
The `lipydomics.plotting` module contains several functions for generating plots from the data stored in a `Dataset` 
instance. All plotting  functions have a similar interface, taking a reference to the `Dataset` instance, a list of 
group names, a directory name to save the image under, and an optional boolean indicating whether to use raw or 
normalized intensities. The list of group names and raw/normalized data indicator are used to reference `Dataset.stats`
and retrieve the relevant data for plotting.


### Barplot Feature by Group
The `barplot_feature_bygroup` function generates barplots of average intensities (and standard deviations) for features 
from a defined set of groups. Feature identifiers (m/z, retention time, and CCS) are provided as input, along with 
search tolerances for each, and plots are generated for all matching features in the `Dataset`. This function returns 
a boolean indicating whether at least one matching feature was found (and plotted). If no features are found, no plots
are generated. 

```python
from lipydomics.plotting import barplot_feature_bygroup

# m/z = 234.65678, retention time = 2.34, CCS = 123.4
feature = (234.5678, 2.34, 123.4)
# tight search tolerance
tol = (0.01, 0.1, 1.0)
found_feat = barplot_feature_bygroup(dset, ['A', 'B', 'C'], 'analysis/features/', feature, tolerance=tol)
if found_feat:
    print('found matching feature(s), generated plot(s)')
else:
    print('failed to find matching features, no plots generated')
```


### Batch Barplot Features by Group
The `batch_barplot_feature_bygroup` function can generate bar plots for multiple features at once using the same groups
and search tolerance. The call signature is the same as the normal `barplot_feature_bygroup` function, except that
instead of taking a tuple defining an individual feature to plot, it takes the path to a .csv file defining all of the 
features to look for and try to plot. 

```python
from lipydomics.plotting import batch_barplot_feature_bygroup

# tight search tolerance
tol = (0.01, 0.1, 1.0)
batch_barplot_feature_bygroup(dset, ['A', 'B', 'C'], 'analysis/features/', 'plot_these_features.csv', tolerance=tol)
```

*Where `plot_these_features.csv` has the following structure:*

| m/z | rt | ccs |
|:---:|:---:|:---:|
| 234.5678 | 3.45 | 123.4 |
| ... | ... | ... |


### Scatter PCA Projections by Group
Generate a scatter plot of PCA projections for a specified set of groups. Requires the corresponding statistic to
already be present in `Dataset.stats`.

```python
from lipydomics.plotting import scatter_pca3_projections_bygroup

# groups A, B, C, and D, normalized data
scatter_pca3_projections_bygroup(dset, ['A', 'B', 'C', 'D'], 'analysis/', normed=True)
```


### Scatter PLS-DA Projections by Group
Generate a scatter plot of PLS-DA projections for two specified set of groups. Requires the corresponding statistic to
already be present in `Dataset.stats`.

```python
from lipydomics.plotting import scatter_plsda_projections_bygroup

# groups wt vs. ko, normalized data
scatter_plsda_projections_bygroup(dset, ['wt', 'ko'], 'analysis/wt_ko/', normed=True)
```


### S-Plot 
Generates an S-Plot using the x-loadings and Pearson correlation coefficients for all features, computed using the same
two groups. Requires the corresponding statistics to already be present in `Dataset.stats`.

```python
from lipydomics.plotting import splot_plsda_pcorr_bygroup

# groups wt vs. ko, raw data
splot_plsda_pcorr_bygroup(dset, ['wt', 'ko'], 'analysis/wt_ko/')
```


### Scatter PLS-RA Projections by Group
Generate a scatter plot of the PLS-RA projections for a specified set of groups. Requires the corresponding statistic to
already be present in `Dataset.stats`.


### Heatmap of Log2(fold-change) by Lipid Class
Generate a heatmap of Log2(fold-change) of two groups for a specified lipid class. Requires the corresponding statistic 
(`LOG2FC`) to be present in `Dataset.stats` and lipid identifications to be made prior. This function returns a boolean
indicating whether the specified lipid class was found in the lipid identifications. No plot is generated if the lipid
class was not found. 

```python
from lipydomics.plotting import heatmap_lipid_class_log2fc

# Lipid Class: PG
# Groups: Grog, Scanlan
# Save in Dir: ./analysis/plots/
# Data: raw
found_lipids = heatmap_lipid_class_log2fc('PG', dset, ['Grog', 'Scanlan'], 'analysis/plots/')
if found_lipids:
    print('Found PGs in lipid identifications, generated a plot')
else:
    print('Did not find any PGs in lipid identifications, no plot generated')
```


## Identification
The `lipydomics.identification` module allows lipid features to be identified on the basis of their m/z, retention time,
and CCS. Lipid identifications are produced by comparison against a database of experimentally observed lipids, as well
as a database of theoretical values. Different levels of identification can be specified as follows:

| identification level | description |
| :---: | :--- |
| `theo_mz` | match only on theoretical m/z |
| `theo_mz_rt` | match on theoretical m/z and retention time |
| `theo_mz_ccs` | match on theoretical m/z and CCS |
| `theo_mz_rt_ccs` | match on theoretical m/z, retention time, and CCS |
| `meas_mz_ccs` | match on measured m/z and CCS |
| `meas_mz_rt_ccs` | match on measured m/z, retention time, and CCS |
| `any` | start at the highest level (`meas_mz_rt_ccs`) then work down until an identification can be made |

_Example:_
```python
from lipydomics.identification import add_feature_ids

tol = (0.02, 0.2, 2.0)
# identify features at the highest level possible
add_feature_ids(dset, tol, level='any')
```
The resulting lipid identifications are stored in the `Dataset.feat_ids` instance variable as lists of putative IDs for 
each feature. The `Dataset.feat_id_levels` instance variable holds the identification level for each feature, and the 
`Dataset.feat_id_scores` instance variable holds the scores for each putative lipid ID. Multiple calls to 
`add_feature_ids` with different parameters will override results from previous calls. 

**IMPORTANT:** *CCS tolerance is in percent NOT absolute units*


### Retention Time Calibration
All of the retention times (measured or theoretical) in the lipid database correspond to a reference HILIC method 
(Hines, _et al. J. Lipid Res._ **58**, 2017). The `lipydomics.identification.rt_calibration` module allows comparison 
between retention times measured on other (HILIC) methods via the `add_rt_calibration` method:

```python
from lipydomics.identification.rt_calibration import add_rt_calibration

# lipid calibrants, measured and reference retention times
lipids = ['FFA(18:0)', 'PG(36:0)', 'CL(60:0)', 'LysylPG(34:0)']
meas_rt = [0.673, 1.549, 2.843, 3.996]
ref_rt = [0.673, 1.549, 4.393, 7.252]

# apply RT calibration to Dataset
add_rt_calibration(dset, lipids, meas_rt, ref_rt)
```

**Once a retention time calibration has been set up, `add_feature_ids` automatically uses calibrated retention times 
when trying to identify lipids**
