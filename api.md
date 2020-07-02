# API (version 1.6.x)
The following is a brief overview of the package API with usage examples, organized by module.

**Complete documentation for all modules is available in HTML format under `lipydomics/doc/lipydomics.html`.**

__Modules:__
* [Data](#data)
* [Stats](#stats)
* [Plotting](#plotting)
* [Identification](#identification)


## Data
The `lipydomics.data` module contains the `Dataset` object, which organizes all of the relevant data from a lipidomics 
analysis together. This includes raw data, normalized data, group assignments, experimental metadata, computed 
statistics, retention time calibration, and compound identifications. 

### The `Dataset` Object
The `Dataset` object stores all data and metadata relevant to a lipidomics data analysis. This data is accessible through
the instance variables described below.
* `Dataset.csv` (_str_) - name of the .csv file used to initialize the `Dataset`
* `Dataset.esi_mode` (_str_) - used for lipid identification, "pos" for positive mode data, "neg" for negative mode data
* `Dataset.n_features` (_int_) - number of features (_i.e._ rows)
* `Dataset.n_samples` (_int_) - number of samples (_i.e._ columns)
* `Dataset.labels` (_numpy.ndarray_) - feature labels [[mz, rt, ccs], ... ], shape = (n_features, 3)
* `Dataset.intensities` (_numpy.ndarray_) - feature intensities (raw) [[i0, i1, ... ], ... ], shape = (n_features, n_samples)
* `Dataset.normed_intensities` (_numpy.ndarray_) - feature intensities (normalized) [[i0, i1, ... ], ... ], shape = 
(n_features, n_samples), set by the [`Dataset.normalize(...)`](#data-normalization) method
* `Dataset.group_indices` (_dict[str : list(int)]_) - dictionary mapping sample group names to their corresponding indices 
in the `Dataset.intensities` matrix, set by the [`Dataset.assign_groups(...)`](#assigning-groups) or 
[`Dataset.assign_groups_with_replicates(...)`](#assigning-groups) methods
* `Dataset.stats` (_dict[str : numpy.ndarray]_) - dictionary mapping statistic labels to their corresponding data arrays,
entries are added via calls to functions defined in the [`lipydomics.stats`](#stats) module
* `Dataset.feat_ids` (_list(list(str))_) - a list of putative lipid identifications in descending order of confidence 
for each feature (or a single string indicating an unidentified feature), set using the [`lipydomics.identification`](#identification)
module
* `Dataset.feat_id_levels` (_list(str)_) - indicates the confidence level of lipid identifications made for each feature,
set using the [`lipydomics.identification`](#identification) module
* `Dataset.feat_id_scores` (list(list(float))) - confidence scores of all putative lipid identifications for each feature,
set using the [`lipydomics.identification`](#identification) module
* `Dataset.ext_var` (_numpy.ndarray_) - a feature-length array containing an external variable, set when performing 
[PLS-Regression Analysis](#partial-least-squares-regression-analysis)

The current state (group assignments, computed statistics, lipid identifications, _etc._) of a `Dataset` instance can be 
investigated with the informative `__repr__` method build into the object. An example of the output of a call to `print(dset)`:  
```
Dataset(
	csv="/Users/DylanRoss/Documents/GitHub/lipydomics/lipydomics/test/real_data_1.csv",
	esi_mode="neg",
	samples=20,
	features=313,
	identified=102,
	normalized=False,
	rt_calibrated=False,
	ext_var=True,
	group_indices={
		"A": [0, 1, 2, 3]
		"B": [4, 5, 6, 7]
		"C": [8, 9, 10, 11]
		"D": [12, 13, 14, 15]
		"E": [16, 17, 18, 19]
	},
	stats={
		"PCA3_A-B-C-D-E_loadings_raw" (3, 313)
		"PCA3_A-B-C-D-E_projections_raw" (20, 3)
		"PLS-RA_A-B-C-D-E_loadings_raw" (313, 2)
		"PLS-RA_A-B-C-D-E_projections_raw" (20, 2)
	}
)
```

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


### Exporting to Excel Spreadsheet
A `Dataset` instance can be exported to an Excel spreadsheet (`.xlsx`) format for easy viewing of the data, lipid 
identifications, and computed statistics.

```python
dset.export_xlsx('exported_dataset.xlsx')
```


### Dropping Features 
Features may be dropped from a `Dataset` instance according to user-defined criteria. Any feature not meeting the
specified criteria is removed from the `Dataset`, along with all corresponding identifications, statistics, _etc_. Valid
options for criteria are "mintensity", "meantensity", or any feature-length statistic already present in `Dataset.stats`.
For all criteria, one or both of the `upper_bound` and `lower_bound` kwargs must be set to reflect the desired selection.
The "mintensity" and "meantensity" criteria are two different intensity-based filtering techniques:
* "mintensity" - filtering such that at least one sample must have an intensity >= the specified `lower_bound` (the 
`normed` kwarg must be set to a boolean to indicate whether to use normalized or raw intensities)
* "meantensity" - filtering such that mean intensity across all samples must be >= the specified `lower_bound` or <= the 
specified `upper_bound`, or both

For filtering by the value of a computed statistic that contains more than 1 dimension, there is an additional kwarg 
(`axis`) which can be used to select the specific column to use. 

_Example:_ 
```python
# drop features that do not have at least one sample with normalized intensity > 1e3
dset.drop_features('mintensity', lower_bound=1e3, normed=True)

# drop features with average raw intensity < 1e3 or > 1e5
dset.drop_features('meantensity', lower_bound=1e3, upper_bound=1e5, normed=False)

# drop features with an ANOVA p-value > 0.05
dset.drop_features('ANOVA_A-B-C-D_normed', upper_bound=0.05)

# drop features that have a loading coefficient < 500 for principal component 2
dset.drop_features('PCA3_A-B-C-D_loadings_normed', lower_bound=500, axis=1)
```

**IMPORTANT:** *Dropping features cannot be undone, it is advisable to save the instance using `Dataset.save_bin(...)` 
prior to calling `Dataset.drop_features(...)` to be able to roll back the changes if needed.*


### Saving to File
A `Dataset` instance can be saved in a serialized binary format, retaining group assignments, normalization, computed
statistics, retention time calibration, and lipid identifications. Internally, this is done using Python's `pickle` 
library. 

```python
dset.save_bin('saved_dataset.pickle')
```


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

### Two-Group p-value
Computes the p-value from a comparison of two groups using a user-specified statistical test. The statistical test may
be a Student's t-test (`students`), Welch's t-test (`welchs`), or Mann-Whitney u-test (`mann-whitney`) depending on the
type of comparison and associated assumtions the user is willing to make about the two populations. 

```python
from lipydomics.stats import add_2group_pvalue

# compute p-values from a Welch's t-test using normalized data
add_2group_pvalue(dset, ['red', 'blue'], 'welchs', normed=True)
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
tol = [0.01, 0.1, 1.0]  # tolerance must be a list
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
| `meas_mz` | match only on measured m/z |
| `theo_mz_ccs` | match on theoretical m/z and CCS |
| `meas_mz_ccs` | match on measured m/z and CCS |
| `theo_mz_rt` | match on theoretical m/z and HILIC retention time |
| `measured_mz_rt` | match on measured m/z and HILIC retention time |
| `theo_mz_rt_ccs` | match on theoretical m/z, HILIC retention time, and CCS |
| `meas_mz_rt_ccs` | match on measured m/z, HILIC retention time, and CCS |
| `any` | start at the highest level (`meas_mz_rt_ccs`) then work down until an identification can be made |

_Example:_
```python
from lipydomics.identification import add_feature_ids

tol = [0.02, 0.2, 2.0]  # tol must be a list
# identify features at the highest level possible
add_feature_ids(dset, tol, level='any')
```
While this example demonstrates the use of the `any` identification level, any other single identification level
may be specified with the `level` kwarg to perform identifications only using a single confidence level. The `level` 
kwarg can also be a list of identification levels, which are tried in the order specified until an identification made
(essentially a custom version of the `any` identification level).

The resulting lipid identifications are stored in the `Dataset.feat_ids` instance variable as lists of putative IDs for 
each feature. The `Dataset.feat_id_levels` instance variable holds the identification level for each feature, and the 
`Dataset.feat_id_scores` instance variable holds the scores for each putative lipid ID. Multiple calls to 
`add_feature_ids` with different parameters will override results from previous calls. 

**IMPORTANT:** *CCS tolerance is in percent NOT absolute units*


### CCS and HILIC retention time prediction
The theoretical lipid database consists of CCS and HILIC retention time values generated using predictive models trianed
on experimental reference values. These predictive models are directly accessible via two convenience functions: 
`predict_ccs` and `predict_rt`. Both functions take as input a lipid (as defined by lipid class, fatty acid sum 
composition, fatty acid modifier, and MS adduct if relevant). Example:

```python
from lipydomics.identification import predict_ccs, predict_rt

# predict CCS for the lipid: PC(p34:3) [M+H]+
ccs = predict_ccs('PC', 34, 3, '[M+H]+', fa_mod='p')

# predict HILIC retention time for the same lipid, this time we do not need the MS adduct
rt = predict_rt('PC', 34, 3, fa_mod='p')
```


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

