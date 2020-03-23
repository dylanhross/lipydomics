"""
    example_standard_analysis.py
    Dylan H. Ross
    2020/03/23

    description:
        This is an example of an analysis script that can apply a standard set of analyses to a dataset. The analysis
        steps in this example are as follows:

            1. load the dataset
            2. assign groups
            3. normalize intensities to an external standard
            4. compute a PCA on the normalized data
            5. generate a plot of the PCA projections
            6. identify lipids
            7. export the dataset into an Excel spreadsheet (to review lipid identifications)
            8. save the Dataset in a binary format that can be loaded again later for further manual analysis
            9. print an overview of the dataset to the console

        This script assumes the following file structure before running:

            example/                                     <-- top-level directory containing the data
                figures/                                 <-- a directory to save generated plot image files into
                example_pos.csv                          <-- the raw lipidomics data
                example_standard_analysis.py             <-- this script
                example_weights.txt                      <-- a text file with the external weights for each sample

        After running the analysis script, the file structure will look like:

            example/                                     <-- top-level directory containing the data
                figures/                                 <-- a directory to save generated plot image files into
                    PCA3_wt-ko1-ko2-ko12_projections_normed.png <-- a plot of the PCA projections
                example_pos.csv                          <-- the raw lipidomics data
                example_pos.pickle                       <-- the binary file that can be loaded for further analysis
                example_pos.xlsx                         <-- exported spreadsheet
                example_standard_analysis.py             <-- this script
                example_weights.txt                      <-- a text file with the external weights for each sample

        This script can be modified to handle input/output in a flexible fashion (e.g. using argparse or sys.argv), but
        that is beyond the scope of this example.
"""


from numpy import genfromtxt as gft

from lipydomics.data import Dataset
from lipydomics.stats import add_pca3
from lipydomics.plotting import scatter_pca3_projections_bygroup
from lipydomics.identification import add_feature_ids


# Step 1. load the dataset
dset = Dataset('example_pos.csv', esi_mode='pos')

# Step 2. assign groups
# this set has 4 groups and 3 replicates each
groups = ['wt', 'ko1', 'ko2', 'ko12']
dset.assign_groups_with_replicates(groups, 3)

# Step 3. normalize intensities to an external standard
# load the weights into a numpy array from the text file
weights = gft('weights.txt')
dset.normalize(weights)

# Step 4. compute a PCA on the normalized data
add_pca3(dset, groups, normed=True)

# Step 5. generate a plot of the PCA projections
scatter_pca3_projections_bygroup(dset, groups, 'figures/', normed=True)

# Step 6. identify lipids
# we will use the 'any' level in order to get lipid identifications from measured or theoretical data
# this example will assume different retention times, so we will also include the use_rt=False parameter
# we will use a mass tolerance of 0.02 Da
# the retention time is ignored anyways so we will set the tolerance to 1 min
# we will use the standard 3 % CCS tolerance
tol = [0.02, 1.0, 3.0]
add_feature_ids(dset, tol, level='any', use_rt=False)

# Step 7. export the dataset into an Excel spreadsheet (to review lipid identifications)
dset.export_xlsx('example_pos.xlsx')

# Step 8. save the Dataset in a binary format that can be loaded again later for further manual analysis
dset.save_bin('example_pos.pickle')

# Step 9. print an overview of the dataset to the console
print(dset)
