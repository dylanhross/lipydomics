"""
    lipydomics/data.py
    Dylan H. Ross
    2019/02/03

    description:
        TODO
"""


import numpy as np


class Dataset:
    """
Dataset
    description:
        A class for encapsulating a lipidomics dataset, with functions for loading and storing the data
"""

    def __init__(self, dataset_csv, skip_header=1, esi_mode=None):
        """
Dataset.__init__
    description:
        Loads a dataset from .csv file

        Uses numpy.genfromtxt to load in a raw dataset, stored in .csv format. Stores two numpy.ndarrays containing
        the labels (self.labels):
            [[mz, rt, dt], ... ] shape = (n_features, 3)
        and raw intensities (self.intensities):
            [[i0, i1, ... ], ... ] shape = (n_features, n_samples)
        Expects the input data file to be a .csv structured with the following columns:
            mz, rt, dt, i_1, i_2, ..., i_n
        By default, 1 header line is skipped
    parameters:
        dataset_csv (str) -- filename of raw dataset in .csv format
        [skip_header (int)] -- number of header lines in the file to skip when loading data, passed on to 
                                np.genfromtxt [optional, default=1]
        [esi_mode (str)] -- electrospray ionization mode, can be used for identification, either 'neg', 'pos' or None
                            for unspecified [optional, default=None]
"""
        mz, rt, dt, *intensities = np.genfromtxt(dataset_csv, delimiter=',', unpack=True, skip_header=skip_header)
        self.labels, self.intensities = np.array([mz, rt, dt]).T, np.array(intensities).T
        # store the number of features and samples in convenient
        self.n_features, self.n_samples = self.intensities.shape
        # group assignments and normalization can be done at some point
        self.group_indices = None
        self.normed_intensities = None
        # statistics can be added in later
        self.stats = {}
        # ESI mode used for identification purposes
        self.esi_mode = esi_mode

    def assign_groups(self, group_indices):
        """
Dataset.assign_groups
    description:
        Assign group identities to a subset or all of the sample columns. Group identities are assigned using a 
        dict mapping the group name to indices of the sample columns that are members of that group, e.g.:
            intensities = [[i0, i1, i2, i3, i4, i5], ... ]
            group_indices = {"A": [0, 1, 2], "B": [3, 4, 5]} 
            --> first three samples belong to group "A" and last three belong to group "B"
            OR group indices = {"wt": [0, 2], "ko": [1, 5]} 
            --> only a subset are labeled
        This method may be called multiple times to esablish multiple sub-groupings if necessary, duplicate group names 
        will override previous entries in this usage
        Group indices are stored in the self.group_indices instance variable
        WARNING: validity of group indices is not checked here so if for instance the value of one of the indices is 
                 larger than n_samples - 1, an error will occur later if you try to get the data for that group
    parameters:
        group_indices (dict(str:list(int))) -- dictionary of group names mapped to the indices of sample columns 
                                                belonging to each group
"""
        # initialize a new dictionary if groups havent been assigned before
        if not self.group_indices:
            self.group_indices = {}
        # add all of the user-defined groupings
        for group_name in group_indices:
            self.group_indices[group_name] = group_indices[group_name]

    def get_data_bygroup(self, group_names, normed=False):
        """
Dataset.get_data_bygroup
    description:
        Returns a subset (or multiple subsets) of the dataset as defined by group names. The group indices must be 
        defined ahead of time. Normalized data can be fetched if normalize() was called ahead of time.
    parameters:
        group_names (str OR list(str)) -- a single group name or list of group names to retrieve intensity data for
        [normed (bool)] -- if True, return the normalized intensities (requires that self.normalize(...) be called 
                            ahead of time), otherwise return the raw intensities [optional, default=False]
    returns:
        (np.ndarray) -- array of intensities for sample columns belonging to the specified group, shape = 
                        (n_features, group_size). multiple arrays will be returned if multiple group names were 
                        provided as input
"""
        if not self.group_indices:
            e = "Dataset: get_data_bygroup: groups not defined"
            raise ValueError(e)
        if normed and self.normed_intensities is None:
            e = "Dataset: get_data_bygroup: normed set to True but self.normed_intensities not defined"
            raise ValueError(e)
        if type(group_names) == str:
            if group_names not in self.group_indices:
                e = "Dataset: get_data_bygroup: group name '{}' not defined in self.groups".format(group_names)
                raise ValueError(e)
            if normed:
                return self.normed_intensities[:, self.group_indices[group_names]]
            else:
                return self.intensities[:, self.group_indices[group_names]]
        elif type(group_names) == list:
            for name in group_names:
                if name not in self.group_indices:
                    e = "Dataset: get_data_bygroup: group name '{}' not defined in self.groups".format(name)
                    raise ValueError(e)
            if normed:
                return (self.normed_intensities[:, self.group_indices[name]] for name in group_names)
            else:
                return (self.intensities[:, self.group_indices[name]] for name in group_names)
        else:
            e = "Dataset: get_data_bygroup: group_names should be either str or list "
            e += "(got type: {})".format(type(group_names))
            raise TypeError(e)

    def normalize(self, norm_weights):
        """
Dataset.normalize
    description:
        Normalize sample intensities according to a user-supplied array of weights, stores normalized intensities in
        self.normed_intensities
        norm_weights must be a numpy.array of shape: (n_samples,)
    parameters:
        norm_weights (np.ndarray) -- array of weights to use to normalize each sample, must be of shape (n_samples, )
"""
        if norm_weights.shape != (self.n_samples,):
            e = "Dataset: normalize: norm_weights has shape: "
            e += "{} but should have shape: ({},)".format(norm_weights.shape, self.n_samples)
            raise ValueError(e)
        # apply normalization
        self.normed_intensities = np.multiply(self.intensities, norm_weights)

