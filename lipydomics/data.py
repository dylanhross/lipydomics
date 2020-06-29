"""
    lipydomics/data.py
    Dylan H. Ross
    2019/02/03

    description:
        TODO
"""


import numpy as np
import pickle
from csv import reader
from pandas import DataFrame, concat, ExcelWriter

from lipydomics.util import abbreviate_sheet


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
            [[mz, rt, ccs], ... ] shape = (n_features, 3)
        and raw intensities (self.intensities):
            [[i0, i1, ... ], ... ] shape = (n_features, n_samples)
        Expects the input data file to be a .csv structured with the following columns:
            mz, rt, ccs, i_1, i_2, ..., i_n
        By default, 1 header line is skipped
    parameters:
        dataset_csv (str) -- filename of raw dataset in .csv format
        [skip_header (int)] -- number of header lines in the file to skip when loading data, passed on to 
                                np.genfromtxt [optional, default=1]
        [esi_mode (str)] -- electrospray ionization mode, can be used for identification, either 'neg', 'pos' or None
                            for unspecified [optional, default=None]
"""
        # store the name of the .csv file 
        self.csv = dataset_csv
        mz, rt, ccs, *intensities = np.genfromtxt(dataset_csv, delimiter=',', unpack=True, skip_header=skip_header)
        self.labels, self.intensities = np.array([mz, rt, ccs]).T, np.array(intensities).T
        # identifications can be added later
        self.feat_ids, self.feat_id_levels, self.feat_id_scores = None, None, None
        # store the number of features and samples in convenient instance variables
        self.n_features, self.n_samples = self.intensities.shape
        # group assignments and normalization can be done at some point
        self.group_indices = None
        self.normed_intensities = None
        # retention time calibration can be added at some point
        self.rt_calibration = None
        # statistics can be added in later
        self.stats = {}
        # ESI mode used for identification purposes
        self.esi_mode = esi_mode
        # store an external variable (for regression)
        self.ext_var = None


    @staticmethod
    def load_bin(bin_path):
        """
Dataset.load_bin
    description:
        Loads a Dataset instance from a serialized binary file (.pickle)

        This is a static method, so rather than using a Dataset instance, the uninitialized class is used e.g. :
            from lipydomics.data import Dataset
            dset = Dataset.load_bin('saved_dataset.pickle')
    parameters:
        bin_path (str) -- path to the serialized binary file
    returns:
        (Dataset) -- the dataset instance
"""
        with open(bin_path, 'rb') as pf:
            return pickle.load(pf)


    def save_bin(self, bin_path):
        """
Dataset.save_bin
    description:
        Saves this Dataset instance into a serialized binary format (.pickle) that can be reloaded again later
    parameters:
        bin_path (str) -- path to save the serialized binary file under
"""
        with open(bin_path, 'wb') as pf:
            pickle.dump(self, pf)


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


    def assign_groups_with_replicates(self, groups, n_replicates):
        """
Dataset.assign_groups_with_replicates
    description:
        Assign group identities to sample columns using a list of group names and number of replicates in each group.
        This method should ideally be applied to label all samples and the number of groups * number of replicates
        should be equal to the number of sample columns, BUT this is not enforced. A dictionary reflecting the group
        assignments (in order of the group names provided and with the specified number of replicates in each group) is
        produced and passed on to Dataset.assign_groups to actually make the group assignments
    parameters:
        groups (list(str)) -- list of group names
        n_replicates (int) -- number of samples in each group
"""
        group_indices = {}
        i = 0
        for group in groups:
            group_indices[group] = []
            for j in range(n_replicates):
                group_indices[group].append(i)
                i += 1
        self.assign_groups(group_indices)


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
        if type(norm_weights) is not np.ndarray:
            e = 'Dataset: normalize: norm_weights should be a numpy ndarray (type: {})'
            raise TypeError(e.format(type(norm_weights)))
        if norm_weights.shape != (self.n_samples,):
            e = "Dataset: normalize: norm_weights has shape: "
            e += "{} but should have shape: ({},)".format(norm_weights.shape, self.n_samples)
            raise ValueError(e)
        # apply normalization
        self.normed_intensities = np.multiply(self.intensities, norm_weights)


    def __repr__(self):
        """
Dataset.__repr__
    description:
        Prodces a string representation of this Dataset instance that provides information on the data contained and
        the stats analyses that have been performed.
    returns:
        (str) -- string representation of this Dataset instance
"""
        s = 'Dataset(\n'
        s += '\tcsv="{}",\n'.format(self.csv)
        s += '\tesi_mode="{}",\n'.format(self.esi_mode)
        s += '\tsamples={},\n'.format(self.n_samples)
        s += '\tfeatures={},\n'.format(self.n_features)
        s += '\tidentified={},\n'.format(len([_ for _ in self.feat_id_levels if _]) if self.feat_ids else 'False')
        s += '\tnormalized={},\n'.format(self.normed_intensities is not None)
        s += '\trt_calibrated={},\n'.format(self.rt_calibration is not None)
        s += '\text_var={},\n'.format(self.ext_var is not None)
        if self.group_indices is not None:
            s += '\tgroup_indices={\n'
            for name in self.group_indices:
                s += '\t\t"{}": {}\n'.format(name, self.group_indices[name])
            s += '\t},\n'
        else:
            s += '\tgroup_indices=None,\n'
        if not self.stats:
            s += '\tstats={}\n'
        else:    
            s += '\tstats={\n'
            for stat in self.stats:
                x = self.stats[stat]
                x = np.array(x) if type(x) is not np.ndarray else x
                s += '\t\t"{}" {}\n'.format(stat, x.shape)
            s += '\t}\n'
        s += ')'
        return s


    def select_feature_data(self, in_csv, out_csv, tolerance=(0.01, 0.1, 3.0)):
        """
Dataset.select_feature_data
    description:
        Searches features using values from an input .csv file and adds corresponding data for matching features to
        an output .csv file.
        * m/z tolerance is in Da, rt tolerance is in min, CCS tolerance is in percent *
    parameters:
        in_csv (str) -- input file path, .csv format
        out_csv (str) -- output file path, .csv format
        [tolerance (tuple(float))] -- tolerance to use for m/z, rt, and ccs search, CCS tolerance is a percentage not
                                        an absolute tolerance [optional, default=(0.01, 0.1, 3.)]
    returns:
        (bool) -- success (True), failure (no features found, False)
"""
        mzt, rtt, ccst = tolerance
        to_export = []
        with open(in_csv, 'r') as inf:
            next(inf)  # expects  a header row
            rdr = reader(inf)
            for mz, rt, ccs in rdr:  # iterate through query values
                mz, rt, ccs = float(mz), float(rt), float(ccs)
                for i in range(self.n_features): # iterate through all of the feature labels
                    mz_f, rt_f, ccs_f = self.labels[i]
                    if abs(mz - mz_f) <= mzt and abs(rt - rt_f) <= rtt and (100. * (abs(ccs - ccs_f) / ccs)) <= ccst:
                        # matched feature
                        line = [mz_f, rt_f, ccs_f] + self.intensities[i].tolist()
                        # check for normalized data
                        if self.normed_intensities is not None:
                            line += self.normed_intensities[i].tolist()
                        to_export.append(line)
        # check if any features were found
        if to_export != []:
            # determine the header(s)
            top_line = 'm/z,retention time,CCS,raw intensities' + ''.join([',' for _ in range(self.n_samples - 1)])
            if self.normed_intensities is not None:
                top_line += ',normalized intensities' + ''.join([',' for _ in range(self.n_samples - 1)])
            top_line += '\n'
            # write the output file
            with open(out_csv, 'w') as outf:
                outf.write(top_line)
                for line in to_export:
                    outf.write(','.join(str(_) for _ in line) + '\n')
            return True
        else:
            return False


    def export_xlsx(self, xlsx_path):
        """
Dataset.export_xlsx
    description:
        Exports this Dataset instance into an Excel Spreadsheet
    parameters:
        xlsx_path (str) -- path to save the spreadsheet under
"""
        # create a pandas DataFrame
        label_df = DataFrame(self.labels)
        int_df = DataFrame(self.intensities)
        df = concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
        writer = ExcelWriter(xlsx_path, engine='xlsxwriter')
        new_df = df
        columns = ['' for x in df.columns.values]
        if self.group_indices is not None:
            for group in self.group_indices:
                for i in self.group_indices[group]:
                    columns[i + 3] = group if columns[i + 3] == "" else columns[i + 3] + "/" + group
        columns = ['mz', 'rt', 'ccs'] + columns[3:]
        new_df.columns = columns
        new_df.to_excel(writer, sheet_name='Data')
        for key in self.stats:
            stats_df = DataFrame(self.stats[key])
            if "PCA3" in key and "loadings" in key:
                stats_df = stats_df.transpose()
            key = abbreviate_sheet(key) if len(key) > 31 else key
            stats_df.to_excel(writer, sheet_name=key)
        m = 0
        feat_dict = {}
        if self.feat_ids:
            for feat in self.feat_ids:
                if type(feat) is list:
                    m = max(m, len(feat))
        for i in range(0, m):
            s = []
            for feat in self.feat_ids:
                if type(feat) is list:
                    try:
                        s.append(feat[i])
                    except:
                        s.append("")
                else:
                    if i == 0:
                        s.append(feat)
                    else:
                        s.append("")
            feat_dict[i] = s
        cal_df = None
        if self.rt_calibration is not None:
            cal_rts = [self.rt_calibration.get_calibrated_rt(rt) for mz, rt, ccs in self.labels]
            cal_df = DataFrame(cal_rts)
            cal_df.columns = ['calibrated_ rt']
            cal_df.to_excel(writer, sheet_name="rt_calibration")
        if self.normed_intensities is not None:
            norm_df = DataFrame(self.normed_intensities)
            norm_df = concat([label_df, norm_df], axis=1, ignore_index=True, sort=False)
            norm_df.to_excel(writer, sheet_name="Normalized Intensities")
        if feat_dict:
            iden_df = DataFrame(feat_dict)
            level_df = DataFrame(self.feat_id_levels)
            if self.rt_calibration is not None:
                identification_df = concat([label_df, cal_df, level_df, iden_df], axis=1, ignore_index=True, sort=False)
                id_columns = ["putative_id" for x in identification_df.columns.values]
                id_columns = ['mz', 'rt', 'ccs', 'calibrated_rt', 'id_level', 'putative_id'] + id_columns[6:]
            else:
                identification_df = concat([label_df, level_df, iden_df], axis=1, ignore_index=True, sort=False)
                id_columns = ["putative_id" for x in identification_df.columns.values]
                id_columns = ['mz', 'rt', 'ccs', 'id_level', 'putative_id'] + id_columns[5:]
            identification_df.columns = id_columns
            identification_df.to_excel(writer, sheet_name='Identifications')
        writer.save()


    def drop_features_helper(self, target):
        """
Dataset.drop_features_helper
    description:
        A helper function that updates all relevant data arrays (i.e. any feature-length arrays) in the Dataset,
        dropping all of the features set to False in the target array
    parameters:
        target (np.array(bool)) -- a boolean array of shape (self.n_features) indicating which features (rows) should
                                    be kept in the Dataset (indicated by True)
"""
        # ensure target has the proper shape, i.e. (self.n_features,)
        if target.shape != (self.n_features,):
            m = 'Dataset: drop_features_helper: target array has shape {} but expected shape ({},)'
            raise ValueError(m.format(target.shape, self.n_features))

        # first drop features from the intensities matrix (and normalized intensities if present)
        self.intensities = self.intensities[target]
        if self.normed_intensities is not None:
            self.normed_intensities = self.normed_intensities[target]

        # drop features from the labels
        self.labels = self.labels[target]

        # drop identifications (if present)
        # these are lists with some mixed types, so need to iterate through them together with target
        if self.feat_ids is not None:
            feat_ids, feat_id_levels, feat_id_scores = [], [], []
            for fid, fidl, fids, keep in zip(self.feat_ids, self.feat_id_levels, self.feat_id_scores, target):
                if keep:
                    feat_ids.append(fid)
                    feat_id_levels.append(fidl)
                    feat_id_scores.append(fids)
            self.feat_ids = feat_ids
            self.feat_id_levels = feat_id_levels
            self.feat_id_scores = feat_id_scores

        # drop features from stats entries, we can assume they are all already numpy.ndarrays
        for stat_label in self.stats:
            # check that the statistic has column data
            stat_data = self.stats[stat_label]
            if self.n_features in stat_data.shape:
                if stat_data.shape[0] != self.n_features:
                    # need to transpose before indexing, then transpose back to preserve the shape
                    stat_data = stat_data.T[target].T
                else:
                    stat_data = stat_data[target]
                # update self.stats with the dropped feature data
                self.stats[stat_label] = stat_data

        # finally update counts of samples and features
        self.n_features, self.n_samples = self.intensities.shape


    def drop_features(self, criteria, lower_bound=None, upper_bound=None, normed=None, axis=None):
        """
Dataset.drop_features
    description:
        Uses user-specified criteria to selectively remove features from the Dataset. Valid criteria are:
            * "mintensity" - filtering such that at least one sample must have an intensity >= the specified lower_bound
                (the normed kwarg must be set to a boolean to indicate whether to use normalized or raw intensities)
            * "meantensity" - filtering such that mean intensity across all samples must be >= the specified lower_bound
                or <= the specified upper_bound, or both
                (the normed kwarg must be set to a boolean to indicate whether to use normalized or raw intensities)
            * any column-data (i.e. feature-length) statistic defined in Dataset.stats - The entry must be column data,
                meaning that at least one dimension must be Dataset.n_features in length. If the statistic has more than
                1 dimension, the axis kwarg must be set to the index of the column to use. Either lower_bound,
                upper_bound, or both must be set.
    parameters:
        criteria (str) -- specifies the criteria to filter on, which is used in conjunction with the upper_bound,
                            lower_bound, or both kwargs. Valid options include 'mintensity' and 'meantensity' for
                            intensity-based filtering, or any feature-length statistic already defined in Dataset.stats
        [lower_bound (None or float)] -- [optional, default=None]
        [upper_bound (None or float)] -- [optional, default=None]
        [normed (None or bool)] -- when relevant, whether to use normalized intensity data [optional, default=None]
        [axis (None or int)] -- when filtering on a statistic with more than one column (e.g. loadings from 3 component
                                PCA), the axis kwarg indicates which column to use (e.g. axis=0 for PC1 in loadings from
                                a 3 compounent PCA) [optional, default=None]
"""
        # first ensure valid parameters
        if criteria not in ['mintensity', 'meantensity']:
            stats_ok = True
            if criteria not in self.stats.keys():
                stats_ok = False
                m2 = 'criteria not defined in Dataset.stats'
            elif self.stats[criteria].shape != (self.n_features,) and axis is None:
                stats_ok = False
                m2 = 'the statistic has multiple columns but axis was not specified'
            if not stats_ok:
                m = 'Dataset: drop_features: invalid criteria "{}" ({})'.format(criteria, m2)
                raise ValueError(m)
        else:
            # the criteria is "mintensity" or "meantensity" so normed must be set to a boolean
            if normed is None:
                m = 'Dataset: drop_features: criteria is "mintensity" or "meantensity" so normed kwarg must be set'
                raise ValueError(m)
            if criteria == 'mintensity' and lower_bound is None:
                # lower_bound must be set for mintensity
                m = 'Dataset: drop_features: criteria "mintensity" requires lower_bound to be set'
                raise ValueError(m)
        if lower_bound is None and upper_bound is None:
            # at least one of lower_bound or upper_bound must be set, regardless of the criteria being used
            m = 'Dataset: drop_features: at least one of lower_bound or upper_bound kwargs must be set'
            raise ValueError(m)
        elif lower_bound is not None and upper_bound is not None:
            # if both bounds are set, upper must be > lower
            if lower_bound >= upper_bound:
                m = 'Dataset: drop_features: lower_bound must be < upper_bound'
                raise ValueError(m)
        if normed and self.normed_intensities is None:
            # if normed is set, then normalized intensities must be available
            m = 'Dataset: drop_features: normed set to True, but no normalization has been performed'
            raise ValueError(m)

        # mintensity
        if criteria == 'mintensity':
            # fetch intensities
            i = self.normed_intensities if normed else self.intensities
            # max is used here because AT LEAST ONE FEATURE must be >= to lower_bound NOT ALL FEATURES
            target = np.max(i, axis=1) >= lower_bound
            self.drop_features_helper(target)
        # meantensity
        elif criteria == 'meantensity':
            # fetch intensities
            i = self.normed_intensities if normed else self.intensities
            # set the bounds
            lb = lower_bound if lower_bound is not None else 0
            ub = upper_bound if upper_bound is not None else np.inf
            target = np.logical_and((np.mean(i, axis=1) > lb), (np.mean(i, axis=1) < ub))
            self.drop_features_helper(target)
        # stats
        else:
            # fetch the desired statistic
            stat = self.stats[criteria]
            # set the bounds
            lb = lower_bound if lower_bound is not None else 0
            ub = upper_bound if upper_bound is not None else np.inf
            # deal with multiple axes in statistic
            if stat.shape != (self.n_features,):
                # determine orientation of the statistic array, it may need to be transposed before indexing
                if stat.shape[0] == self.n_features:
                    stat = stat.T
                # index to get a single column
                stat = stat[axis]
            target = np.logical_and((stat > lb), (stat < ub))
            self.drop_features_helper(target)
