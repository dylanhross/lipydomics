"""
    lipydomics/identification/__init__.py
    Dylan H. Ross
    2019/10/04

    description:
        Module for performing identification of individual lipid features
"""


from numpy import array


def add_identifications(dataset):
    """
add_identifications
    description:
        Goes through the list of features (mz, rt, ccs) in Dataset.labels and provides an identification for each, 
        storing the identifications in Dataset.identifications. If a feature is unable to be identified using the
        specified criteria, then it is given a generic fixed-width feature name based on mz, rt, and CCS:
            'UNK_{:09.4f}_{:05.2f}_{:06.2f}'.format(mz, rt, ccs)
        
        * if identifications have already been made, subsequent calls to this function will override previous results *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
"""
    ids = []
    for mz, rt, ccs in dataset.labels:

        # try other identification strategies ...

        feat_id = 'UNK_{:09.4f}_{:05.2f}_{:06.2f}'.format(mz, rt, ccs)
        ids.append(feat_id)

    # apply the identifications to the Dataset
    dataset.identifications = array(ids)

