"""
    encoder_params.py
    Dylan H. Ross
    2020/03/06

    description:
        defines all of the lipid classes, FA modifiers, and MS adducts to encode explicitly for CCS and RT prediction
        models
"""


# CCS prediction
ccs_lipid_classes = [
    'CE', 'Cer', 'DG', 'DGDG', 'FA', 'GlcCer', 'LPC', 'LPE', 'LPG', 'LPS', 'MGDG', 'PA', 'PC', 'PE',
    'PG', 'PI', 'PS', 'SM', 'TG',  # these all have ~100 or more samples in the measured database
    'LysylPG', 'AcylPG', 'AcylPE'  # manually add in some bacterial specific classes
]
ccs_fa_mods = ['d', 'o', 'p']
ccs_ms_adducts = [
    '[M-2H]2-', '[M-H]-', '[M+2K]2+', '[M+2Na-H]+', '[M+CH3COO]-', '[M+Cl]-', '[M+H-H2O]+', '[M+H]+', '[M+HCOO]-',
    '[M+Na]+', '[M+NH4]+'
]


# Retention Time prediction
rt_lipid_classes = [
    'PG', 'DGDG', 'PE', 'LPE', 'PE', 'PC', 'CL', 'PI', 'PA', 'Cer', 'AcylPG', 'LysylPG', 'GlcCer', 'MGDG', 'LPG',
    'AcylPE', 'SM', 'LPC', 'PS', 'DG', 'GlcADG', 'AlaPG', 'PIP'
]
rt_fa_mods = ['d', 'p']
