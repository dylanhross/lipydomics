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
    'TG', 'PC', 'PE', 'DG', 'PG', 'SM', 'PI', 'PS', 'LPE', 'SM', 'CE', 'LPC', 'PA', 'PE', 'PS', 'Cer', 'GlcCer', 'LPG',
    'DGDG', 'CL', 'LysylPG', 'AcylPG', 'LPI', 'LPS'
]
ccs_fa_mods = ['d', 'o', 'p']
ccs_ms_adducts = [
    '[M+NH4]+', '[M+H]+', '[M+HCOO]-', '[M-H]-', '[M+Na]+', '[M-CH3]-', '[M+CH3COO]-', '[M+2Na-H]+', '[M-2H]2-',
    '[M+Cl]-', '[M+2K]2+', '[M+H-H2O]+'
]


# Retention Time prediction
rt_lipid_classes = [
    'PG', 'DGDG', 'PE', 'LPE', 'PE', 'PC', 'CL', 'PI', 'PA', 'Cer', 'AcylPG', 'LysylPG', 'GlcCer', 'MGDG', 'LPG',
    'AcylPE', 'SM', 'LPC', 'PS', 'DG', 'GlcADG', 'AlaPG', 'PIP'
]
rt_fa_mods = [
    'd', 'p'
]
