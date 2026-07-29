import pandas as pd
## USDF
my_path =  '/sdf/data/rubin/user/kexcel/'
my_plotspath = my_path + 'plots/'
# ------------------------------------
# DP2
repo = "dp2_prep"
collection = ['LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage1',
                  'LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage2',
                  'LSSTCam/runs/DRP/DP2/v30_0_6_rc1/DM-53881/stage3',
                  'LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage4']
skymap = 'lsst_cells_v2'
INCOLS = [
        'coord_ra',
        'coord_dec',
        'objectId',
        'detect_isIsolated',
        'refExtendedness',
        'tract', 'patch'
    ]
INCOLS += ['griz_model_extendedness']
bands = 'grizy'
for band in bands:
    INCOLS += [
        f'{band}_psfFlux',
        f'{band}_cModelFlux',
        f'{band}_cModelFluxErr',
        f'{band}_psfFluxErr',
        f'{band}_extendedness',
        f'{band}_extendedness_flag',
        f'{band}_psfFlux_flag'
    ]
    INCOLS += [f'{band}_model_extendedness']

# Euclid Q1
tile_to_coord = pd.read_csv('euclid_q1_tile_to_coord.csv')

# DP1
