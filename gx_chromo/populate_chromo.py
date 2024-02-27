import numpy as np
from scipy.io import readsav
from pathlib import Path

def populate_chromo(model_mask):
    s = model_mask.shape

    paths = ['fontenla_b_v3.sav', 'fontenla_d_v3.sav', 'fontenla_f_v3.sav', 'fontenla_h_v3.sav', 'fontenla_p_v3.sav', 'fontenla_r_v3.sav', 'fontenla_s_v3.sav', 'eduard_v3.sav']

    chromo_models = [readsav(Path("gx_chromo/sav_precomputed/") / p) for p in paths]
    ch_len = len(chromo_models[0]["temp"])

    params = ("temp","nne","np","nh","nhi","dh","h")
    chromo = {k: np.zeros((s[0], s[1], ch_len), dtype=np.float32) for k in params}

    for i in range(len(chromo_models)):
        mask_i1 = (model_mask == (i+1))
        for k in params:
            chromo[k][mask_i1, :] = chromo_models[i][k]

    return chromo