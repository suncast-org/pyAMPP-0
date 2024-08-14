import numpy as np
import pdb

def chromo_model(chromo, pbox, refmaps):
    csize = np.shape(chromo["nh"])
    dx = pbox["dr"][0]
    dy = pbox["dr"][1]
    dz = pbox["dr"][2]

    h = np.stack([np.reshape(chromo["dh"][..., 0:int(np.max(chromo["dh"][-2, -1]) + 1], (csize[1], csize[2], int(np.max(chromo["dh"][-2, -1]) + 1)), order='F'),
                              np.sum(chromo["dh"][..., 1:], axis=-1).squeeze() / 696000 / dz])

    idx = np.where((chromo["nne"] != 0) & (chromo["temp"] != 0))
    idx_idz = np.indices(idx)[0]

    z = h[..., np.newaxis, :] / dz
    c0size = np.shape(z)

    tr_cube = np.zeros((c0size[1], c0size[2], c0size[3]))
    tr_cube[idx] = chromo["dh"][idx][..., np.newaxis] / 696000 / dz

    izmax = np.max(idx_idz[:, :, -1]) + 1
    z = z[..., np.newaxis, :].astype(np.int32)
    z = np.where(z >= izmax, np.ceil(np.max(z)).reshape((c0size[1], c0size[2], 1)) + izmax - 1, z)

    csize = np.shape(z)
    deltaz = dz * (np.max(z) - np.min(z)) / csize[2]

    bcube = np.zeros((csize[0], csize[1], csize[2], 3))

    tr = np.sum(tr_cube, axis=-1)

    dh = chromo["dh"][..., 0:izmax] / 696000
    total_dh = np.sum(dh, axis=-1)
    hmax = (z[:, :, izmax] * dz).reshape((csize[0], csize[1], 1))
    dh[..., izmax] = hmax - total_dh

    for i in range(csize[0]):
        for j in range(csize[1]):
            for k in range(csize[2]):
                bcube[i, j, k, :] = interp2d(z[i, j, :], h[i, j, :], kind='linear')(z[i, j, np.newaxis].T)[0]

    t = idx.reshape((-1,))
    t = chromo["temp"][idx]
    n = idx.reshape((-1,))
    n = chromo["nne"][idx]
    nh = np.squeeze(chromo["nh"][idx])
    nhi = np.squeeze(chromo["nhi"][idx])
    np = np.squeeze(chromo["np"][idx])

    box = {
        "bcube": bcube,
        "refmaps": refmaps,
        "dr": [dx, dy, deltaz],
        "idx": idx,
        "n": n,
        "t": t,
        "n_htot": nh,
        "n_hi": nhi,
        "n_p": np,
        "dz": dh,
        "tr": tr,
        "chromo_layers": csize[2]
    }

    return box
