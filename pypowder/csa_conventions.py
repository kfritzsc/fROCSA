import numpy as np


def fromhaeberlen(diso, delta, eta):
    """
    Converts from the haeberlen convention to the PAS CSA tensor.

    Parameters
    ----------
    diso : float
         diso == 1/3(d11 + d22 + d33)
    delta : flaot
         dzz - diso
    eta : float
       (dyy - dxx)/delta (0 <= eta <= 1)

    Returns
    ----------
    tensor
    """
    if delta > 0.0:
        d11 = diso + delta
        d22 = diso - delta * (1 - eta) / 2
        d33 = diso - delta * (1 + eta) / 2
        return np.diag([d11, d22, d33])

    # delta < 0.0
    d11 = diso - delta * (1 + eta) / 2
    d22 = diso - delta * (1 - eta) / 2
    d33 = diso + delta
    return np.diag([d11, d22, d33])


def tohaeberlen(tensor, **kwargs):
    """
    Converts from the PAS CSA tensor to haeberlen convention.

    Parameters
    ----------
    tensor : 1x3 or 3x3 diagonal slicable of floats
        (d11, d22, d33)
    diso : float
        If for whatever reason diso != 1/3(d11 + d22 + d33).

    Returns
    ----------
    diso, delta, eta
    """
    shape = np.shape(tensor)
    if shape == (3,):
        pass
    elif shape == (3, 3) and isdiagonal(tensor):
        tensor = np.diag(tensor)
    else:
        raise ValueError('Tensor should have shape (3,) or be diag')

    # Get/Calculate delta isotropic (diso).
    try:
        diso = kwargs['diso']
    except KeyError:
        diso = np.sum(tensor, axis=-1) / 3.

    # |dzz-diso| > |dxx-diso| > |dyy-diso|
    dif = np.abs(tensor[:3] - diso)
    order = np.argsort(dif, axis=-1)
    dxx, dyy, dzz = (tensor[order[1]],
                     tensor[order[0]],
                     tensor[order[-1]])

    # Calculate delta and diso
    delta = dzz - diso
    eta = (dyy - dxx)/delta
    return diso, delta, eta


def fromherzfeld_berger(diso, span, skew):
    """
    Converts from the haeberlen convention to the PAS CSA tensor.

    Parameters
    ----------
    diso : float
         diso == 1/3(d11 + d22 + d33)
    span : flaot
         d33-d11
    skew : float
       3(d22 - diso)/span (-1 <= skew <= 1)

    Returns
    ----------
    tensor
    """
    d22 = diso + (span * skew)/3
    d33 = (3 * diso - d22 - span)/2
    d11 = 3 * diso - d22 - d33
    return np.diag([d11, d22, d33])


def toherzfeld_berger(tensor, **kwargs):
    """
    Converts from the PAS CSA tensor to herzfeld-berger convention.

    Parameters
    ----------
    tensor = (d11, d22, d33): array of floats.
    diso : float
        If for whatever reason diso != 1/3(d11 + d22 + d33).

    Returns
    ----------
    diso, span and skew
    """
    shape = np.shape(tensor)
    if shape == (3,):
        pass
    elif shape == (3, 3) and isdiagonal(tensor):
        tensor = np.diag(tensor)
    else:
        raise ValueError('Should have shape (3,) or diag.')

    # Get/Calculate delta isotropic (diso).
    try:
        diso = kwargs['diso']
    except KeyError:
        diso = np.sum(tensor, axis=-1) / 3.

    # Calculate span and skew
    span = tensor[2] - tensor[0]
    skew = (3 * (tensor[1] - diso)) / span
    return diso, span, skew


def isdiagonal(matrix):
    """
    True if matrix is diagonal or a zero matrix. Otherwise False.

    Paramaters
    ----------
    matrix : ndarray
    """
    shape = np.shape(matrix)
    if shape[0] != shape[1]:
        return False
    matrix -= np.diag(matrix)
    if not matrix.sum() <= (np.size(matrix) * np.spacing(1)):
        return False
    return True