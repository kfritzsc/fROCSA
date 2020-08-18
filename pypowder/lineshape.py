import numpy as np
import scipy.special


def powder_isotropic(omega, pas):
    """
    Frequency domain calculation over an isotropic powder for CSA tensor.

    The expressions are evaluated using complete elliptic integrals of the first kind.

    Parameters
    ----------
    omega : array
        frequency
    pas : array
        [omega_11, omega_22, omega_33]

    Returns
    -------
    intensity(omega) array
        for an isotropic powder
    """
    # Ensure  \omega_11 < omega_{22} \le omega_{33}
    pas = sorted(pas)

    # \omega_{11} < \omega > \omega_{33}
    intensity = np.zeros_like(omega)

    # \omega_{33} \ge \omega > \omega_{22}
    where = ((pas[2] >= omega) & (omega > pas[1]))
    m = (((pas[1] - pas[0]) * (pas[2] - omega[where])) /
         ((pas[2] - pas[1]) * (omega[where] - pas[0])))
    a = np.pi * np.sqrt((omega[where] - pas[0]) * (pas[2] - pas[1]))
    k = scipy.special.ellipk(m)
    intensity[where] = k / a

    # \omega_{22} > \omega \ge \omega_{11}
    where = ((pas[1] > omega) & (omega >= pas[0]))
    m = (((omega[where] - pas[0]) * (pas[2] - pas[1])) /
         ((pas[2] - omega[where]) * (pas[1] - pas[0])))
    k = scipy.special.ellipk(m)
    a = np.pi * np.sqrt((pas[2] - omega[where]) * (pas[1] - pas[0]))
    intensity[where] = k / a
    return intensity


# Line shape simulations.
def sim_gauss_fwhm(x, x0, fwhm):
    return np.exp(-(x - x0) ** 2 * 4 * np.log(2) / (fwhm ** 2))


def sim_lorentz_fwhm(x, x0, fwhm):
    return (0.5 * fwhm) ** 2 / ((0.5 * fwhm) ** 2 + (x - x0) ** 2)


def lorentz_kernel(x, fwhm, kernel_length=10):
    """
    :param: fwhm : Full-Width-Half-Max of Lorentzian.
    """
    step = abs(x[1] - x[0])
    if step < fwhm:
        limit = round(kernel_length * fwhm) * step
        kernel_x = np.arange(-limit, limit + step, step)

        return sim_lorentz_fwhm(kernel_x, 0, fwhm)
    else:
        raise ValueError("FWHM can't be < step size of input.")


# Generate a kernel from the line-shapes
def gauss_kernel(x, fwhm, kernel_length=10):
    """
    :param: fwhm : Full-Width-Half-Max of Gaussian.
    """
    step = abs(x[1] - x[0])
    if step < fwhm:
        limit = round(kernel_length * fwhm) * step
        kernel_x = np.arange(-limit, limit + step, step)
        return sim_gauss_fwhm(kernel_x, 0, fwhm)
    else:
        raise ValueError("FWHM can't be < step size of input.")


def filter1d(x, y, filter_type='gauss', **kwargs):
    """
    Apply a convolution of y with a line-shape function defined by a
    FWHM in units of x.

    Parameters
    ----------
    :param x: array of length n, position
    :param y: array of length n, intensity
    :param filter_type: 'lorentz', 'gauss', 'voigt' or 'pvoigt' filter-type
    :param fwhm: Full-Width-Half-Max of Lorentzian or Gaussian
    :param fwhm_g:Full-Width-Half-Max of Gaussian part of Voigt or
           psudo-Voigt
    :param fwhm_l: Full-Width-Half-Max of Lorentzian part of Voigt or
           psudo-Voigt
    :param kernel_length: int to indicate ~ kernel length / 2. A
           reasonable default is set for each function but you do
           you.

    Returns
    --------
    :return array of length n

    """

    filter_dispatch = {'gauss': gauss_kernel,
                       'lorentz': lorentz_kernel}

    filter_function = filter_dispatch[filter_type]
    kernel = filter_function(x, **kwargs)
    kernel /= kernel.sum()
    return np.convolve(y, kernel, mode='same')
