import numpy as np
from pypowder.csa_conventions import fromhaeberlen
from pypowder.lineshape import powder_isotropic, filter1d


def frocsa_slice(x, diso, delta, eta, fwhm=2, scaling_factor=0.272):
    """
    Scaled CSA in an isotropic powder calculated in the frequency domain and
    convoluted with a Gaussian function.

    The scaling factor should be set from https://doi.org/10.1101/2020.07.02.184770

    """
    pas = fromhaeberlen(diso, delta * scaling_factor, eta)

    y = powder_isotropic(x, np.diag(pas))

    y = filter1d(x, y, 'gauss', fwhm=fwhm)

    return y / y.max()


def _get_data_slice(path, diso, section_width=0.6, verbose=False):
    # Read data

    dic, data = ng.bruker.read_pdata(path)
    data /= data.max()
    udic = ng.bruker.guess_udic(dic, data)
    uc0 = ng.fileiobase.uc_from_udic(udic, dim=0)
    ppm0 = uc0.ppm_scale()
    uc1 = ng.fileiobase.uc_from_udic(udic, dim=1)
    ppm1 = uc1.ppm_scale()

    # Determine noise
    noise = np.std(data[uc0(170, 'ppm'):uc0(160, 'ppm'), uc1(110, 'ppm'):uc1(100, 'ppm')])

    # Plot the 2D data
    levels = 6 * noise * 1.5 ** np.arange(1, 12)
    limits = ((220, 150), (190, 182))

    # Plot the 2D
    aspect = (limits[0][0] - limits[0][1]) / (limits[1][0] - limits[1][1])
    plt.contour(ppm1, ppm0, data, levels)
    plt.xlim(limits[1])
    plt.ylim(limits[0])
    plt.axvline(diso, linestyle=':', color='black')
    plt.show()

    # Cut out a section of the data
    section_ind = slice(uc1(diso + section_width / 2, 'ppm'), uc1(diso - section_width / 2, 'ppm'))
    section = data[:, section_ind]

    # Sum of the section to make a slice
    sliced = np.sum(section, axis=1)
    sliced /= sliced.max()

    return ppm0, sliced


if __name__ == "__main__":
    from lmfit import Model
    import matplotlib.pyplot as plt
    import nmrglue as ng
    import corner
    import os


    # Example data
    data_dir = os.path.join(os.path.dirname(__file__),
        'example_data', '600Nitro', 'fROCSA_His_CO', '84', 'pdata', '1')

    # Read in the included example data
    diso = 186.3
    section_width = 0.6
    ppm0, sliced = _get_data_slice(data_dir, diso=diso)
    offset = ppm0 - 186.3

    # fROCSA model
    # Estimate and set limits
    frocsa_model = Model(frocsa_slice)
    params = frocsa_model.make_params()
    params['diso'].set(value=0, vary=False)

    ##
    # Important:
    #
    # This data was collected with "1/2ROCSA(0.467, 0.033)-B"
    #
    ##
    params['scaling_factor'].set(value=0.177, vary=False)

    params['delta'].set(value=-65, min=-100, max=-20)
    params['eta'].set(value=0.5, min=0, max=1)
    params['fwhm'].set(value=2.0, min=1.0, max=3.0, vary=True)

    # Fit the data
    quick_result = frocsa_model.fit(sliced, x=offset, **params, nan_policy='omit', method = 'nadler')

    delta_q = quick_result.values['delta']
    params['delta'].set(value=delta_q, min=delta_q-20, max=delta_q+20)
    params['eta'].set(value=quick_result.values['eta'], min=0, max=1)

    # Plot Fit
    plt.figure(figsize=(3.25, 3.25/1.618))
    plt.plot(offset, sliced, color='black')
    plt.plot(offset, quick_result.best_fit, color='red')
    plt.minorticks_on()
    plt.xlabel('ppm')
    plt.xlim(offset[0], offset[-1])
    plt.show()

    print('Quick fit')
    print(quick_result.fit_report())

    # Estimate the uncertainties
    run_emcee = True
    if run_emcee:
        emcee_kws = dict(steps=5000, burn=500, thin=25, is_weighted=True, progress=True)
        result = frocsa_model.fit(sliced, x=offset, **params,
            method='emcee', fit_kws=emcee_kws, nan_policy='omit')

        # Plot Emcee Fit
        plt.figure(figsize=(3.25, 3.25 / 1.618))
        plt.plot(offset, sliced, color='black')
        plt.plot(offset, quick_result.best_fit, color='red')
        plt.minorticks_on()
        plt.xlabel('ppm')
        plt.xlim(offset[0], offset[-1])
        plt.show()

        # Plot the distribution of fit values
        emcee_plot = corner.corner(result.flatchain, labels=result.var_names,
            show_titles=True, title_kwargs={"fontsize": 12})
        plt.show()

        print('Emcee fit')
        print(result.fit_report())
