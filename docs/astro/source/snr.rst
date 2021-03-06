Supernova Remnant Models
========================

Plot the evolution of radius of the SNR:


 .. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.units import Quantity
    from gammapy.astro.source import SNR, SNRTrueloveMcKee

    snr_models = [SNR, SNRTrueloveMcKee]
    densities = Quantity([1, 0.1], 'cm^-3')
    linestyles = ['-', '--']
    colors = ['b', 'g']
    t = Quantity(np.logspace(0, 5, 100), 'yr')
    
    for color, density in zip(colors, densities):
        for linestyle, snr_model in zip(linestyles, snr_models):
            snr = snr_model(n_ISM=density)
            plt.plot(t.value, snr.radius(t).to('pc').value,
                     label=snr.__class__.__name__ + ' (n_ISM = {0})'.format(density.value),
                     linestyle=linestyle, color=color)
    plt.xlabel('time [years]')
    plt.ylabel('radius [pc]')
    plt.legend(loc=4)
    plt.loglog()
    plt.show()
    

Plot the evolution of the flux of the SNR above 1 TeV and at 1 kpc distance:

 .. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.units import Quantity
    from gammapy.astro.source import SNR, SNRTrueloveMcKee

    densities = Quantity([1, 0.1], 'cm^-3')
    colors = ['b', 'g']

    t = Quantity(np.logspace(0, 5, 100), 'yr')
    
    for color, density in zip(colors, densities):
        snr = SNR(n_ISM=density)
        F = snr.luminosity_tev(t) / (4 * np.pi * Quantity(1, 'kpc') ** 2)
        plt.plot(t.value, F.to('ph s^-1 cm^-2').value, color=color, label='n_ISM = {0}'.format(density.value))
        plt.vlines(snr.sedov_taylor_begin.to('yr').value, 1E-13, 1E-10, linestyle='--', color=color)
        plt.vlines(snr.sedov_taylor_end.to('yr').value, 1E-13, 1E-10, linestyle='--', color=color)
    plt.xlim(1E2, 1E5)
    plt.ylim(1E-13, 1E-10)
    plt.xlabel('time [years]')
    plt.ylabel('flux @ 1kpc [ph s^-1 cm ^-2]')
    plt.legend(loc=4)
    plt.loglog()
    plt.show()