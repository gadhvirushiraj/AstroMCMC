import emcee
import corner
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from multiprocessing import Pool, cpu_count
from norm_algo import norm_df
from log_prob import ln_posterior, ln_likelihood
from predict import param_interpol

def get_inputs():
    
    obs_file_path = input('\nEnter the path to the observed spectrum: ')

    params = input('\nEnter the parameter names: ').split()
    ndim = len(params)

    print('')
    param_range = []
    truth_val = []
    for i in range(ndim):
        param_range.append(list(map(float, input(f'Enter the range for parameter {params[i]}: ').split())))
        truth_val.append(float(input(f'Enter the true value for parameter {params[i]}: ')))
        print('')

    nwalkers = int(input('Enter the number of walkers: '))
    nsteps = int(input('Enter the number of steps: '))

    wave_min = int(input('\nEnter the minimum wavelength (open end): '))
    wave_max = int(input('Enter the maximum wavelength (open end): '))

    is_grid = int(input('\nDo you want to use a model grid(0) or parameter interpolation(1): '))

    print('\nHow do you want to initialize the walkers?')
    print('1. Self define initial parameters')
    print('2. Uniformly distributed in the parameter space')
    print('3. Use the MLE (takes time to find)')
    choice = int(input('Enter your choice: '))

    if choice == 1:
        initial_params = list(map(float, input('\nEnter the initial parameters (seperated by space): ').split()))
        spread = list(map(float, input('Enter the spread around the initial parameters (seperated by space): ').split()))
        pos = initial_params + spread * np.random.randn(nwalkers, ndim)
    elif choice == 2:
        param_range = np.array(param_range)
        spread = list(map(float, input('\nEnter the spread around the initial parameters (seperated by space): ').split()))
        pos = np.random.uniform(low=param_range[:,0], high=param_range[:,1], size=(nwalkers, ndim))
    elif choice == 3:
        spread = list(map(float, input('\nEnter the spread around the initial parameters (seperated by space): ').split()))
        pos = []

    return params, param_range, ndim, nwalkers, nsteps, wave_min, wave_max, spread, pos, truth_val, is_grid, obs_file_path


if __name__ == '__main__':

    params, param_range, ndim, nwalkers, nsteps, wave_min, wave_max, spread, pos, truth_val, is_grid, obs_file_path = get_inputs()

    # read observed spectrum file
    df_obs = pd.read_table(obs_file_path, delim_whitespace=True, header=None)
    df_obs.columns = ['wave', 'flux']

    # trim the observed spectrum in range of interest
    df_obs = df_obs[(df_obs['wave'] > wave_min) & (df_obs['wave'] < wave_max)]

    # telluric regions
    gaprange = [[8200,8390]]
    telluric_ranges = [[6860, 6960],[7550, 7750],[8200, 8430],[8930,9000]]
    telluric_ranges += gaprange

    # trim the observed spectrum
    for i in telluric_ranges:
        df_obs = df_obs[(df_obs['wave'] < i[0]) | (df_obs['wave'] > i[1])]

    # reset the index and normalize the flux
    df_obs = df_obs.reset_index(drop=True)
    df_obs = norm_df(df_obs)

    # calculate the error
    SNR = 32
    error_eso = df_obs['flux'] / SNR

    if len(pos) == 0:
        np.random.seed(42)
        nll = lambda *args: -ln_likelihood(*args)
        initial = np.array(truth_val) + 1e-1 * np.random.randn(3)
        soln = minimize(nll, initial, args=(df_obs['wave'], df_obs['flux'], error_eso, is_grid))
        teff_init, logg_init, m_init = soln.x

        pos = soln.x + [50, 0.2, 0.2] * np.random.randn(nwalkers, ndim)

        print("\nMaximum likelihood estimates:")
        print("Teff = {0:.3f}".format(teff_init))
        print("logg = {0:.3f}".format(logg_init))
        print("m = {0:.3f}".format(m_init))

    # cpu count
    ncpu = cpu_count()
    print('\nNumber of CPUs availble for multiprocessing:', ncpu)

    # warning message
    print("\nWarning -> If you have a data file with the same name as the default name (data.h5) in the same directory, it will be overwritten.\nPlease move or delete the old data file.")
    if input("type 'run' when ready: ").strip() == 'run':
        pass
    else:
        exit()

    # Set up the backend
    filename = "data.h5"
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    print('\nRunning MCMC...')
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_posterior, args=(df_obs['wave'], df_obs['flux'], error_eso, param_range, is_grid), pool=pool, backend=backend, moves=[(emcee.moves.DEMove(), 0.8),(emcee.moves.DESnookerMove(), 0.2),],)
        sampler.run_mcmc(pos, nsteps, progress=True, store=True)

    print('MCMC finished.')

    # acceptance fraction
    print('\nAcceptance fraction:', np.mean(sampler.acceptance_fraction),'\n')

    # autocorrelation time
    print('\nAutocorrelation analysics:', sampler.get_autocorr_time(quiet=True))

    # convergence plot
    samples = sampler.get_chain()
    fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
    fig.suptitle('Trace Plot')

    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], 'k', alpha=0.3)
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel('Step')
    for ax, label in zip(axes, params):
        ax.set_ylabel(label)

    plt.savefig('trace_plot.png')
    plt.show()
    
    # print final parameter values
    print('\nlast parameter values:')
    print(sampler.flatchain[-1],'\n')

    # make corner plot
    fig = corner.corner(sampler.flatchain,
                        labels=params,
                        show_titles=True,
                        range=param_range,
                        quantiles=[0.16, 0.5, 0.84],
                        title_kwargs={"fontsize": 12})
    
    plt.savefig('corner.png')
    plt.show()

    # print the parameter values and errors
    final_theta = []
    error = []
    for i in range(ndim):
        mcmc = np.percentile(sampler.get_chain(flat = True)[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(params[i], '--->',mcmc[1], '  err  ', -q[0],', ', q[1])
        final_theta.append(mcmc[1])
        error.append([q[0],q[1]])

    # make best fit graph
    fig = plt.figure(figsize=(15,5))
    best_fit = param_interpol(final_theta[0], final_theta[1], final_theta[2], df_obs['wave'])
    plt.plot(df_obs['wave'], df_obs['flux'], label="Observation")
    plt.plot(df_obs['wave'], best_fit, label= "Best Fit")
    plt.legend()
    plt.title('Best Fit')
    plt.savefig('best_fit.png')
    plt.show()

    # save best fit data to txt file
    best_fit_df = pd.DataFrame({'wave':df_obs['wave'], 'flux':best_fit})
    best_fit_df.to_csv('best_fit.txt', sep='\t', index=False)

    # best fit with the error band
    fig = plt.figure(figsize=(15,5))
    best_fit = param_interpol(final_theta[0], final_theta[1], final_theta[2], df_obs['wave'])
    lower_error = param_interpol(final_theta[0] - error[0][0], final_theta[1] - error[1][0], final_theta[2] - error[2][0], df_obs['wave'])
    upper_error = param_interpol(final_theta[0] + error[0][1], final_theta[1] + error[1][1], final_theta[2] + error[2][1], df_obs['wave'])
    
    plt.fill_between(df_obs['wave'], lower_error, upper_error, alpha=0.5, label='Error Band')
    plt.plot(df_obs['wave'], df_obs['flux'], label="Observation")
    plt.plot(df_obs['wave'], best_fit, label= "Best Fit")

    plt.legend()
    plt.title('Best Fit with Error Band')
    plt.savefig('best_fit_error.png')
    plt.show()

    # SNR error band and best fit
    fig = plt.figure(figsize=(15,5))
    error_best_fit = best_fit / SNR
    plt.fill_between(df_obs['wave'], best_fit - error_best_fit, best_fit + error_best_fit, alpha=0.5, label='Error Band (SNR = 32)')
    plt.plot(df_obs['wave'], best_fit, label= "Best Fit")
    plt.legend()
    plt.title('Best Fit with Error Band')
    plt.savefig('best_fit_snr.png')
    plt.show()