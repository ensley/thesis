import os
import sys
import contextlib
import subprocess
import argparse



print('\n####################################')
print('#            PARAMETERS            #')
print('####################################\n')
family = input('--> Covariance family (currently supported: matern, dampedcos): ').lower()
if family == 'matern':
	nu = float(input('--> Shape parameter (nu): '))
	rho = float(input('--> Range parameter (rho): '))
	sigma = float(input('--> Scale parameter (sigma): '))
	nugget = float(input('--> Nugget: '))
	params = [nu, rho, sigma, nugget]
elif family == 'dampedcos':
	lam = float(input('--> Lambda parameter (>= 1): '))
	theta = float(input('--> Scale parameter (theta): '))
	sigma = float(input('--> Variance parameter (sigma): '))
	nugget = float(input('--> Nugget: '))
	params = [lam, theta, sigma, nugget]
else:
	print('unknown covariance family\n')
	sys.exit(0)


print('\n####################################')
if family == 'matern':
	print('Matern covariance function with parameters\n\tnu = ' + str(nu) + '\n\trho = ' + str(rho) + '\n\tsigma = ' + str(sigma) + '.')
elif family == 'dampedcos':
	print('Damped cosine covariance function with parameters\n\tlambda = ' + str(lam) + '\n\ttheta = ' + str(theta) + '\n\tsigma = ' + str(sigma) + '.')
print('Nugget of ' + str(nugget) + '.')
print('####################################\n')

defaults = input('--> Do you want to use the default settings? (y,n): ')

if defaults == 'y':
	nobs = 400
	ngrid = 50
	knots = 15
	B = 50000
	N = 100
	print('\n####################################')
	print('Defaults set.\nNumber of non-gridded observations: 400\nSize of grid for gridded observations: 50x50\nNumber of knots: 15\nNumber of Monte Carlo samples: 50000\nBatch size: 100')
	print('####################################\n')
else:
	nobs = int(input('--> Number of non-gridded observations to simulate (400 recommended): '))
	ngrid = int(input('--> Size of grid for gridded observations to simulate (50 recommended): '))
	knots = int(input('--> Number of knots (15 recommended): '))
	B = int(input('--> Number of Monte Carlo samples (50000 recommended): '))
	N = int(input('--> Number of data sets to generate: '))


seed = input('--> Do you want to set a random seed? (y,n): ')


with open('data.log', 'w') as f:
	f.write('The data in this folder was generated with the following properties:\n\n')
	if seed.lower() == 'y':
		f.write('Random seed set\n')
	f.write('{0} covariance function with parameters\n'.format(family))
	if family == 'matern':
		f.write('\tnu = {0}\n'.format(nu))
		f.write('\trho = {0}\n'.format(rho))
		f.write('\tsigma = {0}\n'.format(sigma))
		f.write('\tnugget = {0}\n'.format(nugget))
	elif family == 'dampedcos':
		f.write('\tlam = {0}\n'.format(lam))
		f.write('\ttheta = {0}\n'.format(theta))
		f.write('\tsigma = {0}\n'.format(sigma))
		f.write('\tnugget = {0}\n'.format(nugget))
	f.write('MCMC parameters are as follows:\n')
	f.write('\t{0} non-gridded observations simulated\n'.format(nobs))
	f.write('\t{0} x {0} grid of prediction points simulated\n'.format(ngrid))
	f.write('\t{0} knots\n'.format(knots))
	f.write('\t{0} Monte Carlo samples\n'.format(B))
	f.write('\t{0} data sets generated\n'.format(N))


for i in range(1, N):
	subprocess.run(['Rscript', 'generate_data.R', str(i), family, str(nobs), str(ngrid), str(knots), str(B), seed] + [str(p) for p in params])

