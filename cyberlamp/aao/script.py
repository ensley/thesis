import os
import sys
import contextlib
import subprocess
import argparse
import smtplib
from datetime import datetime


def sendEmail(starttime, usr, pw, fromaddr, toaddr):
    """
    Sends an email message through Gmail once the script is completed.  
    Developed to be used with AWS so that instances can be terminated 
    once a long job is done. Only works for those with Gmail accounts.
    
    starttime : a datetime() object for when to start run time clock
    usr : the Gmail username, as a string
    psw : the Gmail password, as a string 
    
    fromaddr : the email address the message will be from, as a string
    
    toaddr : a email address, or a list of addresses, to send the 
             message to
    """

    # calculate runtime
    runtime = datetime.now() - starttime

    # initialize SMTP server
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.ehlo()
    server.starttls()
    server.login(usr, pw)

    # send email
    send_date = datetime.strftime(datetime.now(), '%Y-%m-%d')
    subject = 'AWS Process Complete'
    body = 'The process has completed.'
    email_text = """\
    Date: %s
    From: %s
    To: %s
    Subject: %s

    %s

    %s
    """ % (send_date, fromaddr, toaddr, subject, body, str(runtime))

    server.sendmail(fromaddr, toaddr, email_text)
    server.close()

    print('Email sent!')


parser = argparse.ArgumentParser()
parser.add_argument('outputpath', help='path to output folder')
args = parser.parse_args()

filepath = os.path.dirname(args.outputpath)
os.makedirs(filepath, exist_ok=True)

# with contextlib.suppress(FileNotFoundError):
#	os.remove(os.path.join(filepath, 'allbetas.csv'))
#	os.remove(os.path.join(filepath, 'errbounds.csv'))

try:
	os.remove(os.path.join(filepath, 'allbetas.csv'))
except FileNotFoundError:
	pass

try:
	os.remove(os.path.join(filepath, 'errbounds.csv'))
except FileNotFoundError:
	pass


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
print('Output will be stored at: ' + filepath)
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
	nbatch = 100
	print('\n####################################')
	print('Defaults set.\nNumber of non-gridded observations: 400\nSize of grid for gridded observations: 50x50\nNumber of knots: 15\nNumber of Monte Carlo samples: 50000\nBatch size: 100')
	print('####################################\n')
else:
	nobs = int(input('--> Number of non-gridded observations to simulate (400 recommended): '))
	ngrid = int(input('--> Size of grid for gridded observations to simulate (50 recommended): '))
	knots = int(input('--> Number of knots (15 recommended): '))
	B = int(input('--> Number of Monte Carlo samples (50000 recommended): '))
	nbatch = int(input('--> Batch size (100 recommended): '))

niter = int(input('--> Iterations: '))
N = niter // nbatch

print('\n####################################')
print('Program will run in ' + str(N) + ' batches.')
print('####################################\n')

seed = input('--> Do you want to set a random seed? (y,n): ')

starttime = datetime.now()

subprocess.call(['Rscript', '02-createdata.R', filepath, family, str(nobs), str(ngrid), str(knots), str(B), seed] + [str(p) for p in params])
subprocess.call(['Rscript', '03-firstchunk.R', filepath, str(nbatch)])

for i in range(1, N):
	subprocess.call(['Rscript', '04-chunk.R', filepath, str(nbatch), str(i)])

subprocess.call(['Rscript', '05-processoutput.R', filepath])

gmail_user = 'johnensley17@gmail.com'
gmail_password = 'put password here'

sendEmail(starttime, gmail_user, gmail_password, 'johnensley17@gmail.com', 'johnensley17@gmail.com')
