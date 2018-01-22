import os
import sys
import contextlib
import subprocess
import argparse
import smtplib
from datetime import datetime


def sendEmail(starttime, usr, pw, fromaddr, toaddr, i, idx_start, idx_end):
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
    body = 'The process number ' + str(i) + ' (' + str(i - idx_start + 1) + ' of ' + str(idx_end - idx_start + 1) + ') has completed.'
    email_text = """\
    Date: %s
    From: %s
    To: %s
    Subject: %s

    %s

    Elapsed time: %s
    """ % (send_date, fromaddr, toaddr, subject, body, str(runtime))

    server.sendmail(fromaddr, toaddr, email_text)
    server.close()

    print('Email sent!')

# sendEmail(datetime.now(), 'johnensley17@gmail.com', 'put password here', 'johnensley17@gmail.com', 'johnensley17@gmail.com', 12, 10, 19)

parser = argparse.ArgumentParser()
parser.add_argument('iterations', help='number of iterations for each dataset')
parser.add_argument('start', help='index of first dataset to analyze (1-99)')
parser.add_argument('end', help='index of last dataset to analyze (1-99)')
parser.add_argument('inputpath', help='path to data folder')
parser.add_argument('outputpath', help='path to output folder')
args = parser.parse_args()

fileinpath = os.path.dirname(args.inputpath)
fileoutpath = os.path.dirname(args.outputpath)
iterations = int(args.iterations)
idx_start = int(args.start)
idx_end = int(args.end)
os.makedirs(fileoutpath, exist_ok=True)

N = iterations // 100

gmail_user = 'johnensley17@gmail.com'
gmail_password = 'put password here'

for i in range(idx_start, idx_end+1):
    curoutpath = os.path.join(fileoutpath, 'dataset' + str(i).zfill(3))
    os.makedirs(curpath)
    starttime = datetime.now()
    subprocess.call(['Rscript', '03-firstchunk.R', fileinpath, curoutpath, str(100)])
    for j in range(1, N+1):
        subprocess.call(['Rscript', '04-chunk.R', curoutpath, str(100), str(j)])
    subprocess.call(['Rscript', '05-processoutput.R', curoutpath])
    sendEmail(starttime, gmail_user, gmail_password, 'johnensley17@gmail.com', 'johnensley17@gmail.com', i, idx_start, idx_end)