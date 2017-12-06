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


parser = argparse.ArgumentParser()
parser.add_argument('iterations', help='number of iterations')
parser.add_argument('inputpath', help='path to data folder')
parser.add_argument('outputpath', help='path to output folder')
args = parser.parse_args()

filepath = os.path.dirname(args.inputpath)
fileoutpath = os.path.dirname(args.outputpath)
iterations = int(args.iterations)
os.makedirs(fileoutpath, exist_ok=True)

N = iterations // 100

gmail_user = 'johnensley17@gmail.com'
gmail_password = 'KiphGuzv0@Jehi'


filename = 'init_args.rds'


starttime = datetime.now()
subprocess.call(['Rscript', '02-createdata.R', filepath, str(15), str(50000)])
subprocess.call(['Rscript', '03-fitdamped.R', filepath, filename, fileoutpath, str(iterations)])

sendEmail(starttime, gmail_user, gmail_password, 'johnensley17@gmail.com', 'johnensley17@gmail.com', 1, 1, 1)
