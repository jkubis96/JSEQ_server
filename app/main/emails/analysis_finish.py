from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
import smtplib
import os
import sys


sys1=sys.argv[1]
sys2=sys.argv[2]
sys3=sys.argv[3]
to=sys.argv[4]
id=sys.argv[5]
project_name=sys.argv[6]



smtp = smtplib.SMTP('smtp.gmail.com', 587)
smtp.ehlo()
smtp.starttls()
smtp.login(str(sys1), str(sys2))


msg = MIMEMultipart()
msg['Subject'] = 'JBSDA - analysis finished'
text = 'Your analysis ' + str(project_name) + ' in JBioSystem - Platforma Database & Analytic has finished. Your validation code: \n<strong>' + str(id) + '</strong> \nNow go to the JBSDA platform and check the results in Validate tab\nGreetings JBioSystem'
msg.attach(MIMEText(text))



to = [str(to)]
smtp.sendmail(from_addr=str(sys3),
              to_addrs=to, msg=msg.as_string())
smtp.quit()