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
acc_name=sys.argv[5]


smtp = smtplib.SMTP('smtp.gmail.com', 587)
smtp.ehlo()
smtp.starttls()
smtp.login(str(sys1), str(sys2))


msg = MIMEMultipart()
msg['Subject'] = 'JBSDA - account creation'
text = 'Hi ' + str(acc_name) + '!\nYour account on the JBioSystem - Database & Analytic platform has been created. From now on, you can freely start your own analyzes and/or explore the data available in the database.\nGreetings JBioSystem'
msg.attach(MIMEText(text))



to = [str(to)]
smtp.sendmail(from_addr=str(sys3),
              to_addrs=to, msg=msg.as_string())
smtp.quit()