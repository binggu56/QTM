#/usr/bin/python

import os
from subprocess import call 



NCMDS = 2

message = str(input('type your message: '))

command1 = 'git commit -m '+message
push = 'git push origin master'
cmd_list = [command1, push]

for cmd in cmd_list:
    call(cmd, shell=True)

