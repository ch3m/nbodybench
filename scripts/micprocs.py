#!/usr/bin/env python

from os import system
import sys
import subprocess
import curses
import re

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = "\033[1m"


def get_param(prompt_string):
     screen.clear()
     screen.border(0)
     screen.addstr(2, 2, prompt_string)
     screen.refresh()
     input = screen.getstr(10, 10, 60)
     return input

def execute_cmd(cmd_string):
     system("clear")
     a = system(cmd_string)
     print ""
     if a == 0:
          print "Command executed correctly"
     else:
          print "Command terminated with error"
     raw_input("Press enter")
     print ""

def get_phi_procs(cpus_load):
  p = subprocess.Popen('/opt/intel/mic/bin/micsmc -c', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  for line in p.stdout.readlines():
#    if re.search("^[[:blank:]]*Core \#", line):
    if re.search("^\s*Core \#", line):
      m = re.findall(r'\d+.\d+', line)
      cpus_load.append(m)
  retval = p.wait()

def print_cmd_cpus_load(cpus_load):
  total=100
  i = 1
  size_bar = 20
  for cpu in cpus_load:
    str = '\rcpu[{0}]'.format(i)
    num_cats = int( float(cpu[0]) * size_bar / total)
    str += '\t[{0:{width}}] {1:>2}%'.format('#' * num_cats , num_cats,width=size_bar)
    num_cats = int( float(cpu[1]) * size_bar / total)
    str += '\t[{0:{width}}] {1:>2}%'.format('#' * num_cats , num_cats,width=size_bar)
    num_cats = int( float(cpu[2]) * size_bar / total)
    str += '\t[{0:{width}}] {1:>2}%'.format('#' * num_cats , num_cats,width=size_bar)
    i += 1
    print str

cpus_load=[]
get_phi_procs(cpus_load)
print_cmd_cpus_load(cpus_load)

x = 0
#while x != ord('4'):
#  print '\r[{0:40}] {1:>2}%'.format('#' * int(progress * 40 /total), progress)
  
#     screen = curses.initscr()
#
#     screen.clear()
#     screen.border(0)
#     progress = 60
#     total = 100
#     screen.addstr(0, 0, '\r[{0:40}]{1:>2}%'.format('#' * int(progress * 40 /total), progress))
#     screen.addstr(4, 4, "1 - Add a user")
#     screen.addstr(5, 4, "2 - Restart Apache")
#     screen.addstr(6, 4, "3 - Show disk space")
#     screen.addstr(7, 4, "4 - Exit")
#     screen.refresh()
#
#     x = screen.getch()
#
#     if x == ord('1'):
#          username = get_param("Enter the username")
#          homedir = get_param("Enter the home directory, eg /home/nate")
#          groups = get_param("Enter comma-separated groups, eg adm,dialout,cdrom")
#          shell = get_param("Enter the shell, eg /bin/bash:")
#          curses.endwin()
#          execute_cmd("useradd -d " + homedir + " -g 1000 -G " + groups + " -m -s " + shell + " " + username)
#     if x == ord('2'):
#          curses.endwin()
#          execute_cmd("apachectl restart")
#     if x == ord('3'):
#          curses.endwin()
#          execute_cmd("df -h")
#
#curses.endwin()
try:
    sys.stdout.close()
except:
    pass
try:
    sys.stderr.close()
except:
    pass
