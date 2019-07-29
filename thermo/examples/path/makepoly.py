#!/usr/bin/env python
import sys
import re
from sys import argv
from qtc import iotools as io
from qtc import tctools as tc


if (len(argv) < 2):
    prompt = 'Specify MESS input or output file\n'
    userInp = raw_input(prompt)
else:
    userInp = argv[1]
if (len(argv) < 3):
    prompt = 'Specify stoichiometry:\n'
    stoich = raw_input(prompt)
else:
    stoich = argv[2]
if (len(argv) < 4):
    prompt = 'Specify 0 K heat of formation:\n'
    Hf0k = raw_input(prompt)
else:
    Hf0k = argv[3]
if 'Partition function (log) and its derivatives:' in io.read_file(userInp):
    messOut = io.read_file(userInp)
else:
    tc.run_pf('messpf', userInp)
    if '.inp' in userInp:
         messOut = io.read_file(userInp.replace('.inp','.dat'))
    else:
         messOut = io.read_file(userInp + '.dat')

io.write_file(messOut, 'pf.dat')
inp = tc.get_thermp_input(stoich,float(Hf0k))
tc.run_thermp(inp,'thermp.dat','pf.dat','thermp')

print('Writing thermochemical file {0}.\n'.format('thermp.out'))
#lines = io.read_file('thermp.out')
#Hf298 = ' h298 final\s*([\d,\-,\.]*)'
#Hf298 = re.findall(Hf298,lines)[-1]
#tc.run_pac99(stoich,'pac99')
#c97file = stoich + '.c97'
#if io.check_file(c97file):
#    c97text  = io.read_file(c97file)
#    las, has, msg = tc.get_coefficients(c97text)
#    chemkinfile = stoich + '.ckin'
#    print('Writing chemkin file {0}.\n'.format(chemkinfile))
#    chemininput = tc.write_chemkin_file(stoich, '', float(Hf0k), float(Hf298), stoich, 0, las, has, chemkinfile)
#else:
#     print('No {} produced, check thermp.out and new.groups'.format(c97file))


