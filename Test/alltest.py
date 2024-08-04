import socket
import sys
import subprocess as sp 
import os
import time

def run_test(testname):
    dumtest = tempdir+testname
    fout = tempdir+testname+'.out'

    #remove if it exist
    if os.path.isdir(dumtest):
        sp.run(['rm', '-rf', dumtest])

    #copy entire directory over
    sp.run(['cp', '-r', testname, dumtest])

    os.chdir(dumtest)

    with open(fout,'w') as f:
        result = sp.run(['python3', runNLcom], stdout=f, stderr=f)

    os.chdir(testdir) #change back

    if result.returncode!=0:
        print('[run_test]: {:15s} has returned returncode = {:3d}'.format(testname.upper(), result.returncode))
        return False
    else:
        return True

if len(sys.argv)>=2:
    testnameL = sys.argv[1:]
else:
    #add stuff here...
    testnameL = ['simple', 'YoudinShu', 'simplesatellites']

testdir = os.getcwd()
runNLcom = testdir+'/../main/runcpd.py'

tempdir = os.path.expanduser('~')+'/Temp/NLtests/'
if not os.path.isdir(tempdir):
    sp.run(['mkdir', tempdir])


for testname in testnameL:
    t0 = time.time()
    print('[dotest.py]:proceeding with test: {:15s}'.format(testname.upper()))
    success = run_test(testname)
    t1 = time.time()
    if success:
        print('[dotest.py]:finished with         {:15s} after {:8.2f} sec'.format(testname.upper(), t1-t0))
    else:
        print('[dotest.py]failure: in test       {:15s} \033[91mFAILED \033[0m'.format(testname.upper()))

print('[dotest.py]:Please inspect {:25s} for detailed output'.format(tempdir))
