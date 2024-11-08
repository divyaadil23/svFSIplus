import os
import time

#os.system("cd /Users/divya/svFSIplus/tests/cases/struct/LV_HolzapfelOgden_TimeCompare ")
os.chdir('/Users/divya/svFSIplus/tests/cases/struct/LV_HolzapfelOgden_TimeCompare')
t0 = time.time()
os.system("/Users/divya/svFSIplus/build/svFSIplus-build/bin/svfsiplus svFSIplus.xml")
tf = time.time()

total = tf - t0

print(f'\tTime taken: {(total)*10**3:.03f} ms')