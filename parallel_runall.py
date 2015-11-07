import os
import sys

for arrayIndex in range(int(sys.argv[1])):
	print 'python -u Outrigger_Mapmaking.py ' + str(arrayIndex) + ' > ./Results/Logs/log_' + str(arrayIndex) + '.txt &'
	os.system('python -u Outrigger_Mapmaking.py ' + str(arrayIndex) + ' > ./Results/Logs/log_' + str(arrayIndex) + '.txt &')


# import os
# import numpy as np

# comps = np.array(['14', '13', '12', '11', '10', '09', '08', '07', '06', '05', '04', '03'])

# for arrayIndex in range(12):
# 	#print 'ssh eor-' + comps[0] +  '/bin/bash << cd ~/Outrigger_Project; ./single_array.sh ' + str(arrayIndex)
# 	os.system('ssh -t eor-' + comps[0] +  '"/bin/bash -c "cd ~/Outrigger_Project; ./single_array.sh ' + str(arrayIndex) + '""')
# 	comps = np.roll(comps, -1)