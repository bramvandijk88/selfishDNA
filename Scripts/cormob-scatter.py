import numpy as np
import matplotlib.pyplot as plt
import math
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--corenr")
args = parser.parse_args()
#print(args)
n = args.corenr
Core = "data/coregene_mobilities/Coregene"+n+".dat"




def main():



    # Toxin = "data/toxin_mobilities/Toxin"+n+".dat"
	titelcore = "Core " + n
	parsefile(Core, c = 'goldenrod', titel = titelcore)
	plt.show()


def parsefile(f, c, titel):

	#go over every line (=timepoint) of a toxmobdatafile
	times = []
	data = []
	#timedatatuples =[]
	pattern_times= "^\d+"
	pattern_data= "\d\.\d+\w*\-*\d*"
	#pattern_split= "\s*,"
	file = open (f,"r")

#Gather data from file
	for line in file:

		times_now = [int(s) for s in re.findall(pattern_times, line)]
		times.append(times_now[0])
		data_now =  [float(i) for i in re.findall(pattern_data, line)]
		data.append(data_now)
		#timedatatuples.append((times_now[0],data_now)) 					# saves data in form:[(time (i), [mob(i), mob(i), (mob(i)])]

		#split option (does not work yet!):
		#split

		# split_data = (pattern_split, line)
		# times.append([int(s) for s in split_data[0]])
		# data.append([float(i) for i in split_data[1:]])
		# print "split_data",split_data



	file.close



	#print timedatatuples

	#Plotten
	x = times
	y = data

	plt.figure()
	for xe, ye in zip(x,y):
			plt.scatter([xe] * len(ye), ye, color = c, s=1,alpha=0.1)
			plt.axis([0,times[-1], -0.05, 1.05])
			plt.xlabel('Time')
			plt.ylabel('Mobility')
			plt.title(titel)




	print f





main()
