import numpy as np
import matplotlib.pyplot as plt
import math
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--toxinnr")
args = parser.parse_args()
#print(args)
n = args.toxinnr
Toxin = "data/toxin_mobilities/Toxin"+n+".dat"
Antitoxin = "data/antitoxin_mobilities/Antitoxin"+n+".dat"



def main():



    # Toxin = "data/toxin_mobilities/Toxin"+n+".dat"
	titeltox = "Toxin " + n
	parsefile(Toxin, c = 'firebrick', titel = titeltox)
	titelantitox = "Antitoxin " + n
    # Antitoxin = "data/antitoxin_mobilities/Antitoxin"+n+".dat"
	parsefile(Antitoxin, c = 'navy', titel = titelantitox)
	plt.show()

	# for i in range(20):
	# 	parsefile("Toxin%i.dat")
	#




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
			plt.scatter([xe] * len(ye), ye, color = c, s=1,alpha=0.05)
			plt.axis([0,times[-1], -0.05, 1.05])
			plt.xlabel('Time')
			plt.ylabel('Mobility')
			plt.title(titel)

	# fig, ax = plt.subplots()
	# for e in timedatatuples:
	# 	if times > 200000:
	# 		break
	# 	t, mobs = e
	# 	for m in mobs:
	# 		ax.scatter(t, m)

    # ax.scatter(x,y)

	# ax.legend()
	# ax.grid(True)
	#plt.plot(times, data)


	print f





main()



#plot time[i]
#for i in range (#nr_of_lines)
#x = times
#y = data  # in list form

#plt.plot(x, y, ro)

#plt.show()

#plt.savefig("Toxin" ++ i ++".png")
#plt.close()
