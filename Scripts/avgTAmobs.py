import numpy as np
import matplotlib.pyplot as plt
import math
import re
import argparse
import matplotlib.cm as cmx
import matplotlib.colors as colors
from datetime import datetime


# parser = argparse.ArgumentParser()
# parser.add_argument("-n", "--toxinnr")
# args = parser.parse_args()
#print(args)
# n = args.toxinnr
import_file = "log.out"

startTime = datetime.now()

colormap = "spectral_r"


def main():




	print import_file
	parsefile(import_file)
	print "total runtime",datetime.now() - startTime
	plt.show()







def parsefile(f):


	times = []
	data_tox = []
	data_ant = []

	file = open (f,"r")


	pattern_data= "\d\.*\d*\w*\-*\d*|-nan"

	for line in file:
		if len(line) > 10:
			init = line[0:10]
			if init == "Summary at":
				times.append(int(re.findall("\d+", line)[0]))
			elif init == "Avgmobtox:":
				# #REGEX findall option
				# data_tox_now = [float(i) for i in re.findall(pattern_data,line)]
				# data_tox.append(data_tox_now)
				#SPLIT OPTION
				data_tox_now =  line.split()[1:]
				data_tox_float = [float(i) for i in data_tox_now]
				data_tox.append(data_tox_float)

			elif init == "Avgmobant:":
				# #REGEX findall option
				# data_ant_now = [float(i) for i in re.findall(pattern_data,line)]
				# data_ant.append(data_ant_now)
				# # SPLIT option
				data_ant_now = line.split()[1:]
				data_ant_float = [float(i) for i in data_ant_now]
				data_ant.append(data_ant_float)








	file.close

	datareadtime_now = datetime.now()
	print "time now",datareadtime_now
	datareadtime = datareadtime_now - startTime
	print "datareadtime",datareadtime

	#Plotten
	x = times
	y = data_tox
	z = data_ant
	N = len(data_tox[0])


	cmap = plt.cm.get_cmap(colormap, N+1)
	fig = plt.figure()
	ax0 = fig.add_subplot(211)
	for i in range(N):
		ely= [el[i] for el in y]
		plt.scatter(x, ely, c=cmap(i), marker = "1", s =2 )
	plt.axis([0,times[-1], -0.05, 1.05])
	plt.xlabel('Time')
	plt.ylabel('Mobility')
	plt.title("AVG toxin mobilities")

	ax1 = fig.add_subplot(212, sharex=ax0)
	for i in range(N-1):
		elz = [el[i] for el in z]
		plt.scatter(x, elz, c=cmap(i), marker= "1", s = 2)
	plt.axis([0,times[-1], -0.05, 1.05])
	plt.xlabel('Time')
	plt.ylabel('Mobility')
	plt.title("AVG antitoxin mobilities")

	plttime = datetime.now() - datareadtime_now
	print"plottime",plttime

	## ATTEMPT TO PLOT LINE
	# plt.subplot(313)
	# for i in range(N-1):
	# 			elz = [el[i] for el in z]
	# 			plt.plot(x, elz, c=cmap(i))
	# plt.axis([0,times[-1], -0.05, 1.05])
	# plt.xlabel('Time')
	# plt.ylabel('Mobility')
	# plt.title("AVG antitoxin mobilities")







	print f
#COLOR stuff
def get_cmap(N):
    # '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    # RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap = colormap)
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color



main()
