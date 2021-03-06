#!/usr/bin/env python
#
#################################################################################
#
# name: evaluate_Ic_mult.py
# date: 10 June 18
# auth: Brandon Sorbom, Zach Hartwig, Erica Salazar
# mail: bsorbom@psfc.mit.edu, hartwig@psfc.mit.edu, erica@psfc.mit.edu
#
# desc: Python Script that fits a standard power law fit of the I-V
#       curves in order to determine the critical current and n-value
#       using a 1 uV/cm threshold.
#
# dpnd: pandas, sys, os, numpy, lmfit, peakutils, matplotlib
#
# 2run: python fit_iv_curve.py <options> (use '-h' to see details)
#
#################################################################################

import pandas as pd
import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit, Model
import peakutils
import csv as csv
from scipy.interpolate import interp1d
import xlsxwriter
import pdb # debugging tool


def shotList():
	date=raw_input("Input date in YYYYMMDD format \n")
	options=raw_input("manual(m) or list(l) \n")
	if options == "m":
		shot_nums = raw_input("Input shot numbers (separated by spaces): \n")
		shot_nums = shot_nums.split()
	else:
		start=raw_input("start shot \n")
		stop=raw_input("stop shot \n")
		shot_nums=[]
		for i in range(int(start),int(stop)+1):
			shot_nums.append(str(i))

	shot_list=[]
	for i in range(0,len(shot_nums)):
		shot_list.append(date+"%03d" % (int(shot_nums[i]),)) # https://stackoverflow.com/questions/134934/display-number-with-leading-zeros
	return shot_list

########
# Main #
########


def main(command_line=True):

	# Set the top level SPARTA data directory path
	# if not 'SPARTA_DATA'in os.environ:
	# 	print "\nError! The SPARTA_DATA environmental variable must be set to use"
	# 	print   "       the 'create_tree.py' script. Please run the 'setup.py'"
	# 	print   "       script included in this repository.\n"
	# 	return -42
	# else:
	# 	data_dir = os.environ['SPARTA_DATA']

	# Set command line options
	if command_line:
		parser = argparse.ArgumentParser()
		parser.add_argument('-s','--shot',
				nargs=1,
				required=True,
				type=int,
				help='PSFC shot number in YYYYMMDDSSS format (Y=year, M=month, D=day, S=shot)')
		parser.add_argument('-l','--length',
				nargs=1,
				required=False,
				type=float,
				help='Override voltage tap length in data file [cm]')
		parser.add_argument('-np','--no-plot',
				required=False,
				action='store_true',
				help='Prevent plotting the results of the I-V critical current fit')
		parser.add_argument('-t','--thresh',
				nargs=1,
				required=False,
				type=float,
				help='Set the minimum current threshold for valid I-V data [A] (Default: 10 A)')
		parser.add_argument('-b','--baseline-range',
				nargs=2,
				required=False,
				type=float,
				help='Set the range for baseline calculation [A] (Default: 10 - 100 A)')

		args = parser.parse_args()

		# Exract command line options

		shot_number = args.shot[0]
		tap_length = args.length[0] if args.length else -1.
		current_thresh = args.thresh[0] if args.thresh else 10.
		baseline_min = args.baseline_range[0] if args.baseline_range else 10
		baseline_max = args.baseline_range[1] if args.baseline_range else 100
		plot_result = False if args.no_plot else True
		use_calibration = None
	else:
		shot_number = int
		current_thresh = float
		baseline_min = float
		baseline_max = float
		reel_id = float
		#shot_number = int(raw_input("Input shot number:\n"))
		#reel_id = str(raw_input("Input reel_id:\n"))
		#shot_list = raw_input("Input shot list (separated by spaces): \n")
		shot_list=shotList()
		current_thresh = 10.
		baseline_min = 10
		baseline_max = 100
		plot_result = True# False
		use_calibration = False

	# Execute I-V curve

	plot_data(shot_list, baseline_min, baseline_max, current_thresh, plot_result)
	#Ic,n = fit_data(shot_number,
	#	baseline_min,
	#	baseline_max,
	#	current_thresh,
	#	plot_result)
	#get_tapestar(reel_id,
	#	Ic,
	#	use_calibration)

#######################
# Import and fit data #
#######################


def fit_data(shot_number, baseline_min, baseline_max, current_thresh, plot_result):

	# Determine full path to the shot data file
	shot_name = str(shot_number)
	#data_dir = os.environ['SPARTA_DATA']
	data_file_path = shot_name + ".csv"
	# Extract the file header information into line strings

	header_lines = 10
	header_info = []

	f = open(data_file_path)
	for i in range(header_lines):
		header_info.append(f.readline())
		header_info[i] = header_info[i][2:-1]
	f.close()

	# Extract relevant header values from the strings. Note that
	# unless the user has specified the tap length on the command
	# line the tap length from the file will be used
	# (i.e. signaled by a value of tap_length<0)

	shot_date = header_info[0].split()[1]
	operator = header_info[1].split()[1:]
	sample_name=''
	for i in range(0,len(header_info[4].split()[1:])):
		sample_name=sample_name+header_info[4].split()[1:][i]+' '
	#sample_name = header_info[4].split()[1:][0]
	#tap_length = float(header_info[3].split()[4])
	tap_length = float(header_info[3].split()[1]) #<<if using new pyplate data
	print(tap_length)

	# Load the data from the data file using pandas

	#data = pd.read_csv(data_file_path, header=5)
	data = pd.read_csv(data_file_path, header=10) #<<if using new pyplate data


	# Remove extraneous data columns

	#data.drop(['DATE'], axis=1, inplace = True)
	#data.drop(['TIME'], axis=1, inplace = True)
	#data.drop(['Status'], axis=1, inplace = True)

	# Preserve only data above a user-specified threshold (default: 10 A)
	data_conditioned = (data[data["Shunt [A]"] > current_thresh])#<<if using new pyplate data
	#data_conditioned = (data[data["Current [A]"] > current_thresh])

	# Truncate data to calculate a user-specified baseline (default 10 A - 100 A)
	baseline_conditioned = (data_conditioned[data_conditioned["Shunt [A]"] > baseline_min])
	baseline_conditioned = (baseline_conditioned[baseline_conditioned["Shunt [A]"] < baseline_max])

	# If there is remaining data, process it
	if not data_conditioned.empty:

		# Define functional fit for V/I data
		voltages1 = data_conditioned["Tap [uV]"]
		currents1 = data_conditioned["Shunt [A]"]
		# removes last data point in data (because the last data point is incorrect due to pyplate error)
		voltages = voltages1[0:-1]
		currents = currents1[0:-1]

		# Calculate the voltage floor as the average of baseline array [uV]
		V_floor = np.mean(baseline_conditioned["Tap [uV]"])

		# Set the critical voltage [uV]
		Vc = tap_length

		# Define the fit functions
		def fit_fc(params, x, data):

			# Specify the free parameters
			Ic = params['Ic'].value
			n = params['n'].value

			# Specify the fit formula
			model = Vc*(x/Ic)**n + V_floor

			# Return the quantity that is being minimized in the fit
			return (model - data)

		# Create initial guesses for the free parameters

		params = Parameters()
		params.add('Ic', value=125, min=1, max=1000)
		#params.add('n', value=15, min=0, max=100)
		params.add('n', value=20, min=0, max=100)

		# Perform the fit and report on the results
		result = minimize(fit_fc, params, args=(currents, voltages))
		#result = minimize(fit_fc, params, args=(data_conditioned["Current [A]"], data_conditioned["QD Voltage [uV]"]))
		report_fit(result.params)

		# Get the fit results

		Ic = result.params['Ic']
		n = result.params['n']
		Ic_error = result.params['Ic'].stderr
		n_error = result.params['n'].stderr
		########################################## adding for plotting
		# Synthesize an I-V curve from the fit results

		fit_xmin = np.min(currents)
		fit_xmax = np.max(currents) * 1.05

		Ic_fit_array = np.linspace(fit_xmin, fit_xmax, num='100', endpoint='true')
		V_fit_array = Vc * (Ic_fit_array / Ic) ** n + V_floor
		Vc_array = np.full(100, (Vc + V_floor))

		if plot_result:
			print"##################plotting!!!############"
			font = {'family': 'sans-serif',
					'weight': 'normal',
					'size': 15}
			plt.rc('font', **font)

			plot_y_min = np.min(voltages)
			plot_y_max = Vc * 2 + V_floor

			fig, ax = plt.subplots()

			plt.plot(currents, voltages, 'o', markeredgecolor='grey', markeredgewidth=1, markersize=8, mfc='none')
			plt.plot(Ic_fit_array, V_fit_array, '-', color='blue', linewidth=3)
			plt.plot(Ic_fit_array, Vc_array, '--', color='blue')
			plt.plot(np.array([Ic, Ic]), np.array([plot_y_min, plot_y_max]), '--', color='blue')

			ax.set_xlabel('Current (A)')
			ax.set_ylabel('Voltage (uV)')
			ax.set_title('Shot %s Ic fit' % shot_name)

			plt.ylim(-3, 20)

			# plt.ylim(plot_y_min, plot_y_max)

			# Add Tc value +/- error box to plot
			props = dict(boxstyle='round', facecolor='white', alpha=0.5)
			ax.text(0.10, 0.90, 'Sample: %s \n $I_c$ = %.2f +/- %.2f A \n n = %.2f +/- %.2f' % (
				sample_name, Ic, Ic_error, n, n_error), transform=ax.transAxes,
					verticalalignment='top', horizontalalignment='left', bbox=props)

			#plt.show()
			plt.savefig(shot_name+"_re-fit.png")


	#************************ adding for plotting


	else:
		print "Empty dataframe!"
	return [Ic ,n ,Ic_error ,n_error , currents, voltages, V_floor, Vc, sample_name]


########
		# end of fit_data

#def plot_data(currents, voltages, V_floor, n, Vc, Ic, plot_result, shot_name):
def plot_data(shot_list, baseline_min, baseline_max, current_thresh, plot_result):

	#shot_list = shot_list.split()
	shot_list = map(int, shot_list)
	n_shots = int(len(shot_list))
	Ic_array = [0]*n_shots
	Ic_error_array = [0]*n_shots
	Vc_array = [0]*n_shots
	V_floor_array = [0]*n_shots
	n_array = [0]*n_shots
	n_error_array = [0]*n_shots
	currents_array = [0]*n_shots
	voltages_array = [0]*n_shots
	sample_name_array = [0]*n_shots
	Iminmax_array = np.zeros([n_shots,2])
	Vminmax_array = np.zeros([n_shots, 2])
	i = 0
	for shot_number in shot_list:
		# params_array = [Ic_temp, n_temp, Ic_error_temp, n_error_temp, currents_temp, voltages_temp, V_floor_temp, Vc_temp]
		params_array = fit_data(shot_number, baseline_min, baseline_max, current_thresh, plot_result)
		Ic_array[i] = params_array[0]/1
		n_array[i] = params_array[1]/1
		Ic_error_array[i] = params_array[2]
		n_error_array[i] = params_array[3]
		currents_array[i] = params_array[4]
		voltages_array[i] = params_array[5]
		V_floor_array[i] = params_array[6]
		Vc_array[i] = params_array[7]
		sample_name_array[i] = params_array[8]
		Iminmax_array[i,0] = np.min(currents_array[i])
		Iminmax_array[i,1] = np.max(currents_array[i])
		Vminmax_array[i,:] = [np.min(voltages_array[i]), np.max(voltages_array[i])]
		i += 1

	# Synthesize an I-V curve from the fit results
	Iminmax_array = np.asarray(Iminmax_array)
	fit_xmin = np.min(Iminmax_array[:,0])
	fit_xmax = np.max(Iminmax_array[:,1]) * 1.05
	Ic_fit_array = np.linspace(fit_xmin,fit_xmax,num='100',endpoint='true')
	V_fit_array = np.zeros([n_shots,int(len(Ic_fit_array))]) #Vc_array*(Ic_fit_array/Ic_array)**n_array + V_floor_array
	#Vc_array2 = np.full(100,(Vc_array+V_floor_array))
	for i in range(0,n_shots):
		V_fit_array[i] = Vc_array[i]*(Ic_fit_array/Ic_array[i])**n_array[i] + V_floor_array[i]

	# Write to Excel
	workbk = xlsxwriter.Workbook('tape_test_xlsx.xlsx')
	worksht = workbk.add_worksheet()
	bold = workbk.add_format({'bold': True})
	worksht.write('A1','HTS tape results', bold)
	worksht.write('A3', 'Sample Name')
	worksht.write('B3', 'Shot number')
	worksht.write('C3', 'Ic [A]')
	worksht.write('D3', 'Ic error [+/- A]')
	worksht.write('E3', 'n')
	worksht.write('F3', 'n error [+/-]')


	row = 3
	col = 0
	for j in range(0,n_shots):
		worksht.write(row, col, sample_name_array[j])
		worksht.write(row, col+1, shot_list[j])
		worksht.write(row, col+2, Ic_array[j])
		worksht.write(row, col+3, Ic_error_array[j])
		worksht.write(row, col+4, n_array[j])
		worksht.write(row, col+5, n_error_array[j])
		print "j = /n %d" % j
		print "shot %d" % shot_list[j]
		row += 1

	workbk.close()

	#pdb.set_trace()
	# Plot data for comparison
	#font = {'family': 'sans-serif','weight': 'normal', 'size': 15}
	#plt.rc('font', **font)
	plot_y_min = np.min(Vminmax_array[:, 0]) - 3.
	plot_y_max = np.max(Vminmax_array[:, 1]) + 2
	plt.figure(n_shots+1)
	plt.plot(Ic_array,Vc_array,'o')
	plt.xlabel('Current [A]')
	plt.ylabel('Voltage [uV]')
#	plt.show()




	## Output the results to the command line

	#print "\n"
	#print "************************   I-V fit results   ************************\n"
	#print "Shot date   : %s" % shot_date
	#print "Shot number : %d" % shot_number
	#rint "Operator    : %s" % operator
	#print "Sample      : %s" % sample_name
	#print "Tap length  : %.2f cm" % tap_length
	#print "I thresh.   : %.2f A" % current_thresh
	#print ""
	#print "Ic = %2.2f +/- %2.2f A" % (Ic,Ic_error)
	#print " n = %2.2f +/- %2.2f" % (n,n_error)
	#print "Vc = %2.2f uV" % Vc
	#print "Vf = %2.2f uV" % V_floor
	#print "Vt = %2.2f uV" % (Vc+V_floor)
	#print ""
	#print "*********************************************************************\n"






######################
# Read TAPESTAR data #
######################

def get_tapestar(reel_id, Ic_meas, use_calibration):

	# Set the file name for the user-selected TAPESTAR reel

	if(reel_id == '104'):
		file_name = "SuperPower_K19039_June2017/M4-426-5-0104-1239to1349-Ic4M-sys2-25mT-Ic-x.csv"
	elif(reel_id == '508'):
		file_name = "SuperPower_K19039_June2017/M4-426-5-0508-1239to1349-Ic4M-sys3-14mT-Ic-x.csv"
	elif (reel_id == '2c'):
		file_name = "SuperPower_K19191_April2018/2c_M3-1318-4 0104 1769-1829m.csv"
	else:
		print "\nError! Unrecognized HTS reel ID!\n"
		return

	data_dir = os.environ['SPARTA_DATA']
	data_file_path = data_dir + "/tapestar/" + file_name

	# Read in the TAPESTAR data

	raw_data = []
	reader = csv.reader(open(data_file_path,'rb'),delimiter=',')
	for index,row in enumerate(reader):
		raw_data.append([float(row[0]),float(row[1])])
	tapestar_data = zip(*np.array(raw_data))

	# Position along the reel [cm]
	x = tapestar_data[0]

	# Critical current [A]
	y = tapestar_data[1]

	# Create a function to be used for interpolation
	f = interp1d(x, y)

	# Set the start and end position of the tape. Note that our
	# present SuperPower reel is numerically reversed
	# (i.e. "start" position is higher than "end" position)

	tape_start = input('Input starting position (in cm): ')
	tape_length = input('Input sample length (in cm): ')
	tape_end = tape_start - tape_length

	# Set array size using TapeStar resolution of 1 mm
	resolution_factor = 10
	Ic_array_size = tape_length * resolution_factor

	# Create arrays to hold position and Ic
	position_array = []
	Ic_array = []

	# Fill an array with 1mm steps along the tape
	for i in range(0, int(Ic_array_size)):
		position_array.append(tape_end + i * 1./resolution_factor)

	# Use the 1D interp function to populate the critical currents
	print "position array = " + str(position_array)
	tapestar_calibration = np.array([1.0000, 1.0895, 1.1022])
	tapestar_index = 0
	#pdb.set_trace()
	if use_calibration:
		if reel_id == 104: tapestar_index = 1
		elif reel_id == 508: tapestar_index = 2

	for position in position_array:
		#try:
		Ic_array.append(f(position) / tapestar_calibration[tapestar_index])
		print "*****************************************"
		print "NORMAL RUN OF FOR LOOP"
		print "position=" + str(position)
		print "f(position)=" + str(f(position))
		#except ValueError:
		#	print "*****************************************"
		#	print "ERROR"
		#	print "position=" + str(position)
		#	print "x=" + str(x)
		#	print "y=" + str(y)

	Ic_min = np.min(Ic_array)
	Ic_max = np.max(Ic_array)
	Ic_avg = np.mean(Ic_array)
	Ic_std = np.std(Ic_array)

	calibration_string = "No"
	if use_calibration: calibration_string = "Yes"

	print ""
	print "************************   TAPESTAR results   ************************\n"
	print "TAPESTAR Calibrated: %s" % (calibration_string)
	print "Calibration factor : %0.4f" % (tapestar_calibration[tapestar_index])
	print ""
	print "Results for tape range %f to %f:" % (tape_end, tape_start)
	print "Minimum Ic: %.1f A" % (Ic_min)
	print "Maximum Ic: %.1f A" % (Ic_max)
	print "Average Ic: %.1f +/- %.1f A (%.1f%%)" % (Ic_avg, Ic_std, (Ic_std/Ic_avg*100))
	print ""
	print "**********************************************************************\n"

	# Create x-axis as relative position along tape
	position_plot = []
	for i in position_array:
		position_plot.append(i - tape_end)

	x_min = 0
	x_max = len(position_array) / resolution_factor

	y_min = Ic_min*0.95
	y_max = Ic_max*1.05

	font = {'family' : 'sans-serif',
		'weight' : 'normal',
		'size'   : 15}
	plt.rc('font', **font)

	f2 = plt.figure(2)
	plt.plot(position_plot, Ic_array, 'o', markeredgecolor='grey', markeredgewidth=1, markersize=8, mfc='none')
	plt.plot(np.array([x_min, x_max]), np.array([Ic_avg, Ic_avg]), '-', color='red', linewidth=2)
	plt.plot(np.array([x_min, x_max]), np.array([Ic_meas, Ic_meas]), '-', color='blue', linewidth=2)

	plt.axhspan(Ic_avg-Ic_std, Ic_avg+Ic_std, alpha=0.25, color='red')

	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)

	title_string = "Start pos. (%s): %0.1f cm" % (file_name[0:13], tape_end)

	plt.title(title_string)
	plt.xlabel("Position on tape (cm) ")
	plt.ylabel("Critical current [A]")

	plt.legend(['TAPESTAR', 'Ic_average', 'Ic_measured'], loc='best')

	f2.show()

	raw_input()


###############################
# Create Python main function #
###############################

if __name__ == "__main__":
	main(command_line=False)
