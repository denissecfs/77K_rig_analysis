#!/usr/bin/env python
#
#################################################################################
#
# name: calculate_IC.py
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
import lmfit
from lmfit import minimize, Parameters, Parameter, report_fit, Model
import peakutils
import csv as csv
from scipy.interpolate import interp1d
import xlsxwriter
import pdb # debugging tool


########
# Main #
########


def main(command_line=True):
	"""
	Depending on whether this script is run the command line or not, a set of instructions will execute
	Parameters: Boolean command line
	Returns: Depend on whether executed in terminal or not. shot_list, current_thresh, baseline_min, baseline_max, plot_result
	Parameters that function plot data will use. 
	"""
	
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
		
		return shot_number, tap_length, current_thresh, baseline_min, baseline_max, plot_result 
	
	else:
		shot_number = int
		current_thresh = float
		baseline_min = float
		baseline_max = float
		#Unless tapestar data is being analyzed the next four lines should be commented out. 
		#reel_id = float 
		#shot_number = int(input("Input shot number:\n"))
		#reel_id = str(input("Input reel_id:\n"))
		#shot_list = input("Input shot list (separated by spaces): \n")
		shot_list=shotList()
		current_thresh = 10.
		baseline_min = 10
		baseline_max = 100
		plot_result = True # False
		use_calibration = False
		
		return shot_list, current_thresh, baseline_min, baseline_max, plot_result

def shotList():
	date=input("Input date in YYYYMMDD format \n")
	options=input("manual(m) or list(l) \n")
	
	if options == "m":
		shot_nums = input("Input shot numbers (separated by spaces): \n")
		shot_nums = shot_nums.split()
	else:
		start=input("start shot \n")
		stop=input("stop shot \n")
		shot_nums=[]
		for i in range(int(start),int(stop)+1):
			shot_nums.append(str(i))
	
	shot_list=[]
	for i in range(0,len(shot_nums)):
		shot_list.append(date+"%03d" % (int(shot_nums[i]),)) # https://stackoverflow.com/questions/134934/display-number-with-leading-zeros
	
	return shot_list


#######################
# Import and fit data #
#######################


def fit_data(shot_number, baseline_min, baseline_max, current_thresh, plot_result):
	# Determine full path to the shot data file
	shot_name = str(shot_number)
	data_dir = os.getcwd()
	data_folder = input("Name of Data Folder")
	data_date_batch = input("Date of Data Batch")
	data_file_path = data_dir + "\\" + data_folder + "\\"+ data_date_batch + "\\" + shot_name + ".csv"
	print(data_file_path)
	# Extract the file header information into line strings
	header_lines = 5
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
	# Load the data from the data file using pandas

	data = pd.read_csv(data_file_path, header=5)
	#data = pd.read_csv(data_file_path, header=10) #<<if using new pyplate data
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
		
		
		Ic_fit_array = np.linspace(fit_xmin, fit_xmax, num=100, endpoint='true')
		V_fit_array = Vc * (Ic_fit_array / Ic) ** n + V_floor
		Vc_array = np.full(100, (Vc + V_floor))

		if plot_result:
			print("##################plotting!!!############")
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
		print("Empty dataframe!")
		return None
	return [Ic ,n ,Ic_error ,n_error , currents, voltages, V_floor, Vc, sample_name]


########
		# end of fit_data

#def plot_data(currents, voltages, V_floor, n, Vc, Ic, plot_result, shot_name):
def plot_data(shot_list, baseline_min, baseline_max, current_thresh, plot_result):
	shot_list = list(map(int, shot_list))
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
	Ic_fit_array = np.linspace(fit_xmin,fit_xmax,num=100,endpoint='true')
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
		print("j = /n %d" % j)
		print("shot %d" % shot_list[j])
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




#	## Output the results to the command line

#	#print( "\n")
#	#print( "************************   I-V fit results   ************************\n")
#	#print("Shot date   : %s" % shot_date)
#	#print("Shot number : %d" % shot_number)
#	#print( "Operator    : %s" % operator)
#	#print("Sample      : %s" % sample_name)
#	#print( "Tap length  : %.2f cm" % tap_length)
#	#print( "I thresh.   : %.2f A" % current_thresh)
#	#print( "")
#	#print( "Ic = %2.2f +/- %2.2f A" % (Ic,Ic_error))
#	#print( " n = %2.2f +/- %2.2f" % (n,n_error))
#	#print( "Vc = %2.2f uV" % Vc)
#	#print( "Vf = %2.2f uV" % V_floor)
#	#print( "Vt = %2.2f uV" % (Vc+V_floor))
#	#print( "")
#	#print( "*********************************************************************\n")



###############################
# Create Python main function #
###############################

if __name__ == "__main__":
	shot_list, current_thresh, baseline_min, baseline_max, plot_result = main(command_line=False) #these variables are only for FALSE
    
    

#Execute I-V curve
plot_data(shot_list, baseline_min, baseline_max, current_thresh, plot_result)
#Ic,n = fit_data(shot_number,baseline_min,baseline_max, current_thresh, plot_result)




