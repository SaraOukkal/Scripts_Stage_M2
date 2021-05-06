#!/usr/bin/python3

#Importer packages: 
import argparse 
import numpy as np
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def number(input_1,input_2, out):
	count1=open(input_1, "r") 
	sp1=count1.readlines()
	count2=open(input_2, "r") 
	sp2=count2.readlines()
	
	dist=[]
	list1=[]
	list2=[]
	
	for l in sp1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				num=0
				num=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				print(num)
				list1.append(num)
	
	for l in sp2:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				num=0
				num=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				print(num)
				list2.append(num)
	
	plt.figure(figsize=(10,10))	
	plt.plot(dist,list1, color='darkviolet')
	plt.plot(dist,list2, color='violet')
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.title("Number of sites around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Number of sites", fontsize=16)
	plt.savefig(out)
	plt.clf()
	
	
		
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-sp1', '--input_1', type=str, help='Path to mutations count', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")
	parser.add_argument('-sp2', '--input_2', type=str, help='Path to mutations count', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")				

	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/Mut_plots/")			

	
	args = parser.parse_args()
	number(args.input_1,args.input_2, args.output_dir)

if "__main__" == __name__:
	main()
