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


def distri_bar(intervals):
	f=open(intervals,"r")
	inter=f.readlines()
	sizes=[]
	
	for l in inter:
		line=l.strip().split("\t")
		st=int(line[1])
		end=int(line[2])	
		size=end-st
		sizes.append(size)
		
	return sizes
	
def distri_inter_bar(intervals):
	f=open(intervals,"r")
	inter=f.readlines()
	sizes=[]
	
	for l in inter:
		line=l.strip().split("\t")
		st=int(line[2])
		end=int(line[3])	
		size=end-st
		sizes.append(size)
		
	return sizes

def plot (distribution1, distribution2, filename, output_dir):
	
	plt.figure(figsize=(10,10))	
	plt.hist(distribution1, bins=20, color="darkviolet", alpha=0.5)
	plt.hist(distribution2, bins=20, color="violet", alpha=0.5)
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.title("Intervals size distribution", fontsize=20)
	plt.xlabel("Size", fontsize=16)
	plt.ylabel("Number", fontsize=16)
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)	
	plt.clf()	
	
				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier .bed des barrières chez le chimpanzé: 
	parser.add_argument('-bar1', '--input_bar1', type=str, help='Path to bar  intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/H_intervals.bed")
	parser.add_argument('-bar2', '--input_bar2', type=str, help='Path to bar intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/H_intervals.bed")
	
	parser.add_argument('-inter1', '--input_inter1', type=str, help='Path to interbar intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/H_intervals.bed")
	parser.add_argument('-inter2', '--input_inter2', type=str, help='Path to interbar intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/H_intervals.bed")
	
	parser.add_argument('-out', '--output', type=str, help='Path to output directory', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step4_results/Plots/")
	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	distribution1= distri_bar(args.input_bar1)
	distribution2= distri_bar(args.input_bar2)
	plot(distribution1, distribution2, "Distri_NIEBs.png", args.output)
	
	distribution3= distri_inter_bar(args.input_inter1)
	distribution4= distri_inter_bar(args.input_inter2)
	plot(distribution3, distribution4,"Distri_inter_NIEBs.png", args.output)

if "__main__" == __name__:
	main()

