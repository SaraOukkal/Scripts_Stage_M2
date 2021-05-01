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
		if size < 10000:
			sizes.append(size)
		
	return sizes
	
def plot(distribution1, distribution2, filename, output_dir):
	
	plt.figure(figsize=(10,10))	
	plt.hist(distribution1, bins=20, color="limegreen", alpha=0.3, ec="darkgreen")
	plt.hist(distribution2, bins=20, color="dodgerblue", alpha=0.3, ec="darkblue")
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
	
def inter_bar_1000(intervals1, intervals2, output_dir):	
	f1=open(intervals1,"r")
	inter1=f1.readlines()
	f2=open(intervals2,"r")
	inter2=f2.readlines()
	
	types=["less", "more"]
	less1=0
	more1=0
	size1=[]
	less2=0
	more2=0
	size2=[]
	
	for l in inter1:
		line=l.strip().split("\t")
		st=int(line[2])
		end=int(line[3])	
		sz=end-st
		if sz < 1000:
			less1+=1
		else:
			more1+=1
			
	size1.append(less1)
	size1.append(more1)
	
	for l in inter2:
		line=l.strip().split("\t")
		st=int(line[2])
		end=int(line[3])	
		sz=end-st
		if sz < 1000:
			less2+=1
		else:
			more2+=1
			
	size2.append(less2)
	size2.append(more2)
	
	plt.figure(figsize=(10,10))	
	plt.bar(types,size1, color="limegreen", alpha=0.3, ec="darkgreen")
	plt.bar(types,size2, color="dodgerblue", alpha=0.3, ec="darkblue")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.title("Size of InterNIEBs", fontsize=20)
	plt.ylabel("Number", fontsize=16)
	filename="More_less.png"
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
	
	inter_bar_1000(args.input_inter1, args.input_inter2, args.output)

if "__main__" == __name__:
	main()

