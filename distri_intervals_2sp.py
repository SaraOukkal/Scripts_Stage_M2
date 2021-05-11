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
import seaborn as sns

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
	tot=len(inter)
	sizes={}
	dens=[]
	length=[]
	liste=[]
	
	for l in inter:
		line=l.strip().split("\t")
		st=int(line[2])
		end=int(line[3])	
		size=end-st
		if size < 2000:
			if size not in sizes.keys():
				sizes[size]=0
			sizes[size]+=1
		
	for k in sorted(sizes.keys()):
		density=sizes[k]/tot
		dens.append(density)
		length.append(sizes[k])
	
	liste.append(length)
	liste.append(dens)
	
	return liste
	
def plot(distribution1, distribution2, filename, output_dir):
	
	plt.figure(figsize=(10,10))	
	plt.plot(distribution1[0],distribution1[1], color='firebrick')
	plt.plot(distribution2[0],distribution2[1], color='tomato')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(117, color='black', alpha=0.5)
	plt.axvline(270, color='black', alpha=0.5)
	plt.axvline(423, color='black', alpha=0.5)
	plt.axvline(576, color='black', alpha=0.5)
	plt.axvline(729, color='black', alpha=0.5)
	plt.axvline(882, color='black', alpha=0.5)
	plt.title("Intervals size distribution", fontsize=20)
	plt.xlabel("Size", fontsize=16)
	plt.ylabel("Density", fontsize=16)
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)	
	plt.clf()	
	
def plot_inter(distribution1, distribution2, filename, output_dir):
	
	plt.figure(figsize=(20,10))	
	sns.kdeplot(distribution1, color="black")
	sns.kdeplot(distribution2, color="grey")
	plt.xlim(0,1500)
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.title("Intervals size distribution", fontsize=20)
	plt.xlabel("Size", fontsize=16)
	plt.ylabel("density", fontsize=16)
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
	All1=0
	size1=[]
	less2=0
	more2=0
	All2=0
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
		All1+=1
	
	less_per1=(less1/All1)*100	
	more_per1=(more1/All1)*100		
	size1.append(less_per1)
	size1.append(more_per1)
	
	for l in inter2:
		line=l.strip().split("\t")
		st=int(line[2])
		end=int(line[3])	
		sz=end-st
		if sz < 1000:
			less2+=1
		else:
			more2+=1
		All2+=1
		
	less_per2=(less1/All2)*100	
	more_per2=(more1/All2)*100		
	size2.append(less_per2)
	size2.append(more_per2)
	
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
	plot_inter(distribution3, distribution4,"Distri_inter_NIEBs.png", args.output)
	
	inter_bar_1000(args.input_inter1, args.input_inter2, args.output)

if "__main__" == __name__:
	main()

