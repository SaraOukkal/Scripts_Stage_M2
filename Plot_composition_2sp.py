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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
	
def ancestral_sub_sites(input_AB1, input_AB2, output_dir): 
	"""
	plot du nombre de sites par distance aux barrières par type de base
	"""
	Gen_count=open(input_AB2, "r")
	AB_count=open(input_AB1, "r")
	Gen=Gen_count.readlines()
	AB=AB_count.readlines()

	types=["A","C","G","T"]
	AB_List={}
	Gen_List={}
	
	dist=[]
	
	for sub in types: 
		AB_List[sub]=[]
		Gen_List[sub]=[]
		
	for l in AB:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				for t in AB_List: 
					AB_List[t].append(int(line[num]))
					num+=1
		
	for a in Gen:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				for t in AB_List: 
					Gen_List[t].append(int(line[num]))
					num+=1
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,AB_List["A"], color='dodgerblue')
	plt.plot(dist,AB_List["T"], color='limegreen')
	plt.plot(dist,Gen_List["A"], color='dodgerblue', linestyle='dashed')
	plt.plot(dist,Gen_List["T"], color='limegreen', linestyle='dashed')
	plt.plot(dist,AB_List["G"], color='orange')
	plt.plot(dist,AB_List["C"], color='tomato')
	plt.plot(dist,Gen_List["G"], color='orange', linestyle='dashed')
	plt.plot(dist,Gen_List["C"], color='tomato', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=15)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of A/T/G/C sites at NIEB loci", fontsize=20)
	plt.xlabel("Distance from NIEB borders", fontsize=20)
	plt.ylabel("Number of sites", fontsize=16)
	filename= "ATGC.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	

				

def ancestral_sub_sites_percent(input_AB1, input_AB2, output_dir): 
	"""
	Plot du pourcentage de chaque type de base par distance aux barrières
	"""
	
	Gen_count=open(input_AB2, "r")
	AB_count=open(input_AB1, "r")
	Gen=Gen_count.readlines()
	AB=AB_count.readlines()

	types=["A","C","G","T"]
	AB_List={}
	Gen_List={}
	
	dist=[]
	
	for sub in types: 
		AB_List[sub]=[]
		Gen_List[sub]=[]
		
	for l in AB:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4]) 
				dist.append(int(line[0]))
				for t in AB_List: 
					AB_List[t].append((int(line[num])/tot)*100)
					num+=1
		
	for a in Gen:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4]) 
				for t in AB_List: 
					Gen_List[t].append((int(line[num])/tot)*100)
					num+=1
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,AB_List["A"], color='dodgerblue')
	plt.plot(dist,AB_List["T"], color='limegreen')
	plt.plot(dist,Gen_List["A"], color='dodgerblue', linestyle='dashed')
	plt.plot(dist,Gen_List["T"], color='limegreen', linestyle='dashed')
	plt.plot(dist,AB_List["G"], color='orange')
	plt.plot(dist,AB_List["C"], color='tomato')
	plt.plot(dist,Gen_List["G"], color='orange', linestyle='dashed')
	plt.plot(dist,Gen_List["C"], color='tomato', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Percentages of A/T/G/C around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEB borders", fontsize=16)
	plt.ylabel("Percentages of A/TG/C sites", fontsize=16)
	filename= "ATGC_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	

def GC_content(input_AB1, input_AB2, output_dir):
	"""
	Plot du pourcentage GC
	"""
	Gen_count=open(input_AB2, "r")
	AB_count=open(input_AB1, "r")
	Gen=Gen_count.readlines()
	AB=AB_count.readlines()
	
	dist=[]
	GC_AB=[]
	GC_Gen=[]

	for l in AB:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				C=int(line[2])
				G=int(line[3]) 
				GC=((C + G)/tot)*100
				dist.append(int(line[0]))
				GC_AB.append(GC)
	
	for a in Gen:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				C=int(line[2])
				G=int(line[3]) 
				GC=((C + G)/tot)*100
				GC_Gen.append(GC)
	
	plt.figure(figsize=(10,10))			
	plt.plot(dist,GC_AB, color='blue')
	plt.plot(dist,GC_Gen, color='red', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axhline(40.40, color='red', alpha=0.5)#%GC moyen du génome
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("GC % around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("GC percentage", fontsize=16)
	filename= "GC_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)	
	plt.clf()
		
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichier des bases: 
	parser.add_argument('-AB1', '--input_AB1', type=str, help='Path to ancestral bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	parser.add_argument('-AB2', '--input_AB2', type=str, help='Path to bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
				
	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output directory', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/AB_plots/")			

	
	args = parser.parse_args()
	
	ancestral_sub_sites(args.input_AB1, args.input_AB2, args.output_dir)
	ancestral_sub_sites_percent(args.input_AB1, args.input_AB2, args.output_dir)
	GC_content(args.input_AB1, args.input_AB2, args.output_dir)
	
if "__main__" == __name__:
	main()

