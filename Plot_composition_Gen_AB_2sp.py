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

def norm_sites(input_AB1, input_Gen1,input_AB2, input_Gen2, output_dir): 
	"""
	Normalisation du nombre d'AB sur le nombre de bases du génome
	"""
	Gen_count1=open(input_Gen1, "r")
	AB_count1=open(input_AB1, "r") 
	Gen1=Gen_count1.readlines()
	AB1=AB_count1.readlines()
	
	Gen_count2=open(input_Gen1, "r")
	AB_count2=open(input_AB1, "r") 
	Gen2=Gen_count2.readlines()
	AB2=AB_count2.readlines()

	
	dist=[]
	NormBase1=[]
	NormBase2=[]
	
	for i in range(len(AB1)):
			if not AB1[i].startswith("d"): #Ignore le header du fichier
				line1=AB1[i].strip().split("\t")
				distance1=int(line1[0])
				if distance1 >= -50 and distance1 <= 500 :
					if not Gen1[i].startswith("d"): #Ignore le header du fichier
							line2=Gen1[i].strip().split("\t")
							distance2=int(line2[0])
							
							if distance1 == distance2:
								dist.append(int(line1[0]))
								num_AB=int(line1[1]) + int(line1[2]) + int(line1[3]) + int(line1[4])
								num_Gen=int(line2[1]) + int(line2[2]) + int(line2[3]) + int(line2[4])
								norm=num_AB/num_Gen
								NormBase1.append(norm)
								
	for i in range(len(AB2)):
			if not AB2[i].startswith("d"): #Ignore le header du fichier
				line1=AB2[i].strip().split("\t")
				distance1=int(line1[0])
				if distance1 >= -50 and distance1 <= 500 :
					if not Gen2[i].startswith("d"): #Ignore le header du fichier
							line2=Gen2[i].strip().split("\t")
							distance2=int(line2[0])
							
							if distance1 == distance2:
								num_AB=int(line1[1]) + int(line1[2]) + int(line1[3]) + int(line1[4])
								num_Gen=int(line2[1]) + int(line2[2]) + int(line2[3]) + int(line2[4])
								norm=num_AB/num_Gen
								NormBase2.append(norm)
								
	plt.figure(figsize=(10,10))
	plt.plot(dist,NormBase1, color='darkblue')
	plt.plot(dist,NormBase2, color='blue')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axhline(0.7952966, color='red', alpha=0.5) 
	plt.title("Proportion of Bases with known ancestral state around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Ancestral Bases / All bases", fontsize=16)
	filename= "AB_norm_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	

def GC_content(input_AB1, input_Gen1,input_AB2, input_Gen2, output_dir):
	"""
	Plot du pourcentage GC
	"""
	Gen_count1=open(input_Gen1, "r")
	AB_count1=open(input_AB1, "r") 
	Gen1=Gen_count1.readlines()
	AB1=AB_count1.readlines()
	
	Gen_count2=open(input_Gen1, "r")
	AB_count2=open(input_AB1, "r") 
	Gen2=Gen_count2.readlines()
	AB2=AB_count2.readlines()
	
	dist=[]
	GC_AB1=[]
	GC_Gen1=[]
	GC_AB2=[]
	GC_Gen2=[]
	
	
	for l in AB1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				C=int(line[2])
				G=int(line[3]) 
				GC=((C + G)/tot)*100
				dist.append(int(line[0]))
				GC_AB1.append(GC)
	
	for a in Gen1:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				C=int(line[2])
				G=int(line[3]) 
				GC=((C + G)/tot)*100
				GC_Gen1.append(GC)
	
	for l in AB2:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				C=int(line[2])
				G=int(line[3]) 
				GC=((C + G)/tot)*100
				GC_AB2.append(GC)
	
	for a in Gen2:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				C=int(line[2])
				G=int(line[3]) 
				GC=((C + G)/tot)*100
				GC_Gen2.append(GC)
	
	plt.figure(figsize=(10,10))			
	plt.plot(dist,GC_AB1, color='orangered')
	plt.plot(dist,GC_Gen1, color='orangered', linestyle='dashed')
	plt.plot(dist,GC_AB2, color='lightcoral')
	plt.plot(dist,GC_Gen2, color='lightcoral', linestyle='dashed')
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
	parser.add_argument('-Gen1', '--input_Gen1', type=str, help='Path to bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	parser.add_argument('-AB2', '--input_AB2', type=str, help='Path to ancestral bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	parser.add_argument('-Gen2', '--input_Gen2', type=str, help='Path to bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output directory', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/AB_plots/")			

	
	args = parser.parse_args()

	norm_sites(args.input_AB1, args.input_Gen1,args.input_AB2, args.input_Gen2, args.output_dir)
	GC_content(args.input_AB1, args.input_Gen1,args.input_AB2, args.input_Gen2, args.output_dir)
	
if "__main__" == __name__:
	main()

