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


def ancestral_sites(input_AB, input_Gen, output_dir): 
	"""
	Plot du nombre de sites d'AB et Gen par distance aux barrières
	"""
	Gen_count=open(input_Gen, "r")
	AB_count=open(input_AB, "r") 
	Gen=Gen_count.readlines()
	AB=AB_count.readlines()
	
	dist=[]
	AncestralBase=[]
	GenomeBase=[]
	
	for l in AB:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				num_AB=0
				num_AB=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				print(num_AB)
				AncestralBase.append(num_AB)
	
	for a in Gen:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				num_B=0
				num_B=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				print(num_B)
				GenomeBase.append(num_B)			
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,AncestralBase, color='royalblue')
	plt.plot(dist,GenomeBase, color='royalblue', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=13)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axhline(1727712, color='red', alpha=0.5) #Nombre de bords de NIEBs dans interNIEBs >= 1000pb
	plt.title("Number of sites around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Ancestral Bases", fontsize=16)
	filename= "AB_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	

def norm_sites(input_AB, input_Gen, output_dir): 
	"""
	Normalisation du nombre d'AB sur le nombre de bases du génome
	"""
	
	Gen_count=open(input_Gen, "r")
	AB_count=open(input_AB, "r")
	Gen=Gen_count.readlines()
	AB=AB_count.readlines()

	
	dist=[]
	NormBase=[]
	
	for i in range(len(AB)):
			if not AB[i].startswith("d"): #Ignore le header du fichier
				line1=AB[i].strip().split("\t")
				distance1=int(line1[0])
				if distance1 >= -50 and distance1 <= 500 :
					if not Gen[i].startswith("d"): #Ignore le header du fichier
							line2=Gen[i].strip().split("\t")
							distance2=int(line2[0])
							
							if distance1 == distance2:
								dist.append(int(line1[0]))
								num_AB=int(line1[1]) + int(line1[2]) + int(line1[3]) + int(line1[4])
								num_Gen=int(line2[1]) + int(line2[2]) + int(line2[3]) + int(line2[4])
								norm=num_AB/num_Gen
								NormBase.append(norm)
								
	plt.figure(figsize=(10,10))
	plt.plot(dist,NormBase, color='royalblue')
	plt.axes().minorticks_on()
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

	
def ancestral_sub_sites(input_AB, input_Gen, output_dir): 
	"""
	plot du nombre de sites par distance aux barrières par type de base
	"""
	Gen_count=open(input_Gen, "r")
	AB_count=open(input_AB, "r")
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
	plt.plot(dist,AB_List["T"], color='mediumaquamarine')
	plt.plot(dist,Gen_List["A"], color='dodgerblue', linestyle='dashed')
	plt.plot(dist,Gen_List["T"], color='mediumaquamarine', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=15)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of A/T around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=20)
	plt.ylabel("Number of A/T", fontsize=16)
	filename= "A_T.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,AB_List["G"], color='orange')
	plt.plot(dist,AB_List["C"], color='tomato')
	plt.plot(dist,Gen_List["G"], color='orange', linestyle='dashed')
	plt.plot(dist,Gen_List["C"], color='tomato', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=15)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of G/C around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Number of G/C", fontsize=16)
	filename= "G_C.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
				

def ancestral_sub_sites_percent(input_AB, input_Gen, output_dir): 
	"""
	Plot du pourcentage de chaque type de base par distance aux barrières
	"""
	
	Gen_count=open(input_Gen, "r")
	AB_count=open(input_AB, "r")
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
	plt.plot(dist,AB_List["T"], color='mediumaquamarine')
	plt.plot(dist,Gen_List["A"], color='dodgerblue', linestyle='dashed')
	plt.plot(dist,Gen_List["T"], color='mediumaquamarine', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Percentages of A/T around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Percentages of A/T", fontsize=16)
	filename= "A_T_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,AB_List["G"], color='orange')
	plt.plot(dist,AB_List["C"], color='tomato')
	plt.plot(dist,Gen_List["G"], color='orange', linestyle='dashed')
	plt.plot(dist,Gen_List["C"], color='tomato', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Percentages of G/C around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Percentages of G/C", fontsize=16)
	filename= "G_C_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	


def GC_content(input_AB, input_Gen, output_dir):
	"""
	Plot du pourcentage GC
	"""
	Gen_count=open(input_Gen, "r")
	AB_count=open(input_AB, "r")
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
	plt.plot(dist,GC_AB, color='orangered')
	plt.plot(dist,GC_Gen, color='orangered', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axhline(37.9, color='red', alpha=0.5)#%GC moyen du génome
	plt.axvline(0, color='black', alpha=0.5)
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
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	parser.add_argument('-Gen', '--input_Gen', type=str, help='Path to bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
				
	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output directory', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/AB_plots/")			

	
	args = parser.parse_args()
	
	ancestral_sites(args.input_AB, args.input_Gen, args.output_dir)
	norm_sites(args.input_AB, args.input_Gen, args.output_dir)
	ancestral_sub_sites(args.input_AB, args.input_Gen, args.output_dir)
	ancestral_sub_sites_percent(args.input_AB, args.input_Gen, args.output_dir)
	GC_content(args.input_AB, args.input_Gen, args.output_dir)
	
if "__main__" == __name__:
	main()

