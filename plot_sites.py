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

def ancestral_sites(bar_AB_count, output_dir): #Plot du nombre de sites par distance aux barrières
	AB_count=open(bar_AB_count, "r") #Fichier avec nombe de bases par type en fonction de la distance aux barrières
	AB=AB_count.readlines()
	print("Loaded files AB")
	
	dist=[]
	AncestralBase=[]
	
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

	plt.figure(figsize=(10,10))
	plt.plot(dist,AncestralBase, color='royalblue')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=13)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of sites in poly nt around NIEBs")
	plt.xlabel("Distance from NIEBs")
	plt.ylabel("Poly bases")
	filename= "Poly_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
def ancestral_sub_sites(input_AB, output_dir): #plot du nombre de sites par distace aux barrières par type de base
	AB_count=open(input_AB, "r")
	AB=AB_count.readlines()

	types=["A","C","G","T"]
	AB_List={}
	
	dist=[]
	
	for sub in types: 
		AB_List[sub]=[]
		
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

	plt.figure(figsize=(10,10))
	plt.plot(dist,AB_List["A"], color='dodgerblue', linestyle='dashed')
	plt.plot(dist,AB_List["T"], color='mediumaquamarine', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=15)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of A/T in poly A/T around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=20)
	plt.ylabel("Number of A/T", fontsize=16)
	filename= "A_T_poly.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,AB_List["G"], color='orange', linestyle='dashed')
	plt.plot(dist,AB_List["C"], color='tomato', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=15)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of G/C in poly G/C around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Number of G/C", fontsize=16)
	filename= "G_C_poly.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
				

def ancestral_sub_sites_percent(bar_AB_count, output_dir): #plot du pourcentage de bases par distace aux barrières
	AB_count=open(bar_AB_count, "r")
	AB=AB_count.readlines()
	types=["A","C","G","T"]
	liste={}
	print("Loaded files AB sub percent")
	
	dist=[]
	
	for sub in types: 
		liste[sub]=[]
		
	#print(liste)

	for l in AB:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4]) 
				dist.append(int(line[0]))
				for t in liste: 
					liste[t].append((int(line[num])/tot)*100)
					num+=1
			
	plt.figure(figsize=(10,10))
	plt.plot(dist,Liste["A"], color='dodgerblue', linestyle='dashed')
	plt.plot(dist,Liste["T"], color='mediumaquamarine', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Percentages of A/T in poly A/T around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Percentages of A/T", fontsize=16)
	filename= "A_T__poly_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,Liste["G"], color='orange', linestyle='dashed')
	plt.plot(dist,Liste["C"], color='tomato', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Percentages of G/C in poly G/C around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Percentages of G/C", fontsize=16)
	filename= "G_C_poly_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
		

def GC_content(bar_AB_count, output_directory):
	AB_count=open(bar_AB_count, "r")
	AB=AB_count.readlines()
	print("Loaded files GC")
	
	dist=[]
	GC_per=[]

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
				GC_per.append(GC)

	plt.figure(figsize=(10,10))		
	plt.plot(dist,GC_per, color='orangered', linestyle='dashed')
	plt.axes().minorticks_on()
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("GC % around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("GC percentage", fontsize=16)
	filename= "GC_poly_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)	
	plt.clf()

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichier des bases: 
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")			
	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output directory', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/AB_plots/")			

	
	args = parser.parse_args()
	
	ancestral_sites(args.input_AB, args.output_dir)
	ancestral_sub_sites(args.input_AB, args.output_dir)
	ancestral_sub_sites_percent(args.input_AB, args.output_dir)
	GC_content(args.input_AB, args.output_dir)
	
if "__main__" == __name__:
	main()

