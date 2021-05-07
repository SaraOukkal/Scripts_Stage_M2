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


def ancestral_sites(input_1, input_2, output_dir): 
	"""
	Plot du nombre de sites d'AB et Gen par distance aux barrières
	"""
	sp1_count=open(input_1, "r")
	sp2_count=open(input_2, "r") 
	sp1=sp1_count.readlines()
	sp2=sp2_count.readlines()
	
	dist=[]
	sp1bases=[]
	sp2bases=[]
	
	for l in sp1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				num_AB=0
				num_AB=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				sp1bases.append(num_AB)
	
	for a in sp2:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				num_B=0
				num_B=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])
				sp2bases.append(num_B)			
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,sp1bases, color='darkblue')
	plt.plot(dist,sp2bases, color='dodgerblue')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=13)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of homopolysites around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("homopolybases", fontsize=16)
	filename= "numpoly_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
def ancestral_sub_sites(input_1, input_2, output_dir): 
	"""
	plot du nombre de sites par distance aux barrières par type de base
	"""
	sp1_count=open(input_1, "r")
	sp2_count=open(input_2, "r") 
	sp1=sp1_count.readlines()
	sp2=sp2_count.readlines()

	types=["A","C","G","T"]
	sp1_List={}
	sp2_List={}
	
	dist=[]
	
	for sub in types: 
		sp1_List[sub]=[]
		sp2_List[sub]=[]
		
	for l in sp1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				for t in sp1_List: 
					sp1_List[t].append(int(line[num]))
					num+=1
		
	for a in sp2:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				for t in sp2_List: 
					sp2_List[t].append(int(line[num]))
					num+=1
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,sp1_List["A"], color='darkblue')
	plt.plot(dist,sp1_List["T"], color='darkgreen')
	plt.plot(dist,sp2_List["A"], color='dodgerblue')
	plt.plot(dist,sp2_List["T"], color='limegreen')
	plt.plot(dist,sp1_List["G"], color='chocolate')
	plt.plot(dist,sp1_List["C"], color='firebrick')
	plt.plot(dist,sp2_List["G"], color='orange')
	plt.plot(dist,sp2_List["C"], color='tomato')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=15)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Number of homopolynucleotides around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Number of A/T/G/C", fontsize=16)
	filename= "poly_ATCG.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
				

def ancestral_sub_sites_percent(input_1, input_2, output_dir): 
	"""
	Plot du pourcentage de chaque type de base par distance aux barrières
	"""
	
	sp1_count=open(input_1, "r")
	sp2_count=open(input_2, "r") 
	sp1=sp1_count.readlines()
	sp2=sp2_count.readlines()

	types=["A","C","G","T"]
	sp1_List={}
	sp2_List={}
	
	dist=[]
	
	for sub in types: 
		sp1_List[sub]=[]
		sp2_List[sub]=[]
		
	for l in sp1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4]) 
				dist.append(int(line[0]))
				for t in sp1_List: 
					sp1_List[t].append((int(line[num])/tot)*100)
					num+=1
		
	for a in sp2:
		if not a.startswith("d"): #Ignore le header du fichier
			line=a.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				tot=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4]) 
				for t in sp2_List: 
					sp2_List[t].append((int(line[num])/tot)*100)
					num+=1
	
	plt.figure(figsize=(10,10))
	plt.plot(dist,sp1_List["A"], color='darkblue')
	plt.plot(dist,sp1_List["T"], color='darkgreen')
	plt.plot(dist,sp2_List["A"], color='dodgerblue')
	plt.plot(dist,sp2_List["T"], color='limegreen')
	plt.plot(dist,sp1_List["G"], color='chocolate')
	plt.plot(dist,sp1_List["C"], color='firebrick')
	plt.plot(dist,sp2_List["G"], color='orange')
	plt.plot(dist,sp2_List["C"], color='tomato')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Percentages of homopolynucleotides around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Percentages of A/T/G/C", fontsize=16)
	filename= "poly_ATCG_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	

def Poly_on_Gen(poly_1, poly_2, Gen_1, Gen_2, output_dir):
	poly1=open(poly_1, "r")
	poly2=open(poly_2, "r") 
	Gen1=open(Gen_1, "r")
	Gen2= open(Gen_2, "r")
	
	types=["A","C","G","T"]
	sp1_List={}
	sp2_List={}
	
	dist=[]
	
	for sub in types: 
		sp1_List[sub]=[]
		sp2_List[sub]=[]
		
	for l in poly1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				for a in Gen1: 
					line2=a.strip().split("\t")
					distance2=int(line2[0])
					if distance == distance2:
						dist.append(distance)
						tot= int(line2[1]) + int(line2[2]) + int(line2[3]) + int(line2[4])
						for t in sp1_List: 
							sp1_List[t].append((int(line[num])/tot)*100)
							num+=1
						break
		
	for l in poly2:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				for a in Gen2: 
					line2=a.strip().split("\t")
					distance2=int(line2[0])
					if distance == distance2:
						dist.append(distance)
						tot= int(line2[1]) + int(line2[2]) + int(line2[3]) + int(line2[4])
						for t in sp2_List: 
							sp2_List[t].append((int(line[num])/tot)*100)
							num+=1
						break
						
	plt.figure(figsize=(10,10))
	plt.plot(dist,sp1_List["A"], color='darkblue')
	plt.plot(dist,sp1_List["T"], color='darkgreen')
	plt.plot(dist,sp2_List["A"], color='dodgerblue')
	plt.plot(dist,sp2_List["T"], color='limegreen')
	plt.plot(dist,sp1_List["G"], color='chocolate')
	plt.plot(dist,sp1_List["C"], color='firebrick')
	plt.plot(dist,sp2_List["G"], color='orange')
	plt.plot(dist,sp2_List["C"], color='tomato')
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Homopolynucleotides coverage (%) around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Homopolynucletotides coverage (%)", fontsize=16)
	filename= "poly_ATCG_cov.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()	
		
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichier des bases: 
	parser.add_argument('-sp1', '--input_1', type=str, help='Path to polybases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	parser.add_argument('-sp2', '--input_2', type=str, help='Path to polybases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	
	parser.add_argument('-Gen1', '--input_Gen1', type=str, help='Path to bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	parser.add_argument('-Gen2', '--input_Gen2', type=str, help='Path to bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")
	
				
	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output directory', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/AB_plots/")			

	
	args = parser.parse_args()
	
	ancestral_sites(args.input_1, args.input_2, args.output_dir)
	ancestral_sub_sites(args.input_1, args.input_2, args.output_dir)
	ancestral_sub_sites_percent(args.input_1, args.input_2, args.output_dir)
	Poly_on_Gen(args.input_1, args.input_2, args.input_Gen1, args.input_Gen2, args.output_dir)
	
if "__main__" == __name__:
	main()

