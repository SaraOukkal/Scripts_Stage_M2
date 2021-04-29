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

def mut_number(input_Mut, output_dir):
	Mut_count=open(input_Mut, "r") 
	Mut=Mut_count.readlines()
	
	dist=[]
	Mutations=[]
	
	for l in Mut:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				num=0
				num=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])+ int(line[5])+ int(line[6])+ int(line[7])+ int(line[8])+ int(line[9])+ int(line[10])+ int(line[11])+ int(line[12])
				print(num)
				Mutations.append(num)
	
	plt.figure(figsize=(10,10))	
	plt.plot(dist,Mutations, color='violet')
	plt.axvline(0, color='black', alpha=0.5)
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.title("Number of mutations around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutations", fontsize=16)
	filename= "Mut_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()


def mut_rate_plot(mut_rate, output_dir):
	data=open(mut_rate, "r")
	mut_R=data.readlines()
	print("Loaded files")
	
	dist=[]
	MR=[]
	
	for l in mut_R:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
				MR.append(float(line_MR[1])*100)
	
	Smooth_MR=Lissage(MR)

	plt.figure(figsize=(10,10))	
	plt.plot(dist,Smooth_MR, color="darkviolet")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Mutation rates around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "Mut_rate.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	

def mut_rates_plots(sub_mut_rates, output_dir):
	
	data=open(sub_mut_rates, "r")
	MR=data.readlines()
	types=["MR_AT","MR_AC","MR_AG","MR_CT","MR_CA","MR_CG","MR_GT","MR_GA","MR_GC","MR_TA","MR_TC","MR_TG"]
	liste={}
	
	MR_dist=[]
	
	for sub in types: 
		liste[sub]=[]
		

	for l in MR:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			num=1
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				MR_dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
			
				for t in liste: 
					liste[t].append(float(line_MR[num])*100) #Calcule le taux de mutations en pourcentage 
					num+=1
			
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,liste["MR_AT"], color="blue")
	plt.plot(MR_dist,liste["MR_TA"], color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("AT_TA", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "AT_TA.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,liste["MR_AC"], color="blue")
	plt.plot(MR_dist,liste["MR_TG"], color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("AC_TG", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "AC_TG.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.plot(MR_dist,liste["MR_CG"], color="red")
	plt.plot(MR_dist,liste["MR_GC"], color="orange")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("CG_GC", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "CG_GC.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,liste["MR_CT"], color="red")
	plt.plot(MR_dist,liste["MR_GA"], color="orange")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("CT_GA", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "CT_GA.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,liste["MR_AG"], color="blue")
	plt.plot(MR_dist,liste["MR_TC"], color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("AG_TC", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "AG_TC.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()

	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,liste["MR_GT"], color="orange")
	plt.plot(MR_dist,liste["MR_CA"], color="red")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("GT_CA", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "GT_CA.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()

def mut_rate_base(sub_mut_rates, output_dir):
	data=open(sub_mut_rates, "r")
	MR=data.readlines()

	bases=["A","C","G","T"]
	liste={}
	
	MR_dist=[]
	
	for sub in bases: 
		liste[sub]=[]

	for l in MR:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				MR_dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
				
				liste["A"].append((float(line_MR[1])+float(line_MR[2])+float(line_MR[3]))*100) #Calcule le taux de mutations de A en pourcentage 
				liste["C"].append((float(line_MR[4])+float(line_MR[5])+float(line_MR[6]))*100) 
				liste["G"].append((float(line_MR[7])+float(line_MR[8])+float(line_MR[9]))*100) 
				liste["T"].append((float(line_MR[10])+float(line_MR[11])+float(line_MR[12]))*100) 
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,liste["A"], color="blue")
	plt.plot(MR_dist,liste["C"], color="red")
	plt.plot(MR_dist,liste["T"], color="green")
	plt.plot(MR_dist,liste["G"], color="orange")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.title("Mutation percentage of each nucleotide type around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "ACGT_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
def Lissage(Liste):
	
	Smooth=[]
	
	for i in range(len(Liste):
		if i >=10 and i<=len(Liste)-10:
			Sum=Liste[i-5] + Liste[i-4]+ Liste[i-3]+ Liste[i-2]+ Liste[i-1]+ Liste[i]+ Liste[i+1]+ Liste[i+2]+ Liste[i+3]+ Liste[i+4]+ Liste[i+5]
			Mean=Sum/10
			Smooth.append(Mean)
			
		else: 
			Smooth.append(Liste[i])
	
	return Smooth

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-MR', '--input_MR', type=str, help='Path to mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")	
	parser.add_argument('-MRS', '--input_MRS', type=str, help='Path to sub mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")	
	parser.add_argument('-Mut', '--input_Mut', type=str, help='Path to mutations count', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")				

	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/Mut_plots/")			

	
	args = parser.parse_args()
	mut_number(args.input_Mut, args.output_dir)
	mut_rate_plot(args.input_MR, args.output_dir)
	mut_rates_plots(args.input_MRS, args.output_dir)
	mut_rate_base(args.input_MRS, args.output_dir)
	

	
if "__main__" == __name__:
	main()


