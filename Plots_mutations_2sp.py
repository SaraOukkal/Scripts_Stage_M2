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

def mut_number(input_Mut1,input_Mut2, output_dir):
	Mut_count1=open(input_Mut1, "r") 
	Mut1=Mut_count1.readlines()
	Mut_count2=open(input_Mut2, "r") 
	Mut2=Mut_count2.readlines()
	
	dist=[]
	Mutations1=[]
	Mutations2=[]
	
	for l in Mut1:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
				num=0
				num=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])+ int(line[5])+ int(line[6])+ int(line[7])+ int(line[8])+ int(line[9])+ int(line[10])+ int(line[11])+ int(line[12])
				print(num)
				Mutations1.append(num)
	
	for l in Mut2:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				num=0
				num=int(line[1]) + int(line[2]) + int(line[3]) + int(line[4])+ int(line[5])+ int(line[6])+ int(line[7])+ int(line[8])+ int(line[9])+ int(line[10])+ int(line[11])+ int(line[12])
				print(num)
				Mutations2.append(num)
	
	plt.figure(figsize=(10,10))	
	plt.plot(dist,Mutations1, color='darkviolet')
	plt.plot(dist,Mutations2, color='violet')
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
	plt.title("Number of mutations around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutations", fontsize=16)
	filename= "Mut_plot.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()


def mut_rate_plot(mut_rate1, mut_rate2, output_dir):
	data1=open(mut_rate1, "r")
	mut_R1=data1.readlines()
	data2=open(mut_rate2, "r")
	mut_R2=data2.readlines()
	
	dist=[]
	MR1=[]
	MR2=[]
	
	for l in mut_R1:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
				MR1.append(float(line_MR[1])*100)
				
	for l in mut_R2:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				MR2.append(float(line_MR[1])*100)

	Smooth_MR1=Lissage(MR1)
	Smooth_MR2=Lissage(MR2)

	plt.figure(figsize=(10,10))	
	plt.plot(dist,Smooth_MR1, color="darkviolet")
	plt.plot(dist,Smooth_MR2, color="violet")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("Mutation rates around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "Mut_rate.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	

def mut_rates_plots(sub_mut_rates1, sub_mut_rates2, output_dir):
	
	data1=open(sub_mut_rates1, "r")
	MR1=data1.readlines()
	data2=open(sub_mut_rates2, "r")
	MR2=data2.readlines()
	
	types=["MR_AT","MR_AC","MR_AG","MR_CT","MR_CA","MR_CG","MR_GT","MR_GA","MR_GC","MR_TA","MR_TC","MR_TG"]
	liste1={}
	liste2={}
	
	MR_dist=[]
	
	for sub in types: 
		liste1[sub]=[]
		liste2[sub]=[]

	for l in MR1:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			num=1
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				MR_dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
				for t in liste: 
					liste1[t].append(float(line_MR[num])*100) #Calcule le taux de mutations en pourcentage 
					num+=1
					
	for l in MR2:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			num=1
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				for t in liste: 
					liste2[t].append(float(line_MR[num])*100) #Calcule le taux de mutations en pourcentage 
					num+=1
	
	Smooth_AT1=Lissage(liste1["MR_AT"])
	Smooth_TA1=Lissage(liste1["MR_TA"])
	Smooth_AT2=Lissage(liste2["MR_AT"])
	Smooth_TA2=Lissage(liste2["MR_TA"])
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_AT1, color="darkblue")
	plt.plot(MR_dist,Smooth_TA1, color="darkgreen")
	plt.plot(MR_dist,Smooth_AT2, color="blue")
	plt.plot(MR_dist,Smooth_TA2, color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("AT_TA", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "AT_TA.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()

	Smooth_AC1=Lissage(liste1["MR_AC"])
	Smooth_TG1=Lissage(liste1["MR_TG"])
	Smooth_AC2=Lissage(liste2["MR_AC"])
	Smooth_TG2=Lissage(liste2["MR_TG"])
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_AC1, color="darkblue")
	plt.plot(MR_dist,Smooth_TG1, color="darkgreen")
	plt.plot(MR_dist,Smooth_AC2, color="blue")
	plt.plot(MR_dist,Smooth_TG2, color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("AC_TG", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "AC_TG.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	Smooth_CG1=Lissage(liste1["MR_CG"])
	Smooth_GC1=Lissage(liste1["MR_GC"])
	Smooth_CG2=Lissage(liste2["MR_CG"])
	Smooth_GC2=Lissage(liste2["MR_GC"])
	
	plt.plot(MR_dist,Smooth_CG1, color="firebrick")
	plt.plot(MR_dist,Smooth_GC1, color="darkorange")
	plt.plot(MR_dist,Smooth_CG2, color="red")
	plt.plot(MR_dist,Smooth_GC2, color="orange")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("CG_GC", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "CG_GC.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	Smooth_CT1=Lissage(liste1["MR_CT"])
	Smooth_GA1=Lissage(liste1["MR_GA"])
	Smooth_CT2=Lissage(liste2["MR_CT"])
	Smooth_GA2=Lissage(liste2["MR_GA"])
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_CT1, color="firebrick")
	plt.plot(MR_dist,Smooth_GA1, color="darkorange")
	plt.plot(MR_dist,Smooth_CT2, color="red")
	plt.plot(MR_dist,Smooth_GA2, color="orange")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("CT_GA", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "CT_GA.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	Smooth_AG1=Lissage(liste1["MR_AG"])
	Smooth_TC1=Lissage(liste1["MR_TC"])
	Smooth_AG2=Lissage(liste2["MR_AG"])
	Smooth_TC2=Lissage(liste2["MR_TC"])
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_AG1, color="darkblue")
	plt.plot(MR_dist,Smooth_TC1, color="darkgreen")
	plt.plot(MR_dist,Smooth_AG2, color="blue")
	plt.plot(MR_dist,Smooth_TC2, color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("AG_TC", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "AG_TC.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	Smooth_GT1=Lissage(liste1["MR_GT"])
	Smooth_CA1=Lissage(liste1["MR_CA"])
	Smooth_GT2=Lissage(liste2["MR_GT"])
	Smooth_CA2=Lissage(liste2["MR_CA"])

	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_GT1, color="darkorange")
	plt.plot(MR_dist,Smooth_CA1, color="firebrick")
	plt.plot(MR_dist,Smooth_GT2, color="orange")
	plt.plot(MR_dist,Smooth_CA2, color="red")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("GT_CA", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "GT_CA.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()

def mut_rate_base(sub_mut_rates1, sub_mut_rates2, output_dir):
	
	data1=open(sub_mut_rates1, "r")
	MR1=data1.readlines()
	data2=open(sub_mut_rates2, "r")
	MR2=data2.readlines()

	bases=["A","C","G","T"]
	liste1={}
	liste2={}
	
	MR_dist=[]
	
	for sub in bases: 
		liste1[sub]=[]
		liste2[sub]=[]

	for l in MR1:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				MR_dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
				
				liste1["A"].append((float(line_MR[1])+float(line_MR[2])+float(line_MR[3]))*100) #Calcule le taux de mutations de A en pourcentage 
				liste1["C"].append((float(line_MR[4])+float(line_MR[5])+float(line_MR[6]))*100) 
				liste1["G"].append((float(line_MR[7])+float(line_MR[8])+float(line_MR[9]))*100) 
				liste1["T"].append((float(line_MR[10])+float(line_MR[11])+float(line_MR[12]))*100) 
				
	for l in MR2:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			distance=int(line_MR[0])
			if distance >= -50 and distance <= 500 :
				liste2["A"].append((float(line_MR[1])+float(line_MR[2])+float(line_MR[3]))*100) #Calcule le taux de mutations de A en pourcentage 
				liste2["C"].append((float(line_MR[4])+float(line_MR[5])+float(line_MR[6]))*100) 
				liste2["G"].append((float(line_MR[7])+float(line_MR[8])+float(line_MR[9]))*100) 
				liste2["T"].append((float(line_MR[10])+float(line_MR[11])+float(line_MR[12]))*100) 
	
	
	Smooth_A1=Lissage(liste1["A"])
	Smooth_T1=Lissage(liste1["T"])
	Smooth_C1=Lissage(liste1["C"])
	Smooth_G1=Lissage(liste1["G"])
	
	Smooth_A2=Lissage(liste2["A"])
	Smooth_T2=Lissage(liste2["T"])
	Smooth_C2=Lissage(liste2["C"])
	Smooth_G2=Lissage(liste2["G"])
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_A1, color="darkblue")
	plt.plot(MR_dist,Smooth_A2, color="blue")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("A mutation percentage around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "A_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_T1, color="darkgreen")
	plt.plot(MR_dist,Smooth_T2, color="green")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("T mutation percentage around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "T_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_C1, color="firebrick")
	plt.plot(MR_dist,Smooth_C2, color="red")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("C mutation percentage around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "C_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
	plt.figure(figsize=(10,10))	
	plt.plot(MR_dist,Smooth_G1, color="darkorange")
	plt.plot(MR_dist,Smooth_G2, color="orange")
	plt.axes().minorticks_on()
	plt.axes().tick_params(axis='both', which='major', direction='in', length= 8, width=2)
	plt.axes().tick_params(axis='both', which='minor', direction='in', length= 4, width=1.5)
	plt.axes().xaxis.set_major_locator(MultipleLocator(100))
	plt.axes().xaxis.set_minor_locator(MultipleLocator(50))
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.axvline(0, color='black', alpha=0.5)
	plt.axvline(137, color='black', alpha=0.2)
	plt.axvline(280, color='black', alpha=0.2)
	plt.title("G mutation percentage around NIEBs", fontsize=20)
	plt.xlabel("Distance from NIEBs", fontsize=16)
	plt.ylabel("Mutation rate (%)", fontsize=16)
	filename= "G_per.png"
	filepath=os.path.join(output_dir, filename)
	plt.savefig(filepath)
	plt.clf()
	
def Lissage(Liste):
	
	Smooth=[]
	
	for i in range(len(Liste)):
		if i == 0:
			Sum=Liste[i] + Liste[i+1]+ Liste[i+2]+ Liste[i+3]+ Liste[i+4]+ Liste[i+5]
			Mean=Sum/6
			Smooth.append(Mean)
		elif i == 1:
			Sum=Liste[i-1] + Liste[i] + Liste[i+1]+ Liste[i+2]+ Liste[i+3]+ Liste[i+4]+ Liste[i+5]
			Mean=Sum/7
			Smooth.append(Mean)
		elif i == 2:
			Sum=Liste[i-2] +Liste[i-1] + Liste[i] + Liste[i+1]+ Liste[i+2]+ Liste[i+3]+ Liste[i+4]+ Liste[i+5]
			Mean=Sum/8
			Smooth.append(Mean)
		elif i == 3:
			Sum=Liste[i-3] +Liste[i-2] +Liste[i-1] + Liste[i] + Liste[i+1]+ Liste[i+2]+ Liste[i+3]+ Liste[i+4]+ Liste[i+5]
			Mean=Sum/9
			Smooth.append(Mean)
		elif i >3 and i<len(Liste)-5:
			Sum=Liste[i-4]+ Liste[i-3]+ Liste[i-2]+ Liste[i-1]+ Liste[i]+ Liste[i+1]+ Liste[i+2]+ Liste[i+3]+ Liste[i+4]+ Liste[i+5]
			Mean=Sum/10
			Smooth.append(Mean)
		else: 
			Sum=Liste[i-5] + Liste[i-4]+ Liste[i-3]+ Liste[i-2]+ Liste[i-1]+ Liste[i]
			Mean=Sum/6
			Smooth.append(Mean)
	
	return Smooth

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-MR1', '--input_MR1', type=str, help='Path to mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")	
	parser.add_argument('-MRS1', '--input_MRS1', type=str, help='Path to sub mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")	
	parser.add_argument('-Mut1', '--input_Mut1', type=str, help='Path to mutations count', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")		
	
	parser.add_argument('-MR2', '--input_MR2', type=str, help='Path to mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")	
	parser.add_argument('-MRS2', '--input_MRS2', type=str, help='Path to sub mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")	
	parser.add_argument('-Mut2', '--input_Mut2', type=str, help='Path to mutations count', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")				
	



	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/Mut_plots/")			

	
	args = parser.parse_args()
	mut_number(args.input_Mut1,args.input_Mut2, args.output_dir)
	mut_rate_plot(args.input_MR1,args.input_MR2, args.output_dir)
	mut_rates_plots(args.input_MRS1,args.input_MRS2, args.output_dir)
	mut_rate_base(args.input_MRS1,args.input_MRS2, args.output_dir)
	
if "__main__" == __name__:
	main()


