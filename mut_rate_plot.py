#!/usr/bin/python3

#Importer packages: 
import argparse 
import numpy as np
import matplotlib.pyplot as plt
import os

def mut_rate_plot(mut_rate, output):
	data=open(mut_rate, "r")
	mut_R=data.readlines()
	print("Loaded files")
	
	dist=[]
	MR=[]
	
	for l in mut_R:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
			MR.append(float(line_MR[1])*100)

	plt.plot(dist,MR)
	plt.title("Mutation rates around NIEBs")
	plt.xlabel("Distance from NIEBs")
	plt.ylabel("Mutation rate (%)")
	plt.savefig(output)

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-MR', '--input_MR', type=str, help='Path to mutation rates', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")					

	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/MR_plot.png")			

	
	args = parser.parse_args()
	
	mut_rate_plot(args.input_MR, args.output)
	
	
if "__main__" == __name__:
	main()


