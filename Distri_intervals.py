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


def distri (intervals):
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

def plot (distribution, output):
	plt.hist(distribution, bins=20, color="salmon")
	plt.title("Intervals size distribution")
	plt.xticks(np.arange(25))
	plt.xlabel("Size")
	plt.ylabel("Number of intervals")

	plt.savefig(output)	
	plt.clf()
		
	
				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier .bed des barrières chez le chimpanzé: 
	parser.add_argument('-inter', '--input_inter', type=str, help='Path to intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/H_intervals.bed")
	parser.add_argument('-out', '--output', type=str, help='Path to output plot', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step4_results/Plots/distri.png")
	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	distribution= distri(args.input_inter)
	plot(distribution, args.output)
	

if "__main__" == __name__:
	main()

