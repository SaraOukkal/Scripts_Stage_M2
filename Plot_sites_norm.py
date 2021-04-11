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


def ancestral_sites(bar_AB_count, bar_genome_count, output): #Plot du nombre de sites par distance aux barrières
	Gen_count=open(bar_genome_count, "r")
	AB_count=open(bar_AB_count, "r") #Fichier avec nombe de bases par type en fonction de la distance aux barrières
	Gen=Gen_count.readlines()
	AB=AB_count.readlines()
	print("Loaded files")
	
	dist=[]
	AncestralBase=[]
	GenBase=[]
	
	for l in AB::
			if not l.startswith("d"): #Ignore le header du fichier
				line1=l.strip().split("\t")
				distance1=int(line1[0])
				if distance1 <= 500 : 
					for a in Gen:
						if not a.startswith("d"): #Ignore le header du fichier
							line2=a.strip().split("\t")
							distance2=int(line2[0])
							
							if distance1 == distance2:
								dist.append(int(line1[0]))
								num_AB=int(line1[1]) + int(line1[2]) + int(line1[3]) + int(line1[4])
								num_Gen=int(line2[1]) + int(line2[2]) + int(line2[3]) + int(line2[4])
								norm=num_AB/num_Gen
								AncestralBase.append(norm)

	plt.plot(dist,AncestralBase, color='#0000cc')
	plt.axvline(0, color='red', linewidth=2, label='NIEBs borders')
	plt.axvspan(-50, 0, zorder=1, alpha=0.1, color='#cc0000', label='Inside NIEBs')
	plt.axvspan(0, 500, zorder=1, alpha=0.1, color='#00cccc', label='Inter NIEBs')
	plt.title("Number of AB sites normalized on all sites around NIEBs")
	plt.xlabel("Distance from NIEBs")
	plt.ylabel("Ancestral Bases/ All bases")
	plt.savefig(output)
	plt.clf()	

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")	
	parser.add_argument('-AB', '--input_Gen', type=str, help='Path to all genome bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step5_results/Genome_count_bar.txt")			
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output plot', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/AB_norm_plot.png")			

	
	args = parser.parse_args()
	
	ancestral_sites(args.input_AB, args.input_Gen, args.output)
	
if "__main__" == __name__:
	main()

