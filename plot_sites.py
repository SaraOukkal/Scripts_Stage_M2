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


def ancestral_sites(bar_AB_count, output_directory): #Plot du nombre de sites par distance aux barrières
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

	plt.plot(dist,AncestralBase, color='#0000cc')
	plt.axvline(0, color='red', linewidth=2, alpha=0.5)
	plt.axhline(1741560, color='green', linewidth=2, linestyle='dashed', alpha=0.5) #Nombre de bords de NIEBs dans interNIEBs >= 1000pb
	plt.axvspan(min(dist), 0, zorder=1, alpha=0.1, color='#cc0000', label='Inside NIEBs')
	plt.axvspan(0, 500, zorder=1, alpha=0.1, color='#00cccc', label='Inter NIEBs')
	plt.title("Number of ancestral sites around NIEBs")
	plt.xlabel("Distance from NIEBs")
	plt.ylabel("Ancestral Bases")
	filename= "AB_plot.png"
	filepath=os.path.join(output_directory, filename)
	plt.savefig(filepath)
	plt.clf()
	
def ancestral_sub_sites(bar_AB_count, output_directory): #plot du nombre de sites par distace aux barrières par type de base
	AB_count=open(bar_AB_count, "r")
	AB=AB_count.readlines()
	types=["A","C","G","T"]
	liste={}
	print("Loaded files AB sub")
	
	dist=[]
	
	for sub in types: 
		liste[sub]=[]
		
	for l in AB:
		if not l.startswith("d"): #Ignore le header du fichier
			line=l.strip().split("\t")
			num=1
			distance=int(line[0])
			if distance >= -50 and distance <= 500 :
				dist.append(int(line[0]))
			
				for t in liste: 
					liste[t].append(int(line[num]))
					num+=1
			
	for i in types:
				#print(liste[i])
				plt.plot(dist,liste[i], color='#0000cc')
				plt.axvline(0, color='red', linewidth=2, label='NIEBs borders')
				plt.axvspan(min(dist), 0, zorder=1, alpha=0.1, color='#cc0000', label='Inside NIEBs')
				plt.axvspan(0, 500, zorder=1, alpha=0.1, color='#00cccc', label='Inter NIEBs')
				plt.title("Number of ancestral base sites per type and per distance to NIEBs")
				plt.xlabel("Distance from NIEBs")
				plt.ylabel("Ancestal Bases")
				filename= "%s.png" % i
				filepath=os.path.join(output_directory, filename)
				plt.savefig(filepath)
				plt.clf()	
				

def ancestral_sub_sites_percent(bar_AB_count, output_directory): #plot du pourcentage de bases par distace aux barrières
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
			
	for i in types:
				#print(liste[i])
				plt.plot(dist,liste[i], color='#0000cc')
				plt.axvline(0, color='red', linewidth=2, label='NIEBs borders')
				plt.axvspan(min(dist), 0, zorder=1, alpha=0.1, color='#cc0000', label='Inside NIEBs')
				plt.axvspan(0, 500, zorder=1, alpha=0.1, color='#00cccc', label='Inter NIEBs')
				plt.title("Percentage of bases per distance to NIEBs")
				plt.xlabel("Distance from NIEBs")
				plt.ylabel("Bases percentage")
				filename= "%spercent.png" % i
				filepath=os.path.join(output_directory, filename)
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

			
	plt.plot(dist,GC_per, color='#0000cc')
	plt.axhline(37.9, color='green', linewidth=2, linestyle='dashed', alpha=0.5)#%GC moyen du génome
	plt.axvline(0, color='red', linewidth=2, label='NIEBs borders')
	plt.axvspan(min(dist), 0, zorder=1, alpha=0.1, color='#cc0000', label='Inside NIEBs')
	plt.axvspan(0, 500, zorder=1, alpha=0.1, color='#00cccc', label='Inter NIEBs')
	plt.title("GC % around NIEBs")
	plt.xlabel("Distance from NIEBs")
	plt.ylabel("GC percentage")
	filename= "GC_plot.png"
	filepath=os.path.join(output_directory, filename)
	plt.savefig(filepath)	
	plt.clf()
		

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
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

