#!/usr/bin/python3

#Importer packages: 
import argparse 

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
	
	pyplot.hist(distribution)
	plt.title("Intervals size distribution")
	plt.xlabel("Size (in nt)")
	plt.ylabel("Number of intervals")

	plt.savefig(output)	
	plt.clf()
		
	
				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier .bed des barrières chez le chimpanzé: 
	parser.add_argument('-inter', '--input_inter', type=str, help='Path to intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/H_intervals.bed")
	parser.add_argument('-out', '--output', type=str, help='Path to output plot', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step4_results/Plots/")
	
	#Lancer fonctions: 
	distribution= distri(args.input_inter)
	plot(distribution, args.output)
	

if "__main__" == __name__:
	main()
	args = parser.parse_args()
	
