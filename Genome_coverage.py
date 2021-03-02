#!/usr/bin/python3

#Importer packages: 
import argparse 

def coverage (intervals):
	f=open(intervals,"r")
	inter=f.readlines()
	count=0
	
	for l in inter:
		line=l.strip().split("\t")
		st=int(line[1])
		end=int(line[2])	
		size=end-st
		count+=size
		
	print(count)


				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier .bed des barrières chez le chimpanzé: 
	parser.add_argument('-inter', '--input_inter', type=str, help='Path to intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Chimp_Step3_results/H_intervals.bed")
	
	args = parser.parse_args()
	
	#Lancer fonction: 
	coverage(args.input_inter)
	

if "__main__" == __name__:
	main()
