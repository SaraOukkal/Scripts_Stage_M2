#!/usr/bin/python3

"""
Genère le fichier .bed à partir d'un aligement Homme/espèce pour réaliser le bedtools intersect
"""

#Importer packages:
import subprocess 
import argparse 

#Charger fichiers dans un dictionnaire: 

def load_files(input_dir,output_file):
	print("Start loading files")
	lof = subprocess.check_output(["ls", input_dir]).decode('utf-8').strip().split("\n")
	dico={}
	for filename in lof: 
		chrom=[]
		chromH=filename.replace(".dat","")
		filepath=input_dir+filename
		data=open(filepath,"r")
		data2=data.readlines()
		for a in range(len(data2)):
			line=data2[a].strip().split("\t")
			ali={}
			ali["chromH"]=chromH	
			ali["st_H"]=line[0]
			ali["end_H"]=line[1]
			chrom.append(ali)	
		dico[chromH]=chrom
		
	f=open(output_file,'w')
	for key in sorted(dico.keys()): 
		for ali in dico[key]:
			f.write("{}\t{}\t{}\n".format(ali["chromH"],ali["st_H"],ali["end_H"]))
	f.close()
	print("end")
	


def main():
	parser = argparse.ArgumentParser()
	
	#Input directory :
	parser.add_argument('-i', '--input_dir', type=str, help='Path to input directory chr files ', default ="/home/soukkal/Bureau/Projet/20200616_hg38_to_panTro5/")
	
	#Output file: 
	parser.add_argument('-o', '--output_file', type=str, help='Path to output file', default = "/home/soukkal/Bureau/Projet/Step1_results/20200728_HC.bed")
	
	args = parser.parse_args()
	
	#Vérifie que le pathway vers un dossier se termine par / sinon ajoute un /:
	if not args.input_dir.endswith("/"): 
		args.input_dir+="/"
	
	load_files(args.input_dir,args.output_file)
		
if "__main__" == __name__:
	main()

#Commandes terminal: 
##python3 script_create_bed.py 
##python3 script_create_bed.py -i /home/soukkal/Bureau/Projet/20200724_hg38_to_gorGor4/ -o /home/soukkal/Bureau/Projet/Step1_results/20200728_HG.bed
##python3 script_create_bed.py -i /home/soukkal/Bureau/Projet/20200724_hg38_to_ponAbe2/ -o /home/soukkal/Bureau/Projet/Step1_results/20200728_HO.bed

