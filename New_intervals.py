#!/usr/bin/python3

"""
Génère les alignements communs aux 4 espèces en prenant l'Homme comme référence
Necessite d'avoir fait en amont un bedtools intersect
"""

#Importer packages: 
import subprocess
import argparse 
import time


def change_intervals(bedpath,input_dir1,input_dir2,input_dir3,output_file):
	"""
	Genère puis écrit dans un fichier output les nouveaux intervalles pour les trois espèces.
	"""
	
	f=open(bedpath, "r")	#ouvre le fichier output de bedtools intersect
	out=open(output_file,'w')
	done_chrom=[]
	sp1=[] #Crée les listes pour pouvoir les supprimer à la ligne 1 
	sp2=[]
	sp3=[]
	for line in f : 
		split_line=line.strip().split("\t")
		chrom_Bed=split_line[0]	#détermine le chromosome Humain
		st_Bed=split_line[1]	#Début de l'intervalle pour l'Homme
		end_Bed=split_line[2]	#Fin de l'intervalle pour l'Homme	

		if chrom_Bed not in done_chrom:
			
			del sp1 #Supprime les listes générées précédement 
			del sp2
			del sp3
			
			sp1= load_file(input_dir1, chrom_Bed)
			sp2= load_file(input_dir2, chrom_Bed)
			sp3= load_file(input_dir3, chrom_Bed)
			
			print(chrom_Bed)
			done_chrom.append(chrom_Bed) #si le chromosome n'est pas déja dans le dictionnaire cela crée une liste vide dont la clé est le chromosome
			index1=0 #réinitialise l'index 
			index2=0 #réinitialise l'index 
			index3=0 #réinitialise l'index 

		res1=calculate_intervals(sp1,st_Bed,end_Bed,index1) 
		new_interval=res1[0]	
		index1=res1[1]
	
		res2=calculate_intervals(sp2,st_Bed,end_Bed,index2) 
		new_interval_2=res2[0]	
		index2=res2[1]
		
		new_interval.append(new_interval_2[2])
		new_interval.append(new_interval_2[3])
		new_interval.append(new_interval_2[4])
		new_interval.append(new_interval_2[5])
		
		res3=calculate_intervals(sp3,st_Bed,end_Bed,index3) 
		new_interval_3=res3[0]	
		index3=res3[1]
		
		new_interval.append(new_interval_3[2])
		new_interval.append(new_interval_3[3])
		new_interval.append(new_interval_3[4])
		new_interval.append(new_interval_3[5])
		#~ print(new_interval)	
		
		out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom_Bed,new_interval[0],new_interval[1],new_interval[2],new_interval[3],new_interval[4],new_interval[5],new_interval[6],new_interval[7],new_interval[8],new_interval[9],new_interval[10],new_interval[11],new_interval[12],new_interval[13]))
		
		

def load_file(input_dir,chrom_Bed):
	"""
	Charge le fichier d'alignements pour une espèce correspondant au chromosome traité.
	"""
	
	sp=[]
	filepath=input_dir+chrom_Bed+".dat"
	fil=open(filepath, "r")
	data=fil.readlines()
	for a in range(len(data)):
		line=data[a].strip().split("\t")
		ali=[]
		ali.append(line[0])	#start Human 0
		ali.append(line[1])	#end Human 1
		ali.append(line[2])	#chromosome specie 2
		ali.append(line[3])	#start specie 3
		ali.append(line[4])	#end specie 4
		ali.append(line[5])	#specie strand 5
		ali.append(line[6])	#chain id 6
		sp.append(ali)	#Ajoute ali à la liste sp1 
	
	fil.close()
	return sp

	

def calculate_intervals(chro,st_Bed,end_Bed,index):
	"""	
	Calcule les nouveaux intervalles des espèces
	"""
	for i in range (index,len(chro),1):	
		if int(st_Bed) >= int(chro[i][0]) and int(end_Bed) <= int(chro[i][1]):
			index=i
			interval=[]
			diff_st=int(st_Bed)-int(chro[i][0]) #Calcule la différence entre le start H avant et après bedtools intersect
			diff_end=int(chro[i][1])-int(end_Bed)	#Calcule la différence entre le end H avant et après bedtools intersect
			
			if chro[i][5]=="+": #Si c'est le brin +
				new_st_sp=int(chro[i][3])+int(diff_st)	#Calcule le start espèce à l'aide de la différence start
				new_end_sp=int(chro[i][4])-int(diff_end)	#Calcule le end espèce à l'aide de la différence end
				
			elif chro[i][5]=="-": #Si c'est le brin -
				new_st_sp=int(chro[i][3])+int(diff_end)
				new_end_sp=int(chro[i][4])-int(diff_st)

			interval.append(st_Bed)	
			interval.append(end_Bed)	
			interval.append(chro[i][2])	#Chromosome specie
			interval.append(new_st_sp)
			interval.append(new_end_sp)
			interval.append(chro[i][5]) #specie strand
			
			break					
	
	return [interval, index] #retourne une liste contenant le nouvel l'intervalle et l'index mis à jour 



def main(): 
	parser = argparse.ArgumentParser()
	
	#Dossiers pour les 3 espèces:
	parser.add_argument('-sp1', '--input_dir1', type=str, help='Path to sp1 input directory', default ="/home/soukkal/Bureau/Projet/20200616_hg38_to_panTro5/")
	parser.add_argument('-sp2', '--input_dir2', type=str, help='Path to sp2 input directory', default ="/home/soukkal/Bureau/Projet/20200724_hg38_to_gorGor4/")
	parser.add_argument('-sp3', '--input_dir3', type=str, help='Path to sp3 input directory', default ="/home/soukkal/Bureau/Projet/20200724_hg38_to_ponAbe2/")
	
	#Fichier de sortie de bedtools intersect:
	parser.add_argument('-bed', '--input_file', type=str, help='Path to bed file ', default ="/home/soukkal/Bureau/Projet/Step1_results/HC_HG_HO_intersect.bed")
	
	#Fichier output nouveaux alignements: 
	parser.add_argument('-o', '--output_file', type=str, help='Path to output file ', default ="/home/soukkal/Bureau/Projet/Step1_results/new_intervals.txt")
	
	args = parser.parse_args()
	
	#Vérifie que le pathway vers un dossier se termine par / sinon ajoute un /:
	if not args.input_dir1.endswith("/"): 
		args.input_dir1+="/"
	if not args.input_dir2.endswith("/"):
		args.input_dir2+="/"	
	if not args.input_dir3.endswith("/"):
		args.input_dir3+="/"	
	
	new_dict=change_intervals(args.input_file,args.input_dir1,args.input_dir2,args.input_dir3,args.output_file)


	
if "__main__" == __name__:
	main()
