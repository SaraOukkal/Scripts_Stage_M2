#!/usr/bin/python3

#Importer packages: 
import argparse 

def select_bar(input_bar,output_file):
	
	out=open(output_file,"w") 
	f=open(input_bar,"r")
	
	data=f.readlines()
	
	pos=[] #Crée la liste pos qui permettra de savoir si une barrière a déja été écrite dans l'output (car inter barrières des deux cotés d'une barrière, pour éviter la redondance
	done_chrom=[] #Crée la liste done_chrom qui permet de suivre l'avancée du programme
	
	for l in data: 
		line=l.strip().split('\t')
		chrom_br=line[0] #Chromosome
		st_br1=int(line[1]) #Start première barrière
		end_br1=int(line[2]) #End première barrière
		st_br2=int(line[3]) #Start Deuxième barrière
		end_br2=int(line[4]) #End deuxième barrière
		inter_br=st_br2-end_br1 #Calcule l'inter barrière
		
		if chrom_br not in done_chrom:
			pos=[]
			done_chrom.append(chrom_br)
			print(chrom_br)
		
		if inter_br >= 1000: #Les cas qu'on veut garder 
			
			if st_br1 in pos: #Si la première barrière a déja été rencontrée (donc déja dans le fichier output) 
				out.write("{}\t{}\t{}\n".format(chrom_br,st_br2,end_br2)) #Ecrit dans l'output uniquement la deuxième barrière
				pos.append(st_br2) #ajoute à Pos la deuxième barrière pour qu'elle soit catégorisée comme déja rencontrée

#Je n'ajoute que la deuxième barrière à pos à chaque fois, car vu que le fichier avance dans le génome et donc c'est uniquement la deuxième barrière du couple qui risque d'être retrouvée en première barrière d'un autre couple. 
				
			else: #Si la première barrière n'a pas déja été rencontrée 
				out.write("{}\t{}\t{}\n".format(chrom_br,st_br1,end_br1))
				out.write("{}\t{}\t{}\n".format(chrom_br,st_br2,end_br2))
				pos.append(st_br2) 
			
			
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	parser.add_argument('-bar', '--input_bar', type=str, help="/home/saraoukkal/Documents/Stage_M1/Barriers/panTro5/20200228_interSmallNFR.dat")					

	#fichier output:				
	parser.add_argument('-o', '--output_file', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step2_results/selected_bar.txt")					
	
	args = parser.parse_args()
	
	select_bar(args.input_bar, args.output_file)


if "__main__" == __name__:
	main()
