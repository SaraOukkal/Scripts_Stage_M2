#!/usr/bin/python3

#Importer packages: 
import argparse 

def select_bar(input_bar, input_interbar,output_file):
	
	out=open(output_file,"w") 
	f1=open(input_bar,"r")
	f2=open(input_interbar, "r")
	
	bar=f1.readlines()
	interbar=f2.readlines()
	
	done_chrom=[] #Crée la liste done_chrom qui permet de suivre l'avancée du programme
	
	dico={}
	
	for l in bar: #Barrières sans alu
			line=l.strip().split('\t')
			chrom=line[0]
			st_br=int(line[1])
			if chrom not in dico.keys():
				dico[chrom]=[]
			dico[chrom].append(st_br)
	
	for a in interbar: 
		line=a.strip().split('\t')
		chrom_br=line[0] #Chromosome
		st_br1=int(line[1]) #Start première barrière
		end_br1=int(line[2]) #End première barrière
		st_br2=int(line[3]) #Start Deuxième barrière
		end_br2=int(line[4]) #End deuxième barrière
		
		if chrom_br not in done_chrom:
			done_chrom.append(chrom_br)
			print(chrom_br)
			
		if st_br1 in dico[chrom_br] and st_br2 in dico[chrom_br]:
			out.write("{}\t{}\t{}\t{}\t{}\n".format(chrom_br,st_br1,end_br1,st_br2,end_br2))

			
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	parser.add_argument('-bar', '--input_bar', type=str, help="/home/saraoukkal/Documents/Stage_M1/Barriers/panTro5/20200228_interSmallNFR.dat")
	parser.add_argument('-interbar', '--input_interbar', type=str, help="/home/saraoukkal/Documents/Stage_M1/Barriers/panTro5/20200228_interSmallNFR.dat")				

	#fichier output:				
	parser.add_argument('-out', '--output_file', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step2_results/selected_inter_bar.txt")					
	
	args = parser.parse_args()
	
	select_bar(args.input_bar, args.input_interbar, args.output_file)


if "__main__" == __name__:
	main()
