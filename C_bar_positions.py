#!/usr/bin/python3

#Importer packages: 
import argparse 


def barriers_pos(input_inter,input_bar,output_file):
	"""
	Ecrit dans un fichier les intervalles du chimpanzé présentant une barrière et la position des barrières dans ces intervalles
	"""
	out=open(output_file,"w") 
	f=open(input_inter,"r") #Ouvre le fichier contenant les intervalles du chimpanzé
	data=f.readlines()
	br_file=open(input_bar,"r")  #Ouvre le fichier contenant les positions des barrières présentes dans les intervalles du chimpanzé
	done_chrom=[] #Crée une liste qui stock le nom des chromosomes traités
	index=0
	for line in br_file: 
		split_line=line.strip().split('\t')
		chrom_br=split_line[0]
		st_br=int(split_line[1])
		end_br=int(split_line[2])
		
		if chrom_br not in done_chrom:
			done_chrom.append(chrom_br)
			print(chrom_br)
			index=0 #réinitialise l'index pour un nouveau chromosome

		res=calculate_pos(chrom_br,st_br,end_br,data,index)
		new_st_br=int(res[0])
		new_end_br=int(res[1])
		C_st=int(res[2])
		C_end=int(res[3])
		index=res[4]

		out.write("{}\t{}\t{}\t{}\t{}\n".format(chrom_br,C_st,C_end,new_st_br,new_end_br))

def calculate_pos(chrom_br,st_br,end_br,data,index):
	"""
	Calcule les positions des barrières dans les intervalles du chimpanzé
	"""
	
	for i in range(index,len(data),1):		
			split_ln=data[i].strip().split('\t')
			chrom=split_ln[0]
			st=int(split_ln[1])
			end=int(split_ln[2])
			
			if chrom_br == chrom:
				if st_br >= st and end_br <= end:
					index=i
					new_st_br=st_br-st
					new_end_br=end_br-st
					C_st=st
					C_end=end
					#~ print(st)
					#~ print(st_br)
					#~ print(new_st_br)
					#~ print(end)
					#~ print(end_br)					
					break

	return [new_st_br, new_end_br,C_st, C_end,index]
					
								
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	parser.add_argument('-inter', '--input_inter', type=str, help='Path to chimp intervals', default ="/home/soukkal/Bureau/Projet/Step2_results/C_intervals_sorted.bed")					
	parser.add_argument('-bar', '--input_bar', type=str, help='Path to barriers intervals', default ="/home/soukkal/Bureau/Projet/Step2_results/C_barriers_in_inter_sorted.bed")									
	
	#fichier output:				
	parser.add_argument('-o', '--output_file', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step2_results/new_bar_pos.txt")					
	
	args = parser.parse_args()
	
	barriers_pos(args.input_inter, args.input_bar, args.output_file)


if "__main__" == __name__:
	main()
