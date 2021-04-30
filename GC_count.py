
#!/usr/bin/python3

#Importer packages: 
import argparse 

def GC_count(input_AB):
	"""
	Calcule le %GC
	"""
	AB=open(input_AB,"r")
	done_chrom=[] 
	c=0 #compteur de bases
	Actual_GC_count=0
	Ancestral_GC_count=0
	Base_count=0

	for l in AB: 
		if c%1000000 == 0 : 
			print(c) #affiche C toutes les 1.000.000 nt 
			
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		nuc=line[2].upper() #Nucléotide à l'état actuel
			if nuc == "G" or nuc == "C":
				Actual_GC_count+=1
				
		nuc_EA=line[3].upper() #nucléotide à l'état ancestral
			if nuc_EA == "G" or nuc_EA == "C":
				Ancestral_GC_count+=1
		
		Base_count+=1
				
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom) #imprime le nouveau chromosome pris en charge 
				
		c += 1 #mets à jour le compteur de bases totales 
	
	print("Number of bases", Base_count)
	Actual_GC_per=(Actual_GC_count/Base_count)*100
	print("Actual GC%", Actual_GC_per)
	Ancestral_GC_per=(Ancestral_GC_count/Base_count)*100
	print("Ancestral GC%", Ancestral_GC_per)
	
				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:				
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to mutation positions and nuc ancestral state', default ="/home/soukkal/Bureau/Projet/Step3_results/AB_chimp.txt")			

	args = parser.parse_args()
	
	#Lancer fonctions:
	GC_count(args.input_AB)

if "__main__" == __name__:
	main()
