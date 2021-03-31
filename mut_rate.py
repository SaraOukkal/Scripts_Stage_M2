#!/usr/bin/python3

#Importer packages: 
import argparse 

def mut_rate(bar_mut_count, bar_AB_count, output):
	
	out=open(output,"w")
	out.write("{}\t{}\n".format("dist","mut_rate"))#Ecrit un header au fichier 

	mut_count=open(bar_mut_count, "r")
	AB_count=open(bar_AB_count, "r")
	
	mut=mut_count.readlines()
	AB=AB_count.readlines()
	
	for l in mut: #Parcours le fichier contenant le nombre de substitutions par type et par distance des NIEBs
		if not l.startswith("d"): #Ignore le header du fichier

			line_mut=l.strip().split("\t")
			mut_dist=line_mut[0] #Charge les données des mutations dans les variables associées 
			mut_AT=int(line_mut[1])
			mut_AC=int(line_mut[2])
			mut_AG=int(line_mut[3])
			mut_CT=int(line_mut[4])
			mut_CA=int(line_mut[5])
			mut_CG=int(line_mut[6])
			mut_GT=int(line_mut[7])
			mut_GA=int(line_mut[8])
			mut_GC=int(line_mut[9])
			mut_TA=int(line_mut[10])
			mut_TC=int(line_mut[11])
			mut_TG=int(line_mut[12])
			tot_mut=mut_AT + mut_AC + mut_AG + mut_CT + mut_CA + mut_CG + mut_GT + mut_GA + mut_GC + mut_TA + mut_TC + mut_TG #Somme des mutations
			
			for i in AB: #Parcours le fichier contenant le nombre de bases ancestrales par type et par distance des NIEBs
				if not i.startswith("d"): #Ignore le header du fichier 
					line_AB=i.strip().split("\t")
					AB_dist=line_AB[0] #Charge les données des bases ancestrales dans les variables associées 
					AB_A=int(line_AB[1])
					AB_C=int(line_AB[2])
					AB_G=int(line_AB[3])
					AB_T=int(line_AB[4])
					tot_AB=AB_A + AB_C + AB_G + AB_T #Somme des bases ancestrales
					
					if mut_dist == AB_dist: #A chaque même distance des barrières
						
						MR=[]
						MR.append(mut_dist)
						MR.append(tot_mut/tot_AB) #Calcul du taux de mutation 
						
						#print(MR)
						out.write("{}\t{}\n".format(MR[0],MR[1]))
					
					else: 
						continue
						

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-mut', '--input_mut', type=str, help='Path to chimp mutations count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/mut_count_bar.txt")					
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases count around barriers', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/AB_count_bar.txt")					
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rate.txt")			

	
	args = parser.parse_args()
	
	print("start counting mutation rates")
	mut_rate(args.input_mut, args.input_AB, args.output)
	
	
if "__main__" == __name__:
	main()

