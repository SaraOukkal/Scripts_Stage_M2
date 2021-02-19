
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


def mut_rates_plots(MR, output_directory):
	
	data=open(MR, "r")
	MR=data.readlines()
	types=["MR_AT","MR_AC","MR_AG","MR_CT","MR_CA","MR_CG","MR_GT","MR_GA","MR_GC","MR_TA","MR_TC","MR_TG"]
	liste={}
	print("Loaded files")
	
	MR_dist=[]
	
	for sub in types: 
		liste[sub]=[]
		
	print(liste)

	for l in MR:
		if not l.startswith("d"): #Ignore le header du fichier
			line_MR=l.strip().split("\t")
			num=1
			
			MR_dist.append(float(line_MR[0])) #Charge les données des mutations dans les variables associées
			
			for t in liste: 
				liste[t].append(float(line_MR[num])*100) #Calcule le taux de mutations en pourcentage 
				num+=1
			
	for i in types:
				print(liste[i])
				plt.plot(MR_dist,liste[i])
				plt.title(i)
				plt.xlabel("Distance from NIEBs")
				plt.ylabel("Mutation rate (%)")
				filename= "%s.png" % i
				filepath=os.path.join(output_directory, filename)
				plt.savefig(filepath)
				plt.clf()


def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-MR', '--input_MR', type=str, help='Path to mutation rates per substitutions type', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/chimp_mut_rates.txt")					

	#fichier de sortie: 
	parser.add_argument('-out', '--output_dir', type=str, help='Path to output directory', default ="/home/saraoukkal/Documents/Stage_M1/Chimp_Step4_results/Plots/")			

	
	args = parser.parse_args()
	
	mut_rates_plots(args.input_MR, args.output_dir)
	
	
if "__main__" == __name__:
	main()

