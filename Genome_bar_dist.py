
#!/usr/bin/python3

#Importer packages: 
import argparse 
from collections import Counter


def load_bar(input_bar):
	"""
	Charge le fichier .bed des barrières dans un dictionnaire
	"""
	
	barriers ={} 
	f=open(input_bar,"r")
	bars=f.readlines()	
	
	for l in bars: 
		line=l.strip().split("\t")
	
		chrom=line[0]

		if chrom not in barriers.keys():
			barriers[chrom]=[]
		
		bar={}
		bar["st1"]=int(line[1])
		bar["end1"]=int(line[2])
		bar["st2"]=int(line[3])
		bar["end2"]=int(line[4])
			
		barriers[chrom].append(bar)
		
	
		
	return barriers #dictionnaire de chromosomes, chaque chromosome est une liste contenant des dictionnaires pour chaque barrière, un dictionnaire de barrière contient les positions start et end
	
	
	
	
def base_bar(barriers, input_AB):
	"""
	Charge le fichier des bases non mutées et les place en fonction des barrières
	"""
	dico={}	
	AB=open(input_AB,"r")
	done_chrom=[] 
	c = 0 #compteur de bases

	for l in AB: 
		if c%100000 == 0 : 
			print(c) #affiche C toutes les 100000 nt 
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #nucléotide 
					
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom) #imprime le nouveau chromosome pris en charge 
			index=0 #réinitialise l'index 
			
		for i in range(index,len(barriers[chrom]),1):
			st1=barriers[chrom][i]["st1"]	#start de la barrière x 
			end1=barriers[chrom][i]["end1"]	#end de la barrière x 
			st2=barriers[chrom][i]["st2"]
			end2=barriers[chrom][i]["end2"]
			
			mid_bar1=end1-((end1-st1)/2) #Milieu de première barrière
			mid_bar2=st2+((end2-st2)/2)	#Milieu de deuxième barrière
			mid_inter_bar=end1+((st2-end1)/2) #Milieu de l'interbarrière
			
			if pos < mid_bar1: #Ignorer les bases avant la premièe barrière
				index=i
				break
			
			elif pos >= mid_bar1 and pos <= mid_bar2: #Bases à prendre en compte
				if pos <= mid_inter_bar: #Si autour du bord de barrière 1
					dist=pos-end1
				else: #Si autour du bord de barrière 2
					dist=st2-pos 
				
				if dist not in dico.keys(): #Si cette distance n'a pas encore été croisée on l'ajoute au dictionnaire 
					dico[dist]=Counter()	
				
				dico[dist][nuc]+=1	
				index=i #mets à jour l'index 
				break
				
		#Si la base est après la deuxième barrière on continue de parcourir les barrières
				
		c += 1 #mets à jour le compteur de bases totales 
					
	return dico															


		
						
def write_out(base_count,output_file):	
	"""
	Ecrit dans le fichier output le compte de chaque type de bases par position par rapport aux barrières.
	"""
	
	out=open(output_file,"w")
	
	out.write("{}\t{}\t{}\t{}\t{}\n".format("dist","A","C","G","T"))#Ecrit un header au fichier 
	
	for dist in sorted(base_count.keys()):
		out.write("{}\t{}\t{}\t{}\t{}\n".format(str(dist),str(base_count[dist]["A"]),str(base_count[dist]["C"]),str(base_count[dist]["G"]),str(base_count[dist]["T"])))


				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier des inter barrières: 
	parser.add_argument('-bar', '--input_bar', type=str, help='Path to barriers intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Barriers/hg38/20200228_interSmallNFR.dat")
	##fichier des bases:				
	parser.add_argument('-B', '--input_B', type=str, help='Path to bases and positions', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step4_results/AB_human.txt")			

	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step5_results/Ali_count_bar.txt")			

	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	print("loading barriers file")
	barriers=load_bar(args.input_bar)
	print("counting bases around barriers")
	base_count=base_bar(barriers, args.input_B)
	print("writing counts in the output file")
	write_out(base_count,args.output)
	print("done")
	

if "__main__" == __name__:
	main()
