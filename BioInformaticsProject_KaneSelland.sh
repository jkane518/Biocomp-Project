#Emily Selland and John Kane
#the goal of this project is to generate a list of the organisms that can make methane and are resistant to pH from a proteome of 50 organisms
#we organized our directory to that we didn't have to use relative paths. This means moving hmmbuild, hmmsearch, and muscle all to the working directory
#we also moved the reference sequences for each gene into their own directories, and put those in our working directory
#same with the proteome directory (moved into our working directory)

#first step: we have reference sequences for the mcrA genes responsible for methane production and the HSP genes associated with pH resistance
#there are 15 sequences of each mcrA and hsp genes, so we want to make an hmm profile for each with aligned sequences so that we can search our proteome 
#file for any matches to both mcrA and hsp reference sequences

#building the mcrA hmm profile: first the reference sequences are not in directories by gene so first we make those directories and moved the fasta files
#of each gene into those directories (using mkdir and mv), within the reference sequence directory
#First: muscle takes 1 file input to align sequences, so first we need to get all of the reference gene sequences into one file
#usage: cat gene1files >> gene1_all
#perform this for all genes that you will utilize to search from the proteome and then move the all files into one work space, where you also have
#muscle and hmmer, so that the need for too many paths is eliminated
#-------------------------------------------------------------------------------------------------------- script starts here, explanation above
#Script Usage: bash BioInformaticsProject_KaneSelland.sh referencegene1 referencegene2 proteomedirectory
#for the project: bash BioInformaticsProject_KaneSelland.sh mcrAgene_all hsp70gene_all proteomes

./muscle3.8.31_i86linux64 -in $1 -out $1_results
./muscle3.8.31_i86linux64 -in $2 -out $2_results

#this output, called $1_results or $2_results is the fuzzy hmm profile of 15 aligned sequences that all code for the gene in question
#we need to search through the 50 (or however many) proteomes for the first reference gene, in this case mcrA, in order to cull any that do not have the gene
#will use hmmbuild and hmmsearch for this
./hmmbuild $1_build $1_results
./hmmbuild $2_build $2_results

for proteomeX in $3/*
do
echo $proteomeX 
./hmmsearch --tblout $1_results $1_build $proteomeX | grep -E -A1000 "Alignments for each domain" | grep -E -B1000 "Internal pipeline statistics summary" | grep -E "WP_[0-9]{9}" | wc -l  
./hmmsearch --tblout $2_results $2_build $proteomeX | grep -E -A1000 "Alignments for each domain" | grep -E -B1000 "Internal pipeline statistics summary" | grep -E "WP_[0-9]{9}" | wc -l 
done >> GenesCounts_Proteomes.txt 

#the above for loop is doing the following: taking the results of hmm build for both mcrA and hsp70 genes for each proteome in the proteome directory (this could be any number of
#proteomes), and doing an hmmsearch. This hmmsearch output is then being put into a pipeline that is pulling the lines that contain the sequence identifier from the gene (either 
#mcrA or hsp70 or any other gene used in this script that someone wants to look for). Then it is counting the number of genes in question the proteome contains
#the output of this for loop is an organized, but not tabulated, list of the number of mcrA and hsp70 genes in each proteome (or any other gene you wanted to use) and then saved
#in a file called GenesCounts_Proteomes.txt


#then we are making a new file for our table, called GenesCounts_Table.txt
#this file is a table, with a comma as the deliminator between columns, that lists the proteome number, number of first genes desired, then number of second genes desired
echo Proteome Number , $1 Number , $2 Number > GenesCounts_Table.txt
cat GenesCounts_Proteomes.txt | tr "\n" "," | sed -E 's/,[p]{1}/\np/g' | sort -t , -k2n -k3n  >> GenesCounts_Table.txt

#the table will be sorted numerically first by the number of the first desired gene (mcrA in this case) and then the number of second desired gene (hsp70 in this case)
#from this table, you can look at the bottom most results and find the proteomes which have presence and high numbers of the two genes desired
#for the BioComputing Project, we decided a that a good place for this grad student to start in terms of further investigation into these proteomes which are resistant to pH,
#would be to choose the proteomes with >0 (presence of) mcrA genes and >8 hsp70 genes, because we saw a clustering of 8 copies of hsp70 and then a large jump in copy number after

#to do this we used the code below to move all proteomes with more than 8 hsp70 copies (we looked at the table to determine what number we would move)
#this file can be looked at by using nano Candidate_pH_resistant_Proteomes.txt to see our final results

cat GenesCounts_Table.txt | tail -n 8 >> Candidate_pH_resistant_Proteomes.txt

