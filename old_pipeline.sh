###Step 1.
#Pool Pydamage-selected contigs into a single file and include filenames (sample names) into contig ID's. Will differ for you most likely.

for i in $(ls ancient_contigs_damage_0_1/*.ancient_filtered.fa)
do 
	k=$(echo $i | sed -e 's/.*\/\(.*\)\.ancient.*$/\1_/g')
	cat $i | sed -e "s/>/>$k/g" | sed -e 's/ /_/g' >> contigs_pool1.fa
done

for i in $(ls ancient-phage/*.pydamage_filtered.fa) 
do 
	k=$(echo $i | sed -e 's/.*\/\(.*\)\.pydamage.*$/\1_/g')
	cat $i | sed -e "s/>/>$k/g" | sed -e 's/ /_/g' >> contigs_pool2.fa
done

cat contigs_pool1.fa contigs_pool2.fa > all.fa

###Step 2.
#Blast contigs against various databases

nohup blastn -query all.fa -db /data/databases/IMG_VR/IMG_VR_2020-10-12_5.1/IMGVR_all_nucleotides.fna -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' -out all_contigs_IMGVR_outfmt6.txt -num_threads 5 &

#In a similar fashion run it against GVD, GPD, MGV...

#Blast against a recent version of crAss-like phages DB
nohup blastn -query all.fa -db /backup/user_data/andreys/crass/final_db_stephen/all_contigs.fa -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' -out all_contigs_newcrasslike_outfmt6.txt -num_threads 5 &

#HMP ref genomes
nohup blastn -query all.fa -db /data/databases/hmp_ref_genomes/all.fna -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' -out all_contigs_hmprefgenomes_outfmt6.txt -num_threads 5 &

#nt
nohup blastn -query all.fa -db /data/databases/blast/nt/nt -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' -out all_contigs_nt_outfmt6.txt -num_threads 10 &

#Unpublished oral virome DB
nohup blastn -query all.fa -db /backup/user_data/andreys/victoria/virsorter_andrey/oral_viral_contigs.fa -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' -out all_contigs_oralvirome_outfmt6.txt -num_threads 10 &

#Optionally: do separate blast against human sequences only to detect human contigs.
#nohup blastn -query all.fa -db /data/databases/blast/nt/nt -taxids 9606 -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' -out all_humandna_outfmt6.txt -num_threads 5 &

###Step 3.
#Run Virsorter
#My original protocol inluded Virsorter 1, which is now outdated
#nohup sudo run-virsorter --fna all.fa --ncpu 4 --db 2 &
#grep 'VIRSorter_' VIRSorter_global-phage-signal.csv | cut -d ',' -f 1 | sed -e 's/_cov_\(.*\)_/\_cov_\1./g' -e 's/-circular//g' -e 's/VIRSorter_//g' > virsorter_positive.ids

#Instead, install and run Virsorter 2.
conda create -n vs2 -c conda-forge -c bioconda "python>=3.6" scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer==3.3 prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click mamba
conda activate vs2
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e .
~/miniconda3/envs/vs2/bin/virsorter setup -d virsorter2_db -j 4
virsorter run -w all_virsorter2_out -i all.fa -j 10
grep '>' all_virsorter2_out/final-viral-combined.fa | cut -d '|' -f 1 | tr -d '>' > virsorter_positive.ids

###Step 4

#Circular contigs. Instead, this information can also be extracted from Virsorter output.
python ~/scripts/VICA-master/find_circular.py -i all.fa
grep '>' all.fa_circular.fna | tr -d '>' > circular_contigs.ids

### Step 5.
#Import the output into R script and identify viral contigs based on BLAST hits (requires 70% coverage of contig with viral blast hits and length of >1000 nt to deem it viral). This needs to be modified to accommodate for additional DBs.

###Step 6.
#Once the list of viral contig IDs is produced, they can be annotated and taxonomically classified

#Prokka
pullseq -i all.fa -n viral_contigs_1000.ids > viral_contigs_1000.fa

#Prokka does not accept long contig names, so we will have to rename them temporarily

awk '/^>/{print ">seq" ++i; next}{print}' < viral_contigs_1000.fa > viral_contigs_1000_renamed.fa
paste <(grep '>' viral_contigs_1000.fa | tr -d '>') <(grep '>' viral_contigs_1000_renamed.fa | tr -d '>') >  viral_contigs_1000_renamed.txt
mv viral_contigs_1000_renamed.fa ~/prokka_wd/palaeo_viral_contigs_1000_renamed.fa
cd ~/prokka_wd/

sudo docker run -v $PWD:/data staphb/prokka:latest prokka --outdir palaeofaeces --prefix prokka_pvogs_crass --gcode 11 --hmms all_vogs.hmm --proteins viral_proteins_plus_crass.fa palaeo_viral_contigs_1000_renamed.fa

#Demovir
nohup prodigal -i viral_contigs_1000.fa -f gff -o prodigal_vir_contigs1000.gff -a translations_vir_contigs1000.faa -p meta >> prodigal_out.txt &
nohup bash /home/andreys/Demovir/demovir.sh translations_vir_contigs1000.faa 1e-10 10 & #NB!: more strict e-value cut-off used


###Step 7.
#vConTACT2

#Download and deploy miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh #Install into $HOME/conda

#Create and activate Conda environment
conda create --name vContact2_<username> python=3
conda activate vContact2_<username>

#Install dependables
wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
conda install -y -c conda-forge hdf5 pytables pypandoc biopython networkx numpy scipy scikit-learn psutil pyparsing
conda install -y -c pandas==0.25.1
conda install -y -c bioconda mcl blast diamond

#Install vContact2 script
wget https://bitbucket.org/MAVERICLab/vcontact2/get/master.tar.gz
tar -xvzf master.tar.gz
cd MAVERICLab-vcontact2-XXXXXXX
pip install .

#Run vConTACT2 with ProkaryoticViralRefSeq85-Merged (however check if a newer version of viral RefSeq will also work, as taxonomy there will be more up-to-date!)
/home/ross/scripts/proteins_csv.sh translations_vir_contigs1000.faa
nohup vcontact --raw-proteins translations_vir_contigs1000.faa --rel-mode Diamond --proteins-fp proteins.csv --db ProkaryoticViralRefSeq85-Merged --pcs-mode MCL --vcs-mode ClusterONE --c1-bin ~/miniconda3/envs/vContact2_<username>/bin/cluster_one-1.0.jar --output-dir vConTACT2_vir_contigs1000_vrs85 &

