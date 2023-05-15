#!/bin/bash -l

#SBATCH -J phyluce_Aptost_all
#SBATCH -D /home/lgnewton/phyluce_Aptost_all
#SBATCH -o /home/lgnewton/phyluce_Aptost_all/out-%j.txt
#SBATCH -e /home/lgnewton/phyluce_Aptost_all/err-%j.txt
#SBATCH -t 10:00:00
#SBATCH -N 16
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1
#SBATCH -p high2
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lgnewton@ucdavis.edu # Email to which notifications will be sent

set -e
set -u


## to run script, use this command ##
sbatch .sh

## to check job, use this command ##
squeue -u lgnewton

#### load modules ####
module load phyluce

#### activate phyluce environment ####
conda activate phyluce

#### spades has already been run ####

########### steps to get phylogeny 
## skip to match contigs to probes ##
phyluce_assembly_match_contigs_to_probes \ 
--contigs /home/lgnewton/Aptost_scaffolds_all \
--probes /home/lgnewton/UCEProbes_blended_Spider2Kv1+Arachnid11Kv1-WPM2020.fasta \
--output /home/lgnewton/phyluce_Aptost_all/probe_match_UCEs_min65 \
--min-coverage 65 \
--min-identity 65 \
--log-path /home/lgnewton/phyluce_Aptost_all/log

#get_match_counts
phyluce_assembly_get_match_counts \
    --locus-db /home/lgnewton/phyluce_Aptost_all/probe_match_UCEs_min65/probe.matches.sqlite \
    --taxon-list-config /home/lgnewton/phyluce_Aptost_all/datasets.conf \
    --taxon-group 'Aptost_IceClade' \
    --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/Aptost_IceClade.conf \
    --incomplete-matrix \
    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log

#get_fastas_from_match_counts
phyluce_assembly_get_fastas_from_match_counts \
     --contigs /home/lgnewton/Aptost_scaffolds_all \
     --locus-db /home/lgnewton/phyluce_Aptost_all/probe_match_UCEs_min65/probe.matches.sqlite \
     --match-count-output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/Aptost_IceClade.conf \
     --incomplete-matrix /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/Aptost_IceClade.incomplete \
     --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/Aptost_IceClade.fasta \
     --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log

# explode the monolithic FASTA by taxon (you can also do by locus)
phyluce_assembly_explode_get_fastas_file \
    --input /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/Aptost_IceClade.fasta \
    --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/exploded-fastas-by-taxon \
    --by-taxon

# get summary stats on the FASTAS
for i in /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/exploded-fastas-by-taxon/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

# individual loci, mafft aligned
#phyluce_align_seqcap_align \
#    --fasta /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/Aptost_IceClade.fasta \
#    --output  /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus/ \
#    --taxa 62 \
#    --aligner mafft \
#    --ambiguous \
#    --no-trim \
#    --cores 8 \
#    --incomplete-matrix
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log

#phyluce_align_get_align_summary_data \
#       --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus/ \
#       --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log \
#       --cores 2


#phyluce_align_remove_locus_name_from_nexus_lines \
#    --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus/ \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean/ \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log \
#    --cores 1

#min 50 percent
#phyluce_align_get_only_loci_with_min_taxa \
#    --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean/ \
#    --taxa 62 \
#    --percent 0.50 \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa/ \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log \
#    --cores 1

### trim alignments with trimal
#phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
#    --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa/ \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa-trimaldefault/ \
#    --input-format nexus \
#    --output-format nexus \
#    --log-path /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/log \
#    --cores 1

#phyluce_align_convert_one_align_to_another \
#    --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa-trimaldefault/ \
#    --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa-trimaldefault-fasta/ \
#    --input-format nexus \
#    --output-format fasta

#phyluce_align_format_nexus_files_for_raxml \
#       --alignments /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa-trimaldefault/ \
#       --output /home/lgnewton/phyluce_Aptost_all/Aptost_IceClade/mafft-nexus-clean-min-50p-taxa-trimaldefault-raxml/ \
#       --charsets

## run IQTree ##

################## steps to get snps




