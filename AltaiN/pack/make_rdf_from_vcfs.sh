#!/bin/sh
#$ -S /bin/sh
#$ -cwd
export PATH=$PATH:$HOME/program:$HOME/program/tabix-0.2.6:$HOME/program/vcftools_0.1.11/bin
export PERL5LIB=$PERL5LIB:$HOME/program/vcftools_0.1.11/lib/perl5/site_perl

# The second to fifth lines above are for Job script of DDBJ SuperComputer

# required programs:
#  VCFtools (vcftools, vcf-merge)
#  Tabix (bgzip, tabix)
#  gt_majority.pl
#  vcf_for_network.pl
#  maf2vcf.pl
#  vcf2rdf.pl
#  uniq_rdf.pl
#  (grep_rdf.pl ... not yet used)

if [ ${#} -ne 5 ]
then
	echo "Usage: ./make_rdf_from_vcfs.sh <basename_of_VCFs> <Chimpanzee_MAF> <chr_name> <hgref_name> <Chimpanzee_name>"
	echo
	echo "Arguments are interpreted as below."
	echo " - 1000 Genomes VCF file = <basename_of_VCFs>.BAM.vcf"
	echo " - DenisovaPinky VCF file = <basename_of_VCFs>.Deni.recode.vcf"
	echo " - AltaiNeanderthal VCF file = <basename_of_VCFs>.Nea.recode.vcf"
	echo " - <chr_name>, <hgref_name>, <Chimpanzee_name> are those written in <Chimpanzee_MAF>"
	exit 1
fi

basename=${1}
pan_maf=${2}
chr_name=${3}
hgref_name=${4}
pan_name=${5}

bgzip -c ${basename}.BAM.vcf > ${basename}.BAM.vcf.gz
tabix ${basename}.BAM.vcf.gz

gt_majority.pl ${basename}.Deni.recode.vcf > ${basename}.Deni.major.vcf
vcftools --vcf ${basename}.Deni.major.vcf --remove-indels --out ${basename}.Deni.nw --recode --recode-INFO-all
bgzip ${basename}.Deni.nw.recode.vcf
tabix ${basename}.Deni.nw.recode.vcf.gz

gt_majority.pl ${basename}.Nea.recode.vcf > ${basename}.Nea.major.vcf
vcftools --vcf ${basename}.Nea.major.vcf --remove-indels --out ${basename}.Nea.nw --recode --recode-INFO-all
bgzip ${basename}.Nea.nw.recode.vcf
tabix ${basename}.Nea.nw.recode.vcf.gz

vcf-merge -R '0|0' ${basename}.BAM.vcf.gz ${basename}.Deni.nw.recode.vcf.gz ${basename}.Nea.nw.recode.vcf.gz > ${basename}.BAMDeniNea.vcf

vcf_for_network.pl site ${basename}.BAMDeniNea.vcf > ${basename}.BAMDeniNea.nw.vcf
vcf2rdf.pl ${basename}.BAMDeniNea.nw.vcf | grep -ve "-" -e "\." | sed -e 's/DenisL/DENISO/' -e 's/AltaiL/NEANDE/' -e '/DenisR/d' -e '/AltaiR/d' > ${basename}.BAMDeniNea.nw.rdf
uniq_rdf.pl ${basename}.BAMDeniNea.nw.rdf
# grep_rdf.pl in --id H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,DENISO,NEANDE ${basename}.BAM.DeniNea.nw.uniq.rdf > ${basename}.BAMDeniNea.nw.2otu.rdf

maf2vcf.pl ${basename}.BAMDeniNea.nw.vcf ${pan_maf} ${chr_name} ${hgref_name} ${pan_name} | vcf2rdf.pl | grep -ve '[-\.]' | sed -e 's/panTrL/ROOT  /' -e 's/DenisL/DENISO/' -e 's/AltaiL/NEANDE/' -e '/DenisR/d' -e '/AltaiR/d' > ${basename}.BAMDeniNea.nw.Chimp.rdf

uniq_rdf.pl ${basename}.BAMDeniNea.nw.Chimp.rdf
# grep_rdf.pl in --id H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,DENISO,NEANDE,ROOT mcph1_LD08all.BAM.DeniNea.nw.Chimp.uniq.rdf > mcph1_LD08all.BAM.DeniNea.nw.Chimp.2otu.rdf
