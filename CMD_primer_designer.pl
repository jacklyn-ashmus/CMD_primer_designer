#!/usr/bin/perl -w
use strict;
use warnings;

#prompt user for input file
print "Please enter input file name.\n"
#assign input file variables
my $inputfile = <>;
chomp $inputfile;
my $FILE = quotemeta $inputfile;

#open input file

print "Looking for: ";
my $input = <>;
chomp $input;
my $re = quotemeta $input;

#open read file, or die
	open ( $FILE, "<", $FILE )
	or die "Could not open file '$FILE' $!";

#prompt user for restriction enzyme
##MAKE IT MULTIPLE CHOICE????
print "Please enter the NUMBER corresponding to the restriction enzyme to be used: \n 1) AluI 2) BamHI 3) BbvCI \n 4) BglII 5) Cla1 6) Dra1 \n 7) EcoP15I 8) EcoRI 9) EcoRII \n 10) EcoRV 11) HaeIII 12) HgaI \n 13) Hhal 14) HindIII 15) HinFI \n 16) KpnI 17) NotI 18) PstI \n 19) PvuII 20) Sau3AI 21) SacI \n 22) SalI23) SmaI 24) SpeI \n 25) SphI 26) StuI27) TaqI \n 28) XbaI 29) XmaI \n 
";
#read restriction enzyme string
my $resenzinputorig = <>;
chomp $resenzinputorig;
my $renz = quotemeta $resenzinputorig;




#prompt user for number of AA needed, starting from the end of the sequence
print "Please enter the number of amino acids (on the C-terminus) from the FASTA sequence to be included in the target primer.\n"
my $aanum_orig = <>;
chomp $aanum_orig;
my $aanum = quotemeta $aanum_orig;

#prompt user for Alanine substitution sites
print "Please enter Alanine substitution sites counting from the end of the sequence as '# SPACE # SPACE #...'\n"
my $substitution_string_orig = <>;
chomp $substitution_string_orig;
my $substitution_string = quotemeta $substitution_string_orig;
#prompt user for output file name
print "Please enter the output file name."
my $output_orig = <>;
chomp $output_orig;
my $output = quotemeta $output_orig;

#split file by ">"
my @seqs_orig
my $b = 0;
#start loop 
while ( <$FILE> ) 
{
		$seqs_orig[$b] = split( /\>/,$FILE);
		$b ++;
};

#get seq number
my @seq_num;
my $a = 0;
while ( <@seq_orig> )
{
	if ($seq_orig[$a] =~ /\>(\s+)/ ) {
		$seq_num[$a] = $1;
	
	};
	$a++;
}

#get sequence name and sequence proper
$a = 0;
my @seq_name;
my @seq;
while (<@seq_num>) 
{

	if ($seq_num[$a] =~ /\(\s+)\n(\s+)\n/ ) {
		$seq_name[$a] = $1;
		$seq[$a] = $2;
	};
	$a ++;
}

#get reverse list
my @seq_rev = reverse(@seq);

#get specified number of AA
my @targetseqrev;
my @tagetseq;

$a = 0;
foreach $aanum (@seq_rev) 
{
	$targetseqrev[$a] = $seqrev[$a];
	$a ++;
}

@targetseq = reverse(@targetseqrev);





#create single AA array
$a = 0;
$b = 0;
my @targetseq_split;
my @seq_split
while (@targetseq) 
{
	$seq_split[$a] = split( /\s/, $targetseq[$a]);
	( push @{ $targetseq_split[$a] }, $b )
	= $seq_split[$a];
	$a ++;
}


#generate Alanine substitution arrays
###### REVIEW HOW THE aa number IS HANDLED
my $c;
my $f = 0;
my $g;
	#split amino acid reference numbers for those to be substituted
	my @asubs
	while ($substitution_string) {
		$asubs[$f] = split(/\w+/, $substitution_string);
		$f ++;
	}

my @targetseq_asub;
my @targetseqrev_asub;

#asub forward
while (@targetseq_split) {
	while ($targetseq_split[$c] ) {
		if ($g == $asub[$g]) {
			( push @{ $targetseq_asub[$c] }, $g ) = "A";
		}
	$g ++;
	}
$c ++;
}

@targetseqrev_asub = reverse(@targetseq_asub);


#generate array[seq#][#codon]
my @targetseq_codon_e;
my @targetseqrev_codon_e;
my @targetseq_asub_codon_e;
my @targetseqrev_asub_codon_e;

my @targetseq_codon_h;
my @targetseqrev_codon_h;
my @targetseq_asub_codon_h;
my @targetseqrev_asub_codon_h;	



################################################################################################################################
#################get DNA seqs FORWARD for ORIGINAL SEQUENCE
my $d = 0;
my $e = 0;
while ($targetseq_split[0...$a]) 
{
	while ($targetseq_split[$d][0...$b])
	{
		#Alanine - A

		if ($targetseq_split[$d][$e] =~ /\A/ ) {
			#Human: GCC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "GCC";
			#E. coli: GCC
			( push @{ $targetseq_codon_e[$d] } , $e ) = "GCA";
		}
		
		#Arginine - R

		if ($targetseq_split[$d][$e] =~ /\R/ ) {
			#Human: AGA
			( push @{ $targetseq_codon_h[$d] } , $e ) = "AGA";
			#E. coli: AGG
			( push @{ $targetseq_codon_e[$d] } , $e ) = "AGG";
		}
			
		#Asparagine - N

		if ($targetseq_split[$d][$e] =~ /\N/ ) {
			#Human: AAC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "AAC";
			#E. coli: AAT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "AAT";
		}

		#Aspartic Acid - D

		if ($targetseq_split[$d][$e] =~ /\D/ ) {
			#Human: GAC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "GAC";
			#E. coli: GAU
			( push @{ $targetseq_codon_e[$d] } , $e ) = "GAU";
		}
		
		#Cysteine - C

		if ($targetseq_split[$d][$e] =~ /\C/ ) {
			#Human: TGC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "TGC";
			#E. coli: TGT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TGT";
		}
					
		#Glutamate - E

		if ($targetseq_split[$d][$e] =~ /\E/ ) {
			#Human: TGC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "TGC";
			#E. coli: TGT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TGT";
		}

	
		#Glutamine - Q

		if ($targetseq_split[$d][$e] =~ /\Q/ ) {
			#Human: CAG
			( push @{ $targetseq_codon_h[$d] } , $e ) = "CAG";
			#E. coli: CAG
			( push @{ $targetseq_codon_e[$d] } , $e ) = "CAG";
		}
	
		
		#Glycine - G

		if ($targetseq_split[$d][$e] =~ /\G/ ) {
			#Human: GGC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "GGC";
			#E. coli: GGT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "GGT";
		}


		#Histidine - H

		if ($targetseq_split[$d][$e] =~ /\H/ ) {
			#Human: CAC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "CAC";
			#E. coli: CAT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "CAT";
		}
		

		#Isoleucine - I

		if ($targetseq_split[$d][$e] =~ /\I/ ) {
			#Human: ATC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "ATC";
			#E. coli: ATT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "ATT";
		}
		
		
		#Leucine - L

		if ($targetseq_split[$d][$e] =~ /\L/ ) {
			#Human: CTC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "CTC";
			#E. coli: TTA
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TTA";
		}
		

		#Lysine - K

		if ($targetseq_split[$d][$e] =~ /\K/ ) {
			#Human: AAG
			( push @{ $targetseq_codon_h[$d] } , $e ) = "AAG";
			#E. coli: AAA
			( push @{ $targetseq_codon_e[$d] } , $e ) = "AAA";
		}
		

		
		#Methionine - M

		if ($targetseq_split[$d][$e] =~ /\M/ ) {
			#Human: ATG
			( push @{ $targetseq_codon_h[$d] } , $e ) = "ATG";
			#E. coli: ATG
			( push @{ $targetseq_codon_e[$d] } , $e ) = "ATG";
		}
		
		#Phenyalanine - F

		if ($targetseq_split[$d][$e] =~ /\F/ ) {
			#Human: TTC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "TTC";
			#E. coli: TTT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TTT";
		}
		
		#Proline - P

		if ($targetseq_split[$d][$e] =~ /\P/ ) {
			#Human: CCC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "CCC";
			#E. coli: CCA
			( push @{ $targetseq_codon_e[$d] } , $e ) = "CCA";
		}
		
		
		#Serine - S

		if ($targetseq_split[$d][$e] =~ /\S/ ) {
			#Human: AGC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "AGC";
			#E. coli: TCT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TCT";
		}
		
		#Threonine - T

		if ($targetseq_split[$d][$e] =~ /\T/ ) {
			#Human: ACC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "ACC";
			#E. coli: ACT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "ACT";
		}


		
		#Tryptophan - W

		if ($targetseq_split[$d][$e] =~ /\W/ ) {
			#Human: TGG
			( push @{ $targetseq_codon_h[$d] } , $e ) = "TGG";
			#E. coli: TGG
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TGG";
		}
		
		
		#Tyrosine - Y

		if ($targetseq_split[$d][$e] =~ /\Y/ ) {
			#Human: TAC
			( push @{ $targetseq_codon_h[$d] } , $e ) = "TAC";
			#E. coli: TAC
			( push @{ $targetseq_codon_e[$d] } , $e ) = "TAC";
		}


		
		#Valine - V

		if ($targetseq_split[$d][$e] =~ /\V/ ) {
			#Human: GTG
			( push @{ $targetseq_codon_h[$d] } , $e ) = "GTG";
			#E. coli: GTT
			( push @{ $targetseq_codon_e[$d] } , $e ) = "GTT";
		}


	
			
		$e ++;
	}
	$d ++;

}



################################################################################################################################
#################get DNA seqs REVERSE for ORIGINAL SEQUENCE












################################################################################################################################
#################get DNA seqs FORWARD for SUBSTITUTION SEQUENCE













################################################################################################################################
#################get DNA seqs REVERSE for SUBSTITUTION SEQUENCE






################################################################################################################################
#################ENZYME RESTRICTION SITES

1) AluI 2) BamHI 3) BbvCI \n 4) BglII 5) Cla1 6) Dra1 \n 7) EcoP15I 8) EcoRI 9) EcoRII \n 10) EcoRV 11) HaeIII 12) HgaI \n 13) Hhal 14) HindIII 15) HinFI \n 16) KpnI 17) NotI 18) PstI \n 19) PvuII 20) Sau3AI 21) SacI \n 22) SalI23) SmaI 24) SpeI \n 25) SphI 26) StuI27) TaqI \n 28) XbaI 29) XmaI \n 

my $renzseqf;
my $renzseqr;


#AluI 
	#AGCT
	#TCGA
if ($renz = 1) {
	$renzseqf = "AGCT";
	$renzseqr = "TCGA";
}

#BamHI
	#GGATCC
	#CCTAGG
#BbvCI
	#CCTCAGC
	#GGAGTCG
#BglII
	#AGATCT
	#TCTAGA	
#Cla1
	#ATCGAT
	#TAGCTA
#Dra1
	#TTTAAA
	#AAATTTT
#EcoP15I
	#CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN
	#GTCGTCNNNNNNNNNNNNNNNNNNNNNNNNNNN	
#EcoRI
	#GAATTC
	#CTTAAG
#EcoRII
	#CCWGG
	#GGWCC
#EcoRV
	#GATATC
	#CTATAG
if ($renz = 2) {
	$renzseqf = "GGATCC";
	$renzseqr = "CCTAGG";
}

if ($renz = 3) {
	$renzseqf = "CCTCAGC";
	$renzseqr = "GGAGTCG";
}

if ($renz = 4) {
	$renzseqf = "AGATCT";
	$renzseqr = "TCTAGA";
}

if ($renz = 5) {
	$renzseqf = "ATCGAT";
	$renzseqr = "TAGCTA";
}

if ($renz = 6) {
	$renzseqf = "TTTAAA";
	$renzseqr = "AAATTTT";
}

if ($renz = 7) {
	$renzseqf = "CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	$renzseqr = "GTCGTCNNNNNNNNNNNNNNNNNNNNNNNNNNN";
}

if ($renz = 8) {
	$renzseqf = "GAATTC";
	$renzseqr = "CTTAAG";
}

if ($renz = 9) {
	$renzseqf = "CCWGG";
	$renzseqr = "GGWCC";
}

if ($renz = 10) {
	$renzseqf = "GATATC";
	$renzseqr = "CTATAG";
}

if ($renz = 11) {
	$renzseqf = "GGCC";
	$renzseqr = "CCGG";
}
if ($renz = 12) {
	$renzseqf = "GACCGC";
	$renzseqr = "CTGGCG";
}



HaeIII 
	GGCC
	CCGG
HgaI
	GACCGC
	CTGGCG
Hhal
	GCGC
	CGCG
HindIII
	AAGCTT
	TTCGAA
HinFI
	GACTC
	CTGAG
KpnI 
	GGTACC
	CCATGG
NotI
	GCGGCCGC
	CGCCGGCG
PstI
	CTGCAG
	GACGTC
PvuII
	CAGCTG
	GTCGAC
Sau3AI 
	GATC
	CTAG
SacI
	GAGCTC
	TCATGA
SalI
	GTCGAC
	CAGCAG
SmaI 
	CCCGGG
	GGGCCC
SpeI
	ACTAGT
	TGATCA
SphI
	GCATGC
	CGTACG
StuI
	AGGCCT
	TCCGGA
TaqI
	TCGA
	AGCT
XbaI
	TCTAGA
	AGATCT
XmaI 
	CCCGGG
	GGGCCC
	
	



#output
	#total number of sequences
	print "Total Number of FASTA Sequences = scalar(@seqs_orig)\n"
	#for each seq:
	#sequence name
	#human original sequence forward primer
	#human original sequence reverse primer
	#human original primers GC content
	#human original primers melting temp
	#skip line
	
	#human alanine substitution forward primer
	#human alanine substitution reverse primer
	#human alanine substitution primers GC content
	#human alanine substitution primers melting temp
	#skip line
	
	#e coli original sequence forward primer 
	#e coli original sequence reverse primer
	#e coli original primers GC content
	#e coli original primers melting temp
	#skip line
	
	#e coli alanine substitution forward primer
	#e coli alanine substitution reverse primer
	#e coli alanine substitution primers GC content
	#e coli alanine substitution primers melting temp
	#skip 4 lines
	
#close output file