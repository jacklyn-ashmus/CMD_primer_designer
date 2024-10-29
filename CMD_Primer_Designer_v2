#!/usr/bin/perl -w
use strict;
use warnings;

#prompt user for input file
print "Please enter input file name.\n";
#assign input file variables
#my $inputfile = @ARGV;
#my $inputfile = "testfile.txt";
my $inputfile = <>;
chomp $inputfile;
#my $FILE = quotemeta $inputfile;

#open input file

#print "Please enter input file name: \n";
#my $input = <>;
#chomp $input;
#my $re = quotemeta $input;
my $FILE;
#open read file, or die
	open ( $FILE, "<", $inputfile )
	or die "Could not open file $FILE $!";




#prompt user for restriction enzyme
print "Please enter the NUMBER corresponding to the FIRST restriction enzyme to be used: \n 1) AluI 2) BamHI 3) BbvCI \n 4) BglII 5) Cla1 6) Dra1 \n 7) EcoP15I 8) EcoRI 9) EcoRII \n 10) EcoRV 11) HaeIII 12) HgaI \n 13) Hhal 14) HindIII 15) HinFI \n 16) KpnI 17) NotI 18) PstI \n 19) PvuII 20) Sau3AI 21) SacI \n 22) SalI23) SmaI 24) SpeI \n 25) SphI 26) StuI27) TaqI \n 28) XbaI 29) XhoI 30) XmaI \n ";
#read restriction enzyme string
my $resenzinputorig1 = <>;
chomp $resenzinputorig1;
my $renz1 = quotemeta $resenzinputorig1;

print "Please enter the NUMBER corresponding to the SECOND restriction enzyme to be used: \n 1) AluI 2) BamHI 3) BbvCI \n 4) BglII 5) Cla1 6) Dra1 \n 7) EcoP15I 8) EcoRI 9) EcoRII \n 10) EcoRV 11) HaeIII 12) HgaI \n 13) Hhal 14) HindIII 15) HinFI \n 16) KpnI 17) NotI 18) PstI \n 19) PvuII 20) Sau3AI 21) SacI \n 22) SalI23) SmaI 24) SpeI \n 25) SphI 26) StuI27) TaqI \n 28) XbaI 29) XhoI 30) XmaI \n ";
#read restriction enzyme string
my $resenzinputorig2 = <>;
chomp $resenzinputorig2;
my $renz2 = quotemeta $resenzinputorig2;


#prompt user for number of AA needed, starting from the end of the sequence

print "Please enter the number of amino acids (on the C-terminus) from the FASTA sequence to be included in the target primer.\n";
my $aanum_orig = <>;
chomp $aanum_orig;
my $aanum = quotemeta $aanum_orig;



#prompt user for Alanine substitution sites

print "Please enter Alanine substitution sites counting from the end of the sequence as '# SPACE # SPACE #...'\n";
my $substitution_string_orig = <>;
chomp $substitution_string_orig;
my $substitution_string = quotemeta $substitution_string_orig;

#prompt user for output file name

print "Please enter the output file name: \n";
#my $output_orig = join($inputfile, "_out");
my $output_orig = <>;
chomp $output_orig;
my $output = quotemeta $output_orig;
#my $output = $output_orig;




my @seqs_orig;
my @seqs_orig_temp;
my $a = 0;
my $b = 0;
my $i = 0;
my $c = 0;
my $xz = 0;
my $xx = 0;
my $xy = 0;
my $yy = 0;
my $u = 0;

my $match_a = "A";
my $match_r = "R";
my $match_n = "N";
my $match_d = "D";
my $match_c = "C";
my $match_q = "Q";
my $match_e = "E";
my $match_g = "G";
my $match_h = "H";
my $match_i = "I";
my $match_l = "L";
my $match_k = "K";
my $match_m = "M";
my $match_f = "F";
my $match_p = "P";
my $match_s = "S";
my $match_t = "T";
my $match_w = "W";
my $match_y = "Y";
my $match_v = "V";






my @seq_name;
my $seqs_number;
my @seq;
my @linetemp2 = <$FILE> ;
my $line = scalar(@linetemp2);
my @linetemp;
my $seqscal;


foreach $line (@linetemp2) {
	$linetemp[$i] = $linetemp2[$i];
	$i = $i + 1;
}

$i = 0;
$a = 0;
$b = 0;
$yy = 0;
$xz = 0;
$xy = 0;
foreach $line (@linetemp) {
	if ($linetemp[$yy] =~ /^\s*([^\w+])/g ) {
		$seq_name[$xz] = $linetemp[$yy];
		$xz = $xz + 1;
	} else { 
		$seqs_orig[$xy] = $linetemp[$yy];
		$xy = $xy + 1;
	}
	$yy = $yy + 1;
}



my $num_of_seq = @seqs_orig;

my $num_of_seq_adjusted = $num_of_seq - 1;



print "\n\n\n\n 108 ::: NUM OF SEQ =  $num_of_seq \n";
#close file

close $FILE or die;

#make array from input for aa to be substituted
my @asubstemp = split('\ ', $substitution_string);




#trim whitespace from Alanine substitution instructions
$i = 0;
foreach (@asubstemp) {
	$asubstemp[$i] =~ s/\\//g;
	$i = $i + 1;
}



	
	


#trim whitespace from Original Strings
for ($i = 0; $i <= 	$num_of_seq; ) {
	$seqs_orig[$i] =~ s/^\s*(.*?)\s*$/$1/;
	$i = $i + 1;
}



$a = 0;
$b = 0;
$i = 0;
my $d = 0;
my $e = 0;
my @temp_seq_split;
my $temp_scal;
my $asubstemp_scal = @asubstemp;
my $seq = undef;
my $doof = undef;
my @temp_rev;
my @seq_rev;
my @seq_rev_t;
my @targetseqrev;
my @temp_targetseq;
my @targetseq;
my @temp_targetseqrev;
my @asubsseq_temp;
my @asubsseq_temprev;
my @asubs;
my @asubs_rev_t;
my @asubs_rev;
my @asubs_t;
my @asubs_t_temp;
print " SEQ Names: :: \n @seq_name \n";
print " ORIG SEQS :: \n @seqs_orig \n\n\n\n";
#print " ASUBSTEMP :: \n @asubstemp \n";
#print " ASUBSTEMPSCAL :: \n $asubstemp_scal \n";


my @a_targetseq_codon_h;
my @targetseq_codon_e;
my @targetseqrev_codon_h;
my @targetseqrev_codon_e;
my @a_targetasubs_t_codon_h;
my @a_targetasubs_t_codon_e;
my @a_targetasubs_trev_codon_h;
my @a_targetasubs_trev_codon_e;
	my @human_f_string;
	my @human_r_string;
	my @ecoli_f_string;
	my @ecoli_r_string;
	my @a_human_f_string;
	my @a_human_r_string;
	my @a_ecoli_f_string;
	my @a_ecoli_r_string;	
my $asub_scal;






my $ecoli_r_string;
my @human_orig_f;
my $stopf = "TGA";
my $stopr = "ACT";
my @human_orig_r_35;

my @human_orig_r_temp;
my @human_orig_r_temp_rev;
my @human_orig_r_53;
my @human_orig_f_gctemp;
my @human_orig_tm;
my @human_orig_r_gctemp;
my $w = 0;
my 	$gc_hof = 0;
my 	$gc_hor = 0;
my 	$gc_haf = 0;
my 	$gc_har = 0;
my 	$gc_eof = 0;
my 	$gc_eor = 0;
my 	$gc_eaf = 0;
my 	$gc_ear = 0;
my 	$scalar_hof = 0;
my 	$scalar_hor = 0;
my 	$scalar_haf = 0;
my 	$scalar_har = 0;
my 	$scalar_eof = 0;
my 	$scalar_eor = 0;
my 	$scalar_eaf = 0;
my 	$scalar_ear = 0;


my $renzseq1_code;
my $renzseq1_comp;
my $renzseq2_code;
my $renzseq2_comp;


my @human_asub_f;
my @ecoli_orig_f_gctemp;
my @final_human_orig_r_53;
my @final_human_orig_gc;
my @final_human_orig_tm;
my @final_human_asub_f;
my @final_human_asub_r_53;
my @final_human_asub_gc;
my @final_human_asub_tm;
my @final_ecoli_orig_f;
my @final_ecoli_orig_r_53;
my @final_ecoli_orig_gc;
my @final_ecoli_orig_r_temp;
my @final_ecoli_asub_f;
my @final_ecoli_asub_r_53;
my @final_ecoli_asub_gc;
my @final_ecoli_asub_tm;

my @human_asub_r_35;
my @human_asub_r_temp;
my @human_asub_r_temp_rev;
my @human_asub_r_53  ;
my @human_asub_f_gctemp ;
my @human_asub_r_gctemp ;
my @human_asub_tm  ;
my @ecoli_orig_f  ;
my @ecoli_orig_r_35  ;
my @ecoli_orig_r_temp ;
my @ecoli_orig_r_temp_rev  ;

my @ecoli_orig_r_gctemp  ;

my @ecoli_asub_f ;
my @ecoli_asub_r_35 ;
my @ecoli_asub_r_temp ;
my @ecoli_asub_r_temp_rev;

my @ecoli_asub_f_gctemp;
my @ecoli_asub_r_gctemp  ;
my @ecoli_asub_gc ;
my @final_human_orig_f;
my @ecoli_asub_tm;





my $human_orig_f_gctemp_scal;
my $human_orig_r_gctemp_scal;
my $human_asub_f_gctemp_scal;
my $human_asub_r_gctemp_scal;
my $ecoli_orig_f_gctemp_scal;
my $ecoli_orig_r_gctemp_scal;
my $ecoli_asub_f_gctemp_scal;
my $ecoli_asub_r_gctemp_scal;



















##################################################################################################################################################################










for  ($i = 0; $i <= $num_of_seq_adjusted;) {
	@temp_seq_split = split( '', $seqs_orig[$i]);
	@temp_rev = reverse(@temp_seq_split);


	
	
	for ($a = 0; $a <= ($aanum-1); ) { 
			$seq_rev_t[$a] = $temp_rev[$a];
			$asubs_rev_t[$a] = $temp_rev[$a];
		$a = $a + 1;
	}


for ($d =0; $d <= ($aanum-1);) {
	$seq_rev[$d] = $seq_rev_t[$d];
	$asubs_rev[$d] = $asubs_rev_t[$d];
	$d = $d + 1;
	}


	for ($e = 0; $e <= ($aanum-1); ) {
		for ($a = 0; $a <= $asubstemp_scal;) {
			if ((scalar($asubstemp[$a])-1) == $e ) {
			$asubs_rev[$e] = "A";
			}
			$a = $a + 1;
		}
	$e = $e + 1;
	}

	
	
	@asubs_t_temp = reverse(@asubs_rev);

	for ($e = 0; $e <= ($aanum-1); ) {
		$asubs_t[$e] = $asubs_t_temp[$e];
	#	print " asubs_t_temp $e = $asubs_t_temp[$e] \n";
	#	print " asubs_t $e = $asubs_t[$e] \n";
		
		$e = $e + 1;
		}
	

	print " ////////// Final values//////////////// \n";


	@seq = reverse(@seq_rev);

	print " Seq number $i :: \n";

	print "tempseqsplit $i :  @temp_seq_split \n";
	print "temp rev  $i : \n @temp_rev \n";
	print "SEQ REV $i :: @seq_rev \n";
	print " SEQ $i ::  @seq \n";
	print "ASUB TEMP REV $i : @asubs_rev \n"; 
	print " ASUB TEMP $i : @asubs_t \n";
	print " ////////////////////////// \n";
		print " ////////////////////////// \n";
			print " ////////////////////////// \n";
				print " ////////////////////////// \n";
	print "\n\n";
	
	
	

	
	
	

	
	
	
	
	shift @seq;
	#shift @asubs_t;
	######################################################################################################################################################
	
	
	$seqscal = @seq;
	$asub_scal = @asubs_t;
	
	for ( $a = 0; $a <= $seqscal; ) {
		#Alanine - A
	
		if ($seq[$a] =~ /$match_a/ ) {
			#Human: GCC
			 $a_targetseq_codon_h[$a] = "GCC";
			#E. coli: GCC
			 $targetseq_codon_e[$a] = "GCA";

			#Coding: GCC, GCA
			#complement: CGG, CGT
			 $targetseqrev_codon_h[$a] = "CGG";
			 $targetseqrev_codon_e[$a] = "CGT";
		}



		#Arginine - R

		if ($seq[$a] =~ /$match_r / ) {
			#Human: AGA
			 $a_targetseq_codon_h[$a] = "AGA";
			#E. coli: AGG
			 $targetseq_codon_e[$a] = "AGG";

			#Coding: AGA, AGG
			#complement: TCT, TCC
			$targetseqrev_codon_h[$a] = "TCT";
			$targetseqrev_codon_e[$a] = "TCC";
		}


	
		#Asparagine - N

		if ($seq[$a] =~ /$match_n/ ) {
			#Human: AAC
			 $a_targetseq_codon_h[$a] = "AAC";
			#E. coli: AAT
			 $targetseq_codon_e[$a] = "AAT";

			#Coding: AAC, AAT
			#complement: TTG, TTA
			 $targetseqrev_codon_h[$a] = "TTG";
			 $targetseqrev_codon_e[$a] = "TTG";
		}



		#Aspartic Acid - D

		if ($seq[$a] =~ /$match_d/ ) {
			#Human: GAC
			 $a_targetseq_codon_h[$a] = "GAC";
			#E. coli: GAT
			 $targetseq_codon_e[$a] = "GAT";

			#Coding: GAC, GAT
			#complement: CTG, CTA
			 $targetseqrev_codon_h[$a] = "CTG";
			 $targetseqrev_codon_e[$a] = "CTA";
		}



		#Cysteine - C
		if ($seq[$a] =~ /$match_c/ ) {
			#Human: TGC
			 $a_targetseq_codon_h[$a] = "TGC";
			#E. coli: TGT
			 $targetseq_codon_e[$a] = "TGT";
			
			#Coding: TGC, TGT
			#complement: ACG, ACA
			 $targetseqrev_codon_h[$a] = "ACG";
			 $targetseqrev_codon_e[$a] = "ACA";
		}
			

		
		#Glutamate - E

		if ($seq[$a] =~ /$match_e/ ) {
			#Human: GAG
			 $a_targetseq_codon_h[$a] = "GAG";
			#E. coli: GAA
			 $targetseq_codon_e[$a] = "GAA";

			#Coding: GAG, GAA
			#complement: CTC, CTT
			 $targetseqrev_codon_h[$a] = "CTC";
			 $targetseqrev_codon_e[$a] = "CTT";
		}

	

		#Glutamine - Q

		if ($seq[$a] =~ /$match_q/ ) {
			#Human: CAG
			 $a_targetseq_codon_h[$a] = "CAG";
			#E. coli: CAG
			 $targetseq_codon_e[$a] = "CAG";

			#Coding: CAG, CAG
			#complement: GTC, GTC
			 $targetseqrev_codon_h[$a] = "GTC";
			 $targetseqrev_codon_e[$a] = "GTC";
		}
	
		

		#Glycine - G

		if ($seq[$a] =~ /$match_g/ ) {
			#Human: GGC
			 $a_targetseq_codon_h[$a] = "GGC";
			#E. coli: GGT
			 $targetseq_codon_e[$a] = "GGT";

			#Coding: GGC, GGT
			#complement: CCG, CCA
			 $targetseqrev_codon_h[$a] = "CCG";
			 $targetseqrev_codon_e[$a] = "CCA";
		}



		#Histidine - H

		if ($seq[$a] =~ /$match_h/ ) {
			#Human: CAC
			 $a_targetseq_codon_h[$a] = "CAC";
			#E. coli: CAT
			 $targetseq_codon_e[$a] = "CAT";
			#Coding: CAC, CAT
			#complement: GTG, GTA
			 $targetseqrev_codon_h[$a] = "GTG";
			 $targetseqrev_codon_e[$a] = "GTA";
		}
		

		#Isoleucine - I

		if ($seq[$a] =~ /$match_i/ ) {
			#Human: ATC
			 $a_targetseq_codon_h[$a] = "ATC";
			#E. coli: ATT
			 $targetseq_codon_e[$a] = "ATT";
			#Coding: ATC, ATT
			#complement: TAG, TAA
			 $targetseqrev_codon_h[$a] = "TAG";
			 $targetseqrev_codon_e[$a] = "TAA";
		}

		
		#Leucine - L

		if ($seq[$a] =~ /$match_l/ ) {
			#Human: CTC
			 $a_targetseq_codon_h[$a] = "CTC";
			#E. coli: TTA
			 $targetseq_codon_e[$a] = "TTA";
			#Coding: CTC, TTA
			#complement: GAG, AAT
			 $targetseqrev_codon_h[$a] = "GAG";
			 $targetseqrev_codon_e[$a] = "AAT";
		}
		

		#Lysine - K

		if ($seq[$a] =~ /$match_k/ ) {
			#Human: AAG
			 $a_targetseq_codon_h[$a] = "AAG";
			#E. coli: AAA
			 $targetseq_codon_e[$a] = "AAA";

			#Coding: AAG, AAA
			#complement: TTC, TTT
			 $targetseqrev_codon_h[$a] = "TTC";
			 $targetseqrev_codon_e[$a] = "TTT";
		}


	
		#Methionine - M

		if ($seq[$a] =~ /$match_m/ ) {
			#Human: ATG
			 $a_targetseq_codon_h[$a] = "ATG";
			#E. coli: ATG
			 $targetseq_codon_e[$a] = "ATG";

			#Coding: ATG, ATG
			#complement: TAC, TAC
			 $targetseqrev_codon_h[$a] = "TAC";
			 $targetseqrev_codon_e[$a] = "TAC";
		}
		


		#Phenyalanine - F

		if ($seq[$a] =~ /$match_f/ ) {
			#Human: TTC
			 $a_targetseq_codon_h[$a] = "TTC";
			#E. coli: TTT
			 $targetseq_codon_e[$a] = "TTT";

			#Coding: TTC, TTT
			#complement: AAG, AAA
			 $targetseqrev_codon_h[$a] = "TTC";
			 $targetseqrev_codon_e[$a] = "TTT";
		}

		#Proline - P
		if ($seq[$a] =~ /$match_p/ ) {
			#Human: CCC
			 $a_targetseq_codon_h[$a] = "CCC";
			#E. coli: CCA
			 $targetseq_codon_e[$a] = "CCA";

			#Coding: CCC, CCA
			#complement: GGG, GGT
			 $targetseqrev_codon_h[$a] = "GGG";
			 $targetseqrev_codon_e[$a] = "GGT";
		}

		
		
		#Serine - S

		if ($seq[$a] =~ /$match_s/ ) {
			#Human: AGC
			 $a_targetseq_codon_h[$a] = "AGC";
			#E. coli: TCT
			 $targetseq_codon_e[$a] = "TCT";

			#Coding: AGC, TCT
			#complement: TCG, AGA
			 $targetseqrev_codon_h[$a] = "TCG";
			 $targetseqrev_codon_e[$a] = "AGA";
		}


		
		#Threonine - T

		if ($seq[$a] =~ /$match_t/ ) {
			#Human: ACC
			 $a_targetseq_codon_h[$a] = "ACC";
			#E. coli: ACT
			 $targetseq_codon_e[$a] = "ACT";

			#Coding: ACC, ACT
			#complement: TGG, TGC
			 $targetseqrev_codon_h[$a] = "TGG";
			 $targetseqrev_codon_e[$a] = "TGC";
		}


		
		#Tryptophan - W

		if ($seq[$a] =~ /$match_w/ ) {
			#Human: TGG
			 $a_targetseq_codon_h[$a] = "TGG";
			#E. coli: TGG
			 $targetseq_codon_e[$a] = "TGG";

			#Coding: TGG, TGG
			#complement: ACC, ACC
			 $targetseqrev_codon_h[$a] = "TGG";
			 $targetseqrev_codon_e[$a] = "TGG";
		}

		
		
		#Tyrosine - Y
		if ($seq[$a] =~ /$match_y/ ) {
			#Human: TAC
			 $a_targetseq_codon_h[$a] = "TAC";
			#E. coli: TAC
			 $targetseq_codon_e[$a] = "TAC";
			#Coding: TAC, TAC
			#complement: ATG, ATG
			$targetseqrev_codon_h[$a] = "TAC";
			$targetseqrev_codon_e[$a] = "TAC";
		}


		
		
		
		#Valine - V
		if ($seq[$a] =~ /$match_v/ ) {
			#Human: GTG
			$a_targetseq_codon_h[$a] = "GTG";
			#E. coli: GTT
			$targetseq_codon_e[$a]= "GTT";

			#Coding: GTG, GTT
			#complement: CAC, CAA
			$targetseqrev_codon_h[$a] = "CAC";
			$targetseqrev_codon_e[$a] = "CAA";
		}	
		
		
		$a = $a + 1;
	}
	
	
	print "\n Seq $i - Human Target Seq Forward ::\n";
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$a_targetseq_codon_h[$e] ";
		$e = $e + 1;
	}
	print "\n Seq $i - E coli Target Seq Forward ::\n";
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$targetseq_codon_e[$e] ";
		$e = $e + 1;
	}
	print "\n Seq $i - Human Target Seq Reverse ::\n";
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$targetseqrev_codon_h[$e] ";
		$e = $e + 1;
	}
	print "\n Seq $i - E coli Target Seq Reverse ::\n";
	
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$targetseqrev_codon_e[$e] ";
		$e = $e + 1;
	}
	print "\n";
	
	
	
		for ( $a = 0; $a <= $asub_scal; ) {
		#Alanine - A
		if ($asubs_t[$a] =~ /$match_a/ ) {
			#Human: GCC
			 $a_targetasubs_t_codon_h[$a] = "GCC";
			#E. coli: GCC
			 $a_targetasubs_t_codon_e[$a] = "GCA";

			#Coding: GCC, GCA
			#complement: CGG, CGT
			 $a_targetasubs_trev_codon_h[$a] = "CGG";
			 $a_targetasubs_trev_codon_e[$a] = "CGT";
		}



		#Arginine - R

		if ($asubs_t[$a] =~ /$match_r / ) {
			#Human: AGA
			 $a_targetasubs_t_codon_h[$a] = "AGA";
			#E. coli: AGG
			 $a_targetasubs_t_codon_e[$a] = "AGG";

			#Coding: AGA, AGG
			#complement: TCT, TCC
			$a_targetasubs_trev_codon_h[$a] = "TCT";
			$a_targetasubs_trev_codon_e[$a] = "TCC";
		}


	
		#Asparagine - N

		if ($asubs_t[$a] =~ /$match_n/ ) {
			#Human: AAC
			 $a_targetasubs_t_codon_h[$a] = "AAC";
			#E. coli: AAT
			 $a_targetasubs_t_codon_e[$a] = "AAT";

			#Coding: AAC, AAT
			#complement: TTG, TTA
			 $a_targetasubs_trev_codon_h[$a] = "TTG";
			 $a_targetasubs_trev_codon_e[$a] = "TTG";
		}



		#Aspartic Acid - D

		if ($asubs_t[$a] =~ /$match_d/ ) {
			#Human: GAC
			 $a_targetasubs_t_codon_h[$a] = "GAC";
			#E. coli: GAT
			 $a_targetasubs_t_codon_e[$a] = "GAT";

			#Coding: GAC, GAT
			#complement: CTG, CTA
			 $a_targetasubs_trev_codon_h[$a] = "CTG";
			 $a_targetasubs_trev_codon_e[$a] = "CTA";
		}
		


		#Cysteine - C
		if ($asubs_t[$a] =~ /$match_c/ ) {
			#Human: TGC
			 $a_targetasubs_t_codon_h[$a] = "TGC";
			#E. coli: TGT
			 $a_targetasubs_t_codon_e[$a] = "TGT";
			
			#Coding: TGC, TGT
			#complement: ACG, ACA
			 $a_targetasubs_trev_codon_h[$a] = "ACG";
			 $a_targetasubs_trev_codon_e[$a] = "ACA";
		}
			

		
		#Glutamate - E

		if ($asubs_t[$a] =~ /$match_e/ ) {
			#Human: GAG
			 $a_targetasubs_t_codon_h[$a] = "GAG";
			#E. coli: GAA
			 $a_targetasubs_t_codon_e[$a] = "GAA";

			#Coding: GAG, GAA
			#complement: CTC, CTT
			 $a_targetasubs_trev_codon_h[$a] = "CTC";
			 $a_targetasubs_trev_codon_e[$a] = "CTT";
		}

	

		#Glutamine - Q

		if ($asubs_t[$a] =~ /$match_q/ ) {
			#Human: CAG
			 $a_targetasubs_t_codon_h[$a] = "CAG";
			#E. coli: CAG
			 $a_targetasubs_t_codon_e[$a] = "CAG";

			#Coding: CAG, CAG
			#complement: GTC, GTC
			 $a_targetasubs_trev_codon_h[$a] = "GTC";
			 $a_targetasubs_trev_codon_e[$a] = "GTC";
		}
	
		

		#Glycine - G

		if ($asubs_t[$a] =~ /$match_g/ ) {
			#Human: GGC
			 $a_targetasubs_t_codon_h[$a] = "GGC";
			#E. coli: GGT
			 $a_targetasubs_t_codon_e[$a] = "GGT";

			#Coding: GGC, GGT
			#complement: CCG, CCA
			 $a_targetasubs_trev_codon_h[$a] = "CCG";
			 $a_targetasubs_trev_codon_e[$a] = "CCA";
		}



		#Histidine - H

		if ($asubs_t[$a] =~ /$match_h/ ) {
			#Human: CAC
			 $a_targetasubs_t_codon_h[$a] = "CAC";
			#E. coli: CAT
			 $a_targetasubs_t_codon_e[$a] = "CAT";
			#Coding: CAC, CAT
			#complement: GTG, GTA
			 $a_targetasubs_trev_codon_h[$a] = "GTG";
			 $a_targetasubs_trev_codon_e[$a] = "GTA";
		}
		

		#Isoleucine - I

		if ($asubs_t[$a] =~ /$match_i/ ) {
			#Human: ATC
			 $a_targetasubs_t_codon_h[$a] = "ATC";
			#E. coli: ATT
			 $a_targetasubs_t_codon_e[$a] = "ATT";
			#Coding: ATC, ATT
			#complement: TAG, TAA
			 $a_targetasubs_trev_codon_h[$a] = "TAG";
			 $a_targetasubs_trev_codon_e[$a] = "TAA";
		}

		
		#Leucine - L

		if ($asubs_t[$a] =~ /$match_l/ ) {
			#Human: CTC
			 $a_targetasubs_t_codon_h[$a] = "CTC";
			#E. coli: TTA
			 $a_targetasubs_t_codon_e[$a] = "TTA";
			#Coding: CTC, TTA
			#complement: GAG, AAT
			 $a_targetasubs_trev_codon_h[$a] = "GAG";
			 $a_targetasubs_trev_codon_e[$a] = "AAT";
		}
		

		#Lysine - K

		if ($asubs_t[$a] =~ /$match_k/ ) {
			#Human: AAG
			 $a_targetasubs_t_codon_h[$a] = "AAG";
			#E. coli: AAA
			 $a_targetasubs_t_codon_e[$a] = "AAA";

			#Coding: AAG, AAA
			#complement: TTC, TTT
			 $a_targetasubs_trev_codon_h[$a] = "TTC";
			 $a_targetasubs_trev_codon_e[$a] = "TTT";
		}


	
		#Methionine - M

		if ($asubs_t[$a] =~ /$match_m/ ) {
			#Human: ATG
			 $a_targetasubs_t_codon_h[$a] = "ATG";
			#E. coli: ATG
			 $a_targetasubs_t_codon_e[$a] = "ATG";

			#Coding: ATG, ATG
			#complement: TAC, TAC
			 $a_targetasubs_trev_codon_h[$a] = "TAC";
			 $a_targetasubs_trev_codon_e[$a] = "TAC";
		}
		


		#Phenyalanine - F

		if ($asubs_t[$a] =~ /$match_f/ ) {
			#Human: TTC
			 $a_targetasubs_t_codon_h[$a] = "TTC";
			#E. coli: TTT
			 $a_targetasubs_t_codon_e[$a] = "TTT";

			#Coding: TTC, TTT
			#complement: AAG, AAA
			 $a_targetasubs_trev_codon_h[$a] = "TTC";
			 $a_targetasubs_trev_codon_e[$a] = "TTT";
		}

		#Proline - P
		if ($asubs_t[$a] =~ /$match_p/ ) {
			#Human: CCC
			 $a_targetasubs_t_codon_h[$a] = "CCC";
			#E. coli: CCA
			 $a_targetasubs_t_codon_e[$a] = "CCA";

			#Coding: CCC, CCA
			#complement: GGG, GGT
			 $a_targetasubs_trev_codon_h[$a] = "GGG";
			 $a_targetasubs_trev_codon_e[$a] = "GGT";
		}

		
		
		#Serine - S

		if ($asubs_t[$a] =~ /$match_s/ ) {
			#Human: AGC
			 $a_targetasubs_t_codon_h[$a] = "AGC";
			#E. coli: TCT
			 $a_targetasubs_t_codon_e[$a] = "TCT";

			#Coding: AGC, TCT
			#complement: TCG, AGA
			 $a_targetasubs_trev_codon_h[$a] = "TCG";
			 $a_targetasubs_trev_codon_e[$a] = "AGA";
		}


		
		#Threonine - T

		if ($asubs_t[$a] =~ /$match_t/ ) {
			#Human: ACC
			 $a_targetasubs_t_codon_h[$a] = "ACC";
			#E. coli: ACT
			 $a_targetasubs_t_codon_e[$a] = "ACT";

			#Coding: ACC, ACT
			#complement: TGG, TGC
			 $a_targetasubs_trev_codon_h[$a] = "TGG";
			 $a_targetasubs_trev_codon_e[$a] = "TGC";
		}


		
		#Tryptophan - W

		if ($asubs_t[$a] =~ /$match_w/ ) {
			#Human: TGG
			 $a_targetasubs_t_codon_h[$a] = "TGG";
			#E. coli: TGG
			 $a_targetasubs_t_codon_e[$a] = "TGG";

			#Coding: TGG, TGG
			#complement: ACC, ACC
			 $a_targetasubs_trev_codon_h[$a] = "TGG";
			 $a_targetasubs_trev_codon_e[$a] = "TGG";
		}

		
		
		#Tyrosine - Y
		if ($asubs_t[$a] =~ /$match_y/ ) {
			#Human: TAC
			 $a_targetasubs_t_codon_h[$a] = "TAC";
			#E. coli: TAC
			 $a_targetasubs_t_codon_e[$a] = "TAC";
			#Coding: TAC, TAC
			#complement: ATG, ATG
			$a_targetasubs_trev_codon_h[$a] = "TAC";
			$a_targetasubs_trev_codon_e[$a] = "TAC";
		}


		
		
		
		#Valine - V
		if ($asubs_t[$a] =~ /$match_v/ ) {
			#Human: GTG
			$a_targetasubs_t_codon_h[$a] = "GTG";
			#E. coli: GTT
			$a_targetasubs_t_codon_e[$a]= "GTT";

			#Coding: GTG, GTT
			#complement: CAC, CAA
			$a_targetasubs_trev_codon_h[$a] = "CAC";
			$a_targetasubs_trev_codon_e[$a] = "CAA";
		}	
		
		
		$a = $a + 1;
	}
	
	
	
	################################################################################################################################
#################ENZYME RESTRICTION SITE 1


#AluI 
	#AGCT
	#TCGA
if ($renz1 == 1) {
	$renzseq1_code = "AGCT";
	$renzseq1_comp = "TCGA";
}

#BamHI
	#GGATCC
	#CCTAGG
if ($renz1 == 2) {
	$renzseq1_code = "GGATCC";
	$renzseq1_comp = "CCTAGG";
}

#BbvCI
	#CCTCAGC
	#GGAGTCG

if ($renz1 == 3) {
	$renzseq1_code = "CCTCAGC";
	$renzseq1_comp = "GGAGTCG";
}

#BglII
	#AGATCT
	#TCTAGA	


if ($renz1 == 4) {
	$renzseq1_code = "AGATCT";
	$renzseq1_comp = "TCTAGA";
}

#Cla1
	#ATCGAT
	#TAGCTA


if ($renz1 == 5) {
	$renzseq1_code = "ATCGAT";
	$renzseq1_comp = "TAGCTA";
}

#Dra1
	#TTTAAA
	#AAATTT
if ($renz1 == 6) {
	$renzseq1_code = "TTTAAA";
	$renzseq1_comp = "AAATTTT";
}

#EcoP15I
	#CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN
	#GTCGTCNNNNNNNNNNNNNNNNNNNNNNNNNNN	
if ($renz1 == 7) {
	$renzseq1_code = "CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	$renzseq1_comp = "GTCGTCNNNNNNNNNNNNNNNNNNNNNNNNNNN";
}

#EcoRI
	#GAATTC
	#CTTAAG
if ($renz1 == 8) {
	$renzseq1_code = "GAATTC";
	$renzseq1_comp = "CTTAAG";
}

#EcoRII
	#CCWGG
	#GGWCC
if ($renz1 == 9) {
	$renzseq1_code = "CCWGG";
	$renzseq1_comp = "GGWCC";
}

#EcoRV
	#GATATC
	#CTATAG
if ($renz1 == 10) {
	$renzseq1_code = "GATATC";
	$renzseq1_comp = "CTATAG";
}

#HaeIII 
	#GGCC
	#CCGG
if ($renz1 == 11) {
	$renzseq1_code = "GGCC";
	$renzseq1_comp = "CCGG";
}

#HgaI
	#GACCGC
	#CTGGCG
if ($renz1 == 12) {
	$renzseq1_code = "GACCGC";
	$renzseq1_comp = "CTGGCG";
}

#Hhal
	#GCGC
	#CGCG
if ($renz1 == 13) {
	$renzseq1_code = "GCGC";
	$renzseq1_comp = "CGCG";
}

#HindIII
	#AAGCTT
	#TTCGAA
if ($renz1 == 14) {
	$renzseq1_code = "AAGCTT";
	$renzseq1_comp = "TTCGAA";
}

#HinFI
	#GACTC
	#CTGAG
if ($renz1 == 15) {
	$renzseq1_code = "GACTC";
	$renzseq1_comp = "CTGAG";
}

#KpnI 
	#GGTACC
	#CCATGG
if ($renz1 == 16) {
	$renzseq1_code = "GGTACC";
	$renzseq1_comp = "CCATGG";
}

#NotI
	#GCGGCCGC
	#CGCCGGCG
if ($renz1 == 17) {
	$renzseq1_code = "GCGGCCGC";
	$renzseq1_comp = "CGCCGGCG";
}

#PstI
	#CTGCAG
	#GACGTC
if ($renz1 == 18) {
	$renzseq1_code = "CTGCAG";
	$renzseq1_comp = "GACGTC";
}


#PvuII
	#CAGCTG
	#GTCGAC
if ($renz1 == 19) {
	$renzseq1_code = "CAGCTG";
	$renzseq1_comp = "GTCGAC";
}

#Sau3AI 
	#GATC
	#CTAG
if ($renz1 == 20) {
	$renzseq1_code = "GATC";
	$renzseq1_comp = "CTAG";
}

#SacI
	#GAGCTC
	#TCATGA
if ($renz1 == 21) {
	$renzseq1_code = "GAGCTC";
	$renzseq1_comp = "TCATGA";
}

#SalI
	#GTCGAC
	#CAGCAG
if ($renz1 == 22) {
	$renzseq1_code = "GTCGAC";
	$renzseq1_comp = "CAGCAG";
}

#SmaI 
	#CCCGGG
	#GGGCCC
if ($renz1 == 23) {
	$renzseq1_code = "CCCGGG";
	$renzseq1_comp = "GGGCCC";
}

#SpeI
	#ACTAGT
	#TGATCA
if ($renz1 == 24) {
	$renzseq1_code = "ACTAGT";
	$renzseq1_comp = "TGATCA";
}

#SphI
	#GCATGC
	#CGTACG
if ($renz1 == 25) {
	$renzseq1_code = "GCATGC";
	$renzseq1_comp = "CGTACG";
}

#StuI
	#AGGCCT
	#TCCGGA
if ($renz1 == 26) {
	$renzseq1_code = "AGGCCT";
	$renzseq1_comp = "TCCGGA";
}

#TaqI
	#TCGA
	#AGCT
if ($renz1 == 27) {
	$renzseq1_code = "TCGA";
	$renzseq1_comp = "AGCT";
}

#XbaI
	#TCTAGA
	#AGATCT

if ($renz1 == 28) {
	$renzseq1_code = "TCTAGA";
	$renzseq1_comp = "AGATCT";
}

#XhoI 
	#CTCGAG
	#GAGCTC
if ($renz1 == 29) {
	$renzseq1_code = "CTCGAG";
	$renzseq1_comp = "GAGCTC";
}
	
#XmaI 
	#CCCGGG
	#GGGCCC
if ($renz1 == 30) {
	$renzseq1_code = "CCCGGG";
	$renzseq1_comp = "GGGCCC";
}








################################################################################################################################
#################ENZYME RESTRICTION SITE 2



#AluI 
	#AGCT
	#TCGA
if ($renz2 == 1) {
	$renzseq2_code = "AGCT";
	$renzseq2_comp = "TCGA";
}

#BamHI
	#GGATCC
	#CCTAGG
if ($renz2 == 2) {
	$renzseq2_code = "GGATCC";
	$renzseq2_comp = "CCTAGG";
}
#BbvCI
	#CCTCAGC
	#GGAGTCG

if ($renz2 == 3) {
	$renzseq2_code = "CCTCAGC";
	$renzseq2_comp = "GGAGTCG";
}

#BglII
	#AGATCT
	#TCTAGA	


if ($renz2 == 4) {
	$renzseq2_code = "AGATCT";
	$renzseq2_comp = "TCTAGA";
}

#Cla1
	#ATCGAT
	#TAGCTA


if ($renz2 == 5) {
	$renzseq2_code = "ATCGAT";
	$renzseq2_comp = "TAGCTA";
}

#Dra1
	#TTTAAA
	#AAATTT
if ($renz2 == 6) {
	$renzseq2_code = "TTTAAA";
	$renzseq2_comp = "AAATTTT";
}

#EcoP15I
	#CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN
	#GTCGTCNNNNNNNNNNNNNNNNNNNNNNNNNNN	
if ($renz2 == 7) {
	$renzseq2_code = "CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	$renzseq2_comp = "GTCGTCNNNNNNNNNNNNNNNNNNNNNNNNNNN";
}

#EcoRI
	#GAATTC
	#CTTAAG
if ($renz2 == 8) {
	$renzseq2_code = "GAATTC";
	$renzseq2_comp = "CTTAAG";
}

#EcoRII
	#CCWGG
	#GGWCC
if ($renz2 == 9) {
	$renzseq2_code = "CCWGG";
	$renzseq2_comp = "GGWCC";

}

#EcoRV
	#GATATC
	#CTATAG
if ($renz2 == 10) {
	$renzseq2_code = "GATATC";
	$renzseq2_comp = "CTATAG";
}

#HaeIII 
	#GGCC
	#CCGG
if ($renz2 == 11) {
	$renzseq2_code = "GGCC";
	$renzseq2_comp = "CCGG";

}

#HgaI
	#GACCGC
	#CTGGCG
if ($renz2 == 12) {
	$renzseq2_code = "GACCGC";
	$renzseq2_comp = "CTGGCG";
}

#Hhal
	#GCGC
	#CGCG
if ($renz2 == 13) {
	$renzseq2_code = "GCGC";
	$renzseq2_comp = "CGCG";
}

#HindIII
	#AAGCTT
	#TTCGAA
if ($renz2 == 14) {
	$renzseq2_code = "AAGCTT";
	$renzseq2_comp = "TTCGAA";
}

#HinFI
	#GACTC
	#CTGAG
if ($renz2 == 15) {
	$renzseq2_code = "GACTC";
	$renzseq2_comp = "CTGAG";
}

#KpnI 
	#GGTACC
	#CCATGG
if ($renz2 == 16) {
	$renzseq2_code = "GGTACC";
	$renzseq2_comp = "CCATGG";

}

#NotI
	#GCGGCCGC
	#CGCCGGCG
if ($renz2 == 17) {
	$renzseq2_code = "GCGGCCGC";
	$renzseq2_comp = "CGCCGGCG";
}

#PstI
	#CTGCAG
	#GACGTC
if ($renz2 == 18) {
	$renzseq2_code = "CTGCAG";
	$renzseq2_comp = "GACGTC";
}


#PvuII
	#CAGCTG
	#GTCGAC
if ($renz2 == 19) {
	$renzseq2_code = "CAGCTG";
	$renzseq2_comp = "GTCGAC";
}

#Sau3AI 
	#GATC
	#CTAG
if ($renz2 == 20) {
	$renzseq2_code = "GATC";
	$renzseq2_comp = "CTAG";
}

#SacI
	#GAGCTC
	#TCATGA
if ($renz2 == 21) {
	$renzseq2_code = "GAGCTC";
	$renzseq2_comp = "TCATGA";
}

#SalI
	#GTCGAC
	#CAGCAG
if ($renz2 == 22) {
	$renzseq2_code = "GTCGAC";
	$renzseq2_comp = "CAGCAG";
}

#SmaI 
	#CCCGGG
	#GGGCCC
if ($renz2 == 23) {
	$renzseq2_code = "CCCGGG";
	$renzseq2_comp = "GGGCCC";
}

#SpeI
	#ACTAGT
	#TGATCA
if ($renz2 == 24) {
	$renzseq2_code = "ACTAGT";
	$renzseq2_comp = "TGATCA";
}

#SphI
	#GCATGC
	#CGTACG
if ($renz2 == 25) {
	$renzseq2_code = "GCATGC";
	$renzseq2_comp = "CGTACG";
}

#StuI
	#AGGCCT
	#TCCGGA
if ($renz2 == 26) {
	$renzseq2_code = "AGGCCT";
	$renzseq2_comp = "TCCGGA";
}

#TaqI
	#TCGA
	#AGCT
if ($renz2 == 27) {
	$renzseq2_code = "TCGA";
	$renzseq2_comp = "AGCT";
}

#XbaI
	#TCTAGA
	#AGATCT

if ($renz2 == 28) {
	$renzseq2_code = "TCTAGA";
	$renzseq2_comp = "AGATCT";
}

#XhoI 
	#CTCGAG
	#GAGCTC
if ($renz2 == 29) {
	$renzseq2_code = "CTCGAG";
	$renzseq2_comp = "GAGCTC";
}
	
#XmaI 
	#CCCGGG
	#GGGCCC
if ($renz2 == 30) {
	$renzseq2_code = "CCCGGG";
	$renzseq2_comp = "GGGCCC";
}



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	




	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	print "\n OrigSeq $i : ";
	print "$seqs_orig[$i] \n";
	print "tempseqsplit $i :  @temp_seq_split \n";
	print "temp rev  $i : \n @temp_rev \n";
	print "SEQ REV $i :: @seq_rev \n";
	print " SEQ $i ::  @seq \n";
	print "ASUB TEMP REV $i : @asubs_rev \n"; 
	

	for ($e = 0; $e <= ($aanum-1); ) {
print "	$asubs_t[$e]";
$e = $e + 1;
}
	print " ////////////////////////// \n";
		print " ////////////////////////// \n";
			print " ////////////////////////// \n";
				print " ////////////////////////// \n";
	print "\n\n";
	
	
	
	print "\n Seq $i - ASub Human Target Seq Forward ::\n";
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$a_targetasubs_t_codon_h[$e] ";
		$e = $e + 1;
	}
	print "\n Seq $i - ASub E coli Target Seq Forward ::\n";
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$a_targetasubs_t_codon_e[$e] ";
		$e = $e + 1;
	}
	print "\n Seq $i - ASub Human Target Seq Reverse ::\n";
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$a_targetasubs_trev_codon_h[$e] ";
		$e = $e + 1;
	}
	print "\n Seq $i - ASub E coli Target Seq Reverse ::\n";
	
	for ($e = 0; $e <= ($aanum-1); ) {
		print "$a_targetasubs_trev_codon_e[$e] ";
		$e = $e + 1;
	}
	print "\n";
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	$human_f_string[$i] = join ('', @a_targetseq_codon_h);
	$human_r_string[$i] = join ('', @targetseqrev_codon_h);
	$ecoli_f_string[$i] = join ('', @targetseq_codon_e);
	$ecoli_r_string[$i] = join ('', @targetseqrev_codon_e);

	
	$a_human_f_string[$i] = join ('', @a_targetasubs_t_codon_h);
	$a_human_r_string[$i] = join ('', @a_targetasubs_trev_codon_h);
	$a_ecoli_f_string[$i] = join ('', @a_targetasubs_t_codon_e);
	$a_ecoli_r_string[$i] = join ('', @a_targetasubs_trev_codon_e);	
	
	
	print "\n ResEnz 1 #: $renz1 , ResEnz 1 Code: $renzseq1_code , ResEnz 1 Comp: $renzseq1_comp \n";

print "\n ResEnz 2 #: $renz2 , ResEnz 2 Code: $renzseq2_code , ResEnz 2 Comp: $renzseq2_comp \n";
	

	
	
	print "\n $i Human F: $human_f_string[$i] \n $i Human R: $human_r_string[$i] \n $i Ecoli F:  $ecoli_f_string[$i] \n $i Ecoli R: $ecoli_r_string[$i] \n";
	
	print "\n $i Asub Human F: $a_human_f_string[$i] \n $i Human R: $a_human_r_string[$i] \n $i Ecoli F:  $a_ecoli_f_string[$i] \n $i Ecoli R: $a_ecoli_r_string[$i] \n";
	
	
	
	
	
	
	#human original sequence reverse translation
	# $human_f_string[$i];
	
	#human original sequence reverse translation COMPLEMENT
	#$human_r_string[$i];

	#human alanine substitution sequence reverse translation
	# $a_human_f_string[$i];

	#human alanine substitution sequence reverse translation COMPLEMENT
	# $a_human_r_string[$i];

	#e coli original sequence reverse translation
	# $ecoli_f_string[$i];

	#e coli alanine substitution sequence reverse translation
	#$a_ecoli_f_string[$i];
						


	#e coli original sequence reverse translation COMPLEMENT
	#$ecoli_r_string[$i];   

	#e coli alanine substitution sequence reverse translation COMPLEMENT
     #$a_ecoli_r_string[$i];   


	#human original sequence forward primer
	# resenz cdn cdn cdn cdn cdn cdn cdn cdn cdn STOP rezenz
	$human_orig_f[$i] = join('', $renzseq1_code,$human_f_string[$i],$stopf,$renzseq2_code);	
	
	#human original sequence reverse primer
	$human_orig_r_35[$i] = join('', $renzseq1_comp,$human_r_string[$i],$stopr,$renzseq2_comp);
	
		@human_orig_r_temp = split(//, $human_orig_r_35[$i]);
		@human_orig_r_temp_rev = reverse( @human_orig_r_temp);
		$human_orig_r_53[$i] = join('', @human_orig_r_temp_rev);  
	
	
	
	
	#human original primers GC content, human original primers melting temp
	@human_orig_f_gctemp = split(//, $human_orig_f[$i]);
	
	$human_orig_f_gctemp_scal = @human_orig_f_gctemp;
		$w = 0;
		for ($w = 0; $w <= $human_orig_f_gctemp_scal ;) {
			if ($human_orig_f_gctemp[$w] =~ /$match_g/) {
			
				$gc_hof = $gc_hof + 1;
			}
			
			if ($human_orig_f_gctemp[$w] =~ /$match_c/) {
			
				$gc_hof = $gc_hof + 1;
			}
			$w ++;
		}
		

		$scalar_hof = @human_orig_f_gctemp;
		print " human_orig_f_gctemp $i = @human_orig_f_gctemp \n";
		print " human_orig_f_gctemp_scal $i = $human_orig_f_gctemp_scal \n";
		print "gc_hof $gc_hof \n";
		print "scalar_hof $scalar_hof \n";
		$w = 0;
			
		


	@human_orig_r_gctemp = split(//, $human_orig_r_35[$i]);
		
		$human_orig_r_gctemp_scal = @human_orig_r_gctemp ;
		for ($w = 0; $w <= $human_orig_r_gctemp_scal; ) {
			if ($human_orig_r_gctemp[$w] =~ /$match_g/) {
			
				$gc_hor = $gc_hor + 1;
			}
			
			if ($human_orig_r_gctemp[$w] =~ /$match_c/) {
			
				$gc_hor = $gc_hor + 1;
			}
			$w ++;
		}
		$scalar_hor = @human_orig_f_gctemp;
		print "gc_hor $gc_hor \n";
		print "scalar_hor $scalar_hor \n";
		
		$w = 0;
	
		
		
		#scalar_hof = total number of NTP in human original primer
		#gc_hof = total number of Gs and Cs in human original primer
		
		
		
		print "1863 gc_hof # $i = $gc_hof \n";
		print "1864 gc_hor # $i = $gc_hor \n";
		$final_human_orig_gc[$i] = (((  $gc_hof + $gc_hor  )/2) / ( ($scalar_hor + $scalar_hof ) /2)) * 100 ;
		print " final_human_orig_gc $i  = $final_human_orig_gc[$i] \n ";
		
	
		#Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
		$human_orig_tm[$i] = 81.5 + 16.6 * 1 + (0.41*($final_human_orig_gc[$i])) - (600/(($scalar_hor + $scalar_hof)/2));


	#human alanine substitution forward primer
	$human_asub_f[$i] = join('',$renzseq1_code,$a_human_f_string[$i],$stopf,$renzseq2_code);
	#human alanine substitution reverse primer
	$human_asub_r_35[$i] = join('',$renzseq1_comp,$a_human_r_string[$i],$stopr,$renzseq2_comp);

		@human_asub_r_temp = split(//, $human_asub_r_35[$i]);
		@human_asub_r_temp_rev = reverse( @human_asub_r_temp);
		$human_asub_r_53[$i] = join('', @human_asub_r_temp_rev);  


	#human alanine substitution primers GC content
	@human_asub_f_gctemp = split(//, $human_asub_f[$i]);
		$human_asub_f_gctemp_scal = @human_asub_f_gctemp;
		for ($w = 0; $w <= $human_asub_f_gctemp_scal; ) {
			if ($human_asub_f_gctemp[$w] =~ /$match_g/) {
			
				$gc_haf = $gc_haf + 1;
			}
			
			if ($human_asub_f_gctemp[$w] =~ /$match_c/) {
			
				$gc_haf = $gc_haf + 1;
			}
			$w ++;
		}
		$scalar_haf = @human_asub_f_gctemp;
		
		print "human_asub_f_gctemp_scal $i = $human_asub_f_gctemp_scal \n";
		print "gc_haf $gc_haf \n";
		print "scalar_haf $scalar_haf \n";
		$w = 0;
		
	@human_asub_r_gctemp = split(//, $human_asub_r_35[$i]);
		$human_asub_r_gctemp_scal = @human_asub_r_gctemp;
for ($w = 0; $w <= $human_asub_r_gctemp_scal ;) {
			if ($human_asub_r_gctemp[$w] =~ /$match_g/) {
			
				$gc_har = $gc_har + 1;
			}
			
			if ($human_asub_r_gctemp[$w] =~ /$match_c/) {
			
				$gc_har = $gc_har + 1;
			}
			$w ++;
		}
		$scalar_har = @human_asub_r_gctemp;
		$w = 0;

		
		print "gc_har $i = $gc_har \n";
		print "scalar_har $i = $scalar_har \n";

	#	$final_human_asub_gc[$i] = ((($gc_haf + $gc_har )/2 )/($scalar_har + $scalar_haf)/2)*100 ;
		$final_human_asub_gc[$i] = (((  $gc_haf + $gc_har  )/2) / ( ($scalar_har + $scalar_haf ) /2)) * 100 ;
		
	#human alanine substitution primers melting temp
		#Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
		$human_asub_tm[$i] = 81.5 + 16.6*1 + (0.41*($final_human_asub_gc[$i])) - (600/(($scalar_har + $scalar_haf)/2));

	

	
	#e coli original sequence forward primer 
	$ecoli_orig_f[$i] = join('',$renzseq1_code,$ecoli_f_string[$i],$stopf,$renzseq2_code);
	
	
	#e coli original sequence reverse primer
	$ecoli_orig_r_35[$i] = join('',$renzseq1_comp,$ecoli_r_string[$i],$stopf,$renzseq2_comp);
	
		@ecoli_orig_r_temp = split(//, $ecoli_orig_r_35[$i]);
		@ecoli_orig_r_temp_rev = reverse( @ecoli_orig_r_temp);
		$final_ecoli_orig_r_53[$i] = join('', @ecoli_orig_r_temp_rev);  

	
	#e coli original primers GC content
	
	@ecoli_orig_f_gctemp = split(//, $ecoli_orig_f[$i]);
		$ecoli_orig_f_gctemp_scal = @ecoli_orig_f_gctemp;
for ($w = 0; $w <= $ecoli_orig_f_gctemp_scal ; ) {
			if ($ecoli_orig_f_gctemp[$w] =~ /$match_g/) {
			
				$gc_eof = $gc_eof + 1;
			}
			
			if ($ecoli_orig_f_gctemp[$w] =~ /$match_c/) {
			
				$gc_eof = $gc_eof + 1;
			}
			$w ++;
		}
		
		$scalar_eof = @ecoli_orig_f_gctemp;
		$w = 0;

		
		
	@ecoli_orig_r_gctemp = split(//, $ecoli_orig_r_35[$i]);
		$ecoli_orig_r_gctemp_scal = @ecoli_orig_r_gctemp;
for ($w = 0; $w <= $ecoli_orig_r_gctemp_scal ; ) {
			if ($ecoli_orig_r_gctemp[$w] =~ /$match_g/) {
			
				$gc_eor = $gc_eor + 1;
			}
			
			if ($ecoli_orig_r_gctemp[$w] =~ /$match_c/) {
			
				$gc_eor = $gc_eor + 1;
			}
			$w ++;
		}
		$scalar_eor = @ecoli_orig_f_gctemp;
		$w = 0;

	
	
		# $final_ecoli_orig_gc[$i] = ((($gc_eof + $gc_eor)/2)/($scalar_eor + $scalar_eof)/2)*100;
		
		$final_ecoli_orig_gc[$i] = (((  $gc_eof + $gc_eor  )/2) / ( ($scalar_eor + $scalar_eof ) /2)) * 100 ;
		
	#e coli original primers melting temp	
		#Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
		$ecoli_orig_r_temp[$i] = 81.5 + 16.6 + (0.41*$final_ecoli_orig_gc[$i]) - (600/(($scalar_eor + $scalar_eof)/2));


	#e coli alanine substitution forward primer
	$ecoli_asub_f[$i] = join('',$renzseq1_code,$a_ecoli_f_string[$i],$stopf,$renzseq2_code);
	#e coli alanine substitution reverse primer
	$ecoli_asub_r_35[$i] = join('',$renzseq1_comp,$a_ecoli_r_string[$i],$stopr,$renzseq2_comp);

		@ecoli_asub_r_temp = split(//, $ecoli_asub_r_35[$i]);
		@ecoli_asub_r_temp_rev = reverse( @ecoli_asub_r_temp);
		$final_ecoli_asub_r_53[$i] = join('', @ecoli_asub_r_temp_rev);  
		
		

	#e coli alanine substitution primers GC content
	
	@ecoli_asub_f_gctemp = split(//, $ecoli_asub_f[$i]);
		$ecoli_asub_f_gctemp_scal = @ecoli_asub_f_gctemp;
for ($w = 0; $w <= $ecoli_asub_f_gctemp_scal ;) {
			if ($ecoli_asub_f_gctemp[$w] =~ /$match_g/) {
			
				$gc_eaf = $gc_eaf + 1;
			}
			
			if ($ecoli_asub_f_gctemp[$w] =~ /$match_c/) {
			
				$gc_eaf = $gc_eaf + 1;
			}
			$w ++;
		}
		$scalar_eaf = @ecoli_asub_f_gctemp;
		$w = 0;
		

		
	@ecoli_asub_r_gctemp = split(//, $ecoli_asub_r_35[$i]);
		$ecoli_asub_r_gctemp_scal = @ecoli_asub_r_gctemp;
for ($w = 0; $w <= $ecoli_asub_r_gctemp_scal ; ) {
			if ($final_ecoli_asub_r_53[$w] =~ /$match_g/) {	
				$gc_ear = $gc_ear + 1;
			}
			
			if ($ecoli_asub_r_gctemp[$w] =~ /$match_c/) {

				$gc_ear = $gc_ear + 1;
			}
			$w ++;
		}
		$scalar_ear = @ecoli_asub_f_gctemp;
		$w = 0;
		print " $i ecoli_asub_r_gctemp = @ecoli_asub_r_gctemp \n";
		print " $i ecoli_asub_r_gctemp_scal = $ecoli_asub_r_gctemp_scal \n";
				print "2060 gc_ear $gc_ear \n";
		print "2061 scalar_ear $scalar_ear \n";

	#	$ecoli_asub_gc[$i] = ((($gc_eaf + $gc_ear)/2)/($scalar_ear + $scalar_eaf)/2)*100;
	
	$ecoli_asub_gc[$i] = (((  $gc_eaf + $gc_ear  )/2) / ( ($scalar_ear + $scalar_eaf ) /2)) * 100 ;
	
	
	
	#e coli alanine substitution primers melting temp
		$ecoli_asub_tm[$i] = 81.5 + 16.6*1 + (0.41*($ecoli_asub_gc[$i])) - (600/(($scalar_har + $scalar_haf)/2));
	
			print " 2015 $i a_ecoli_r_string = $a_ecoli_r_string[$i] \n";
		print " 2015 $i renzseq1_comp = $renzseq1_comp \n";
		print " 2015 $i ecoli_asub_r_temp = @ecoli_asub_r_temp \n";
		print " 2015 $i ecoli_asub_r_temp_rev = @ecoli_asub_r_temp_rev \n";
		print " 2015 $i final_ecoli_asub_r_53= $final_ecoli_asub_r_53[$i] \n";
		print "gc_hof $gc_hof \n";
		print "scalar_hof $scalar_hof \n";
		print "gc_hor $gc_hor \n";
		print "scalar_hor $scalar_hor \n";
		print "gc_haf $gc_haf \n";
		print "scalar_haf $scalar_haf \n";
		print "gc_har $gc_har \n";
		print "scalar_har $scalar_har \n";
		print "gc_eof $gc_eof \n";
		print "scalar_eof $scalar_eof \n";
		print "gc_eor $gc_eor \n";
		print "scalar_eor $scalar_eor \n";
		print "gc_eaf $gc_eaf \n";
		print "scalar_eaf $scalar_eaf \n";
		print "gc_ear $gc_ear \n";
		print "scalar_ear $scalar_ear \n";
	
$final_human_orig_f[$i] = $human_orig_f[$i];
$final_human_orig_r_53[$i] = $human_orig_r_53[$i];

$final_human_orig_tm[$i] = $human_orig_tm[$i];
$final_human_asub_f[$i] = $human_asub_f[$i];
$final_human_asub_r_53[$i]= $human_asub_r_53[$i];

$final_human_asub_tm[$i] = $human_asub_tm[$i];
$final_ecoli_orig_f[$i] = $ecoli_orig_f[$i];
$final_ecoli_orig_r_temp[$i] = $ecoli_orig_r_temp[$i];
$final_ecoli_asub_f[$i] = $ecoli_asub_f[$i];
$final_ecoli_asub_gc[$i] = $ecoli_asub_gc[$i];
$final_ecoli_asub_tm[$i] = $ecoli_asub_tm[$i];

	

	
	
	
	
	
	
	
	
	$gc_hof = 0;
	$gc_hor = 0;
	$gc_haf = 0;
	$gc_har = 0;
	$gc_eof = 0;
	$gc_eor = 0;
	$gc_eaf = 0;
	$gc_ear = 0;
	$scalar_hof = 0;
	$scalar_hor = 0;
	$scalar_haf = 0;
	$scalar_har = 0;
	$scalar_eof = 0;
	$scalar_eor = 0;
	$scalar_eaf = 0;
	$scalar_ear = 0;


	
	

print " Seq $i : \n final_human_orig_f $final_human_orig_f[$i] \n  final_human_orig_r_53 $final_human_orig_r_53[$i] \n final_human_orig_gc
$final_human_orig_gc[$i]\n final_human_orig_tm
$final_human_orig_tm[$i] \n final_human_asub_f
$final_human_asub_f[$i] \n final_human_asub_r_53
$final_human_asub_r_53[$i]\n final_human_asub_gc
$final_human_asub_gc[$i]\n final_human_asub_tm
$final_human_asub_tm[$i] \n final_ecoli_orig_f
$final_ecoli_orig_f[$i]\n final_ecoli_orig_r_53
$final_ecoli_orig_r_53[$i] \n final_ecoli_orig_gc
$final_ecoli_orig_gc[$i] \n final_ecoli_orig_r_temp
$final_ecoli_orig_r_temp[$i] \n final_ecoli_asub_f
$final_ecoli_asub_f[$i] \n final_ecoli_asub_r_53
$final_ecoli_asub_r_53[$i]\n final_ecoli_asub_gc
$final_ecoli_asub_gc[$i] \n final_ecoli_asub_tm
$final_ecoli_asub_tm[$i] \n ";

	

#open output file, or die
	open ( my $op, ">>", $output )
	or die "Could not open file '$FILE' $!";

	#total number of sequences
	print $op "Total Number of FASTA Sequences = $seqs_number \n";


	#for each seq:

	#sequence name
	print $op "Sequence Name: $seq_name[$i] \n";

	#human original sequence forward primer
	print $op "\t Human Original Sequence Primers: \n";
	print $op "\t\t Forward Primer: 5' - $final_human_orig_f[$i] - 3' \n";
	#human original sequence reverse primer
	print $op "\t\t Reverse Primer: 5' - $final_human_orig_r_53[$i] - 3' \n";
	#human original primers GC content
	print $op "\t\t\t GC%: $final_human_orig_gc[$i] % \n";
	#human original primers melting temp: 
	print $op "\t\t\t Melting Temperature: $final_human_orig_tm[$i] C \n";
	#skip line
	print $op "\n";
	#human alanine substitution forward primers
	print $op "\t Human Alanine Substitution Sequence Primers: \n";
	print $op "\t\t Forward Primer: 5' - $final_human_asub_f[$i] - 3' \n";
	#human alanine substitution reverse primer
	print $op "\t\t Reverse Primer: 5' - $final_human_asub_r_53[$i] - 3' \n";
	#human alanine substitution primers GC content
	print $op "\t\t\t GC%: $final_human_asub_gc[$i] % \n";
	#human alanine substitution primers melting temp
	print $op "\t\t\t Melting Temperature: $final_human_asub_tm[$i] C \n";
	#skip line
	print $op "\n";
	
	#e coli original sequence forward primer 
	print $op "\t Escherichia coli Original Sequence Primers: \n";
	print $op "\t\t Forward Primer: 5' - $final_ecoli_orig_f[$i] - 3' \n";
	#e coli original sequence reverse primer
	print $op "\t\t Reverse Primer: 5' - $final_ecoli_orig_r_53[$i] - 3' \n";
	#e coli original primers GC content
	print $op "\t\t\t GC%: $final_ecoli_orig_gc[$i] % \n";
	#e coli original primers melting temp
	print $op "\t\t\t Melting Temperature: $final_ecoli_orig_r_temp[$i] C \n";
	#skip line
	print $op "\n";
	#e coli alanine substitution forward primer
	print $op "\t Escherichia coli Alanine Substitution Sequence Primers: \n";
	print $op "\t\t Forward Primer: 5' - $final_ecoli_asub_f[$i] - 3' \n";
	#e coli alanine substitution reverse primer
	print $op "\t\t Reverse Primer: 5' - $final_ecoli_asub_r_53[$i] - 3' \n";
	#e coli alanine substitution primers GC content
	print $op "\t\t\t GC%: $final_ecoli_asub_gc[$i] % \n";
	#e coli alanine substitution primers melting temp
	print $op "\t\t\t Melting Temperature: $final_ecoli_asub_tm[$i] C \n";
	#skip 4 lines
	print $op "\n\n\n\n";	
	
	#close output file
close ($op);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	$i = $i + 1;
}



















	








































