#!/usr/bin/perl -w counts stray elements in a hash
use strict;
#THIS CODE CALCULATES BASIC SUMMARY STATISTICS FROM ISO-SEQ OUTPUT TABLE.
#Input file must be arrange: #chr	start	end	isoform	strand_gmap	biotype	isolen	N_exc	rp_ovlapbp	sgexonclass	score	freads	preads	gene	transcript	strand_matchannot
my(%countermultiNz,%countermulti, %countermultinovel,%counter, %null, %hash, %counterz, %countermultiN,%counterzmulti,
$file, $line, $ID, $score,$strand,$biotype,$fraction, $isoforms,$add, $difference,$filtered, $mono_protein, $mono_lncRNA, $mono_NA, $mono_pseudo, $mono_lncRNA_RefSeq, $mono_TUCP, $mono_lncRNA_Cabili, $mono_ncRNA, $mono_IG_TR,
$multi_protein, $multi_lncRNA, $multi_NA, $multi_pseudo, $multi_lncRNA_RefSeq, $multi_TUCP, $multi_lncRNA_Cabili, $multi_ncRNA, $multi_IG_TR,
$multiN_protein, $multiN_lncRNA, $multiN_NA, $multiN_pseudo, $multiN_lncRNA_RefSeq, $multiN_TUCP, $multiN_lncRNA_Cabili, $multiN_ncRNA, $multiN_IG_TR,
$exon_count,$exon_class,$isolen,$repeat_overlap, $removed,$multiexons, $monoexons,$novelmultiexons, $structure, $structuremultiN, $removed_monoexons_lncRNA,$removed_multiexons_lncRNA,
$total,$seq,$AA, $variance, $position,   $N, $element,$exon_class_0,$exon_class_2,$exon_class_3,$exon_class_4,$exon_class_5,$exon_class_6, $removed_monoexons_protein,$removed_multiexons_protein,
$zero_protein, $zero_lncRNA, $zero_NA, $zero_pseudo, $zero_RefSeq, $zero_TUCP, $zero_Cabili, $zero_ncRNA, $zero_IGTR, $structuremulti,$variancelnc,
$two_protein, $two_lncRNA, $two_NA, $two_pseudo, $two_RefSeq, $two_TUCP, $two_Cabili, $two_ncRNA, $two_IGTR, %counterlnc, %counterlncz,$structurelnc, $structurelncz,%counterpro, %counterproz,$structurepro, $structureproz,
$three_protein, $three_lncRNA, $three_NA, $three_pseudo, $three_RefSeq, $three_TUCP, $three_Cabili, $three_ncRNA, $three_IGTR,
$five_protein, $five_lncRNA, $five_NA, $five_pseudo, $five_RefSeq, $five_TUCP, $five_Cabili, $five_ncRNA, $five_IGTR, $variancez,
$NA_protein, $NA_lncRNA, $NA_NA, $NA_pseudo, $NA_RefSeq, $NA_TUCP, $NA_Cabili, $NA_ncRNA, $NA_IGTR, $removed_monoexons, $removed_multiexons,
$exon2_protein, $exon2_lncRNA, $exon2_NA, $exon2_pseudo, $exon2_RefSeq, $exon2_TUCP, $exon2_Cabili, $exon2_ncRNA, $exon2_IGTR,
$exon3_protein, $exon3_lncRNA, $exon3_NA, $exon3_pseudo, $exon3_RefSeq, $exon3_TUCP, $exon3_Cabili, $exon3_ncRNA, $exon3_IGTR,
$exon4_protein, $exon4_lncRNA, $exon4_NA, $exon4_pseudo, $exon4_RefSeq, $exon4_TUCP, $exon4_Cabili, $exon4_ncRNA, $exon4_IGTR,
$exon5_protein, $exon5_lncRNA, $exon5_NA, $exon5_pseudo, $exon5_RefSeq, $exon5_TUCP, $exon5_Cabili, $exon5_ncRNA, $exon5_IGTR,
$exon6_protein, $exon6_lncRNA, $exon6_NA, $exon6_pseudo, $exon6_RefSeq, $exon6_TUCP, $exon6_Cabili, $exon6_ncRNA, $exon6_IGTR,
$exon7_protein, $exon7_lncRNA, $exon7_NA, $exon7_pseudo, $exon7_RefSeq, $exon7_TUCP, $exon7_Cabili, $exon7_ncRNA, $exon7_IGTR,
$exon8_protein, $exon8_lncRNA, $exon8_NA, $exon8_pseudo, $exon8_RefSeq, $exon8_TUCP, $exon8_Cabili, $exon8_ncRNA, $exon8_IGTR,
$exon9_protein, $exon9_lncRNA, $exon9_NA, $exon9_pseudo, $exon9_RefSeq, $exon9_TUCP, $exon9_Cabili, $exon9_ncRNA, $exon9_IGTR,
$exon10_protein, $exon10_lncRNA, $exon10_NA, $exon10_pseudo, $exon10_RefSeq, $exon10_TUCP, $exon10_Cabili, $exon10_ncRNA, $exon10_IGTR,
$exonmore_protein, $exonmore_lncRNA, $exonmore_NA, $exonmore_pseudo, $exonmore_RefSeq, $exonmore_TUCP, $exonmore_Cabili, $exonmore_ncRNA, $exonmore_IGTR,
$exonN2_protein, $exonN2_lncRNA, $exonN2_NA, $exonN2_pseudo, $exonN2_RefSeq, $exonN2_TUCP, $exonN2_Cabili, $exonN2_ncRNA, $exonN2_IGTR,
$exonN3_protein, $exonN3_lncRNA, $exonN3_NA, $exonN3_pseudo, $exonN3_RefSeq, $exonN3_TUCP, $exonN3_Cabili, $exonN3_ncRNA, $exonN3_IGTR,
$exonN4_protein, $exonN4_lncRNA, $exonN4_NA, $exonN4_pseudo, $exonN4_RefSeq, $exonN4_TUCP, $exonN4_Cabili, $exonN4_ncRNA, $exonN4_IGTR,
$exonN5_protein, $exonN5_lncRNA, $exonN5_NA, $exonN5_pseudo, $exonN5_RefSeq, $exonN5_TUCP, $exonN5_Cabili, $exonN5_ncRNA, $exonN5_IGTR,
$exonN6_protein, $exonN6_lncRNA, $exonN6_NA, $exonN6_pseudo, $exonN6_RefSeq, $exonN6_TUCP, $exonN6_Cabili, $exonN6_ncRNA, $exonN6_IGTR,
$exonN7_protein, $exonN7_lncRNA, $exonN7_NA, $exonN7_pseudo, $exonN7_RefSeq, $exonN7_TUCP, $exonN7_Cabili, $exonN7_ncRNA, $exonN7_IGTR,
$exonN8_protein, $exonN8_lncRNA, $exonN8_NA, $exonN8_pseudo, $exonN8_RefSeq, $exonN8_TUCP, $exonN8_Cabili, $exonN8_ncRNA, $exonN8_IGTR,
$exonN9_protein, $exonN9_lncRNA, $exonN9_NA, $exonN9_pseudo, $exonN9_RefSeq, $exonN9_TUCP, $exonN9_Cabili, $exonN9_ncRNA, $exonN9_IGTR,
$exonN10_protein, $exonN10_lncRNA, $exonN10_NA, $exonN10_pseudo, $exonN10_RefSeq, $exonN10_TUCP, $exonN10_Cabili, $exonN10_ncRNA, $exonN10_IGTR,
$exonNmore_protein, $exonNmore_lncRNA, $exonNmore_NA, $exonNmore_pseudo, $exonNmore_RefSeq, $exonNmore_TUCP, $exonNmore_Cabili, $exonNmore_ncRNA, $exonNmore_IGTR);

open (INFILE, $ARGV[0]) or die "dead\n";
%counter = ();
$total = 0;
$multi_IG_TR =0; $multi_IG_TR=0;
$multiN_IG_TR =0; $three_ncRNA = 0; $exonN2_RefSeq=0; $exonN2_TUCP=0; $exonNmore_Cabili=0; $exonNmore_TUCP=0; $exonNmore_RefSeq=0; $exonN10_TUCP=0; $exonN10_RefSeq=0; $exonN9_TUCP=0; $exonN9_RefSeq=0; $exonN8_TUCP=0; $exonN8_RefSeq=0; $exonN7_Cabili=0; $exonN7_TUCP=0; $exonN7_RefSeq=0; $exonN6_Cabili=0; $exonN6_TUCP=0; $exonN6_RefSeq=0; $exonN5_Cabili=0; $exonN5_TUCP=0; $exonN5_RefSeq=0; $exonN4_Cabili=0; $exonN4_TUCP=0; $exonN4_RefSeq=0; $exonN3_Cabili=0; $exonN3_TUCP=0; $exonN3_RefSeq=0; $exonN2_ncRNA=0; $exonN2_Cabili=0;
$mono_lncRNA_RefSeq=0; $multi_lncRNA_RefSeq=0; $exon4_TUCP=0; $exon4_Cabili=0; $exonmore_ncRNA=0; $exonmore_Cabili=0; $exonmore_TUCP=0; $exonmore_RefSeq=0; $exon10_TUCP=0; $exon10_RefSeq=0; $exon9_TUCP=0; $exon9_RefSeq=0; $exon8_TUCP=0; $exon8_RefSeq=0; $exon7_ncRNA=0; $exon7_Cabili=0; $exon7_TUCP=0; $exon7_RefSeq=0; $exon6_ncRNA=0; $exon6_Cabili=0; $exon5_RefSeq=0; $exon6_TUCP=0; $exon6_RefSeq=0; $exon5_ncRNA=0; $exon5_Cabili=0; $exon5_TUCP=0; $exon4_ncRNA=0; $multiN_lncRNA_RefSeq=0; $mono_TUCP=0; $multi_TUCP=0; $multiN_TUCP=0; $mono_lncRNA_Cabili=0; $multi_lncRNA_Cabili=0; $multiN_lncRNA_Cabili=0; $mono_ncRNA=0; $multi_ncRNA=0; $multiN_ncRNA=0; $zero_RefSeq=0; $two_RefSeq=0; $three_RefSeq=0; $five_RefSeq=0; $zero_TUCP=0; $two_TUCP=0; $three_TUCP=0; $five_TUCP=0; $zero_Cabili=0; $two_Cabili=0; $three_Cabili=0; $five_Cabili=0; $zero_ncRNA=0; $two_ncRNA=0; $exon2_RefSeq=0; $exon2_TUCP=0; $exon2_Cabili=0; $exon2_ncRNA=0; $exon3_RefSeq=0; $exon3_TUCP=0; $exon3_Cabili=0; $exon3_ncRNA=0; $exon4_RefSeq=0;
$zero_IGTR =0; $exonN3_ncRNA=0; $exonN5_ncRNA=0; $exon9_ncRNA=0; $mono_IG_TR=0;
$two_IGTR =0; $exonN2_IGTR=0; $exonN3_IGTR=0; $exonN4_ncRNA=0; $exonN4_IGTR=0; $exonN5_IGTR=0; $exonN6_ncRNA=0; $exonN6_IGTR=0; $exonN7_ncRNA=0; $exonN7_IGTR=0; $exonN8_ncRNA=0; $exonN8_Cabili=0; $exonN8_IGTR=0; $exonN9_Cabili=0; $exonN9_ncRNA=0; $exonN9_IGTR=0; $exonN10_Cabili=0; $exonN10_ncRNA=0; $exonN10_IGTR=0; $exonNmore_ncRNA=0; $exonNmore_IGTR=0;
$three_IGTR=0; $exon2_IGTR=0; $exon3_IGTR=0; $exon4_IGTR=0; $exon5_IGTR=0; $exon6_IGTR=0; $exon7_IGTR=0; $exon8_Cabili=0; $exon8_ncRNA=0; $exon8_IGTR=0; $exon9_Cabili=0; $exon9_IGTR=0; $exon10_Cabili=0; $exon10_ncRNA=0; $exon10_IGTR=0; $exonmore_IGTR=0;
$NA_protein=0; $NA_lncRNA=0; $NA_pseudo=0; $NA_RefSeq=0; $NA_TUCP=0; $NA_Cabili=0; $five_ncRNA=0; $NA_ncRNA=0; $five_IGTR=0; $NA_IGTR=0;

while (<INFILE>)
{
        $line=$_;
        chomp $line;
        
if ($line =~ m/chr/)
{$total++;

	($isoforms) = 		($line =~ m/\S+\s+\S+\s+\S+\s+PB.(\S+)\.\S+/); #print $isoforms; $die;
	($strand) = 		($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
	($biotype) = 		($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
	($isolen) = 		($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
	($exon_count) = 	($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
	($repeat_overlap) = 	($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
	($exon_class) = 	($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
        ($score) = 		($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);

	$add = $repeat_overlap + 0.00001;
	$fraction = $add/$isolen;
	$difference = $isolen - $repeat_overlap;
	#print $difference; $die;
		
	if ($fraction > .20 and $difference < 500)
	{$removed++;
	 
	 	my $out2 = $ARGV[0] . ".removed.txt";# Print out all removed isoforms
		open (OUTFILE, ">>$out2") or die;
		print OUTFILE "$line\n";
		
		if ($exon_count == 1)
		{$removed_monoexons++;
		  if($biotype=~/^protein_coding$/)	{$removed_monoexons_protein++;}
			if($biotype=~/^lncRNA$/ or $biotype=~/^lncRNA_RefSeq$/ or $biotype=~/^lncRNA_Cabili$/)	{$removed_monoexons_lncRNA++;}
		}
		
		else
		{$removed_multiexons++;
		     if($biotype=~/^protein_coding$/)	{$removed_multiexons_protein++;}
		     if($biotype=~/^lncRNA$/ or $biotype=~/^lncRNA_RefSeq$/ or $biotype=~/^lncRNA_Cabili$/)	{$removed_multiexons_lncRNA++;}
		}
		
		
		
		
		
	 }
	
	else
	{
		if ($exon_count == 1)
		{$monoexons++;
		
		if ($biotype=~/^protein_coding$/){$mono_protein++;}
		if ($biotype=~/^lncRNA$/){$mono_lncRNA++;}
		if ($biotype=~/^NA$/){$mono_NA++;}
		if ($biotype=~/^pseudogene$/){$mono_pseudo++;}
		if ($biotype=~/^lncRNA_RefSeq$/){$mono_lncRNA_RefSeq++;}
		if ($biotype=~/^TUCP$/){$mono_TUCP++;}
		if ($biotype=~/^lncRNA_Cabili$/){$mono_lncRNA_Cabili++;}
		if ($biotype=~/^ncRNA$/){$mono_ncRNA++;}
		if ($biotype=~/^IG_TR$/){$mono_IG_TR++;}
		
		if ($exon_class == 0){$exon_class_0++;}
		if ($exon_class == 2){$exon_class_2++;}
		if ($exon_class == 3){$exon_class_3++;}
		if ($exon_class == 4){$exon_class_4++;}
		if ($exon_class == 5){$exon_class_5++;}
		if ($exon_class == 6){$exon_class_6++;}
		
		
		if (defined($counter{$isoforms}))
		{$counter{$isoforms}++;}
		else {$counter{$isoforms}=1;}
		
		my $out = $ARGV[0] . ".monoexons.txt";# Print out all Mono-exons
		open (OUTFILE, ">>$out") or die;
		print OUTFILE "$line\n"; 

		}
		
		elsif($exon_count > 1)
		{$multiexons++;
		 
		if ($biotype=~/^protein_coding$/){$multi_protein++;}
		if ($biotype=~/^lncRNA$/){$multi_lncRNA++;}
		if ($biotype=~/^NA$/){$multi_NA++;}
		if ($biotype=~/^pseudogene$/){$multi_pseudo++;}
		if ($biotype=~/^lncRNA_RefSeq$/){$multi_lncRNA_RefSeq++;}
		if ($biotype=~/^TUCP$/){$multi_TUCP++;}
		if ($biotype=~/^lncRNA_Cabili$/){$multi_lncRNA_Cabili++;}
		if ($biotype=~/^ncRNA$/){$multi_ncRNA++;}
		if ($biotype=~/^IG_TR$/){$multi_IG_TR++;}
		
		if ($biotype=~/^lncRNA$/ or $biotype=~/^lncRNA_RefSeq$/ or $biotype=~/^lncRNA_Cabili$/)
		{		if (defined($counterlnc{$isoforms}))
				{$counterlnc{$isoforms}++;}
				else {$counterlnc{$isoforms}=1;}
		}
		
			 	my $out3 = $ARGV[0] . ".multiexons.txt";# Print out all multi-exon isoforms
				open (OUTFILE, ">>$out3") or die;
				print OUTFILE "$line\n"; 
		 
			if ($score == 2)
			{$novelmultiexons++;
			 
			 		if ($biotype=~/^protein_coding$/){$multiN_protein++;}
					if ($biotype=~/^lncRNA$/){$multiN_lncRNA++;}
					if ($biotype=~/^NA$/){$multiN_NA++;}
					if ($biotype=~/^pseudogene$/){$multiN_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$multiN_lncRNA_RefSeq++;}
					if ($biotype=~/^TUCP$/){$multiN_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$multiN_lncRNA_Cabili++;}
					if ($biotype=~/^ncRNA$/){$multiN_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$multiN_IG_TR++;}
					
					if($exon_count == 2)
		{
					if ($biotype=~/^protein_coding$/){$exonN2_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN2_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN2_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN2_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN2_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN2_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN2_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN2_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN2_IGTR++;}
		}
		
				if($exon_count == 3)
		{
					if ($biotype=~/^protein_coding$/){$exonN3_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN3_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN3_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN3_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN3_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN3_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN3_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN3_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN3_IGTR++;}
		}

				if($exon_count == 4)
		{
					if ($biotype=~/^protein_coding$/){$exonN4_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN4_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN4_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN4_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN4_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN4_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN4_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN4_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN4_IGTR++;}
		}

				if($exon_count == 5)
		{
					if ($biotype=~/^protein_coding$/){$exonN5_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN5_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN5_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN5_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN5_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN5_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN5_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN5_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN5_IGTR++;}
		}

				if($exon_count == 6)
		{
					if ($biotype=~/^protein_coding$/){$exonN6_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN6_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN6_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN6_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN6_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN6_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN6_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN6_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN6_IGTR++;}
		}

				if($exon_count == 7)
		{
					if ($biotype=~/^protein_coding$/){$exonN7_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN7_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN7_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN7_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN7_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN7_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN7_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN7_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN7_IGTR++;}
		}

				if($exon_count == 8)
		{
					if ($biotype=~/^protein_coding$/){$exonN8_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN8_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN8_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN8_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN8_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN8_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN8_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN8_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN8_IGTR++;}
		}

				if($exon_count == 9)
		{
					if ($biotype=~/^protein_coding$/){$exonN9_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN9_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN9_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN9_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN9_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN9_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN9_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN9_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN9_IGTR++;}
		}

				if($exon_count == 10)
		{
					if ($biotype=~/^protein_coding$/){$exonN10_protein++;}
					if ($biotype=~/^lncRNA$/){$exonN10_lncRNA++;}
					if ($biotype=~/^NA$/){$exonN10_NA++;}
					if ($biotype=~/^pseudogene$/){$exonN10_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonN10_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonN10_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonN10_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonN10_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonN10_IGTR++;}
		}

				if($exon_count > 10)
		{
					if ($biotype=~/^protein_coding$/){$exonNmore_protein++;}
					if ($biotype=~/^lncRNA$/){$exonNmore_lncRNA++;}
					if ($biotype=~/^NA$/){$exonNmore_NA++;}
					if ($biotype=~/^pseudogene$/){$exonNmore_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonNmore_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonNmore_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonNmore_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonNmore_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonNmore_IGTR++;}
		}

		
						if (defined($countermultiN{$isoforms}))
						{$countermultiN{$isoforms}++;}
						else {$countermultiN{$isoforms}=1;}

						my $out4 = $ARGV[0] . ".novel_multiexons.txt";# Print out all multi-exon isoforms
						open (OUTFILE, ">>$out4") or die;
						print OUTFILE "$line\n"; 
					
			 }
			
			if ($score == 0)
			{
					if ($biotype=~/^protein_coding$/){$zero_protein++;}
					if ($biotype=~/^lncRNA$/){$zero_lncRNA++;}
					if ($biotype=~/^NA$/){$zero_NA++;}
					if ($biotype=~/^pseudogene$/){$zero_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$zero_RefSeq++;}
					if ($biotype=~/^TUCP$/){$zero_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$zero_Cabili++;}
					if ($biotype=~/^ncRNA$/){$zero_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$zero_IGTR++;}
			}


			if ($score == 2)
			{
					if ($biotype=~/^protein_coding$/){$two_protein++;}
					if ($biotype=~/^lncRNA$/){$two_lncRNA++;}
					if ($biotype=~/^NA$/){$two_NA++;}
					if ($biotype=~/^pseudogene$/){$two_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$two_RefSeq++;}
					if ($biotype=~/^TUCP$/){$two_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$two_Cabili++;}
					if ($biotype=~/^ncRNA$/){$two_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$two_IGTR++;}
			}
			
			if ($score == 3)
			{
					if ($biotype=~/^protein_coding$/){$three_protein++;}
					if ($biotype=~/^lncRNA$/){$three_lncRNA++;}
					if ($biotype=~/^NA$/){$three_NA++;}
					if ($biotype=~/^pseudogene$/){$three_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$three_RefSeq++;}
					if ($biotype=~/^TUCP$/){$three_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$three_Cabili++;}
					if ($biotype=~/^ncRNA$/){$three_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$three_IGTR++;}
				
			}
			if ($score == 5)
			{
					if ($biotype=~/^protein_coding$/){$five_protein++;}
					if ($biotype=~/^lncRNA$/){$five_lncRNA++;}
					if ($biotype=~/^NA$/){$five_NA++;}
					if ($biotype=~/^pseudogene$/){$five_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$five_RefSeq++;}
					if ($biotype=~/^TUCP$/){$five_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$five_Cabili++;}
					if ($biotype=~/^ncRNA$/){$five_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$five_IGTR++;}
			}

			if ($score =~/^NA$/)
			{
					if ($biotype=~/^protein_coding$/){$NA_protein++;}
					if ($biotype=~/^lncRNA$/){$NA_lncRNA++;}
					if ($biotype=~/^NA$/){$NA_NA++;}
					if ($biotype=~/^pseudogene$/){$NA_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$NA_RefSeq++;}
					if ($biotype=~/^TUCP$/){$NA_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$NA_Cabili++;}
					if ($biotype=~/^ncRNA$/){$NA_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$NA_IGTR++;}
			}
			
		if (defined($countermulti{$isoforms}))
		{$countermulti{$isoforms}++;}
		else {$countermulti{$isoforms}=1;}

		}
		
		
		if($exon_count == 2)
		{
					if ($biotype=~/^protein_coding$/){$exon2_protein++;}
					if ($biotype=~/^lncRNA$/){$exon2_lncRNA++;}
					if ($biotype=~/^NA$/){$exon2_NA++;}
					if ($biotype=~/^pseudogene$/){$exon2_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon2_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon2_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon2_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon2_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon2_IGTR++;}
		}
		
				if($exon_count == 3)
		{
					if ($biotype=~/^protein_coding$/){$exon3_protein++;}
					if ($biotype=~/^lncRNA$/){$exon3_lncRNA++;}
					if ($biotype=~/^NA$/){$exon3_NA++;}
					if ($biotype=~/^pseudogene$/){$exon3_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon3_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon3_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon3_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon3_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon3_IGTR++;}
		}

				if($exon_count == 4)
		{
					if ($biotype=~/^protein_coding$/){$exon4_protein++;}
					if ($biotype=~/^lncRNA$/){$exon4_lncRNA++;}
					if ($biotype=~/^NA$/){$exon4_NA++;}
					if ($biotype=~/^pseudogene$/){$exon4_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon4_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon4_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon4_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon4_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon4_IGTR++;}
		}

				if($exon_count == 5)
		{
					if ($biotype=~/^protein_coding$/){$exon5_protein++;}
					if ($biotype=~/^lncRNA$/){$exon5_lncRNA++;}
					if ($biotype=~/^NA$/){$exon5_NA++;}
					if ($biotype=~/^pseudogene$/){$exon5_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon5_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon5_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon5_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon5_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon5_IGTR++;}
		}

				if($exon_count == 6)
		{
					if ($biotype=~/^protein_coding$/){$exon6_protein++;}
					if ($biotype=~/^lncRNA$/){$exon6_lncRNA++;}
					if ($biotype=~/^NA$/){$exon6_NA++;}
					if ($biotype=~/^pseudogene$/){$exon6_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon6_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon6_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon6_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon6_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon6_IGTR++;}
		}

				if($exon_count == 7)
		{
					if ($biotype=~/^protein_coding$/){$exon7_protein++;}
					if ($biotype=~/^lncRNA$/){$exon7_lncRNA++;}
					if ($biotype=~/^NA$/){$exon7_NA++;}
					if ($biotype=~/^pseudogene$/){$exon7_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon7_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon7_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon7_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon7_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon7_IGTR++;}
		}

				if($exon_count == 8)
		{
					if ($biotype=~/^protein_coding$/){$exon8_protein++;}
					if ($biotype=~/^lncRNA$/){$exon8_lncRNA++;}
					if ($biotype=~/^NA$/){$exon8_NA++;}
					if ($biotype=~/^pseudogene$/){$exon8_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon8_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon8_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon8_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon8_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon8_IGTR++;}
		}

				if($exon_count == 9)
		{
					if ($biotype=~/^protein_coding$/){$exon9_protein++;}
					if ($biotype=~/^lncRNA$/){$exon9_lncRNA++;}
					if ($biotype=~/^NA$/){$exon9_NA++;}
					if ($biotype=~/^pseudogene$/){$exon9_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon9_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon9_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon9_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon9_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon9_IGTR++;}
		}

				if($exon_count == 10)
		{
					if ($biotype=~/^protein_coding$/){$exon10_protein++;}
					if ($biotype=~/^lncRNA$/){$exon10_lncRNA++;}
					if ($biotype=~/^NA$/){$exon10_NA++;}
					if ($biotype=~/^pseudogene$/){$exon10_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exon10_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exon10_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exon10_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exon10_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exon10_IGTR++;}
		}
				if($exon_count > 10)
		{
					if ($biotype=~/^protein_coding$/){$exonmore_protein++;}
					if ($biotype=~/^lncRNA$/){$exonmore_lncRNA++;}
					if ($biotype=~/^NA$/){$exonmore_NA++;}
					if ($biotype=~/^pseudogene$/){$exonmore_pseudo++;}
					if ($biotype=~/^lncRNA_RefSeq$/){$exonmore_RefSeq++;}
					if ($biotype=~/^TUCP$/){$exonmore_TUCP++;}
					if ($biotype=~/^lncRNA_Cabili$/){$exonmore_Cabili++;}
					if ($biotype=~/^ncRNA$/){$exonmore_ncRNA++;}
					if ($biotype=~/^IG_TR$/){$exonmore_IGTR++;}
		}
	}
}
}
$filtered=$total - $removed;
print "Total\t$total\nIsoforms passing filter\t$filtered\t>%20 repeats and <500bp free\nMulti-exons\t$multiexons\nMulti-exons (novel)\t$novelmultiexons\nMono-exons\t$monoexons\n\n";

print "Mono-exons filtered out: $removed_monoexons of which $removed_monoexons_protein are protein-coding and $removed_monoexons_lncRNA are lncRNAs\n";
print "Multi-exons filtered out: $removed_multiexons of which $removed_multiexons_protein are protein-coding and  $removed_multiexons_lncRNA are lncRNAs\n\n";

print "\tMono-Exons\tMulti-exons\tMulti-exons (novel)\n";
print "protein_coding\t$mono_protein\t$multi_protein\t$multiN_protein\n";
print "lncRNA\t$mono_lncRNA\t$multi_lncRNA\t$multiN_lncRNA\n";
print "NA\t$mono_NA\t$multi_NA\t$multiN_NA\n";
print "pseudogene\t$mono_pseudo\t$multi_pseudo\t$multiN_pseudo\n";
print "lncRNA_RefSeq\t$mono_lncRNA_RefSeq\t$multi_lncRNA_RefSeq\t$multiN_lncRNA_RefSeq\n";
print "TUCP\t$mono_TUCP\t$multi_TUCP\t$multiN_TUCP\n";
print "lncRNA_Cabili\t$mono_lncRNA_Cabili\t$multi_lncRNA_Cabili\t$multiN_lncRNA_Cabili\n";
print "ncRNA\t$mono_ncRNA\t$multi_ncRNA\t$multiN_ncRNA\n";
print "IG_TR\t$mono_IG_TR\t$multi_IG_TR\t$multiN_IG_TR\n\n";

print "MatchAnnot\t0\t2\t3\t5\tNA\n";
print "protein_coding\t$zero_protein\t$two_protein\t$three_protein\t$five_protein\t$NA_protein\n";
print "lncRNA\t$zero_lncRNA\t$two_lncRNA\t$three_lncRNA\t$five_lncRNA\t$NA_lncRNA\n";
print "NA\t$zero_NA\t$two_NA\t$three_NA\t$five_NA\t$NA_NA\n";
print "pseudogene\t$zero_pseudo\t$two_pseudo\t$three_pseudo\t$five_pseudo\t$NA_pseudo\n";
print "lncRNA_RefSeq\t$zero_RefSeq\t$two_RefSeq\t$three_RefSeq\t$five_RefSeq\t$NA_RefSeq\n";
print "TUCP\t$zero_TUCP\t$two_TUCP\t$three_TUCP\t$five_TUCP\t$NA_TUCP\n";
print "lncRNA_Cabili\t$zero_Cabili\t$two_Cabili\t$three_Cabili\t$five_Cabili\t$NA_Cabili\n";
print "ncRNA\t$zero_ncRNA\t$two_ncRNA\t$three_ncRNA\t$five_ncRNA\t$NA_ncRNA\n";
print "IG_TR\t$zero_IGTR\t$two_IGTR\t$three_IGTR\t$five_IGTR\t$NA_IGTR\n\n";

print "SgExonClass\tMono-exon\n";
print "0\t$exon_class_0\n";
print "2\t$exon_class_2\n";
print "3\t$exon_class_3\n";
print "4\t$exon_class_4\n";
print "5\t$exon_class_5\n";
print "6\t$exon_class_6\n\n";

print "Exon/Gene\tprotein_coding\tlncRNA\tNA\tpseudogene\tlncRNA_RefSeq\tTUCP\tlncRNA_Cabili\tncRNA\tIG_TR\n";
print "2\t$exon2_protein\t$exon2_lncRNA\t$exon2_NA\t$exon2_pseudo\t$exon2_RefSeq\t$exon2_TUCP\t$exon2_Cabili\t$exon2_ncRNA\t$exon2_IGTR\n";
print "3\t$exon3_protein\t$exon3_lncRNA\t$exon3_NA\t$exon3_pseudo\t$exon3_RefSeq\t$exon3_TUCP\t$exon3_Cabili\t$exon3_ncRNA\t$exon3_IGTR\n";
print "4\t$exon4_protein\t$exon4_lncRNA\t$exon4_NA\t$exon4_pseudo\t$exon4_RefSeq\t$exon4_TUCP\t$exon4_Cabili\t$exon4_ncRNA\t$exon4_IGTR\n";
print "5\t$exon5_protein\t$exon5_lncRNA\t$exon5_NA\t$exon5_pseudo\t$exon5_RefSeq\t$exon5_TUCP\t$exon5_Cabili\t$exon5_ncRNA\t$exon5_IGTR\n";
print "6\t$exon6_protein\t$exon6_lncRNA\t$exon6_NA\t$exon6_pseudo\t$exon6_RefSeq\t$exon6_TUCP\t$exon6_Cabili\t$exon6_ncRNA\t$exon6_IGTR\n";
print "7\t$exon7_protein\t$exon7_lncRNA\t$exon7_NA\t$exon7_pseudo\t$exon7_RefSeq\t$exon7_TUCP\t$exon7_Cabili\t$exon7_ncRNA\t$exon7_IGTR\n";
print "8\t$exon8_protein\t$exon8_lncRNA\t$exon8_NA\t$exon8_pseudo\t$exon8_RefSeq\t$exon8_TUCP\t$exon8_Cabili\t$exon8_ncRNA\t$exon8_IGTR\n";
print "9\t$exon9_protein\t$exon9_lncRNA\t$exon9_NA\t$exon9_pseudo\t$exon9_RefSeq\t$exon9_TUCP\t$exon9_Cabili\t$exon9_ncRNA\t$exon9_IGTR\n";
print "10\t$exon10_protein\t$exon10_lncRNA\t$exon10_NA\t$exon10_pseudo\t$exon10_RefSeq\t$exon10_TUCP\t$exon10_Cabili\t$exon10_ncRNA\t$exon10_IGTR\n";
print ">10\t$exonmore_protein\t$exonmore_lncRNA\t$exonmore_NA\t$exonmore_pseudo\t$exonmore_RefSeq\t$exonmore_TUCP\t$exonmore_Cabili\t$exonmore_ncRNA\t$exonmore_IGTR\n\n";


print "Exon/Gene (novel)\tprotein_coding\tlncRNA\tNA\tpseudogene\tlncRNA_RefSeq\tTUCP\tlncRNA_Cabili\tncRNA\tIG_TR\n";
print "2\t$exonN2_protein\t$exonN2_lncRNA\t$exonN2_NA\t$exonN2_pseudo\t$exonN2_RefSeq\t$exonN2_TUCP\t$exonN2_Cabili\t$exonN2_ncRNA\t$exonN2_IGTR\n";
print "3\t$exonN3_protein\t$exonN3_lncRNA\t$exonN3_NA\t$exonN3_pseudo\t$exonN3_RefSeq\t$exonN3_TUCP\t$exonN3_Cabili\t$exonN3_ncRNA\t$exonN3_IGTR\n";
print "4\t$exonN4_protein\t$exonN4_lncRNA\t$exonN4_NA\t$exonN4_pseudo\t$exonN4_RefSeq\t$exonN4_TUCP\t$exonN4_Cabili\t$exonN4_ncRNA\t$exonN4_IGTR\n";
print "5\t$exonN5_protein\t$exonN5_lncRNA\t$exonN5_NA\t$exonN5_pseudo\t$exonN5_RefSeq\t$exonN5_TUCP\t$exonN5_Cabili\t$exonN5_ncRNA\t$exonN5_IGTR\n";
print "6\t$exonN6_protein\t$exonN6_lncRNA\t$exonN6_NA\t$exonN6_pseudo\t$exonN6_RefSeq\t$exonN6_TUCP\t$exonN6_Cabili\t$exonN6_ncRNA\t$exonN6_IGTR\n";
print "7\t$exonN7_protein\t$exonN7_lncRNA\t$exonN7_NA\t$exonN7_pseudo\t$exonN7_RefSeq\t$exonN7_TUCP\t$exonN7_Cabili\t$exonN7_ncRNA\t$exonN7_IGTR\n";
print "8\t$exonN8_protein\t$exonN8_lncRNA\t$exonN8_NA\t$exonN8_pseudo\t$exonN8_RefSeq\t$exonN8_TUCP\t$exonN8_Cabili\t$exonN8_ncRNA\t$exonN8_IGTR\n";
print "9\t$exonN9_protein\t$exonN9_lncRNA\t$exonN9_NA\t$exonN9_pseudo\t$exonN9_RefSeq\t$exonN9_TUCP\t$exonN9_Cabili\t$exonN9_ncRNA\t$exonN9_IGTR\n";
print "10\t$exonN10_protein\t$exonN10_lncRNA\t$exonN10_NA\t$exonN10_pseudo\t$exonN10_RefSeq\t$exonN10_TUCP\t$exonN10_Cabili\t$exonN10_ncRNA\t$exonN10_IGTR\n";
print ">10\t$exonNmore_protein\t$exonNmore_lncRNA\t$exonNmore_NA\t$exonNmore_pseudo\t$exonNmore_RefSeq\t$exonNmore_TUCP\t$exonNmore_Cabili\t$exonNmore_ncRNA\t$exonNmore_IGTR\n\n";


foreach $ID (keys %counter)
{
	$structure=$counter{$ID};
	if (defined($counterz{$structure}))
	{$counterz{$structure}++;}
	else {$counterz{$structure}=1;}
}

print "Iso/Gene\tMono-exons\n";
foreach $structure (keys %counterz)
{
	if($structure == 1) {print "1\t$counterz{$structure}\n";}
	if($structure == 2) {print "2\t$counterz{$structure}\n";}
	if($structure == 3) {print "3\t$counterz{$structure}\n";}
	if($structure == 4) {print "4\t$counterz{$structure}\n";}
	if($structure == 5) {print "5\t$counterz{$structure}\n";}
	if($structure == 6) {print "6\t$counterz{$structure}\n";}
	if($structure == 7) {print "7\t$counterz{$structure}\n";}
	if($structure == 8) {print "8\t$counterz{$structure}\n";}
	if($structure == 9) {print "9\t$counterz{$structure}\n";}
	if($structure == 10) {print "10\t$counterz{$structure}\n";}
	if($structure > 10) {$variance+=$counterz{$structure}};
} print ">10\t$variance\n\n";


foreach my $IDs (keys %countermulti)
{
	my $structuremulti=$countermulti{$IDs};
	if (defined($counterzmulti{$structuremulti}))
	{$counterzmulti{$structuremulti}++;}
	else {$counterzmulti{$structuremulti}=1;}
}

print "Iso/Gene\tMulti-exons\n";
foreach $structuremulti (keys %counterzmulti)
{
	if($structuremulti == 1) {print "1\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 2) {print "2\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 3) {print "3\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 4) {print "4\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 5) {print "5\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 6) {print "6\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 7) {print "7\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 8) {print "8\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 9) {print "9\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti == 10) {print "10\t$counterzmulti{$structuremulti}\n";}
	if($structuremulti > 10) {$variancez+=$counterzmulti{$structuremulti}};
}
print ">10\t$variancez\n\n";


foreach my $IDz (keys %countermultiN)
{
	 $structuremultiN=$countermultiN{$IDz};
	if (defined($countermultiNz{$structuremultiN}))
	{$countermultiNz{$structuremultiN}++;}
	else {$countermultiNz{$structuremultiN}=1;}
}

print "Iso/Gene\tMulti-exons (novel)\n";
foreach $structuremultiN (keys %countermultiNz)
{
	if($structuremultiN == 1) {print "1\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 2) {print "2\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 3) {print "3\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 4) {print "4\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 5) {print "5\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 6) {print "6\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 7) {print "7\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 8) {print "8\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 9) {print "9\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN == 10) {print "10\t$countermultiNz{$structuremultiN}\n";}
	if($structuremultiN > 10) {$variancez+=$countermultiNz{$structuremultiN}};
}
print ">10\t$variancez\n\n";


foreach my $IDlnc (keys %counterlnc)
{
	$structurelnc=$counterlnc{$IDlnc};
	if (defined($counterlncz{$structurelnc}))
	{$counterlncz{$structurelnc}++;}
	else {$counterlncz{$structurelnc}=1;}
}

print "Iso/Gene\tMulti-exons (lncRNAs)\n";
foreach $structurelncz (keys %counterlncz)
{
	if($structurelncz == 1) {print "1\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 2) {print "2\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 3) {print "3\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 4) {print "4\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 5) {print "5\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 6) {print "6\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 7) {print "7\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 8) {print "8\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 9) {print "9\t$counterlncz{$structurelncz}\n";}
	if($structurelncz == 10) {print "10\t$counterlncz{$structurelncz}\n";}
	if($structurelncz > 10) {$variancelnc+=$counterlncz{$structurelncz}};
}
print ">10\t$variancelnc\n";


