use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Getopt::Long;

my $haplotypes_file = "";
# my $haplotypes_pairs_file = shift;
my $prefix = "";
my $samplesheet = "";

GetOptions(
    'haplotypes=s'    => \$haplotypes_file,
    'prefix=s' => \$prefix,
    'samplesheet=s'     => \$samplesheet
) or print "Invalid options\n";

my $selected_file_dqa = "DQA.selected.txt";
my $discarded_file_dqa = "DQA.discarded.txt";
my $summary_file_dqa = "summary.DQA.matrix.txt";
my $selected_file_dqb = "DQB.selected.txt";
my $discarded_file_dqb = "DQB.discarded.txt";
my $summary_file_dqb = "summary.DQB.matrix.txt";
my $selected_file_drb3 = "DRB3.selected.txt";
# my $discarded_file_drb3 = "DRB3.discarded.txt";
my $summary_file_drb3 = "summary.DRB3.matrix.txt";



my %dqa_counts = ();
my %dqb_counts = ();
my %drb3_counts = ();
my %dqa_counts_discarded = ();
my %dqb_counts_discarded = ();
my %drb3_counts_discarded = ();

my $line = 0;
open(P1, "$summary_file_dqa") or die "Cannot open $summary_file_dqa\n";
while(<P1>){
	chomp $_;
	$line++;
	if ($line > 1){
		@words = split("\t", $_);
		$dqa_counts{$words[0]} = $words[27] + $words[28];
		$dqa_counts_discarded{$words[0]} = $words[29];
	}
}
$line = 0;
open(P2, "$summary_file_dqb") or die "Cannot open $summary_file_dqb\n";
while(<P2>){
	chomp $_;
	$line++;
	if ($line > 1){
		@words = split("\t", $_);
		$dqb_counts{$words[0]} = $words[27] + $words[28];
		$dqb_counts_discarded{$words[0]} = $words[29];
	}
}
$line = 0;
open(P3, "$summary_file_drb3") or die "Cannot open $summary_file_drb3\n";
while(<P3>){
	chomp $_;
	$line++;
	if ($line > 1){
		@words = split("\t", $_);
		$drb3_counts{$words[0]} = $words[27] + $words[28];
		$drb3_counts_discarded{$words[0]} = $words[29];
	}
}

open (OUT, ">$prefix.haplotyping_mhcii.log");
open (OUT_HP, ">$prefix.hp_mhcii.txt") or die "Cannot write $prefix.hp\n";
open (OUT_HP_GENES, ">$prefix.hp_genes_mhcii.txt") or die "Cannot write $prefix.hp.genes\n";
open (OUT_NONHP_DQA, ">$prefix.dqa_nonhp.txt") or die "Cannot write $prefix.dqa.nonhp\n";
open (OUT_NONHP_DQB, ">$prefix.dqb_nonhp.txt") or die "Cannot write $prefix.dqb.nonhp\n";
open (OUT_NONHP_DRB3, ">$prefix.drb3_nonhp.txt") or die "Cannot write $prefix.drb3.nonhp\n";

print OUT_HP "Sample_ID\tMHCII_comments_general\tMHCII_comments_specific\tShared alleles\tRead_count_DRB3_selected_for_haplotypes\tRead_count_DQA_selected_for_haplotypes\tRead_count_DQB_selected_for_haplotypes\tNumber of DRB3 alleles\tNumber of DQA alleles\tNumber of DQB alleles\tRead_count_DRB3_discarded_mapped\tRead_count_DQA_discarded_mapped\tRead_count_DQB_discarded_mapped\tDRB3_total_haplotypes\tDQA_total_haplotypes\tDQB_total_haplotypes\tDRB3_unassigned\tDQA_unassigned\tDQB_unassigned\tRead_count_DRB3_in_haplotypeslotypes\tRead_count_DQA_in_haplotypes\tRead_count_DQB_in_haplotypes\tHaplotype1\tHaplotype2\n";
print OUT_HP_GENES "\tHaplotype1 DRB3\tHaplotype1 DQA\tHaplotype1 DQA\tHaplotype1 DQB\tHaplotype1 DQB\t\tHaplotype2 DRB3\tHaplotype2 DQA\tHaplotype2 DQA\tHaplotype2 DQB\tHaplotype2 DQB\n";
print OUT_NONHP_DRB3 "\tUnassinged DRB3\n";
print OUT_NONHP_DQA "\tUnassinged DQA\n";
print OUT_NONHP_DQB "\tUnassinged DQB\n";

my %haplotypes = ();
if (-e $haplotypes_file){
	open (HAP, "$haplotypes_file") or die "Cannot open $haplotypes_file\n";
	while(<HAP>){
		chomp $_;
		my @words = split ("\t", $_);
		my $hp = $words[0];
		# print "$hp";
		if (! exists $haplotypes{$hp}){
			$haplotypes{$hp} = 1;
			if ($words[1] ne ""){
				$drb3_gene{$hp} = $words[1];
				$$hp{$words[1]} = 1;
				# print "\t$words[1]";
			}
			if ($words[2] ne ""){
				$dqa1_gene{$hp} = $words[2];
				$$hp{$words[2]} = 2;
				# print "\t$words[2]";
			}
			if ($words[3] ne ""){
				$dqa2_gene{$hp} = $words[3];
				$$hp{$words[3]} = 3;
				# print "\t$words[3]";
			}
			if ($words[4] ne ""){
				$dqb1_gene{$hp} = $words[4];
				$$hp{$words[4]} = 4;
				# print "\t$words[4]";
			}
			if ($words[5] ne ""){
				$dqb2_gene{$hp} = $words[5];
				$$hp{$words[5]} = 5;
				# print "\t$words[5]";
			}
			# print "\n";
		}
		else{
			print "double $hp\n";
		}
	}
}
else{
	die "Cannot open $haplotypes_file\n";
}

my %selected_dqa = ();
my %selected_dqb = ();
my %selected_drb3 = ();
if (-e $selected_file_dqa){
	open (SELECTED, "$selected_file_dqa");
	while(<SELECTED>){
		chomp $_;
		my @words = split("\t");
		my $sample = $words[0];
		if ($dqa_counts{$sample} > 100){
			$total_dqa_alleles{$sample} = scalar(@words) - 1;
			if (scalar(@words) > 1){
				for (my $i = 1; $i < scalar(@words); $i++){
					if (! exists $selected_dqa{$sample}){
						$selected_dqa{$sample} = $words[$i];
					}
					else{
						$selected_dqa{$sample} = $selected_dqa{$sample}.";".$words[$i];
					}
				}	
			}
		}
		else{
			$total_dqa_alleles{$sample} = 0;
		}
	}
}
if (-e $selected_file_dqb){
	open (SELECTED, "$selected_file_dqb");
	while(<SELECTED>){
		chomp $_;
		my @words = split("\t");
		my $sample = $words[0];
		if ($dqb_counts{$sample} > 100){
			$total_dqb_alleles{$sample} = scalar(@words) - 1;
			if (scalar(@words) > 1){
				for (my $i = 1; $i < scalar(@words); $i++){
					if (! exists $selected_dqb{$sample}){
						$selected_dqb{$sample} = $words[$i];
					}
					else{
						$selected_dqb{$sample} = $selected_dqb{$sample}.";".$words[$i];
					}
				}	
			}
		}
		else{
			$total_dqb_alleles{$sample} = 0;
		}
	}
}
if (-e $selected_file_drb3){
	open (SELECTED, "$selected_file_drb3");
	while(<SELECTED>){
		chomp $_;
		my @words = split("\t");
		my $sample = $words[0];
		if ($drb3_counts{$sample} > 100){
			$total_drb3_alleles{$sample} = scalar(@words) - 1;
			if (scalar(@words) > 1){
				for (my $i = 1; $i < scalar(@words); $i++){
					if (! exists $selected_drb3{$sample}){
						$selected_drb3{$sample} = $words[$i];
					}
					else{
						$selected_drb3{$sample} = $selected_drb3{$sample}.";".$words[$i];
					}
				}
			}	
		}
		else{
			$total_drb3_alleles{$sample} = 0;
		}
	}
}

my %discarded_dqa = ();
my %discarded_dqb = ();
if (-e $discarded_file_dqa){
	open (DISCARDED, "$discarded_file_dqa");
	while(<DISCARDED>){
		chomp $_;
		my @words = split("\t");
		my $sample = $words[0];
		if ($dqa_counts{$sample} > 100){
			if (scalar(@words) > 1){
				for (my $i = 1; $i < scalar(@words); $i++){
					my @info = split(/\|/, $words[$i]);
					my @freq = split("-", $info[1]);
					if ($freq[1] > 1){
						if (! exists $discarded_dqa{$sample}){
							$discarded_dqa{$sample} = $words[$i];
						}
						else{

							$discarded_dqa{$sample} = $discarded_dqa{$sample}.";".$words[$i];
						}
					}
				}	
			}
		}
	}
}
if (-e $discarded_file_dqb){
	open (DISCARDED, "$discarded_file_dqb");
	while(<DISCARDED>){
		chomp $_;
		my @words = split("\t");
		my $sample = $words[0];
		if ($dqb_counts{$sample} > 100){
			if (scalar(@words) > 1){
				for (my $i = 1; $i < scalar(@words); $i++){
					my @info = split(/\|/, $words[$i]);
					my @freq = split("-", $info[1]);
					if ($freq[1] > 1){
						if (! exists $discarded_dqb{$sample}){
							$discarded_dqb{$sample} = $words[$i];
						}
						else{
							$discarded_dqb{$sample} = $discarded_dqb{$sample}.";".$words[$i];
						}
					}
				}	
			}
		}
	}
}

if (-e $samplesheet){
	open (SAMPLES, "$samplesheet") or die "Cannot open $samplesheet\n";
	while(<SAMPLES>){
		chomp $_;
		my @words = split("\t", $_);
		my $sample = $words[0];
		print OUT_HP "$sample";
		print OUT_HP_GENES "$sample";
		print OUT_NONHP_DRB3 "$sample";
		print OUT_NONHP_DQA "$sample";
		print OUT_NONHP_DQB "$sample";
		print OUT "\n**** $sample ****";
		print OUT "\t$dqa_counts{$sample}($total_dqa_alleles{$sample})\t$dqb_counts{$sample}($total_dqb_alleles{$sample})\t$drb3_counts{$sample}($total_drb3_alleles{$sample})\n";
		if ($dqa_counts{$sample} > 100 or $dqb_counts{$sample} > 100 or $drb3_counts{$sample} > 100){
			# print OUT_HP "$sample";
			# print OUT_HP_GENES "$sample";
			# print OUT_NONHP_DQ "$sample";
			# print OUT_NONHP_DRB3 "$sample";
			# print OUT_NONHP_DQA "$sample";
			# print OUT_NONHP_DQB "$sample";
			my $gene = "";
			my %genes = ();
			my %count_dqa = ();
			my %count_dqb = ();
			my %count_drb3 = ();
			my %freq_genes = ();
			my $comment_general = "";
			my $comment_specific = "";
			my $shared_comment = "";
			my %discarded_genes = ();
			if (exists $selected_dqa{$words[0]}){
				my @list = split (";", $selected_dqa{$words[0]});
				foreach my $item (@list){
					my @info = split(/\|/, $item);
					$gene = $info[0];
					if (! exists $genes{$gene}){
						@freq = split("-", $info[1]);
						$freq_genes{$gene} = $freq[0];
						$genes{$gene} = $info[1];
						$count_dqa{$gene} = 0;
						print OUT "selected DQA: $gene|$genes{$gene}\n";
					}
				}
			}
			if (exists $selected_dqb{$words[0]}){
				my @list = split (";", $selected_dqb{$words[0]});
				foreach my $item (@list){
					my @info = split(/\|/, $item);
					$gene = $info[0];
					if (! exists $genes{$gene}){
						@freq = split("-", $info[1]);
						$freq_genes{$gene} = $freq[0];
						$genes{$gene} = $info[1];
						$count_dqb{$gene} = 0;
						print OUT "selected DQB: $gene|$genes{$gene}\n";
					}
				}
			}

			if (exists $selected_drb3{$words[0]}){
				my @list = split (";", $selected_drb3{$words[0]});
				foreach my $item (@list){
					my @info = split(/\|/, $item);
					$gene = $info[0];
					if (! exists $genes{$gene}){
						@freq = split("-", $info[1]);
						$freq_genes{$gene} = $freq[0];
						$genes{$gene} = $info[1];
						$count_drb3{$gene} = 0;
						print OUT "selected DRB3: $gene|$genes{$gene}\n";
					}
				}
			}
			
			my $key_drb3 = keys %count_drb3;
			my $key_dqa = keys %count_dqa;
			my $key_dqb = keys %count_dqb;

			my %possible_hp = ();
			my %checked_genes = ();
			#To find all possible haplotypes in animals
			LOOP1: foreach my $hp (keys %haplotypes){
				foreach my $gene (keys %$hp){
					if (exists $genes{$gene}){
						# print OUT "$hp $gene 1\n";
						$flag = 0;
					}
					else{
						$flag = 1;
						# print OUT "$hp $gene 0\n";
						next LOOP1;
					}
				}
				if ($flag == 0){
					$possible_hp{$hp} = keys %$hp;
					print OUT "Possible $hp\n";
					foreach my $gene (keys %$hp){
						$checked_genes{$gene}++;
					}
				}
			}

			my %final_hp = ();
			my %checked_genes1 = ();
			#To remove haplotypes - one or more alleles ar repeated
			foreach my $hp (sort {$possible_hp{$a} <=> $possible_hp{$b}} keys %possible_hp){
				$flag_not_repeat = 0;
				foreach my $gene (keys %$hp){
					if ($checked_genes{$gene} == 1){
						$flag_not_repeat = 1;
					}
				}
				if ($flag_not_repeat eq 1){
					$final_hp{$hp} = $possible_hp{$hp};
					foreach my $gene (keys %$hp){
						$checked_genes1{$gene} = 1;
					}
				}
				else{
					print OUT "$hp\trepeated";
					foreach my $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
						$checked_genes{$gene}--;
					}
					print OUT "\n";
				}
			}
			
			my $flag_nonhp_dqa = 0;
			my $flag_nonhp_dqb = 0;
			my $flag_nonhp_drb3 = 0;
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %count_dqa){
				if (!exists $checked_genes1{$gene} and $dqa_counts{$sample} > 100){
					$flag_nonhp_dqa = 1;
				}
			}
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %count_dqb){
				if (!exists $checked_genes1{$gene} and $dqb_counts{$sample} > 100){
					$flag_nonhp_dqb = 1;
				}
			}
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %count_drb3){
				if (!exists $checked_genes1{$gene} and $drb3_counts{$sample} > 100){
					$flag_nonhp_drb3 = 1;
				}
			}
			print OUT "Unassigned - DRB3:$flag_nonhp_drb3 DQA:$flag_nonhp_dqa DQB:$flag_nonhp_dqb\n";
			my %checked_rescued_genes = ();
			my %possible_rescued_hp = ();
			my $flag_rescued = 0;
			my %rescued_hp = ();
			if ($flag_nonhp_drb3 eq 1 or $flag_nonhp_dqa eq 1 or $flag_nonhp_dqb eq 1){
				# print "$sample This is running $discarded_dqb{$sample}\n";
				if (exists $discarded_dqa{$sample}){
					my @list = split (";", $discarded_dqa{$sample});
					foreach my $item (@list){
						my @info = split(/\|/, $item);
						$gene = $info[0];
						if (! exists $genes{$gene}){
							@freq = split("-", $info[1]);
							$freq_genes{$gene} = $freq[0];
							$discarded_genes{$gene} = $info[1];
							$genes{$gene} = $info[1];
							print OUT "discarded DQA: $gene|$discarded_genes{$gene}\n";
						}
					}
				}
				if (exists $discarded_dqb{$sample}){
					my @list = split (";", $discarded_dqb{$sample});
					foreach my $item (@list){
						my @info = split(/\|/, $item);
						$gene = $info[0];
						if (! exists $genes{$gene}){
							@freq = split("-", $info[1]);
							$freq_genes{$gene} = $freq[0];
							$discarded_genes{$gene} = $info[1];
							$genes{$gene} = $info[1];
							print OUT "discarded DQB: $gene|$discarded_genes{$gene}\n";
						}
					}
				}
				if (exists $discarded_dqa{$sample} or exists $discarded_dqb{$sample}){
					$flag = 1;
					LOOP1: foreach my $hp (keys %haplotypes){
						if (!exists $final_hp{$hp} and !exists $possible_hp{$hp}){
							foreach my $gene (keys %$hp){
								if (exists $genes{$gene}){
									$flag = 0;
								}
								else{
									$flag = 1;
									next LOOP1;
								}
							}
							if ($flag == 0){
								$possible_rescued_hp{$hp} = keys %$hp;
								print OUT "Possible rescued $hp\n";
								foreach my $gene (keys %$hp){
									$checked_rescued_genes{$gene}++;
								}
							}
						}
					}

					#To remove haplotypes - one or more alleles ar repeated
					foreach my $hp (sort {$possible_rescued_hp{$a} <=> $possible_rescued_hp{$b}} keys %possible_rescued_hp){
						$flag_not_repeat = 0;
						foreach my $gene (keys %$hp){
							if ($checked_rescued_genes{$gene} == 1){
								$flag_not_repeat = 1;
							}
						}
						if ($flag_not_repeat eq 1){
							$flag_rescued = 1;
							$rescued_hp{$hp} = 1;
							$final_hp{$hp} = $possible_hp{$hp};
						}
						else{
							print OUT "$hp\trepeated";
							foreach my $gene (keys %$hp){
								if (exists $discarded_genes{$gene}){
									print OUT "\t$gene|$discarded_genes{$gene}";
								}
								else{
									print OUT "\t$gene|$genes{$gene}";
								}
								$checked_rescued_genes{$gene}--;
							}
							print OUT "\n";
						}
					}
				}
			}

			my $total_unique_dqa = 0;
			my $total_unique_dqb = 0;
			my $total_unique_drb3 = 0;
			my $total_hp_dqa = 0;
			my $total_hp_dqb = 0;
			my $total_hp_drb3 = 0;
			my %checked_genes2 = ();
			foreach my $hp (sort {$final_hp{$a} <=> $final_hp{$b}} keys %final_hp){
				print OUT "final $hp";
				foreach $gene (sort {$$hp{$a} <=> $$hp{$b}} keys %$hp){
					print OUT "\t$gene|$genes{$gene}";
					if ($gene =~ /DRB3/){
						if (!exists $checked_genes2{$gene}){
							$checked_genes2{$gene} = 0;
							$total_unique_drb3 = $total_unique_drb3 + $freq_genes{$gene};
						}
						$total_hp_drb3 = $total_hp_drb3 + $freq_genes{$gene};
					}
					if ($gene =~ /DQA/){
						if (!exists $checked_genes2{$gene}){
							$checked_genes2{$gene} = 0;
							$total_unique_dqa = $total_unique_dqa + $freq_genes{$gene};
						}
						$total_hp_dqa = $total_hp_dqa + $freq_genes{$gene};
					}
					if ($gene =~ /DQB/){
						if (!exists $checked_genes2{$gene}){
							$checked_genes2{$gene} = 0;
							$total_unique_dqb = $total_unique_dqb + $freq_genes{$gene};
						}
						$total_hp_dqb = $total_hp_dqb + $freq_genes{$gene};
					}
					$checked_genes2{$gene}++;
				}
				print OUT "\n";
			}

			my $total_per_dqa = 0;
			my $total_per_dqb = 0;
			my $total_per_drb3 = 0;
			if ($total_hp_dqa > 0){
				$total_per_dqa = sprintf("%.2f", ($total_hp_dqa / $total_unique_dqa) * 100);
			}
			if ($total_hp_dqb > 0){
				$total_per_dqb = sprintf("%.2f", ($total_hp_dqb / $total_unique_dqb) * 100);
			}
			if ($total_hp_drb3 > 0){
				$total_per_drb3 = sprintf("%.2f", ($total_hp_drb3 / $total_unique_drb3) * 100);
			}

			my $shared_genes = "";
			my $shared_flag = 0;
			if (($total_per_dqa > 100) or ($total_per_dqb > 100) or ($total_per_drb3 > 100)){
				foreach $gene (keys %checked_genes){
					if ($checked_genes2{$gene} > 1){
						print OUT "$gene is shared $checked_genes{$gene}\n";
						$shared_flag = 1;
						$shared_genes = $shared_genes."$gene|$genes{$gene},";
					}
				}
				chop $shared_genes;
			}

			my $count_nonhp_dqa = 0;
			my $count_nonhp_dqb = 0;
			my $count_nonhp_drb3 = 0;
			my $per_freq_dqa_nonhp = 0;
			my $per_freq_dqb_nonhp = 0;
			my $per_freq_drb3_nonhp = 0;
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %count_dqa){
				if (!exists $checked_genes2{$gene} and $dqa_counts{$sample} > 100 and !exists $discarded_genes{$gene}){
					print OUT_NONHP_DQA "\t$gene|$genes{$gene}";
					print OUT "Unassigned DQA: $gene|$genes{$gene}\n";
					$count_nonhp_dqa = $count_nonhp_dqa + $freq_genes{$gene};
				}
			}
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %count_dqb){
				if (!exists $checked_genes2{$gene} and $dqb_counts{$sample} > 100 and !exists $discarded_genes{$gene}){
					print OUT_NONHP_DQB "\t$gene|$genes{$gene}";
					print OUT "Unassigned DQB: $gene|$genes{$gene}\n";
					$count_nonhp_dqb = $count_nonhp_dqb + $freq_genes{$gene};
				}
			}
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %count_drb3){
				if (!exists $checked_genes2{$gene} and $drb3_counts{$sample} > 100){
					print OUT_NONHP_DRB3 "\t$gene|$genes{$gene}";
					print OUT "Unassigned DRB3: $gene|$genes{$gene}\n";
					$count_nonhp_drb3 = $count_nonhp_drb3 + $freq_genes{$gene};
				}
			}
			
			
			if ($count_nonhp_dqa > 0 and $dqa_counts{$sample} > 0){
				$per_freq_dqa_nonhp = sprintf("%.2f", ($count_nonhp_dqa / $dqa_counts{$sample}) * 100);
			}
			if ($count_nonhp_dqb > 0 and $dqb_counts{$sample}){
				$per_freq_dqb_nonhp = sprintf("%.2f", ($count_nonhp_dqb / $dqb_counts{$sample}) * 100);
			}
			if ($count_nonhp_drb3 > 0 and $drb3_counts{$sample}){
				$per_freq_drb3_nonhp = sprintf("%.2f", ($count_nonhp_drb3 / $drb3_counts{$sample}) * 100);
			}
			
			if ($dqa_counts{$sample} <= 100 ){
				$comment_general = $comment_general."No DQA,";
			}
			if ($dqb_counts{$sample} <= 100 ){
				$comment_general = $comment_general."No DQB,";
			}
			if ($drb3_counts{$sample} <= 100 ){
				$comment_general = $comment_general."No DRB3,";
			}

			my $count = keys %final_hp;
			if ($count eq 1){
				$comment_general = $comment_general."Homozygous,";
			}
			elsif ($count eq 3){
				$comment_general = $comment_general."3 haplotypes,";
			}
			elsif ($count eq 4){
				$comment_general = $comment_general."4 haplotypes,";
			}
			elsif ($count > 4){
				$comment_general = $comment_general."More than 4 haplotypes,";
			}
			if ($flag_rescued eq 1){
				$comment_general = $comment_general."Rescued,";
			}
			if ($count_nonhp_dqa > 0 or $count_nonhp_dqb > 0 or $count_nonhp_drb3 > 0){
				$comment_general = $comment_general."Unassigned alleles,";
			}
			chop $comment_general;
			chop $comment_specific;
			
			print OUT_HP "\t$comment_general\t$comment_specific\t$shared_genes\t$drb3_counts{$sample}\t$dqa_counts{$sample}\t$dqb_counts{$sample}\t$total_drb3_alleles{$sample}\t$total_dqa_alleles{$sample}\t$total_dqb_alleles{$sample}\t$drb3_counts_discarded{$sample}\t$dqa_counts_discarded{$sample}\t$dqb_counts_discarded{$sample}\t$total_unique_drb3\t$total_unique_dqa\t$total_unique_dqb\t$per_freq_drb3_nonhp\t$per_freq_dqa_nonhp\t$per_freq_dqb_nonhp\t$total_per_drb3\t$total_per_dqa\t$total_per_dqb";
			foreach my $hp (sort { $final_hp{$b} <=> $final_hp{$a} } keys %final_hp){
				print OUT_HP "\t$hp";
				if (exists $rescued_hp{$hp}){
					print OUT_HP "#";
				}
				if (exists $drb3_gene{$hp}){
					my $drb3 = $drb3_gene{$hp};
					print OUT_HP_GENES "\t$drb3|$genes{$drb3}";
				}
				else{
					print OUT_HP_GENES "\t";
				}
				if (exists $dqa1_gene{$hp}){
					my $dqa1 = $dqa1_gene{$hp};
					print OUT_HP_GENES "\t$dqa1|$genes{$dqa1}";
				}
				else{
					print OUT_HP_GENES "\t";
				}
				if (exists $dqa2_gene{$hp}){
					my $dqa2 = $dqa2_gene{$hp};
					print OUT_HP_GENES "\t$dqa2|$genes{$dqa2}";
				}
				else{
					print OUT_HP_GENES "\t";
				}
				if (exists $dqb1_gene{$hp}){
					my $dqb1 = $dqb1_gene{$hp};
					print OUT_HP_GENES "\t$dqb1|$genes{$dqb1}";
				}
				else{
					print OUT_HP_GENES "\t";
				}
				if (exists $dqb2_gene{$hp}){
					my $dqb2 = $dqb2_gene{$hp};
					print OUT_HP_GENES "\t$dqb2|$genes{$dqb2}";
				}
				else{
					print OUT_HP_GENES "\t";
				}
				print OUT_HP_GENES "\t";
			}
		}
		else{
			print OUT_HP "\tNo data\t\t\t$drb3_counts{$sample}\t$dqa_counts{$sample}\t$dqb_counts{$sample}\t$total_drb3_alleles{$sample}\t$total_dqa_alleles{$sample}\t$total_dqb_alleles{$sample}\t$drb3_counts_discarded{$sample}\t$dqa_counts_discarded{$sample}\t$dqb_counts_discarded{$sample}";
		}
		print OUT_HP "\n";	
		print OUT_HP_GENES "\n";
		print OUT_NONHP_DQA "\n";
		print OUT_NONHP_DQB "\n";
		print OUT_NONHP_DRB3 "\n";
	}
}
else{
	die "Cannot open $samplesheet\n";
}
