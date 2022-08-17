use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Getopt::Long;

my $haplotypes_file = "";
my $seleced_file = "";
my $discarded_file = "";
my $summary_file = "";
my $fasta = "";
my $prefix = "";


GetOptions(
    'haplotypes=s'    => \$haplotypes_file,
    'filtered=s'     => \$seleced_file,
    'discarded=s' => \$discarded_file,
    'summary=s'     => \$summary_file,
    'database=s'     => \$fasta,
    'prefix=s' => \$prefix
) or print "Invalid options\n";


my %database = ();
open(FASTA, "$fasta");
while(<FASTA>){
	chomp $_;
	if (/^>(\S+)/){
		$database{$1} = 1;
	}
}

my %primer1_counts = ();
my %primer2_counts = ();
my %allele_info_for1_selected = ();
my %allele_info_for3_selected = ();
my %allele_info_for1_discarded = ();
my %allele_info_for3_discarded = ();

open(P1, "$summary_file") or die "Cannot open $summary_file\n";
while(<P1>){
	chomp $_;
	@words = split("\t", $_);
	$primer1_counts{$words[0]} = $words[4];
	$primer2_counts{$words[0]} = $words[27];
	$allele_info_for1_selected{$words[0]} = $words[12];
	$allele_info_for3_selected{$words[0]} = $words[35];
	$allele_info_for1_discarded{$words[0]} = $words[20];
	$allele_info_for3_discarded{$words[0]} = $words[43];
}

open (OUT, ">$prefix.haplotyping_mhci.log");
open (OUT_HP, ">$prefix.hp_mhci.txt") or die "Cannot write $prefix.hp\n";
open (OUT_HP_GENES, ">$prefix.hp_alleles_mhci.txt") or die "Cannot write $prefix.hp.genes\n";
open (OUT_NONHP, ">$prefix.nonhp_mhci.txt") or die "Cannot write $prefix.nonhp\n";
open (OUT_CONT, ">$prefix.contamination_mhci.txt") or die "Cannot write $prefix.contamination\n";
# open (OUT_BR22, ">$prefix.br22") or die "Cannot write $prefix.br22\n";

print OUT_HP "Sample_ID\tAnimal_ID\tBreed\tLineage\tMHCI_comments_general\tMHCI_comments_specific\tTotal_Read_count_For_1\tTotal_read_count_For_3\tRead_count_For1_selected_for_haplotypes\tRead_count_For3_selected_for_haplotypes\tRead_count_For1_discarded_mapped\tRead_count_For3_discarded_mapped\tRead_count_For1_in_haplotypes\tRead_count_For3_in_haplotypes\tFor1_total_haplotypes\tFor3_total_haplotypes\tFor1_unassigned\tFor3_unassigned\tFor1_contamination_counts\tFor3_contamination_counts\tHaplotype1\tHaplotype1_total_For1_read_counts\tHaplotype1_total_For1_read_counts_%\tHaplotype1_total_For3_read_counts\tHaplotype1_total_For3_read_counts_%\n";
			
my %haplotypes = ();
my $noOfGenes = 0;
if (-e $haplotypes_file){
	open (HAP, "$haplotypes_file") or die "Cannot open $haplotypes_file\n";
	while(<HAP>){
		chomp $_;
		my @words = split ("\t", $_);
		my $hp = $words[0];
		if (! exists $haplotypes{$hp}){
			$haplotypes{$hp} = 0;
		}
		else{
			print "double $hp\n";
		}
		
		if ($noOfGenes < scalar(@words)){
			$noOfGenes = scalar(@words);
		}
		for (my $i = 1; $i < scalar(@words); $i++){
			$$hp{$words[$i]} = 1;
		}
		$higher_expressed{$hp} = $words[1];
	}
}
else{
	die "Cannot open $haplotypes_file\n";
}

my %discarded = ();
if (-e $discarded_file){
	open (DISCARDED, "$discarded_file");
	while(<DISCARDED>){
		chomp $_;
		my @words = split("\t");
		my $sample = $words[0];
		if (scalar(@words) > 1){
			for (my $i = 1; $i < scalar(@words); $i++){
				if (! exists $discarded{$sample}){
					$discarded{$sample} = $words[$i];
				}
				else{
					$discarded{$sample} = $discarded{$sample}.",".$words[$i];
				}
			}	
		}
	}
}

if (-e $seleced_file){
	open (GENES_IN, "$seleced_file") or die "Cannot open $genes_file\n";
	while(<GENES_IN>){
		chomp $_;
		my @words = split("\t", $_);
		my $sample = $words[0];
		if ($allele_info_for1_selected{$sample} > 100 or $allele_info_for3_selected{$sample} > 100){
			my $gene = "";
			my %genes = ();
			my %selected_genes = ();
			my $flag = 1;
			my %toWrite_hp = ();
			my %text = ();
			my %freq_genes = ();
			my %freq_genes1 = ();
			my %freq_genes2 = ();
			my %freq_genes1_per = ();
			my %freq_genes2_per = ();
			my %per_freq1 = ();
			my %per_freq2 = ();
			my %discarded_genes = ();
			my $comment_general = "";
			my $comment_specific = "";
			my %single_id = ();
			my %double_id = ();

			print OUT "\n**** $sample ****";
			print OUT_CONT "$sample";
			print OUT_NONHP "$sample";
			# print OUT_BR22 "$sample";

			print OUT "\t$allele_info_for1_selected{$sample}\t$allele_info_for3_selected{$sample}\n";

			for (my $i = 1; $i < scalar(@words); $i++){
				my @info = split(/\|/, $words[$i]);
				if ($info[0] =~ /\// and !exists $database{$info[0]}){
					@ids = split(/\//, $info[0]);
					$double_id{$info[0]} = 0;
				}
				else{
					@ids = $info[0];
					$single_id{$info[0]} = 0;
				}
				foreach my $gene (@ids){
					if (! exists $genes{$gene}){
						@freq1 = split("-", $info[1]);
						@freq2 = split("-", $info[2]);
						$freq_genes1{$gene} = $freq1[0];
						$freq_genes1_per{$gene} = $freq1[1];
						$freq_genes2{$gene} = $freq2[0];
						$freq_genes2_per{$gene} = $freq2[1];
						$freq_genes{$gene} = $freq1[0] + $freq2[0];
						$genes{$gene} = $info[1]."|".$info[2];
					}
				}
			}
			if (exists $discarded{$words[0]}){
				my @list = split (",", $discarded{$words[0]});
				foreach my $item (@list){
					my @info = split(/\|/, $item);
					if ($info[0] =~ /\// and !exists $database{$info[0]}){
						@ids = split(/\//, $info[0]);
					}
					else{
						@ids = $info[0];
					}
					foreach my $gene (@ids){
						if (! exists $genes{$gene}){
							@freq1 = split("-", $info[1]);
							@freq2 = split("-", $info[2]);
							$freq_genes1{$gene} = $freq1[0];
							$freq_genes1_per{$gene} = $freq1[1];
							$freq_genes2{$gene} = $freq2[0];
							$freq_genes2_per{$gene} = $freq2[1];
							$freq_genes{$gene} = $freq1[0] + $freq2[0];
							$genes{$gene} = $info[1]."|".$info[2];
							$discarded_genes{$gene} = 1;
						}
					}
				}
			}

			my %possible_hp = ();
			#To find all possible haplotypes in animals
			LOOP1: foreach my $hp (keys %haplotypes){
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
					$possible_hp{$hp} = keys %$hp;
				}
			}

			my %possible_hp1 = ();
			my %checked_genes1 = ();
			#To remove mishaplotype with all discarded alleles
			foreach my $hp (sort {$possible_hp{$a} <=> $possible_hp{$b}} keys %possible_hp){
				$flag_discarded = 0;
				foreach my $gene (keys %$hp){
					if (!exists $discarded_genes{$gene}){
						$flag_discarded = 1;
					}
				}
				if ($flag_discarded eq 1){
					foreach my $gene (keys %$hp){
						$checked_genes1{$gene}++;
					}
					$possible_hp1{$hp} = $possible_hp{$hp};
				}
				else{
					print OUT "$hp\tdiscarded_hp";
					foreach my $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";						
					}
					print OUT "\n";
				}
			}

			my %possible_hp2 = ();
			my %checked_genes2 = ();
			#To remove haplotypes - one or more alleles ar repeated
			foreach my $hp (sort {$possible_hp1{$a} <=> $possible_hp1{$b}} keys %possible_hp1){
				$flag_not_repeat = 0;
				foreach my $gene (keys %$hp){
					if ($checked_genes1{$gene} == 1){
						$flag_not_repeat = 1;
					}
				}
				if ($flag_not_repeat eq 1){
					foreach my $gene (keys %$hp){
						$checked_genes2{$gene}++;
					}
					$possible_hp2{$hp} = $possible_hp{$hp};
				}
				else{
					print OUT "$hp\trepeated";
					foreach my $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
						$checked_genes1{$gene}--;
					}
					print OUT "\n";
				}
				
			}

			my %final_hp = ();
			foreach my $hp (sort {$possible_hp2{$a} <=> $possible_hp2{$b}} keys %possible_hp2){
				$flag_discarded = 0;
				foreach my $gene (keys %$hp){
					if ($checked_genes2{$gene} eq 1 and !exists $discarded_genes{$gene}){
						$flag_discarded = 1;
						print OUT "Flag of non discarded $gene\n";
					}
				}
				# print "\n";
				if ($flag_discarded eq 1){
					$final_hp{$hp} = $possible_hp{$hp};
					foreach my $gene (keys %$hp){
						print OUT "$hp $gene|$genes{$gene} $checked_genes2{$gene}\n";
					}
				}
				else{
					print OUT "$hp\tdiscarded_allele";
					foreach my $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
						$checked_genes2{$gene}--;
					}
					print OUT "\n";
				}
			}

			# print "$subtotal_haplotypes1\t$subtotal_haplotypes2\n";
			my $count_gap = 0;

			my $subtotal_haplotypes1 = 0;
			my $subtotal_haplotypes2 = 0;
			my %checked_genes = ();
			my %hp_count1 = ();
			my %hp_count2 = ();

			#To calculate haplotype frequency
			foreach my $hp (sort { $final_hp{$b} <=> $final_hp{a} } keys %final_hp){
				$count1 = 0;
				$count2 = 0;
				$haplotypes{$hp}++;
				foreach $gene (keys %$hp){
					if (!exists $checked_genes{$gene}){
						$subtotal_haplotypes1 = $subtotal_haplotypes1 + $freq_genes1{$gene};
						$subtotal_haplotypes2 = $subtotal_haplotypes2 + $freq_genes2{$gene};
						$checked_genes{$gene} = 1;
					}
					$count1 = $count1 + $freq_genes1{$gene};
					$count2 = $count2 + $freq_genes2{$gene};
				}
				$hp_count1{$hp} = $count1;
				$hp_count2{$hp} = $count2;
				# print "$sample\t$hp\t$count\n";

				# $hp1 = $hp;
				# $hp1 =~ s/[a-z]$//g;
				# if ($haplotypes{$hp1} eq 1){
				# 	open (SCATTER, ">$scatterplot_dir/$hp1.txt");
				# }
				# else{
				# 	open (SCATTER, ">>$scatterplot_dir/$hp1.txt");
				# }

				# foreach $gene (keys %$hp){
				# 	$norm1 = 0;
				# 	$norm2 = 0;
				# 	if ($hp_count1{$hp} > 0){
				# 		$norm1 = ($freq_genes1{$gene} * 100) / $hp_count1{$hp};
				# 	}
				# 	if ($hp_count2{$hp} > 0){
				# 		$norm2 = ($freq_genes2{$gene} * 100) / $hp_count2{$hp};
				# 	}
				# 	print SCATTER "$sample\t$hp1\t$hp\t$gene\t$freq_genes1_per{$gene}\t$freq_genes2_per{$gene}\t$norm1\t$norm2\n";
				# }
				# close (SCATTER);

			}

			my $subtotal_haplotypes1_final = 0;
			my $subtotal_haplotypes2_final = 0;
			%checked_genes = ();
			my $contamination_flag = 0;
			my $contamination_flag1 = 0;
			my $contaminant = "";
			my $contaminant_counts1 = 0;
			my $contaminant_counts2 = 0;
			my $pcr_contamination = 0;
			my $pcr_contaminant = "";
			foreach my $hp (sort { $final_hp{$a} <=> $final_hp{$b} } keys %final_hp){
				$contamination_flag = 0;
				my $per1 = 0;
				my $per2 = 0;
				if ($subtotal_haplotypes1 > 0){
					$per1 = sprintf("%.2f", ($hp_count1{$hp} / $subtotal_haplotypes1) * 100);
				}
				if ($subtotal_haplotypes2 > 0){
					$per2 = sprintf("%.2f", ($hp_count2{$hp} / $subtotal_haplotypes2) * 100);
				}
				if ($per1 < 5 and $per2 < 5){
					print OUT "$hp\tcontamination $per1 $per2";
					foreach $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
					}
					$higher_expressed_allele = $higher_expressed{$hp};
					print OUT "\t($higher_expressed_allele\t$freq_genes1_per{$higher_expressed_allele}\t$freq_genes2_per{$higher_expressed_allele})";
					if ($freq_genes1_per{$higher_expressed_allele} <= 1 and $freq_genes2_per{$higher_expressed_allele} <= 1 ){
						$contamination_flag = 1;
						$contamination_flag1 = 1;
					}
					print OUT "\n";
				}
				elsif ($per2 < 5 and $allele_info_for3_selected{$sample} > 100){
					print OUT "$hp\tcontamination $per1 $per2";
					foreach $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
					}
					$higher_expressed_allele = $higher_expressed{$hp};
					print OUT "\t($higher_expressed_allele\t$freq_genes2_per{$higher_expressed_allele})";
					if ($freq_genes2_per{$higher_expressed_allele} <= 1){
						$contamination_flag = 1;
						$contamination_flag1 = 1;
						$pcr_contamination = 1;
					}
					print OUT "\n";
				}
				elsif ($per1 < 5 and $allele_info_for1_selected{$sample} > 100){
					print OUT "$hp\tcontamination $per1 $per2";
					foreach $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
					}
					$higher_expressed_allele = $higher_expressed{$hp};
					print OUT "\t($higher_expressed_allele\t$freq_genes1_per{$higher_expressed_allele})";
					if ($freq_genes1_per{$higher_expressed_allele} <= 1){
						$contamination_flag = 1;
						$contamination_flag1 = 1;
						$pcr_contamination = 1;
					}
					print OUT "\n";
				}
				if ($contamination_flag eq 1){
					%text = ();
					$discarded_flag = 0;
					foreach $gene (keys %$hp){
						$text{$gene} = $freq_genes{$gene};
						if (! exists $selected_genes{$gene}){
							$contaminant_counts1 = $contaminant_counts1 + $freq_genes1{$gene};
							$contaminant_counts2 = $contaminant_counts2 + $freq_genes2{$gene};
						}
						$selected_genes{$gene}++;
						if (exists $discarded_genes{$gene}){
							$discarded_flag = 1;
						}
					}
					$t = "";
					$count_gap = 0;
					$origin_hp = $hp;
					foreach $gene (sort {$text{$b} <=> $text{$a}} keys %text){
						if (exists $discarded_genes{$gene}){
							$t = $t."$gene#|$genes{$gene}\t";
						}
						else{
							$t = $t."$gene|$genes{$gene}\t";
						}
						$count_gap++;
					}
					chop $t;
					if ($count_gap < $noOfGenes){
						for($i = $count_gap; $i < $noOfGenes; $i++){
							$t = $t."\t";
						}
					}
					if ($discarded_flag == 1){
						$hp = $hp."#";
					}
					if ($pcr_contamination eq 1){
						$pcr_contaminant = $pcr_contaminant."$hp,";
					}
					else{
						$contaminant = $contaminant."$hp,";
					}
					print OUT_CONT "\t$t";
				}
				else{
					print OUT "$hp\tfinal";
					%text = ();
					$discarded_flag = 0;
					foreach $gene (keys %$hp){
						print OUT "\t$gene|$genes{$gene}";
						$text{$gene} = $freq_genes{$gene};
						$selected_genes{$gene}++;
						if (exists $discarded_genes{$gene}){
							$discarded_flag = 1;
						}
						if (!exists $checked_genes{$gene}){
							$subtotal_haplotypes1_final = $subtotal_haplotypes1_final + $freq_genes1{$gene};
							$subtotal_haplotypes2_final = $subtotal_haplotypes2_final+ $freq_genes2{$gene};
							$checked_genes{$gene} = 1;
						}
					}
					print OUT "\n";
					$origin_hp = $hp;
					$t = "";
					$count_gap = 0;
					foreach $gene (sort {$text{$b} <=> $text{$a}} keys %text){
						if (exists $discarded_genes{$gene}){
							$t = $t."$gene#|$genes{$gene}\t";
						}
						else{
							$t = $t."$gene|$genes{$gene}\t";
						}
						$count_gap++;
					}
					chop $t;
					if ($count_gap < $noOfGenes){
						for($i = $count_gap; $i < $noOfGenes; $i++){
							$t = $t."\t";
						}
					}
					if ($discarded_flag == 1){
						$hp = $hp."#";
						$toWrite_hp{$hp} = $t;
					}
					else{
						$toWrite_hp{$hp} = $t;
					}
					$order_towrite{$hp} = $hp_count1{$origin_hp} + $hp_count2{$origin_hp};
				}
			}

			$count1 = 0;
			$count2 = 0;
			my %selected_nonhp = ();
			my $toWriteNonhp = "";
			my $flag_double_id = 0;
			foreach my $gene (sort {$freq_genes{$b} <=> $freq_genes{$a}} keys %genes){
				if (!exists $selected_genes{$gene} and !exists $discarded_genes{$gene} and !exists $selected_nonhp{$gene}){
					if (exists $single_id{$gene}){
						# if ($gene =~ /br22:/){
						# 	print OUT_BR22 "\t$gene|$genes{$gene}";
						# }
						# else{
							print OUT_NONHP "\t$gene|$genes{$gene}";
						# }
						$count1 = $count1 + $freq_genes1{$gene};
						$count2 = $count2 + $freq_genes2{$gene};
						$selected_nonhp{$gene} = 1;
					}
					else{
						foreach my $i (keys %double_id){
							$flag_double_id = 0;
							# print "$sample $i\n";
							@ids = split(/\//, $i);
							$new_gene = "";
							foreach $id (@ids){
								if (exists $selected_genes{$id}){
									$flag_double_id = 1;
								}
								$selected_nonhp{$id} = 1;
								$new_gene = $new_gene."$id/";
							}
							chop $new_gene;
							if ($flag_double_id eq 0){
								# if ($new_gene =~ /br22:/){
								# 	print OUT_BR22 "\t$gene|$genes{$gene}";
								# }
								# else{
									print OUT_NONHP "\t$new_gene|$genes{$gene}";
								# }
								$count1 = $count1 + $freq_genes1{$gene};
								$count2 = $count2 + $freq_genes2{$gene};
								# print "Not hp $new_gene|$genes{$gene}\n";
							}
							else{
								# print "yes hp $gene|$genes{$gene}\n";
							}
						}
					}
				}
			}

			if ($count1 > 0){
				$per_freq1_nonhp = sprintf("%.2f", ($count1 / $primer1_counts{$sample}) * 100);
			}
			else{
				$per_freq1_nonhp = 0;
			}
			if ($count2 > 0){
				$per_freq2_nonhp = sprintf("%.2f", ($count2 / $primer2_counts{$sample}) * 100);
			}
			else{
				$per_freq2_nonhp = 0;
			}
			# print "$sample\n";
			my $total_per1 = 0;
			my $total_per2 = 0;
			my $count_final_hp = 0;
			foreach $hp (sort {$order_towrite{$b} <=> $order_towrite{$a}} keys %toWrite_hp){
				$origin_hp = $hp;
				$origin_hp =~ s/\#//;
				$count_final_hp++;
				if ($subtotal_haplotypes1_final > 0){
					$per_freq1{$origin_hp} = sprintf("%.2f", ($hp_count1{$origin_hp} / $subtotal_haplotypes1_final) * 100);
					$total_per1 = sprintf("%.2f", ($total_per1 + ($hp_count1{$origin_hp} / $subtotal_haplotypes1_final) * 100));
				}
				else{
					$per_freq1{$origin_hp} = 0;
				}
				if ($subtotal_haplotypes2_final > 0){
					$per_freq2{$origin_hp} = sprintf("%.2f", ($hp_count2{$origin_hp} / $subtotal_haplotypes2_final) * 100);
					$total_per2 = sprintf("%.2f", ($total_per2 + ($hp_count2{$origin_hp} / $subtotal_haplotypes2_final) * 100));
				}
				else{
					$per_freq2{$origin_hp} = 0;
				}
			}

			if ($count_final_hp eq 1){
				$comment_general = $comment_general."Homozygous,";
			}
			elsif ($count_final_hp eq 3){
				$comment_general = $comment_general."3 haplotypes,";
			}
			elsif ($count_final_hp eq 4){
				$comment_general = $comment_general."4 haplotypes,";
			}
			elsif ($count_final_hp >4){
				$comment_general = $comment_general."More than 4 haplotypes";
			}

			if ($contamination_flag1 eq 1){
				if ($pcr_contamination eq 1){
					chop $pcr_contaminant;
					$comment_general = $comment_general."PCR contamination removed,";
					$comment_specific = $comment_specific."PCR contaminant: $pcr_contaminant,";
				}
				else{
					chop $contaminant;
					$comment_general = $comment_general."contamination removed,";
					$comment_specific = $comment_specific."contaminant: $contaminant,";
				}
			}
			
			
			if ($allele_info_for1_selected{$sample} <= 100 and $allele_info_for3_selected{$sample} <= 100){
				$comment_general = $comment_general."No data,";
			}
			elsif ($allele_info_for1_selected{$sample} <= 100 ){
				$comment_general = $comment_general."No For1,";
			}
			elsif ($allele_info_for3_selected{$sample} <= 100 ){
				$comment_general = $comment_general."No For3,";
			}
			
			if (($total_per2 > 100) or ($total_per1 > 100)){
				my $shared_genes = "";
				my $flag = 0;
				foreach $gene (keys %selected_genes){
					if ($selected_genes{$gene} > 1){
						$flag = 1;
						$shared_genes = $shared_genes."$gene|$genes{$gene},";
					}
				}
				chop $shared_genes;
				if ($flag eq 1){
					$comment_general = $comment_general."Shared alleles,";
					$comment_specific = $comment_specific."$shared_genes,";
				}
			}

			if ($count1 > 0 or $count2 > 0){
				$comment_general = $comment_general."Unassigned alleles,";
			}

			chop $comment_general;
			chop $comment_specific;
			
			print OUT_HP "$words[0]\t\t\t\t$comment_general\t$comment_specific\t$primer1_counts{$sample}\t$primer2_counts{$sample}\t$allele_info_for1_selected{$sample}\t$allele_info_for3_selected{$sample}\t$allele_info_for1_discarded{$sample}\t$allele_info_for3_discarded{$sample}\t$subtotal_haplotypes1_final\t$subtotal_haplotypes2_final\t$total_per1\t$total_per2\t$per_freq1_nonhp\t$per_freq2_nonhp\t$contaminant_counts1\t$contaminant_counts2";
			print OUT_HP_GENES "$words[0]";

			foreach $hp (sort {$order_towrite{$b} <=> $order_towrite{$a}} keys %toWrite_hp){
				$origin_hp = $hp;
				$origin_hp =~ s/\*//;
				$origin_hp =~ s/\#//;
				print OUT_HP "\t$hp\t$hp_count1{$origin_hp}\t$per_freq1{$origin_hp}\t$hp_count2{$origin_hp}\t$per_freq2{$origin_hp}";
				print OUT_HP_GENES "\t$toWrite_hp{$hp}";

			}
		}
		else{
			print OUT_HP "$words[0]\t\t\t\tNo data\t\t$primer1_counts{$sample}\t$primer2_counts{$sample}\t$allele_info_for1_selected{$sample}\t$allele_info_for3_selected{$sample}\t$allele_info_for1_discarded{$sample}\t$allele_info_for3_discarded{$sample}\t\t\t\t\t\t\t\t";
			print OUT_HP_GENES "$words[0]";
			print OUT_CONT "$words[0]";
			print OUT_NONHP "$words[0]";
		}
		print OUT_HP "\n";
		print OUT_HP_GENES "\n";
		print OUT_CONT "\n";
		print OUT_NONHP "\n";
		# print OUT_BR22 "\n";
	}
}
else{
	die "Cannot open $genes_file\n";
}
