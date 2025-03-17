use Getopt::Long;
use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;

my $samplesheet = "samplesheet.txt";
my $prefix      = "sample";
my $cutoff      = 0;
my $in_fasta    = "fasta/MHCI.fasta";
my $foldchange  = 0;

GetOptions(
    'samplesheet=s' => \$samplesheet,
    'prefix=s'      => \$prefix,
    'cutoff=s'      => \$cutoff,
    'database=s'    => \$in_fasta,
    'fc=i'          => \$foldchange
) or print "Invalid options\n";


my %new         = ();
my %new_samples = ();
my %new_freq1   = ();
my %new_freq2   = ();

my %redundant1 = ();
my %redundant2 = ();
my %for1       = ();
my %for3       = ();

open(IN, "$in_fasta");
while(<IN>) {
	chomp $_;
	if (/^>(.+)/){
		$id = $1;
	}
	else{
		my $sequence = $_;
		if ($id !~ /BoLA-NC/){
			if (length($sequence) eq 318){
				$for3{$id} = substr $sequence, 0, 318;
			}
			elsif (length($sequence) eq 378){
				$for1{$id} = substr $sequence, 0, 378;
			}
			elsif (length($sequence) eq 376){
				$for1{$id} = substr $sequence, 0, 376;
			}
			elsif (length($sequence) eq 407){
				$for1{$id} = substr $sequence, 0, 375;
				$for3{$id} = substr $sequence, 92, 407;
			}
			elsif (length($sequence) eq 408){
				$for1{$id} = substr $sequence, 0, 376;
				$for3{$id} = substr $sequence, 92, 408;
			}
			elsif (length($sequence) eq 410){
				$for1{$id} = substr $sequence, 0, 378;
				$for3{$id} = substr $sequence, 92, 410;
			}
		}
	}
}

foreach my $id1 (keys %for1){
	foreach my $id2 (keys %for1){
		if ($id1 ne $id2 and $for1{$id1} eq $for1{$id2}){
			if (!exists $redundant1{$id1}){
				$redundant1{$id1} = $id2;
			}
			else{
				$redundant1{$id1} = $redundant1{$id1}.",".$id2;
			}
			print "For1 $id1 $id2\n";
		}
	}
}

foreach my $id1 (keys %for3){
	foreach my $id2 (keys %for3){
		if ($id1 ne $id2 and $for3{$id1} eq $for3{$id2}){
			if ($for1{$id1} ne $for1{$id2}){
				if (!exists $redundant2{$id1}){
					$redundant2{$id1} = $id2;
				}
				else{
					$redundant2{$id1} = $redundant2{$id1}.",".$id2;
				}
				print "For3 $id1 $id2\n";
			}
			else{
				# print "$for1{$id1}\t$for1{$id2}\n";
				print "$id1 and $id2 are same in both pcr\n";
			}
		}
	}
}

open(SUMMARY, ">$prefix.summary.MHCI.tsv");
print SUMMARY "sample\tTotal reads\tQuality trimmed\tOverlapping reads\tPrimer1 reads\tPrimer1  Number of clusters\tPrimer1  Singletons\tPrimer1  Reads: size difference\tPrimer1  Reads: Chimeras\tPrimer1  Reads: 1bp variants\tPrimer1  Reads: Low counts\tPrimer1  Reads: ambiguous calls\tPrimer1  Filtered reads\tPrimer1  Clusters: size difference\tPrimer1  Clusters: Chimeras\tPrimer1  Clusters: 1bp variants\tPrimer1  Clusters: Low counts\tPrimer1  Clusters: ambiguous calls\tPrimer1  Filtered clusters\tPrimer1  Reads: Mapped on database\tPrimer1  Reads: Discarded mapped on database\tPrimer1  Reads: Novel\tPrimer1  Reads: Trimmed end\tPrimer1  Clusters: Mapped on database\tPrimer1  Clusters: Discarded mapped on database\tPrimer1  Clusters: Novel\tPrimer1  Clusters: Trimmed end\tPrimer2 reads\tPrimer2 Number of clusters\tPrimer2 Singletons\tPrimer2 Reads: size difference\tPrimer2 Reads: Chimeras\tPrimer2 Reads: 1bp variants\tPrimer2 Reads: Low counts\tPrimer2 Reads: ambiguous calls\tPrimer2 Filtered reads\tPrimer2 Clusters: size difference\tPrimer2 Clusters: Chimeras\tPrimer2 Clusters: 1bp variants\tPrimer2 Clusters: Low counts\tPrimer2 Clusters: ambiguous calls\tPrimer2 Filtered clusters\tPrimer2 Reads: Mapped on database\tPrimer2 Reads: Discarded mapped on database\tPrimer2 Reads: Novel\tPrimer2 Reads: Trimmed end\tPrimer2 Clusters: Mapped on database\tPrimer2 Clusters: Discarded mapped on database\tPrimer2 Clusters: Novel\tPrimer2 Clusters: Trimmed end\n";
#Get in all new sequences...
open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	my @words             = split("\t", $_);
	my $sample            = $words[0];
	$sickle_log           = "results/$sample/$sample.trimming_by_sickle.log";
	my $total_reads       = 0;
	my $trimmed_paired    = 0;
	my $overlapped_paired = 0;

	#To get total read counts and quality trimmed counts
	if (-e $sickle_log){
		open (IN, "$sickle_log");
		LOOP1: while (<IN>){
			chomp $_;
			if (/^FastQ paired records kept: (\d+) \((\d+) pairs/){
				$trimmed_paired = $2;
				$total_reads = $total_reads + $2;
			}
			elsif (/^FastQ single records kept: (\d+) \(from PE1: (\d+), from PE2: (\d+)/){
				$total_reads = $total_reads + $1;
			}
			elsif (/^FastQ paired records discarded: (\d+) \((\d+) pairs/){
				$total_reads = $total_reads + $2;
			}
		}
		close(IN);
	}
	else{
		print "Cannot open $sickle_log\n";
	}

	#To get overlapped counts
	$flash_log = "results/$sample/$sample.overlap_by_flash.log";
	if (-e $flash_log){
		open (IN, "$flash_log");
		LOOP2: while (<IN>){
			chomp $_;
			if (/Combined pairs:\s+(\d+)/){
				$overlapped_paired = $1;
				last LOOP2;
			}
		}
		close(IN);	
	}
	else{
		print "Cannot find $flash_log\n";
	}

	print SUMMARY "$sample\t$total_reads\t$trimmed_paired\t$overlapped_paired";


	foreach $primer ("For1Rev2", "For3Rev1"){
		my $primers_overlap     = 0;
		my $clusters            = 0;
		my $singles             = 0;
		my $lenDiff_reads       = 0;
		my $chimera_reads       = 0;
		my $variants_reads      = 0;
		my $low_reads           = 0;
		my $gap_reads           = 0;
		my $selected_reads      = 0;
		my $lenDiff             = 0;
		my $chimera             = 0;
		my $variants            = 0;
		my $low                 = 0;
		my $gap                 = 0;
		my $selected            = 0;
		my $mapped_reads        = 0;
		my $discarded_reads     = 0;
		my $new_reads           = 0;
		my $splicevariant_reads = 0;
		my $mapped              = 0;
		my $discarded           = 0;
		my $splicevariant       = 0;
		my $new                 = 0;
		
		$sort_primer = "results/$sample/$sample.$primer.sort.stats.tsv";
		if (-e $sort_primer){
			open (IN, "$sort_primer");
			while (<IN>){
				chomp $_;
				@words           = split("\t", $_);
				$primers_overlap = $words[2];
				$clusters        = $words[3];
				$singles         = $words[4];
			}
			close(IN);
		}
		else{
			print "Cannot find $sort_primer\n";
		}

		$filter_log = "results/$sample/$sample.$primer.clusters.stats.tsv";
		if (-e $filter_log){
			open (IN, "$filter_log");
			while (<IN>){
				chomp $_;
				@words          = split("\t", $_);
				$lenDiff_reads  = $words[1];
				$chimera_reads  = $words[2];
				$variants_reads = $words[3];
				$low_reads      = $words[4];
				$gap_reads      = $words[5];
				$selected_reads = $words[6];
				$lenDiff        = $words[7];
				$chimera        = $words[8];
				$variants       = $words[9];
				$low            = $words[10];
				$gap            = $words[11];
				$selected       = $words[12];
			}
			close(IN);
		}
		else{
			print "Cannot find $filter_log\n";
		}

		$mapping_log = "results/$sample/$sample.$primer.clusters.blast.stats.tsv";
		if (-e $mapping_log){
			open (IN, "$mapping_log");
			while (<IN>){
				chomp $_;
				@words               = split("\t", $_);
				$mapped_reads        = $words[2];
				$discarded_reads     = $words[3];
				$new_reads           = $words[4];
				$splicevariant_reads = $words[5];
				$mapped              = $words[6];
				$discarded           = $words[7];
				$splicevariant       = $words[8];
				$new                 = $words[9];
			}
			close(IN);
		}
		else{
			print "Cannot find $mapping_log\n";
		}

		print SUMMARY "\t$primers_overlap\t$clusters\t$singles\t$lenDiff_reads\t$chimera_reads\t$variants_reads\t$low_reads\t$gap_reads\t$selected_reads\t$lenDiff\t$chimera\t$variants\t$low\t$gap\t$selected\t$mapped_reads\t$discarded_reads\t$new_reads\t$splicevariant_reads\t$mapped\t$discarded\t$new\t$splicevariant";
	}
	print SUMMARY "\n";

	$file = "results/$sample/$sample.mhcI.overlapped.fasta";
	if (-s $file){
		open(IN, "$file");
		while(<IN>) {
			chomp $_;
			if (/^>(.+)/){
				$id = $1;
			}
			else{
				my $sequence = $_;
				# >BoLA-3*002:01|2-1.55039-good|3-0.02166-good
				my @split_id = split(/\|/, $id);
				@freq1 = split ("-", $split_id[1]);
				@freq2 = split ("-", $split_id[2]);
				if ($id =~ /-new/){
					$new{$sequence}++;
					$new_samples{$sequence} = $new_samples{$sequence}."$sample($split_id[1], $split_id[2])\t";
					$new_count{$sequence}++;
					$new_freq1{$sequence} = $new_freq1{$sequence}."$freq1[1],";
					$new_freq2{$sequence} = $new_freq2{$sequence}."$freq2[1],";
				}
			}
		}
	}
	else{
		print "To get new sequences: Cannot find overlapping reads or it is empty: $file\n";
	}
}

#Get information of potential new sequences
my $count = 0;
my %new_id = ();
foreach my $seq (sort {$new{$b} <=> $new{$a}} keys %new){
	$count++;
	my $id        = "$prefix"."_new$count";
	$new_id{$seq} = $id;
	chop $new_freq1{$seq};
	chop $new_freq2{$seq};
	@freq1        = split(",", $new_freq1{$seq});
	@freq2        = split(",", $new_freq2{$seq});
	$total_freq1  = 0;
	$total_freq2  = 0;

	foreach $freq (@freq1){
		$total_freq1 = $total_freq1 + $freq;
	}
	
	foreach $freq (@freq2){
		$total_freq2 = $total_freq2 + $freq;
	}

	$ave1{$seq} = sprintf("%.2f", ($total_freq1/scalar(@freq1)));
	$ave2{$seq} = sprintf("%.2f", ($total_freq2/scalar(@freq2)));
}

my %final_new = ();
open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	my @words = split("\t", $_);
	my $sample = $words[0];

	open (OUT_STAT,      ">summary/$sample.mhcI.stat.tsv"     ) or die "Cannot write $work_dir/$sample.mhcI.stat.txt\n";
	open (OUT_SELECTED,  ">summary/$sample.mhcI.selected.tsv" ) or die "Cannot write $work_dir/$sample.mhcI.selected.txt\n";
	open (OUT_DISCARDED, ">summary/$sample.mhcI.discarded.tsv") or die "Cannot write $work_dir/$sample.mhcI.discarded.txt\n";
	open (OUT_AMBI,      ">summary/$sample.mhcI.ambiguous.tsv") or die "Cannot write $work_dir/$sample.mhcI.discarded.txt\n";
	open (OUT_NC,        ">summary/$sample.mhcI.nc.tsv"       ) or die "Cannot write $work_dir/$sample.mhcI.nc.txt\n";

	print OUT_SELECTED "$sample";
	print OUT_DISCARDED "$sample";
	print OUT_NC "$sample";
	print OUT_AMBI "$sample";
	# print OUT_RESCUED "$sample";

	my $count_nc1        = 0;
	my $count_nc2        = 0;
	my $count_selected1  = 0;
	my $count_discarded1 = 0;
	my $count_rescued1   = 0;
	my $count_rescued2   = 0;
	my $count_selected2  = 0;
	my $count_discarded2 = 0;
	my $count_ambiguous1 = 0;
	my $count_ambiguous2 = 0;
	my %flag_redundant   = ();
	my %all;

	$file = "results/$sample/$sample.mhcI.overlapped.fasta";

	print "\n\nRunning $sample...\n";
	if (-s $file){
		open(IN, "$file");
		$count = 0;

		while(<IN>) {
			chomp $_;
			if (/^>(.+)/){
				$id = $1;
			}
			else{
				my $sequence = $_;
				@info  = split(/\|/, $id);
				@freq1 = split("-", $info[1]);
				@freq2 = split("-", $info[2]);
				$count++;
				$order    {$info[0]} = $count;
				$all      {$info[0]} = $id;
				$all_freq1{$info[0]} = $freq1[0];
				$all_freq2{$info[0]} = $freq2[0];
			}
		}

		foreach $gene (sort {$order{$a} <=> $order{$b}} keys %order){
			if (exists $redundant1{$gene}){
				my @paired_genes = split(",", $redundant1{$gene});
				foreach my $paired_gene (@paired_genes){
					if (exists $all{$paired_gene} and $all_freq1{$paired_gene} eq $all_freq1{$gene} and (($all_freq2{$paired_gene}/$all_freq2{$gene}) > 5)){
						$flag_redundant{$gene} = "For1";
						print "$all{$gene} is ambiguous with $all{$paired_gene} in For1Rev2\n";
					}
				}
			}

			if (exists $redundant2{$gene}){
				my @paired_genes = split(",", $redundant2{$gene});
				foreach my $paired_gene (@paired_genes){
					if (exists $all{$paired_gene} and $all_freq2{$paired_gene} eq $all_freq2{$gene} and (($all_freq1{$paired_gene}/$all_freq1{$gene}) > 5)){
						$flag_redundant{$gene} = "For3";
						print "$all{$gene} is ambiguous with $all{$paired_gene} in For3Rev1\n";
					}
				}
			}
		}

		open(IN, "$file");
		while(<IN>) {
			chomp $_;
			if (/^>(.+)/){
				$id = $1;
			}
			else{
				@info     = split(/\|/, $id);
				@freq1    = split("-", $info[1]);
				@freq2    = split("-", $info[2]);
				$freq1[1] = sprintf("%.2f", $freq1[1]);
				$freq2[1] = sprintf("%.2f", $freq2[1]);
				my $sequence = $_;
				
				if ($id !~ /lenDiff/){
					if (! exists $flag_redundant{$info[0]}){
						if ($info[0] =~ /NC/){
							$count_nc1 = $count_nc1 + $freq1[0];
							$count_nc2 = $count_nc2 + $freq2[0];
							print OUT_NC "\t$id";
							print "Found NC $id\n";
						}
						else{
							if ((($freq1[2] eq "good" and $freq2[2] eq "good") and ($freq1[1] > $cutoff or $freq2[1] > $cutoff)) or (($freq1[2] eq "good" and $freq1[1] > $cutoff and $freq2[0] eq 0) or ($freq2[2] eq "good"  and $freq2[1] > $cutoff and $freq1[0] eq 0))) {
								if (exists $new_id{$sequence}){
									print OUT_SELECTED "\t$new_id{$sequence}|$freq1[0]-$freq1[1]|$freq2[0]-$freq2[1]";
									$final_new{$sequence}++;
									print "New $new_id{$sequence}|$freq1[0]-$freq1[1]|$freq2[0]-$freq2[1]\n";
								}
								else{
									print OUT_SELECTED "\t$info[0]|$freq1[0]-$freq1[1]|$freq2[0]-$freq2[1]";
									print "$id\n";
								}
								$count_selected1 = $count_selected1 + $freq1[0];
								$count_selected2 = $count_selected2 + $freq2[0];
							}
							else{
								if ($info[1] =~ /1bpVariant/ and $info[2] =~ /1bpVariant/){
									@info1 = split(":", $info[1]);
									@info2 = split(":", $info[2]);
									if ($info1[1] < $foldchange and $info2[1] < $foldchange and ($freq1[1] > $cutoff or $freq2[1] > $cutoff)){
										if (exists $new_id{$sequence}){
											$final_new{$sequence}++;
											print OUT_SELECTED "\t$new_id{$sequence}|$info[1]|$info[2]";
										}
										else{
											print OUT_SELECTED "\t$info[0]|$info[1]|$info[2]";
										}
										$count_rescued1 = $count_rescued1 + $freq1[0];
										$count_rescued2 = $count_rescued2 + $freq2[0];
									}
									else{
										if (exists $new_id{$sequence}){
											print OUT_DISCARDED "\t$new_id{$sequence}|$info[1]|$info[2]";
										}
										else{
											print OUT_DISCARDED "\t$id";
										}
										$count_discarded1 = $count_discarded1 + $freq1[0];
										$count_discarded2 = $count_discarded2 + $freq2[0];
									}
								}
								elsif ($info[1] =~ /1bpVariant/ and ($info[2] =~ /good/ or $freq2[0] eq 0) and $freq1[1] > $cutoff){
									@info1 = split(":", $info[1]);
									if ($info1[1] < $foldchange){
										if (exists $new_id{$sequence}){
											$final_new{$sequence}++;
											print OUT_SELECTED "\t$new_id{$sequence}|$info[1]|$info[2]";
										}
										else{
											print OUT_SELECTED "\t$info[0]|$info[1]|$info[2]";
										}
										$count_rescued1 = $count_rescued1 + $freq1[0];
										$count_rescued2 = $count_rescued2 + $freq2[0];
									}
									else{
										if (exists $new_id{$sequence}){
											print OUT_DISCARDED "\t$new_id{$sequence}|$info[1]|$info[2]";
										}
										else{
											print OUT_DISCARDED "\t$id";
										}
										$count_discarded1 = $count_discarded1 + $freq1[0];
										$count_discarded2 = $count_discarded2 + $freq2[0];
									}
								}
								elsif ($info[2] =~ /1bpVariant/ and ($info[1] =~ /good/ or $freq1[0] eq 0) and $freq2[1] > $cutoff){
									@info2 = split(":", $info[2]);
									if ($info2[1] < $foldchange){
										if (exists $new_id{$sequence}){
											$final_new{$sequence}++;
											print OUT_SELECTED "\t$new_id{$sequence}|$info[1]|$info[2]";
										}
										else{
											print OUT_SELECTED "\t$info[0]|$info[1]|$info[2]";
										}
										$count_rescued1 = $count_rescued1 + $freq1[0];
										$count_rescued2 = $count_rescued2 + $freq2[0];
									}
									else{
										if (exists $new_id{$sequence}){
											print OUT_DISCARDED "\t$new_id{$sequence}|$info[1]|$info[2]";
										}
										else{
											print OUT_DISCARDED "\t$id";
										}
										$count_discarded1 = $count_discarded1 + $freq1[0];
										$count_discarded2 = $count_discarded2 + $freq2[0];
									}
								}
								else{
									if (exists $new_id{$sequence}){
										print OUT_DISCARDED "\t$new_id{$sequence}|$info[1]|$info[2]";
										$count_discarded1 = $count_discarded1 + $freq1[0];
										$count_discarded2 = $count_discarded2 + $freq2[0];
									}
									else{
										print OUT_DISCARDED "\t$id";
										$count_rescued1 = $count_rescued1 + $freq1[0];
										$count_rescued2 = $count_rescued2 + $freq2[0];
									}									
								}
							}
						}
					}
					else{
						if ($flag_redundant{$info[0]} eq "For1"){
							$count_ambiguous2 = $count_ambiguous2 + $freq2[0];
							print OUT_AMBI "\t$id";
						}
						elsif ($flag_redundant{$info[0]} eq "For3"){
							$count_ambiguous1 = $count_ambiguous1 + $freq1[0];
							print OUT_AMBI "\t$id";
						}
						else{
							print "Problem: $id\n";
						}
					}
				}
			}	
		}
	}				
	else{
		print "To make final table: Cannot find overlapping reads or it is empty: $file\n";
	}

	print OUT_SELECTED  "\n";
	print OUT_DISCARDED "\n";
	print OUT_NC        "\n";
	print OUT_AMBI      "\n";
	print OUT_STAT      "$sample\t$count_selected1\t$count_selected2\t$count_rescued1\t$count_rescued2\t$count_discarded1\t$count_discarded2\t$count_ambiguous1\t$count_ambiguous2\t$count_nc1\t$count_nc2\n";
}

open ("NEW_FASTA", ">fasta/new.MHCI.fasta");
open ("NEW_INFO", ">fasta/new.MHCI.samples.txt");

foreach my $seq (sort {$final_new{$b} <=> $final_new{$a}} keys %final_new){
	print NEW_FASTA ">$new_id{$seq}\n$seq\n";
	print NEW_INFO "$new_id{$seq}\t$new_count{$seq}\t$ave1{$seq}\t$ave2{$seq}\t$new_samples{$seq}\n";
}

