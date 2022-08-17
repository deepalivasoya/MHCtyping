use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Time::Piece;
use Getopt::Long;

my $sample = "sample";
my $work_dir = "results";
my $prefix = "sample";
my $ampliconSize = 0;
my $noOfClusters = 1000;

GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir,
    'primers=s' => \$prefix,
    'ampliconSize=i'     => \$ampliconSize
) or print "Invalid options\n";


my $localtime = localtime;
print "Filtering: starting at $localtime\n\n";

#chimera detection algorithm
sub possible_chimaera {
	my $sequence1 = $_[0];
	my $sequence2 = $_[1];
	my $candidate = $_[2];
	$candidate = uc($candidate);
	$sequence1 = uc($sequence1);
	$sequence2 = uc($sequence2);
	my %base1 = ();
	my %base2 = ();

	for (my $i = 0; $i < length($sequence1); $i++) {
		$base1{$i} = substr($sequence1, $i, 1 );
	}
	for (my $i = 0; $i < length($sequence2); $i++) {
		$base2{$i} = substr($sequence2, $i, 1 );
	}
	$len_diff2 = length($candidate) - length($sequence2);

	my $stop_base = 0;
	my $stop_base1 = 0;

	LOOP1: for (my $i = 0; $i < length($candidate); $i++) {
		my $base = substr($candidate, $i, 1);
		if ($base1{$i} ne $base) {
			$stop_base = $i;
			last LOOP1;
		}
	}
	LOOP2: for (my $i = (length($candidate) - 1); $i >= 0; $i--) {
		my $base = substr($candidate, $i, 1);
		if ($base2{$i - $len_diff2} ne $base) {
			$stop_base1 = $i;
			last LOOP2;
		}
	}
	if ($stop_base1 > 0 and $stop_base > 0 and $stop_base1 < $stop_base){
		return 1;
	}
	else{
		return 0;
	}
}


my $count_variants_reads = 0;
my $count_chimera_reads = 0;
my $count_selected_reads = 0;
my $count_low_reads = 0;
my $count_gap_reads = 0;
my $count_variants = 0;
my $count_chimera = 0;
my $count_selected = 0;
my $count_low = 0;
my $count_gap = 0;
my $reads_len = 0;
my $counts_len = 0;
open (INFO, ">$work_dir/$sample.$prefix.clusters.details.csv") or die "Cannot write $work_dir/$sample.$prefix.clusters.details.csv\n";

print INFO "Cluster\tRead_Counts\tPercentage\tLength\tVariant_type\tDetail\tSequence\n";
my $clusters_file = "$work_dir/$sample.$prefix.clusters.fasta";
if ( -e $clusters_file ){
	my %variants = ();
	my %variants_ids = ();
	my %variants_order = ();
	my %lengths = ();

	open (IN, "$clusters_file");
	LOOP1: while(<IN>){
		chomp $_;
		if (/^>(\S+)/){
			$id = $1;
		}
		else{
			@words = split(/:/, $id);
			$variants{$_} = $words[1];
			$variants_ids{$_} = $id;
			$variants_order{$_} = $words[0];
			$lengths{$_} = length($_);
		}
	}
	$total = keys %variants;
	print "$total clusters to check....\n";
	foreach my $seq (sort {$variants_order{$a} <=> $variants_order{$b}} keys %variants){
	
		my $flag_pcr_error = 0;
		my $flag_pcr_error1 = 0;
		my @to_find_ave_len;
		my $flag_chimera;
		my $mismatch_counts = 0;
		my $mismatch = "";
		my $mismatch_with = "";
		my $fc = 0;
		my $flag_gap = 0;

		my $id = $variants_ids{$seq};
		print "$id: ";
		my @ids = split(":", $id);
		if ($ids[0] <= $noOfClusters){
			my @nucl = split //, $seq;
			LOOPN: for($i = 0; $i < scalar(@nucl); $i++){
				if ($nucl[$i] !~ /[A|T|C|G]/i){
					$flag_gap = 1;
					last LOOPN;
				}
			}

			if ($flag_gap == 1){ #Check if sequence has N 
				print "N in sequence\n";
				$details{$seq} = "N\t-";
				$count_gap_reads = $count_gap_reads + $variants{$seq};
				$count_gap++;
			}
			else{
				if (length($seq) > ($ampliconSize + 9) or length($seq) < ($ampliconSize - 9)){
					$diff = $ampliconSize - length($seq);
					print "Len difference $ampliconSize:$diff\n";
					$details{$seq} = "lenDiff\t$ampliconSize:$diff";
					$reads_len = $reads_len + $variants{$seq};
					$counts_len++;
				}
				else{
					LOOP7: foreach $seq1 (sort {$variants{$b} <=> $variants{$a}} keys %variants){
						if ($variants{$seq1} > $variants{$seq}){
							if (length($seq1) eq length($seq)){
								$fc = sprintf('%.2f', $variants{$seq1} / $variants{$seq});
								$mismatch_with = $variants_ids{$seq1};
								# print "PCR of $mismatch_with\n";
								$mismatch_counts = 0;
								$mismatch = undef;
								my @nucl = split //, $seq;
								my @nucl1 = split //, $seq1;
								for (my $i = 0; $i < scalar(@nucl); $i++){
									if ( $nucl[$i] ne $nucl1[$i] ){
										$mismatch_counts++;
										my $j = $i + 1;
										$mismatch = "$j/$nucl[$i]->$nucl1[$i]";
									}
								}
								if ($mismatch_counts eq 1){
									last LOOP7;
								}
							}
						}
						else{
							last LOOP7;
						}
					}
					# print "$mismatch_counts\n";
					if ($mismatch_counts eq 1){
						$flag_pcr_error = 1;
						$details{$seq} = "1bpVariant\t$mismatch_with".":$fc:$mismatch";
						print "1bpVariant\t$mismatch_with".":$fc:$mismatch\n";
						$count_variants++;
						$count_variants_reads = $count_variants_reads + $variants{$seq};
					}

					if ($flag_pcr_error eq 0){
						$flag_chimera = 0;
						LOOP6: foreach $seq1 (sort {$variants{$b} <=> $variants{$a}} keys %variants){
							if ($seq ne $seq1){
								if ($variants{$seq1} > $variants{$seq}){
									LOOP8: foreach $seq2 (sort {$variants{$b} <=> $variants{$a}} keys %variants){
										if ($seq ne $seq2 and $seq1 ne $seq2){
											if ($variants{$seq2} > $variants{$seq}){
												$check = possible_chimaera($seq1, $seq2, $seq);
												if ($check == 1){
													$flag_chimera = 1;
													$details{$seq} = "Chimera\t$variants_ids{$seq1}-$variants_ids{$seq2}";
													print "chimera: $variants_ids{$seq1}-$variants_ids{$seq2}\n";
													$count_chimera++;
													$count_chimera_reads = $count_chimera_reads + $variants{$seq};
													last LOOP6;
												}
											}
											else{
												last LOOP8;
											}
										}
									}
								}
								else{
									last LOOP6;
								}
							}
						}
						if ($flag_chimera == 0){
							if ($flag_gap eq 1){
								$count_gap_reads = $count_gap_reads + $variants{$seq};
								$count_gap++;
							}
							print "Good\n";
							$details{$seq} = "good\t-";
							$count_selected++;
							$count_selected_reads = $count_selected_reads + $variants{$seq};
						}
					}
				}
			}
			print INFO "$id\t$variants{$seq}\t$ids[2]\t$lengths{$seq}\t$details{$seq}\t$seq\n";
		}
		else{
			$count_low++;
			$count_low_reads = $count_low_reads + $variants{$seq};
			print INFO "$id\t$variants{$seq}\t$ids[2]\t$lengths{$seq}\tremoved\t-\t$seq\n";
		}

	}
}
else{
	print "Cannot read $clusters_file\n";
}

open (STATS, ">$work_dir/$sample.$prefix.clusters.stats.csv") or die "Cannot write $work_dir/$sample.$prefix.clusters.stats.csv\n";
print STATS "$sample\t$reads_len\t$count_chimera_reads\t$count_variants_reads\t$count_low_reads\t$count_gap_reads\t$count_selected_reads\t$counts_len\t$count_chimera\t$count_variants\t$count_low\t$count_gap\t$count_selected\n";
close(STATS);


$localtime = localtime;
print "\n\nFiltering: finished at $localtime\n";









