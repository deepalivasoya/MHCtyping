use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Time::Piece;
use Getopt::Long;

my $sample = "sample";
my $work_dir = "results";
my $prefix = "sample";
my $ampliconSize = 0;
my $cutoff = 0;
my $fc = 0;

GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir,
    'primers=s' => \$prefix,
    'ampliconSize=i'     => \$ampliconSize,
    'cutoff=s'	=> \$cutoff,
    'fc=i'	=> \$fc
) or print "Invalid options\n";

my $localtime = localtime;
print "Filtering: starting at $localtime\n\n";

$blast_details = "$work_dir/$sample.$prefix.clusters.blast.details.tsv";
if ( -e $blast_details ){
	open (IN, "$work_dir/$sample.$prefix.clusters.blast.details.tsv");
	while(<IN>){
		chomp $_;
		@words = split("\t", $_);
		$blast_status{$words[0]} = $words[3];
		if ($words[3] eq "mapped"){
			$blast_ref{$words[0]} = $words[4];
		}
		else{
			$blast_ref{$words[0]} = $words[0];
		}
	}
}
else{
	print "Cannot read $blast_details\n";
}

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


my $count_1bpvariants_reads = 0;
my $count_chimera_reads = 0;
my $count_good_reads = 0;
my $count_gap_reads = 0;
my $count_lendiff_reads = 0;
my $count_mapped_good_reads = 0;
my $count_unmapped_good_reads = 0;
my $count_mapped_chimera_reads = 0;
my $count_unmapped_chimera_reads = 0;
my $count_mapped_1bpvariant_reads = 0;
my $count_unmapped_1bpvariant_reads = 0;


my $count_1bpvariants = 0;
my $count_chimera = 0;
my $count_good = 0;
my $count_gap = 0;
my $count_lendiff = 0;
my $count_mapped_good = 0;
my $count_unmapped_good = 0;
my $count_mapped_chimera = 0;
my $count_unmapped_chimera = 0;
my $count_mapped_1bpvariant = 0;
my $count_unmapped_1bpvariant = 0;

open (INFO, ">$work_dir/$sample.$prefix.clusters.details.tsv") or die "Cannot write $work_dir/$sample.$prefix.clusters.details.tsv\n";
open (FASTA, ">$work_dir/$sample.$prefix.selected.fasta");
my $clusters_file = "$work_dir/$sample.$prefix.clusters.fasta";
if ( -e $clusters_file ){
	my %variants = ();
	my %variants_ids = ();
	my %variants_order = ();
	my %variants_per = ();
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
			$variants_per{$_} = $words[2];
			$lengths{$_} = length($_);
		}
	}
	$total = keys %variants;
	print "$total clusters to check....\n";
	foreach my $seq (sort {$variants_order{$a} <=> $variants_order{$b}} keys %variants){
	
		my $flag_pcr_error = 0;
		my $flag_chimera;
		my $mismatch_counts = 0;
		my $mismatch = "";
		my $mismatch_with = "";
		my $fc = 0;
		my $flag_gap = 0;

		my $id = $variants_ids{$seq};
		print "$id: ";
		my @ids = split(":", $id);

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
				$count_lendiff_reads = $count_lendiff_reads + $variants{$seq};
				$count_lendiff++;
			}
			else{
				LOOP7: foreach $seq1 (sort {$variants{$b} <=> $variants{$a}} keys %variants){
					if ($variants{$seq1} > $variants{$seq}){
						if (length($seq1) eq length($seq)){
							$fc_local = sprintf('%.2f', $variants{$seq1} / $variants{$seq});
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
					$details{$seq} = "1bpVariant\t$mismatch_with".":$fc_local:$mismatch";
					print "1bpVariant\t$mismatch_with".":$fc_local:$mismatch\n";
					if ($blast_status{$id} eq "mapped"){
						print FASTA ">$blast_ref{$id}|$variants{$seq}|$variants_per{$seq}|1bpVariant:$fc_local\n$seq\n";
						$count_mapped_1bpvariant_reads = $count_mapped_1bpvariant_reads + $variants{$seq};
						$count_mapped_1bpvariant++;
					}
					elsif ($variants_per{$seq} >= $cutoff and $fc_local <= $fc){
						print FASTA ">$sample_new-$blast_ref{$id}|$variants{$seq}|$variants_per{$seq}|1bpVariant:$fc_local\n$seq\n";
						$count_unmapped_1bpvariant_reads = $count_unmapped_1bpvariant_reads + $variants{$seq};
						$count_unmapped_1bpvariant++;
					}
					else{
						$count_1bpvariants++;
						$count_1bpvariants_reads = $count_1bpvariants_reads + $variants{$seq};
					}
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
												if ($blast_status{$id} eq "mapped"){
													print FASTA ">$blast_ref{$id}|$variants{$seq}|$variants_per{$seq}|chimera\n$seq\n";
													$count_mapped_chimera_reads = $count_mapped_chimera_reads + $variants{$seq};
													$count_mapped_chimera++;
												}
												elsif ($variants_per{$seq} >= $cutoff and $fc_local <= $fc){
													print FASTA ">$sample_new-$blast_ref{$id}|$variants{$seq}|$variants_per{$seq}|chimera\n$seq\n";
													$count_unmapped_chimera_reads = $count_unmapped_chimera_reads + $variants{$seq};
													$count_unmapped_chimera++;
												}
												else{
													$count_chimera++;
													$count_chimera_reads = $count_chimera_reads + $variants{$seq};
												}
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
						$details{$seq} = "good\t-";
						print "good: $variants_ids{$seq1}-$variants_ids{$seq2}\n";
						if ($blast_status{$id} eq "mapped"){
							print FASTA ">$blast_ref{$id}|$variants{$seq}|$variants_per{$seq}|good\n$seq\n";
							$count_mapped_good_reads = $count_mapped_good_reads + $variants{$seq};
							$count_mapped_good++;
						}
						elsif ($variants_per{$seq} >= $cutoff){
							print FASTA ">$sample_new-$blast_ref{$id}|$variants{$seq}|$variants_per{$seq}|good\n$seq\n";
							$count_unmapped_good_reads = $count_unmapped_good_reads + $variants{$seq};
							$count_unmapped_good++;
						}
						else{
							$count_good++;
							$count_good_reads = $count_good_reads + $variants{$seq};
						}
					}
				}
			}
		}
		
		print INFO "$id\t$variants{$seq}\t$ids[2]\t$lengths{$seq}\t$blast_status{$id}\t$blast_ref{$id}\t$details{$seq}\t$seq\n";
	}
}
else{
	print "Cannot read $clusters_file\n";
}

open (STATS, ">$work_dir/$sample.$prefix.clusters.stats.tsv") or die "Cannot write $work_dir/$sample.$prefix.clusters.stats.tsv\n";
print STATS "$sample\t$count_lendiff_reads\t$count_lendiff";
print STATS "\t$count_gap_reads\t$count_gap";
print STATS "\t$count_chimera_reads\t$count_chimera\t$count_mapped_chimera_reads\t$count_mapped_chimera\t$count_unmapped_chimera_reads\t$count_unmapped_chimera";
print STATS "\t$count_1bpvariants_reads\t$count_1bpvariants\t$count_mapped_1bpvariant_reads\t$count_mapped_1bpvariant\t$count_unmapped_1bpvariant_reads\t$count_unmapped_1bpvariant";
print STATS "\t$count_good_reads\t$count_good\t$count_mapped_good_reads\t$count_mapped_good\t$count_unmapped_good_reads\t$count_unmapped_good\n";
close(STATS);

$localtime = localtime;
print "\n\nFiltering: finished at $localtime\n";









