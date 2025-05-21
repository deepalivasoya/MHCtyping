use Getopt::Long;
use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;

my $samplesheet = "samplesheet.txt";
my $prefix      = "sample";
my $primer      = "";
my $cutoff      = 0;
my $fc          = 0;

GetOptions(
    'samplesheet=s' => \$samplesheet,
    'prefix=s'      => \$prefix,
    'primer=s'      => \$primer,
    'cutoff=i'      => \$cutoff,
    'fc=i'          => \$fc
) or print "Invalid options\n";


my %new         = ();
my %new_samples = ();
my %new_freq1   = ();
my %final_new   = ();

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	my @words  = split("\t", $_);
	my $sample = $words[0];
	chomp $sample;

	my $file = "results/$sample/$sample.$primer.selected.fasta";
	
	if (-s $file){
		open(IN, "$file");

		while(<IN>) {
			chomp $_;
			if (/^>(.+)/){
				$id = $1;
			}
			else{
				my $sequence = $_;
				# >37043M-new;11-0.30598;good
				my @split_id = split(";", $id);
				@freq1 = split ("-", $split_id[1]);
				if ($id =~ /$sample-new/){
					$new{$sequence}++;
					$new_samples{$sequence} = $new_samples{$sequence}."$sample($split_id[1]-$split_id[2])\t";
					$new_count{$sequence}++;
					$new_freq1{$sequence} = $new_freq1{$sequence}."$freq1[1],";
				}
			}
		}
	}
	else{
		print "To get new sequences: Cannot find overlapping reads or it is empty: $file\n";
	}
}


my $count  = 0;
my %new_id = ();

foreach my $seq (sort {$new{$b} <=> $new{$a}} keys %new){
	$count++;
	my $id        = "$prefix.$primer"."_new$count";
	$new_id{$seq} = $id;
	chop $new_freq1{$seq};
	@freq1        = split(",", $new_freq1{$seq});
	$total_freq1  = 0;

	foreach $freq (@freq1){
		$total_freq1 = $total_freq1 + $freq;
	}

	$ave1{$seq} = sprintf("%.2f", ($total_freq1/scalar(@freq1)));
}

open(SUMMARY, ">$prefix.summary.$primer.tsv");
print SUMMARY "Sample\t$primer:Total_reads\t$primer:Trimmed_reads\t$primer:Overlapped_reads\t$primer:Primer_found\t$primer:No_of_clusters\t$primer:singletons\t$primer:lenDiff_reads\t$primer:chimera_reads\t$primer:variants_reads\t$primer:low_reads\t$primer:gap_reads\t$primer:selected_reads\t$primer:lenDiff_clusters\t$primer:chimera_clusters\t$primer:variants_clusters\t$primer:low_clusters\t$primer:gap_clusters\t$primer:selected_clusters\t$primer:mapped_reads\t$primer:discarded_reads\t$primer:new_reads\t$primer:splicevariant_reads\t$primer:mapped_clusters\t$primer:discarded_clusters\t$primer:splicevariant_clusters\t$primer:new_clusters\t$primer:count_selected\t$primer:count_rescued\t$primer:count_discarded\n";
#Get in all new sequences...

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	my @words             = split("\t", $_);
	my $sample            = $words[0];
	chomp $sample;
	$sickle_log           = "results/$sample/$sample.trimming_by_sickle.log";
	my $total_reads       = 0;
	my $trimmed_paired    = 0;
	my $overlapped_paired = 0;

	if (-e $sickle_log){
		open (IN, "$sickle_log");
		LOOP1: while (<IN>){
			chomp $_;
			# FastQ paired records kept: 18604 (9302 pairs)
			# FastQ single records kept: 1158 (from PE1: 622, from PE2: 536)
			# FastQ paired records discarded: 4168 (2084 pairs)
			# FastQ single records discarded: 1158 (from PE1: 536, from PE2: 622)
			if (/^FastQ paired records kept: (\d+) \((\d+) pairs/){
				$trimmed_paired = $2;
				$total_reads    = $total_reads + $2;
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

	$flash_log = "results/$sample/$sample.overlap_by_flash.log";
	if (-e $sickle_log) {
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
	if (-e $sort_primer) {
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

	print SUMMARY "\t$primers_overlap\t$clusters\t$singles\t$lenDiff_reads\t$chimera_reads\t$variants_reads\t$low_reads\t$gap_reads\t$selected_reads\t$lenDiff\t$chimera\t$variants\t$low\t$gap\t$selected\t$mapped_reads\t$discarded_reads\t$new_reads\t$splicevariant_reads\t$mapped\t$discarded\t$splicevariant\t$new\t";

	open (OUT_SELECTED,  ">summary/$sample.$primer.selected.tsv" ) or die "Cannot write $work_dir/$sample.$primer.selected.txt\n";
	open (OUT_DISCARDED, ">summary/$sample.$primer.discarded.tsv") or die "Cannot write $work_dir/$sample.$primer.discarded.txt\n";

	print OUT_SELECTED  "$sample";
	print OUT_DISCARDED "$sample";

	my $count_selected1  = 0;
	my $count_discarded1 = 0;
	my $count_rescued    = 0;

	$file = "results/$sample/$sample.$primer.selected.fasta";
	print "\n\nRunning $sample $primer...\n";
	if (-s $file){
		open(IN, "$file");
		$count = 0;
		while(<IN>) {
			chomp $_;
			if (/^>(.+)/){
				$id = $1;
			}
			else {
				my $sequence = $_;
				@info        = split(/\;/, $id);
				@freq1       = split("-", $info[1]);
				$freq1[1]    = sprintf("%.2f", $freq1[1]);

				if ($info[2] eq "good" and $freq1[1] > $cutoff) {
					if (exists $new_id{$sequence}) {
						print OUT_SELECTED "\t$new_id{$sequence}|$freq1[0]-$freq1[1]";
						$final_new{$sequence}++;
					}
					else {
						print OUT_SELECTED "\t$info[0]|$freq1[0]-$freq1[1]";
					}
					$count_selected1 = $count_selected1 + $freq1[0];
				}
				else {
					if ($info[2] =~ /1bpVariant/){
						@info1 = split(":", $info[2]);
						# print "$sample $id: $info[2] $info1[0] $info1[1] $fc \n";
						if ($info1[1] < $fc and $freq1[1] > $cutoff){
							print "rescued $sample $id: $info[2]\n";
							if (exists $new_id{$sequence}){
								$final_new{$sequence}++;
								print OUT_SELECTED "\t$new_id{$sequence}|$info[1]-$info[2]";
							}
							else {
								print OUT_SELECTED "\t$info[0]|$info[1]-$info[2]";
							}
							$count_rescued = $count_rescued + $freq1[0];
						}
						elsif($info1[1] < 5 and $freq1[1] > 1) {
							if (exists $new_id{$sequence}) {
								print OUT_DISCARDED "\t$new_id{$sequence}|$info[1]-$info[2]";
							}
							else {
								print OUT_DISCARDED "\t$info[0]|$info[1]-$info[2]";
							}
							$count_discarded1 = $count_discarded1 + $freq1[0];
						}
					}
					elsif($info[2] !~ /lenDiff/ and $freq1[1] > 1) {
						if (exists $new_id{$sequence}) {
							print OUT_DISCARDED "\t$new_id{$sequence}|$info[1]-$info[2]";
						}
						else {
							print OUT_DISCARDED "\t$info[0]|$info[1]-$info[2]";
						}
						$count_discarded1 = $count_discarded1 + $freq1[0];
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
	print SUMMARY "$count_selected1\t$count_rescued\t$count_discarded1\n";
}

open ("NEW_FASTA", ">fasta/new.$primer.fasta");
open ("NEW_INFO",  ">fasta/new.$primer.samples.txt");

foreach my $seq (sort {$final_new{$b} <=> $final_new{$a}} keys %final_new){
	print NEW_FASTA ">$new_id{$seq}\n$seq\n";
	print NEW_INFO "$new_id{$seq}\t$new_count{$seq}\t$ave1{$seq}\t$new_samples{$seq}\n";
}
