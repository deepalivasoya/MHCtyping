use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Getopt::Long;

my $sample = "sample";
my $work_dir = "results";
my $prefix = "sample";
my $primers = "primer";

GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir,
    'primers=s' => \$prefix,
    'primer_seq=s'     => \$primers
) or print "Invalid options\n";

my $localtime = localtime;
print "Primer sorting: starting at $localtime\n\n";
print "$sample is searched for $prefix primers using $primers sequences\n\n";

#read and initiate primer sequences - forward and reverse complements
my %forward_primer = ();
my %reverse_primer = ();
my %forward_primer_revcomp = ();
my %reverse_primer_revcomp = ();
my $count_for = 0;
my $count_rev = 0;


open (IN, "$primers") or print "Cannot read/find $primers\n";
while(<IN>){
	chomp $_;
	if ($_ =~ /^>(\S+)/){
		$id = $1;
	}
	else{
		$seq = $_;
		$seq_revcomp = "";
		if ($id =~ /for/i){
			$count_for++;
			$forward_primer{$count_for} = $seq;
			$seq_revcomp = reverse $seq;
			$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
			$forward_primer_revcomp{$count_for} = $seq_revcomp;
		} 
		elsif ($id =~ /rev/i){
			$count_rev++;
			$reverse_primer{$count_rev} = $seq;
			$seq_revcomp = reverse $seq;
			$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
			$reverse_primer_revcomp{$count_rev} = $seq_revcomp;
		}
	}
}

print "** Forward primers: Forward and reverse complement sequences **\n";
foreach my $p (sort {$a <=> $b} keys %forward_primer){
	print "$p $forward_primer{$p} $forward_primer_revcomp{$p}\n";
}

print "** Reverse primers: Forward and reverse complement sequences **\n";
foreach my $p (sort {$a <=> $b} keys %reverse_primer){
	print "$p $reverse_primer{$p} $reverse_primer_revcomp{$p}\n";
}
my %unique_sequences = ();
my %primer_group = ();
my $total_primer_reads_flash = 0;
my $total_primer_reads_single = 0;
my $total_unique_seq = 0;
my $id = "";
my $sequence = "";
my $trimmed_seq = "";
my $flag = 0;

sub check_primer{
	$sequence = $_[0];
	$flag = 0;
	LOOP1: foreach my $primer1 (keys %forward_primer){
		if ($sequence =~ /^(\w+|)$forward_primer{$primer1}(\w+|)/) {
			foreach my $primer2 (keys %reverse_primer_revcomp){
				if ($2 =~ /(\w+|)$reverse_primer_revcomp{$primer2}(\w+|)$/) {
					if (exists $unique_sequences{$1}){
						$unique_sequences{$1} = $unique_sequences{$1} + 1;
					}
					else {
						$unique_sequences{$1} = 1;
					}
					$flag = 1;
					$primer_group{"$sample\t$forward_primer{$primer1}\t$reverse_primer_revcomp{$primer2}\t+"}++;
					$trimmed_seq = $1;
					last LOOP1;
				}
			}
		}
	}
	if ($flag eq 0){
		LOOP2: foreach my $primer1 (keys %reverse_primer){
			if ($sequence =~ /^(\w+|)$reverse_primer{$primer1}(\w+|)/) {
				foreach my $primer2 (keys %forward_primer_revcomp){
					if ($2 =~ /(\w+|)$forward_primer_revcomp{$primer2}(\w+|)$/) {
						my $seq_revcomp = reverse $1;
						$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
						if (exists $unique_sequences{$seq_revcomp}){
							$unique_sequences{$seq_revcomp} = $unique_sequences{$seq_revcomp} + 1;
						} else {
							$unique_sequences{$seq_revcomp} = 1;
						}
						$flag = 1;
						$primer_group{"$sample\t$forward_primer_revcomp{$primer2}\t$reverse_primer{$primer1}\t-"}++;
						$trimmed_seq = $seq_revcomp;
						# print "original: $1\nconverted: $seq_revcomp\n";
						last LOOP2;
					}
				}
			}
		}
	}

	if ($flag eq 1){
		return 1;
	}
	else{
		return 0;
	}
}

my $flash_reads = "$work_dir/$sample.extendedFrags.fastq";
# open (OUT_PRIMERS, ">$work_dir/$sample.$prefix.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.fasta\n";
my %flash_len = ();
my $total_primer_reads = 0;

if ( -e $flash_reads ){

	open(IN, "$flash_reads");
	$line = 0;
	while(<IN>){
		chomp $_;
		if ($line eq 0){
			$id = $_;
			$line = 1;
			next;
		}
		if ($line eq 1){
			$sequence = $_;
			my $flag_primer = check_primer($sequence);
			if ($flag_primer eq 1){
				$total_primer_reads++;
				$len = length($trimmed_seq);
				$flash_len{$len}++;
				print OUT_PRIMERS "$id/$sample/extended/$len\n$trimmed_seq\n";
			}
			$line = 2;
			next;
		}
		if ($line eq 2){
			$line = 3;
			next;
		}
		if ($line eq 3){
			$line = 0;
			next;
		}
	}
}
else{
	print "Cannot read $flash_reads.\n";
}
		
print "Writting lengh histogram..."
open (OUT_PRIMERS_LEN, ">$work_dir/$sample.$prefix.seq.len_hist.csv") or die "Cannot write $work_dir/$sample/$sample.$prefix.seq.len_hist.csv\n";
foreach my $len (sort{$a <=> $b} keys %flash_len){
	print OUT_PRIMERS_LEN "$len\t$flash_len{$len}\n";
}

print "Writting unique variants log $work_dir/$sample.primers.info.txt...\n";
open (INFO, ">$work_dir/$sample.$prefix.primers.info.txt") or die "Cannot write $work_dir/$sample.$prefix.primers.info.txt\n";
print INFO "Sample\tForward primer\tReverse primer\tOrientation\tread counts\n";
foreach $primer (sort {$primer_group{$b} <=> $primer_group{$a}} keys %primer_group){
	print INFO "$primer\t$primer_group{$primer}\n";
}
close (INFO);

print "Clustering sequences...\n";
open (OUT_UNIQUE, ">$work_dir/$sample.$prefix.clusters.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.clusters.fasta\n";
my $count_clusters = 0;
my $count_singles = 0;
my $count = 0;

if($total_primer_reads > 0){
	open (SINGLES, ">$work_dir/$sample.$prefix.singletons.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.singletons.fasta\n";
	
	foreach $seq (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
		$count++;
		my $per = sprintf("%.5f", ($unique_sequences{$seq} / $total_primer_reads) * 100);
		if ($unique_sequences{$seq} eq 1){
			$count_singles++;
			print SINGLES ">$count:$unique_sequences{$seq}:$per\n$seq\n";
		}
		else{
			$count_clusters++;
			print OUT_UNIQUE ">$count:$unique_sequences{$seq}:$per\n$seq\n";
		}
	}

	
	close (SINGLES);  	
}
close (OUT_UNIQUE);

print "Total sequences with $prefix primers in overlapped reads: $total_primer_reads\n";
print "Total cluster: $count_clusters\n";
print "Total singletons: $count_singles\n";

open (STAT, ">", "$work_dir/$sample.$prefix.sort.stats.csv");
print STAT "$sample\t$prefix\t$total_primer_reads\t$count_clusters\t$count_singles\n";
close (STAT);


$localtime = localtime;
print "\n\nPrimer sorting: finished at $localtime\n";





