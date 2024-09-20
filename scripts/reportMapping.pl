use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $sample = "sample";
my $work_dir = "results";
my $prefix = "sample";
my $blast_file = "";
my $filtered_file = "";
my $cluster_fasta = "";
my $ampliconSize = 0;

GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir,
    'primers=s' => \$prefix,
    'ampliconSize=i'     => \$ampliconSize,
    'blast_file=s'     => \$blast_file,
    'cluster_fasta=s'     => \$cluster_fasta
) or print "Invalid options\n";

my $localtime = localtime;
print "Checking blast: starting at $localtime\n\n";

my %counts = ();
my %per = ();

my $id = "";
my %refs = ();
my %identity = ();
# my %error = ();
# my %gap = ();
my %match = ();
my %q_len = ();
my %q_start = ();
my %q_end = ();
my %r_len = ();
my %r_start = ();
my %r_end = ();

# 0 1:4652:24.28102
# 1 BoLA-DQA*021:02
# 2 100.000
# 3 261
# 4 261
# 5 1
# 6 261
# 7 261
# 8 1
# 9 261
# 10 0
# 11 0
# 12 1.60e-139
# 13 483

if (-e $blast_file) {
	open (BLAST, "$blast_file");
	print "Looking for known Alleles from database...\n";
	while (<BLAST>){
		chomp $_;
		my @words = split("\t", $_);
		if ($words[0] ne $id and $words[2] =~ /100.00/ and ($words[3] eq $words[4] or $words[3] eq $words[7])){
			$id = $words[0];
			# print "$id\t$words[2]\t$words[1]\t$words[3]\t$words[4]\n";
			$refs{$id} = $words[1];
			$identity{$id} = $words[2];
		    # $error{$id} = $words[10];
		    # $gap{$id} = $words[11];
		    $match{$id} = $words[3];
		    $q_len{$id} = $words[4];
		    $q_start{$id} = $words[5];
		    $q_end{$id} = $words[6];
		    $r_len{$id} = $words[7];
		    $r_start{$id} = $words[8];
		    $r_end{$id} = $words[9];
		}
		elsif ($identity{$id} eq $words[2] and ($words[3] eq $words[4] or $words[3] eq $words[7])){
			$refs{$id} = $refs{$id}.",".$words[1];
			# print "Double mapping for $id: $words[1]\n";
		}
	}
	close(BLAST);
}
else{
	print "Cannot find $blast_file\n\n";
}
print "\n\n";

open (LOG, ">$work_dir/$sample.$prefix.clusters.blast.details.tsv") or print "Cannot write $work_dir/$sample.$prefix.clusters.blast.details.tsv\n";

my $mapped = 0;
my $mapped_reads = 0;
my $unmapped = 0;
my $unmapped_reads = 0;
my %references = ();

if (-e $blast_file) {
	open(IN, "$cluster_fasta");
	while(<IN>){
		chomp $_;
		if (/^>(\S+)/){
			$id = $1;
			@words = split(":", $id);
			$counts{$id} = $words[1];
			$per{$id} = $words[2];
		}
		else{
			$sequence = $_;
			if (exists $refs{$id}){
				@order_ref = split(",", $refs{$id});
				@ordered_ref = sort { $a cmp $b } @order_ref;
				$to_write_ref = "";
				foreach $ref (@ordered_ref){
					# print "$to_write_ref\n";
					$to_write_ref = $to_write_ref."$ref,";
				}
				chop $to_write_ref;

				if (! exists $references{$to_write_ref}){
					$references{$to_write_ref} = $id;
					print LOG "$id\t$counts{$id}\t$per{$id}\tmapped\t$to_write_ref\t$identity{$id}\t$match{$id}\t$q_start{$id}-$q_end{$id}:$q_len{$id}\t$r_start{$id}-$r_end{$id}:$r_len{$id}\t$sequence\n";
					$mapped_reads = $mapped_reads + $counts{$id};
					$mapped++;
				}
				else{
					print LOG "$id\t$counts{$id}\t$per{$id}\tunmapped\t$to_write_ref\t$identity{$id}\t$match{$id}\t$q_start{$id}-$q_end{$id}:$q_len{$id}\t$r_start{$id}-$r_end{$id}:$r_len{$id}\t$sequence\n";
					$unmapped_reads = $unmapped_reads + $counts{$id};
					$unmapped++;
				}
			}
			else{
				print LOG "$id\t$counts{$id}\t$per{$id}\tunmapped\t-\t-\t-\t-\t-\t$sequence\n";
				$unmapped_reads = $unmapped_reads + $counts{$id};
				$unmapped++;
			}
		}
	}
	close(IN);
}

open (LOG, ">$work_dir/$sample.$prefix.clusters.blast.stats.tsv") or print "Cannot write $work_dir/$sample.$prefix.clusters.blast.details.tsv\n";
print LOG "$sample\t$prefix\t$mapped_reads\t$unmapped_reads\t$mapped\t$unmapped\n";
close(LOG);

print "Mapped: $mapped_reads ($mapped)\n";
print "Unmapped mapped: $unmapped_reads ($unmapped)\n";


