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
my $cutoff = 0;
my $fc = 0;


GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir,
    'primers=s' => \$prefix,
    'blast_file=s'     => \$blast_file,
    'cluster_details=s'	=> \$filtered_file,
    'cluster_fasta=s'     => \$cluster_fasta,
    'cutoff=i'	=> \$cutoff,
    'fc=i'	=> \$fc
) or print "Invalid options\n";

my $localtime = localtime;
print "Checking blast: starting at $localtime\n\n";

my %counts = ();
my %details = ();
my %per = ();

#Initiate cluster sequences...
if (-e $filtered_file) {
	open (IN, "$filtered_file");
	while (<IN>){
		chomp $_;
		my @words = split("\t", $_);
		$details{$words[0]} = $words[4];
		if ($words[4] eq "1bpVariant"){
			@info = split(":", $words[5]);
			$parent{$words[0]} = "$info[0]:$info[1]:$info[2]";
			$details{$words[0]} = $details{$words[0]}.":$info[3]";
		}
		# $len{$words[0]} = $words[3];
		$counts{$words[0]} = $words[1];
		$per{$words[0]} = $words[2];
	}
	close(IN);
}
else{
	print "Cannnot open $filtered_file or it may be empty\n";
}

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

if (-e $blast_file) {
	open (BLAST, "$blast_file");
	print "Looking for known Alleles from database...\n";
	while (<BLAST>){
		chomp $_;
		my @words = split("\t", $_);
		if ($words[0] ne $id and $words[2] =~ /100.00/ and ($words[3] eq $words[4] or $words[3] eq $words[7])){
			$id = $words[0];
			print "$id\t$words[2]\t$words[1]\t$words[3]\t$words[4]\n";
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
			print "Double mapping for $id: $words[1]\n";
		}
	}
	close(BLAST);
}
else{
	print "Cannot find $blast_file\n\n";
}
print "\n\n";

open (LOG, ">$work_dir/$sample.$prefix.clusters.blast.details.csv") or print "Cannot write $work_dir/$sample.$prefix.clusters.blast.details.csv\n";
open (FASTA, ">$work_dir/$sample.$prefix.clusters.selected.fasta") or print "Cannot write $work_dir/$sample.$prefix.clusters.selected.fasta\n";

my $mapped = 0;
my $mapped_reads = 0;
my $new = 0;
my $new_reads = 0;
my $discarded = 0;
my $discarded_reads = 0;
my $splicevariant = 0;
my $splicevariant_reads = 0;
my %references = ();

if (-e $blast_file) {
	open(IN, "$cluster_fasta");
	while(<IN>){
		chomp $_;
		if (/^>(\S+)/){
			$id = $1;
		}
		else{
			$sequence = $_;
			if ($details{$id} !~ /lenDiff/ or $details{$id} eq "N"){
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
						print "$id\tmapped\t$details{$id}\t$to_write_ref\n";
						print LOG "$id\t$counts{$id}\t$per{$id}\t$details{$id}\tmapped\t$to_write_ref\t$identity{$id}\t$match{$id}\t$q_start{$id}-$q_end{$id}:$q_len{$id}\t$r_start{$id}-$r_end{$id}:$r_len{$id}\t$sequence\n";
						print FASTA ">$to_write_ref;$counts{$id}-$per{$id};$details{$id}\n$sequence\n";
						if ($per{$id} > $cutoff and $details{$id} eq "good"){
							$mapped_reads = $mapped_reads + $counts{$id};
							$mapped++;
						}
						else{
							$discarded_reads = $discarded_reads + $counts{$id};
							$discarded++;
						}
					}
					else{
						print "$id\tsplice variant\t$details{$id}\t$to_write_ref\n";
						print "$id is splice variant of $references{$to_write_ref}\n";
						print LOG "$id\t$counts{$id}\t$per{$id}\t$details{$id}\tsplice variant\t$to_write_ref\t$identity{$id}\t$match{$id}\t$q_start{$id}-$q_end{$id}:$q_len{$id}\t$r_start{$id}-$r_end{$id}:$r_len{$id}\t$sequence\n";
						$splicevariant_reads = $splicevariant_reads + $counts{$id};
						$splicevariant++;
					}	
				}
				else{
					if ($per{$id} > $cutoff){
						print "$id\tnew\t$details{$id}\n";
						print LOG "$id\t$counts{$id}\t$per{$id}\t$details{$id}\tnew\t$sample-new-$id\t-\t-\t-\t-\t-\t$sequence\n";
						if ($details{$id} =~ /1bpVariant/){
							my $p = $parent{$id};
							my @de = split(":", $details{$id});
							if ($details{$p} eq "good" and $de[1] <= $fc){
								print FASTA ">$sample-new;$counts{$id}-$per{$id};$details{$id}\n$sequence\n";
								$new_reads = $new_reads + $counts{$id};
								$new++;
							}
						}
						elsif($details{$id} !~ /Chimera/){
							print FASTA ">$sample-new;$counts{$id}-$per{$id};$details{$id}\n$sequence\n";
							$new_reads = $new_reads + $counts{$id};
							$new++;
						}
					}
				}
			}
		}
	}
	close(IN);
}

open (LOG, ">$work_dir/$sample.$prefix.clusters.blast.stats.csv") or print "Cannot write $work_dir/$sample.$prefix.clusters.blast.details.csv\n";
print LOG "$sample\t$prefix\t$mapped_reads\t$discarded_reads\t$new_reads\t$splicevariant_reads\t$mapped\t$discarded\t$splicevariant\t$new\n";
close(LOG);

print "Mapped: $mapped_reads ($mapped)\n";
print "Discarded mapped: $discarded_reads ($discarded)\n";
print "Splicevariant: $splicevariant_reads ($splicevariant)\n";
print "Novel: $new_reads ($new)\n";





