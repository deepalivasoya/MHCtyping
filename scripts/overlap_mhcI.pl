use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Time::Piece;
use Getopt::Long;

my $sample = "sample";
my $primer1_file = "file1";
my $primer2_file = "file2";

GetOptions(
    'sample=s'    => \$sample,
    'primer1=s' => \$primer1_file,
    'primer2=s'     => \$primer2_file
) or print "Invalid options\n";

my $work_dir = "results/$sample";
print "\n\n***** Running $sample *****\n";

my %primer1;
my %primer2;
my $count_primer1 = 0;
my $count_primer2 = 0;
my $flag1 = 0;

if (-s $primer1_file){
	$flag1 = 1;
	open(IN, "$primer1_file");
	while(<IN>){
		chomp $_;
		if (/^>(\S+)/){
			$id = $1;
		}
		else{
			$sequence = $_;
			$count_primer1++;
			$primer1{$sequence} = $id;
			$primer1_order{$sequence} = $count_primer1;
		}
	}
}
else{
	print "Cannot find selected variants fasta file or it is empty: $primer1_file\n";
}

my $flag2 = 0;
if (-s $primer2_file){
	$flag2 = 1;
	open(IN, "$primer2_file");
	while(<IN>){
		chomp $_;
		if (/^>(\S+)/){
			$id = $1;
		}
		else{
			$sequence = $_;
			$count_primer2++;
			$primer2{$sequence} = $id;
			$primer2_order{$sequence} = $count_primer2;
		}
	}
}
else{
	print "Cannot find selected variants fasta file or it is empty: $primer2_file\n";
}

open(MERGED_OUT, ">$work_dir/$sample.mhcI.overlapped.fasta") or die "Cannot write $work_dir/$sample.mhcI.overlapped.fasta\n";
open(MERGED_INFO_OUT, ">$work_dir/$sample.mhcI.overlapped.info.tsv") or die "Cannot write $work_dir/$sample.mhcI.overlapped.info.tsv\n";

my %check1 = ();
my %check2 = ();
my $flag_overlapped = 0;
my $new_id = "";
my $overlapped = "";

if ($flag1 eq 1 and $flag2 eq 1){
	my $count = 0;
	foreach my $seq1 (sort {$primer1_order{$a} <=> $primer1_order{$b}} keys %primer1){
		my $full_id1 = $primer1{$seq1};
		# print "\n$full_id1 is checking in For3Rev1\n";
		my @info1 = split(/\|/, $full_id1);
		my $overlapping1 = substr $seq1, 92;
		foreach my $seq2 (sort {$primer2_order{$a} <=> $primer2_order{$b}} keys %primer2){
			my $full_id2 = $primer2{$seq2};
			$flag_overlapped = 0;
			# print "\n$full_id2 is checking in For3Rev1\n";
			my @info2 = split(/\|/, $full_id2);
			my $overlapping2 = substr $seq2, 0, -32;
			if ($overlapping1 eq $overlapping2){
				$flag_overlapped = 1;
				$overlapped = $seq1.(substr $seq2, -32);
				my $id_check = 0;
				my @ids1 = split(/,/, $info1[0]);
				my @ids2 = split(/,/, $info2[0]);
				# print "@ids1 and @ids2\n";
				foreach $id1 (@ids1){
					foreach $id2 (@ids2){
						if ($id1 eq $id2){
							$new_id = "$id1|$info1[1]-$info1[2]-$info1[3]|$info2[1]-$info2[2]-$info2[3]";
							# print "overlapping $new_id ---- $full_id1 and $full_id2\n";
							$id_check = 1;
						}
					}
				}
				if ($id_check eq 1){
					print MERGED_OUT ">$new_id\n$overlapped\n";
					print MERGED_INFO_OUT "$new_id\t$full_id1\t$full_id2\n";
					$overlapping_primer1{$full_id1} = 1;
					$overlapping_primer2{$full_id2} = 1;
				}
				else{
					$count++;
					$overlapped = $seq1.(substr $seq2, -32);
					$new_id = "$sample-new$count|$info1[1]-$info1[2]-$info1[3]|$info2[1]-$info2[2]-$info2[3]";
					print MERGED_OUT ">$new_id\n$overlapped\n";
					print MERGED_INFO_OUT "$new_id\t$full_id1\t$full_id2\n";
					# print "New overlapping: $new_id\t$full_id1 and $full_id2\n";
					$overlapping_primer1{$full_id1} = 1;
					$overlapping_primer2{$full_id2} = 1;
				}
			}
		}
	}

	foreach my $seq1 (sort {$primer1_order{$a} <=> $primer1_order{$b}} keys %primer1){
		$new_id = "";
		my $full_id1 = $primer1{$seq1};
		if (!exists $overlapping_primer1{$full_id1}){
			my @info1 = split(/\|/, $full_id1);
			my @ids1 = split(",", $info1[0]);
			foreach my $id (@ids1){
				$new_id = $new_id.$id."/";
			}
			chop $new_id;
			print MERGED_OUT ">$new_id|$info1[1]-$info1[2]-$info1[3]|0-0-0\n$seq1\n";
			print MERGED_INFO_OUT "$new_id|$info1[1]-$info1[2]-$info1[3]|0-0-0\t$full_id1\n";
			# print "Not overlapping primer1: $full_id1\n";
		}
	}
	foreach my $seq2 (sort {$primer2_order{$a} <=> $primer2_order{$b}} keys %primer2){
		my $full_id2 = $primer2{$seq2};
		$new_id = "";
		if (!exists $overlapping_primer2{$full_id2}){
			my @info2 = split(/\|/, $full_id2);
			my @ids2 = split(",", $info2[0]);
			foreach my $id (@ids2){
				$new_id = $new_id.$id."/";
			}
			chop $new_id;
			print MERGED_OUT ">$new_id|0-0-0|$info2[1]-$info2[2]-$info2[3]\n$seq2\n";
			print MERGED_INFO_OUT "$new_id|0-0-0|$info2[1]-$info2[2]-$info2[3]\t$full_id2\n";
			# print "Not overlapping primer2: $full_id2\n";
		}
	}
}
elsif ($flag1 eq 1){
	foreach my $seq1 (sort {$primer1_order{$a} <=> $primer1_order{$b}} keys %primer1){
		$new_id = "";
		my $full_id1 = $primer1{$seq1};
		my @info1 = split(/\|/, $full_id1);
		my @ids1 = split(",", $info1[0]);
		foreach my $id (@ids1){
			$new_id = $new_id.$id."/";
		}
		chop $new_id;
		print MERGED_OUT ">$new_id|$info1[1]-$info1[2]-$info1[3]|0-0-0\n$seq1\n";
		print MERGED_INFO_OUT "$new_id|$info1[1]-$info1[2]-$info1[3]|0-0-0\t$full_id1\n";
		# print "Not overlapping primer1: $full_id1\n";
	}
}
elsif ($flag2 eq 1){
	foreach my $seq2 (sort {$primer2_order{$a} <=> $primer2_order{$b}} keys %primer2){
		$new_id = "";
		my $full_id2 = $primer2{$seq2};
		my @info2 = split(/\|/, $full_id2);
		my @ids2 = split(",", $info2[0]);
		foreach my $id (@ids2){
			$new_id = $new_id.$id."/";
		}
		chop $new_id;
		print MERGED_OUT ">$new_id|0-0-0|$info2[1]-$info2[2]-$info2[3]\n$seq2\n";
		print MERGED_INFO_OUT "$new_id|0-0-0|$info2[1]-$info2[2]-$info2[3]\t$full_id2\n";
		# print "Not overlapping primer2: $full_id2\n";
	}
}
else{
	print "$sample does not have $primer1_file and $primer2_file\n\n";
}
