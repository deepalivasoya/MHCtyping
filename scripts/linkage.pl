$in = shift;
#0 			1 				2 				3 		4 		5  			6
# Sample	DRB1_allele_1	DRB1_allele_2	Flock	Breed	Haplotype1	Haplotype2
# BF141	301	308	C	BF	HP:05:01	HP:01:01

open(OUT, ">linkage1.txt");
open(OUT1, ">linkage2.txt");
open(OUT2, ">linkage3.txt");
open(OUT3, ">linkage4.txt");
open(OUT4, ">linkage5.txt");
open(OUT5, ">linkage6.txt");
open(OUT6, ">linkage7.txt");
open(OUT7, ">linkage8.txt");

open(IN, "$in");
while(<IN>){
	chomp $_;
	@words = split("\t", $_);

	print OUT "$words[5]\t$words[1]\n";
	print OUT "$words[6]\t$words[1]\n";
	print OUT "$words[5]\t$words[2]\n";
	print OUT "$words[6]\t$words[2]\n";
	if (!exists $class1{$words[5]}){
		$class1{$words[5]} = 0;
	}
	if (!exists $class1{$words[6]}){
		$class1{$words[6]} = 0;
	}
	if (!exists $class2{$words[1]}){
		$class2{$words[1]} = 0;
	}
	if (!exists $class2{$words[2]}){
		$class2{$words[2]} = 0;
	}
}
close(OUT);

open(IN, "linkage1.txt");

while(<IN>){
	chomp $_;
	@words = split("\t", $_);
	$class2_hp = $words[1];
	$class1_hp = $words[0];
	$$class2_hp{$class1_hp}++;
	$class1{$words[0]}++;
	$class2{$words[1]}++;
}

# $total_class2 = 0;
%total_hp_class2 = ();
foreach $class2_hp (sort {$a<=>$b} keys %class2){
	print OUT1 "\t$class2_hp";
	print OUT2 "\t$class2_hp";
}
print OUT1 "\tTotal\n";
print OUT2 "\tTotal\tMax\n";

foreach $class1_hp (sort {$a<=>$b} keys %class1){
	$max = 0;
	print OUT1 "$class1_hp";
	print OUT2 "$class1_hp";
	foreach $class2_hp (sort {$a<=>$b} keys %class2){
		if (!exists $$class2_hp{$class1_hp}){
			$$class2_hp{$class1_hp} = 0;
		}
		$norm = sprintf("%.2f", (($$class2_hp{$class1_hp} / $class1{$class1_hp}) * 2));
		if ($norm > $max){
			$max = $norm;
			$max_class2 = $class2_hp;
		}
		$$class1_hp{$class2_hp} = $norm;
		print OUT1 "\t$$class2_hp{$class1_hp}";
		print OUT2 "\t$norm";
		$total_hp_class2{$class2_hp} = $total_hp_class2{$class2_hp} + $norm;
	}
	print OUT1 "\t".$class1{$class1_hp}."\n";
	print OUT2 "\t".$class1{$class1_hp}."\t$max\n";
}

print OUT1 "Total";
print OUT2 "Total";
foreach $class2_hp (sort {$a<=>$b} keys %class2){
	print OUT1 "\t$class2{$class2_hp}";
	print OUT2 "\t$total_hp_class2{$class2_hp}";
}
print OUT1 "\t\n";
print OUT2 "\t\n";

foreach $class1_hp (sort {$a<=>$b} keys %class1){
	print OUT3 "$class1_hp";
	foreach $class2_hp ( sort {$$class1_hp{$b} <=> $$class1_hp{$a}} keys %$class1_hp){
		if ($$class1_hp{$class2_hp} > 0){
			print OUT3 "\t$class2_hp|$$class1_hp{$class2_hp}";
		}
	}
	print OUT3 "\n";
}
close(OUT3);
# open(IN, "class2.linkage4.txt");
# while(<IN>){
# 	chomp $_;
# 	@words = split("\t", $_);
# 	# print scalar(@words)."\n";
# 	$hp = $words[0];
# 	for ($i = 1; $i < scalar(@words); $i++){
# 		@info = split(/\|/, $words[$i]);
# 		$$hp{$i} = $info[0];
# 	}
# }
# foreach $class1_hp (sort {$a<=>$b} keys %class1){
# 	foreach $linkage (sort {$a <=> $b} keys %$class1_hp){
# 		print "$class1_hp\t$linkage|$$class1_hp{$linkage}";
# 	}
# 	print "\n";
# }
%final_link = ();
open(IN, "$in");
while(<IN>){
	chomp $_;
	@words = split("\t", $_);
	$class2_hp1 = $words[1];
	$class2_hp2 = $words[2];
	$hp1 = $words[5];
	$hp2 = $words[6];
	# print "$hp1\t";
	$final_hp1 = "";
	$final_hp2 = "";
			
	$link_hp1_drb1 = $$hp1{$class2_hp1};
	$link_hp1_drb2 = $$hp1{$class2_hp2};
	$link_hp2_drb1 = $$hp2{$class2_hp1};
	$link_hp2_drb2 = $$hp2{$class2_hp2};

	print OUT4 "$words[0]\t$words[3]\t$words[4]\t$hp1\t$link_hp1_drb1\t$link_hp1_drb2\t$hp2\t$link_hp2_drb1\t$link_hp2_drb2\t$class2_hp1\t$class2_hp2\t";

	if ($link_hp1_drb1 > $link_hp1_drb2){
		$final_hp1 = $class2_hp1;	
	}
	else{
		$final_hp1 = $class2_hp2;
	}
	
	if ($link_hp2_drb1 > $link_hp2_drb2){
		$final_hp2 = $class2_hp1;	
	}
	else{
		$final_hp2 = $class2_hp2;
	}

	if ($final_hp1 eq $final_hp2){
		print OUT4 "\tQuery";
		if ($final_hp1 eq $class2_hp1){
			if ($link_hp1_drb1 >= $link_hp2_drb1){
				$final_hp2 = $class2_hp2;
			}
			else{
				$final_hp1 = $class2_hp2;
			}
		}
		if ($final_hp1 eq $class2_hp2){
			if ($link_hp1_drb2 >= $link_hp2_drb2){
				$final_hp2 = $class2_hp1;
			}
			else{
				$final_hp1 = $class2_hp1;
			}
		}
	}
	else{
		print OUT4 "\tOk";
	}
	$final_link{$hp1} = $final_link{$hp1}."$final_hp1,";
	$final_link{$hp2} = $final_link{$hp2}."$final_hp2,";
	# print OUT4 "\n";
	print OUT4 "\t$hp1\t$final_hp1\t$$hp1{$final_hp1}\t$hp2\t$final_hp2\t$$hp2{$final_hp2}\n";

	$both = "$hp1 and $final_hp1";
	$farm{$both} = $farm{$both}."$words[3],";
	# print "$both\t$words[3]\n";

	$both = "$hp2 and $final_hp2";
	$farm{$both} = $farm{$both}."$words[3],";
	# print "$both\t$words[3]\n";


}

print OUT6 "HP\tTotal\tLinks\tMHCII\tAnimals\tFlock\n";
print OUT7 "HP\tTotal\tLinks\tMHCII\tAnimals\n";

foreach $hp (keys %final_link){
	my %count = ();
	chop $final_link{$hp};
	my @links = split(",", $final_link{$hp});
	foreach $link (@links){
		$count{$link}++;
	}
	
	my $total_ob = $class1{$hp}/2;
	my $count_link = keys %count;
	print OUT5 "$hp\t$total_ob\t$count_link";
	open (HP_OUT, ">Alluvium/$hp.alluvium.txt");
	print HP_OUT "HP\tTotal\tLinks\tMHCII\tAnimals\tflock\n";
	foreach $link (sort{$count{$b}<=>$count{$a}} keys %count){
		$per = sprintf("%.2f", (($count{$link}/$class1{$hp}) * 200));
		print OUT5 "\t$link\t$$hp{$link}\t$count{$link}\t$per";
		$both = "$hp and $link";
		chop $farm{$both};
		print "$hp\t$farm{$both}\n";
		my @flocks = split(",", $farm{$both});
		my %count_flock = ();
		foreach $br (@flocks){
			$count_flock{$br}++;
		}
		my $to_print = "";
		foreach $br (sort {$count_flock{$b} <=> $count_flock{$a}} keys %count_flock){
			# print "$br\n";
			$to_print = $to_print."$br($count_flock{$br}),";
			print HP_OUT "$hp\t$total_ob\t$count_link\t$link\t$count_flock{$br}\t$br\n";
			print OUT6 "$hp\t$total_ob\t$count_link\t$link\t$count_flock{$br}\t$br\n";
		}
		print OUT7 "$hp\t$total_ob\t$count_link\t$link\t$count{$link}\n";
		
		chop $to_print;
		print OUT5 "\t$to_print";
	}
	print OUT5 "\n";
	close(HP_OUT);
}





