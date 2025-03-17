use Getopt::Long;
use Cwd;

my $samplesheet = "samplesheet.txt";
my $work_dir    = "results";
my $prefix      = "run";
# my $primer1 = "For1Rev2";
# my $primer2 = "For3Rev1";

my $cutoff = 0.2;
my $fc = 3;

GetOptions(
    'samplesheet=s' => \$samplesheet,
    'work_dir=s'    => \$work_dir,
    'prefix=s'      => \$prefix
    # 'primer1=s'     => \$primer1,
    # 'primer2=s'	  => \$primer2,
    # 'cutoff=i'	  => \$cutoff,
    # 'fc=i'	      => \$fc
) or print "Invalid options\n";


print "=> Apply overlap_mhcI.pl on samples ...\n";
open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	$sample = $words[0];
	print "==> Running $sample...\n";
    print "===> CMD: perl scripts/overlap_mhcI.pl \\
        --sample=$sample \\
        --primer1=results/$sample/$sample.For1Rev2.selected.fasta \\
        --primer2=results/$sample/$sample.For3Rev1.selected.fasta\n\n";
	system ("perl scripts/overlap_mhcI.pl --sample=$sample --primer1=results/$sample/$sample.For1Rev2.selected.fasta --primer2=results/$sample/$sample.For1Rev2.selected.fasta");
    print "===> DONE $sample\n\n";
}
print "=> Apply overlap_mhcI.pl on samples ... Done \n\n\n";


print "=> Running makeSummaryOverlappedPrimers.pl";
print "==> CMD: perl scripts/makeSummaryOverlappedPrimers.pl \\
    --samplesheet=$samplesheet \\
    --prefix=$prefix \\
    --cutoff=$cutoff \\
    --database=fasta/MHCI.fasta \\
    --fc=$fc\n\n";
system ("perl scripts/makeSummaryOverlappedPrimers.pl --samplesheet=$samplesheet --prefix=$prefix --cutoff=$cutoff --database=fasta/MHCI.fasta --fc=$fc");
print "==> DONE\n\n";


system ("cat summary/*.mhcI.selected.tsv > $prefix.mhcI.selected.tsv");
system ("cat summary/*.mhcI.discarded.tsv > $prefix.mhcI.discarded.tsv");
system ("cat summary/*.mhcI.nc.tsv > $prefix.mhcI.nc.tsv");

