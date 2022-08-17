use Getopt::Long;
use Cwd;

my $samplesheet = "samplesheet.txt";
my $work_dir = "results";
my $prefix = "run";
# my $primer1 = "For1Rev2";
# my $primer2 = "For3Rev1";

my $cutoff = 0.2;
my $fc = 3;

GetOptions(
    'samplesheet=s'    => \$samplesheet,
    'work_dir=s'     => \$work_dir,
    'prefix=s' => \$prefix
    # 'primer1=s'     => \$primer1,
    # 'primer2=s'	=> \$primer2,
    # 'cutoff=i'	=> \$cutoff,
    # 'fc=i'	=> \$fc
) or print "Invalid options\n";

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	$sample = $words[0];
	print "Running $sample...\n";
	system ("perl scripts/overlap_mhcI.pl --sample=$sample --work_dir=results/$sample --primer1=results/$sample/$sample.For1Rev2.clusters.selected.fasta --primer2=results/$sample/$sample.For1Rev2.clusters.selected.fasta");
}

system ("perl scripts/makeSummaryOverlappedPrimers.pl --samplesheet=$samplesheet --prefix=$prefix --cutoff=$cutoff --database=fasta/MHCI.fasta --fc=$fc");
system ("cat summary/*.mhcI.selected.txt > $prefix.mhcI.selected.txt");
system ("cat summary/*.mhcI.discarded.txt > $prefix.mhcI.discarded.txt");
system ("cat summary/*.mhcI.nc.txt > $prefix.mhcI.nc.txt");

