use Getopt::Long;
use Cwd;

my $samplesheet = "samples.txt";
my $prefix = "run";
my $pcr = "DRB3";
my $cutoff = 2;

GetOptions(
    'samplesheet=s'    => \$samplesheet,
    'prefix=s' => \$prefix,
    'pcr=s' => \$pcr,
    'cutoff=s'	=> \$cutoff,
) or print "Invalid options\n";

print "$samplesheet\n";
system ("perl scripts/makeSummarySinglePrimers.pl --samplesheet=$samplesheet --prefix=$prefix --pcr=$pcr --cutoff=$cutoff");
system ("cat summary/*.$pcr.selected.tsv > $pcr.selected.tsv");
system ("cat summary/*.$pcr.discarded.tsv > $pcr.discarded.tsv");
