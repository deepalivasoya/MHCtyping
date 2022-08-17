use Getopt::Long;
use Cwd;

my $samplesheet = "samplesheet.txt";
my $prefix = "run";
my $primers = "";
my $cutoff = 2;
my $fc = 3;

GetOptions(
    'samplesheet=s'    => \$samplesheet,
    'prefix=s' => \$prefix,
    'primers=s'     => \$primers,
    'cutoff=i'	=> \$cutoff,
    'fc=i'	=> \$fc
) or print "Invalid options\n";


system ("perl scripts/makeSummarySinglePrimers.pl --samplesheet=$samplesheet --prefix=$prefix --primers=$primers --cutoff=$cutoff --fc=$fc");
system ("cat summary/*.$primers.selected.txt > $primers.selected.txt");
system ("cat summary/*.$primers.discarded.txt > $primers.discarded.txt");
