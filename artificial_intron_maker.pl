#artificial_intron_maker.pl
use strict;
use warnings FATAL => "all";
use Getopt::Std;
use Data::Dumper;
use vars qw($opt_g $opt_l $opt_n $opt_p $opt_i $opt_s $opt_k $opt_b $opt_h $opt_m $opt_t);
getopts('g:l:n:p:c:s:k:hm:i:t:');

## Command line and global variables ##

my $usage = "
usage: $0 <kmer file> <all introns>
  -g <float> starting GC frequency of introns [0.35]
  -l <int>   length of intron [250]
  -n <int>   number of generations [200]
  -p <int>   population size [100]
  -i <float> immortal fraction [0.1]
  -s <float> scaling factor for KL distance[2.0]
  -t <float> scaling factor for motif [0.5]
  -b <int>   average intron k-mer length in background (use 1 or 2) [1]
  -m <float> mutation rate [0.01]
  -h         show help
";

if ($opt_h) {die $usage}

my $GC_FREQ     = $opt_g ? $opt_g : 0.35;
my $LENGTH      = $opt_l ? $opt_l : 250;
my $GENERATIONS = $opt_n ? $opt_n : 200;
my $POPULATION  = $opt_p ? $opt_p : 200;
my $IMMORTAL    = $opt_i ? $opt_i : 0.1;
my $SCALE_KL    = $opt_s ? $opt_s : 2;
my $SCALE_MOTIF = $opt_t ? $opt_t : 0.5;
my $BK          = $opt_b ? $opt_b : 1;
my $MUT         = $opt_m ? $opt_m : 0.01;
my $K; # inferred below

die $usage unless @ARGV == 2;
my ($kmer_file, $normal_file) = @ARGV;

## Global Variables : IME and Background ##

my %IME_score;
open (my $fh, $kmer_file) or die;
while (my $line = <$fh>){
	chomp $line;
	my ($kmer, $f1, $f2, $lod) = split(/\s+/, $line);
	$IME_score{$kmer} = $lod;
	$K = length($kmer) unless defined $K;
}
close ($fh);

my %BG_freq;
my %BG_count;
my $introme;
open (my $fh1, $normal_file) or die;
while (<$fh1>) {
	next if /^>/;
	chomp;
	$introme .= $_;
}
my $total = 0;
for (my $i = 0; $i < length($introme) - $BK + 1; $i++) {
	$BG_count{substr($introme, $i, $BK)}++;
	$total++;
}
foreach my $kmer (keys %BG_count) {
	$BG_freq{$kmer} = $BG_count{$kmer} / $total;
}

## Set up initial population of random introns ##

for (my $i = 1; $i <= 20; $i ++){
	my @pop;
	for (my $i = 0; $i < $POPULATION; $i++) {
		push @pop, random_individual($LENGTH, $GC_FREQ);
	}

## Breed population ##
	{
	open (my $h, '>>',  "messing/test.$GENERATIONS.$SCALE_MOTIF.txt") or die "can't make file";
		for (my $g = 0; $g < $GENERATIONS; $g++) {

	# score all individuals in population
		foreach my $ind (@pop) {
			$ind->{ime} = score_ime($ind->{seq}, \%IME_score);
			$ind->{comp} = score_comp($ind->{seq}, \%BG_freq);
			$ind->{motif} = count_motif($ind->{seq});
			$ind->{fitness} = ($ind->{ime}) * ((1 - $ind->{comp}) ** $SCALE_KL) * ($ind->{motif} ** $SCALE_MOTIF);
		}
		
	# sort individuals by fitness
	@pop = sort {$b->{fitness} <=> $a->{fitness}} @pop;
	printf {$h} "%d %.3f %.3f %d %.3f %d\n", $g, $pop[0]{ime}, $pop[0]{comp}, $pop[0]{motif}, $pop[0]{fitness}, $i;
	
	# replace unfit individuals with childen of best parents
	for (my $i = $POPULATION * $IMMORTAL; $i < $POPULATION; $i++) {
		my $p1 = rand($POPULATION * $IMMORTAL);
		my $p2 = rand($POPULATION * $IMMORTAL);
		$pop[$i] = mate($pop[$p1]{seq}, $pop[$p2]{seq});
	}
	
}
	close ($h);
}
}
############################################################################

sub random_individual {
	my ($len, $gc) = @_;
	my $at = 1 - $gc;	
	my $string = 'A' x ($at * 100) . 'C' x ($gc * 100) . 
		'G' x ($gc * 100) . 'T' x ($at * 100);
	my $seq = "";
	for (my $i = 0; $i < $len; $i++) {
		$seq .= substr($string, rand(length($string)), 1);
	}
	return {
		seq => $seq,
		ime => 0,
		comp => 0,
		fitness => 0,
		motif => 0,
	};
}

sub score_ime{
	my ($seq, $ime_score) = @_;
	my @parts;
	my $sum = 0;
	for (my $i = 0; $i < $LENGTH - $K + 1; $i ++ ){
		my $kmer = substr($seq, $i, $K);
		push @parts, $kmer;
	}
	foreach my $kmer (@parts){
		my $logodds = $ime_score->{$kmer};
		$sum += $logodds;
	}
	return $sum;
}

sub score_comp{
	my ($seq, $BG_freq) = @_;
	my $total = 0;
	my $n = 0;
	my (%kld_count, %kld_freq);
	for (my $i = 0; $i < length($seq) - $BK + 1; $i++) {
		$kld_count{substr($seq, $i, $BK)}++;
		$n ++;
	}
	foreach my $kmer (keys %kld_count) {
		$kld_freq{$kmer} = $kld_count{$kmer} / $n;
	}
	foreach my $kmer (keys %kld_freq){
		my $x = $kld_freq{$kmer};
		my $y = $BG_freq->{$kmer};
		my $logodds = $y * log($y/$x);
		my $log2 = $logodds/ log(2);
		$total += $log2;
	}
	return $total;
}

sub mate{
	my ($p1, $p2) = @_;
	my $child = "";
	my @parent = ("$p1", "$p2");
	while (length($child) < $LENGTH){
		my $rp = int(rand(1));
		my $rs = int(rand($LENGTH-$K));
		my $chunk = substr($parent[$rp], $rs, $K);
		$child .= $chunk;
	}
	$child = substr($child, 0, $LENGTH);
	my @mutate = split("", $child);
	my @letters = qw(A T C G);
	for (my $i = 0; $i <= $MUT * ($LENGTH - 1); $i++){
		my $letter = $letters[int(rand(3))];
		splice @mutate, rand($LENGTH), 1, $letter; 
	}
	my $progeny = join ("", @mutate);
	return {
		seq => $progeny,
		ime => 0,
		comp => 0,
		fitness => 0,
		motif => 0,
	};
}

sub count_motif{
	my ($seq) = @_;
	my @match = $seq =~ /TT[ACTG]GAT[TC]TG/g;
	my $count = @match;
	die "count uninitialized in subroutine count_motif" if not defined $count;
	return $count;
}

