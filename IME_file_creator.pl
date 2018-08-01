#KL_distances.pl
use strict; use warnings FATAL => "all";

my ($file, $ilimit, $klimit) = @ARGV;
my @intron;

open (my $fh, $file) or die "Couldn't open";

while (my $line = <$fh>){
	chomp($line);
		if ($line =~ /_i(\d+)_/){
			my $d;
			my ($pl) = $line =~ /\_(\d+)\_/;
			if ($pl <= 400){
				$d = 0;
			}elsif ($pl <= 850){
				$d = 1;
			}
			else{
				$d =2;
			}		
			my $seq = <$fh>;
			chomp($seq);
			push (@{$intron[$d]}, $seq);
		}	
}
close ($fh);

my $directory;
main();

for (my $k = 3; $k <= $klimit; $k ++){
	my $fr1;
	my $fr2;
	for (my $i = 0; $i <= $ilimit; $i ++){
		$fr1 = kmer_hash($intron[$i], $k);
		for (my $j = $i+1; $j <= $ilimit; $j++){
			$fr2 = kmer_hash($intron[$j], $k);
			my ($new1, $new2) = pool($fr1, $fr2);
			my $n1 = kmer_freq($new1);
			my $n2 = kmer_freq($new2);
			my $kld = kld($n1, $n2);
			{
				open (my $fh, '>',  "$directory/report.$k.$i.$j.txt") or die "can't make file";
			
				my %logodds = logodds($n1, $n2);
				my $logodds = \%logodds;
				my @key = keys(%logodds);
				foreach my $kmer (@key){
					print {$fh} "$kmer ${$n1}{$kmer} ${$n2}{$kmer} ${$logodds}{$kmer}\n";
				}
				close ($fh)
			}
		}
	}
}

sub kmer_hash{
	my @seq = @{$_[0]};
	my $l = $_[1];
	my %count;
	foreach my $seq (@seq){
		for (my $i = 0; $i < (length($seq)-$l); $i++){
			my $kmer = substr($seq, $i, $l);
			$count{$kmer} ++ ;
		}
	}
	my $count = \%count;
	return ($count);
}
sub pool{
	my ($p, $q) = @_;
	foreach my $kmer(keys %{$p}){
		if (exists $q->{$kmer}) {
			$q->{$kmer}++;
		}else{
			$q->{$kmer} = 1;
		}
	}
	foreach my $kmer(keys %{$q}){
		if (exists $p->{$kmer}){
			$p->{$kmer}++;
		}else{
			$p->{$kmer} = 1;
		}
	}
	return ($p, $q);
}
sub kmer_freq{
	my $total = 0;
	my ($count) = @_;
	my %frequency;
	foreach my $kmer (keys %{$count}){
		$total += $count->{$kmer};
	}	
	foreach my $kmer (keys %{$count}){
		my $kmer_fr = ($count->{$kmer})/$total;
		$frequency{$kmer} = $kmer_fr;
	}
	my $freq = \%frequency;
	return ($freq);
}

sub kld{
	my $sum = 0;
	my ($fr1, $fr2) = @_;
	foreach my $kmer (keys %{$fr1}){
		my $x = $fr1->{$kmer};
		my $y = $fr2->{$kmer};
		my $part = ($x * log($x/$y));
		my $log = $part/ log(2);
		$sum += $log;	
	}
	return ($sum);
}

sub logodds{
	my ($fr1, $fr2) = @_;
	my $logratio;
	my %ratio;
	foreach my $kmer (keys %{$fr1}){
		my $x = $fr1->{$kmer};
		my $y = $fr2->{$kmer};
		$logratio = log($x/$y);
		$ratio{$kmer} = $logratio;
	}
	return (%ratio);
}

sub main {
    $directory = "report";
   
    unless(mkdir $directory) {
        die "Unable to create $directory\n";
    }
}





