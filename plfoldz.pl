#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2018-10-12 15:52:01 mtw>
#
# Compute z-score for opening energies of intervals
#
# Idea: We are interested in the opening energy of the interval XXX in
# its genomic context
#       ----|augaagaugaXXXaugaagauga|---
# relative to the opening energy in a shuffled context
#       ----|uaaggaaaguXXXgauugaaaga|----

# The flanking regions upstream and downstream of the motif of
# interest are shuffled n times, while the motif itself is
# unchanged. (Mean) opening energies are computed for each suffled
# sequence and a z-score is computed.

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Carp;
use Cwd;
use Getopt::Long;
use Path::Class;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Bio::ViennaNGS::Fasta;
use Bio::ViennaNGS::FeatureIO;
use Bio::ViennaNGS::FeatureChain;
use Bio::ViennaNGS::Util qw(mkdircheck);
use Statistics::Lite;
use Data::UUID;
use Digest::SHA qw(sha1_base64);
use YAML::XS;
use RNA;
use Ushuffle;
use diagnostics;

my ($fastaO,$intervals,$od,$bn,$yn,$outfile,$outfh,$yamlfile,$logfile,$logfh);
my $infile_fa = undef;
my $infile_bed = undef;
my $outdir = undef;
my $wantlog = undef;
my $shuffle = 1000;
my $winlength = 100;
my $plfoldW = 100;
my $plfoldL = 100;
my $R=0.00198717; # R in kcal/mol/K
my $T0=273.15;
my $T=37; # default 37 celsius
my $tempK=$T0+$T;
my $kT=$R*$tempK;
my $dinuc=undef;
my $verbose=undef;
my $nozscore=undef;
my $klet = 0;
my @data = ();
my %params = (); # plfoldz runtim parameters

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless
  GetOptions("f|fa=s"     => \$infile_fa,
	     "b|bed=s"    => \$infile_bed,
	     "o|outdir=s" => \$outdir,
	     "log"        => sub{$wantlog=1},
	     "n=s"        => \$shuffle,
	     "d"          => sub{$dinuc=2},
	     "N|nozscore" => sub{$nozscore=1},
	     "v|verbose"  => sub{$verbose=1},
	     "w|window=s" => \$winlength,
	     "plfoldW=i"  => \$plfoldW,
	     "plfoldL=i"  => \$plfoldL,
	     "man"        => sub{pod2usage(-verbose => 2)},
	     "help|h"     => sub{pod2usage(1)}
	    );

unless (defined ($infile_fa)){
  warn "Please provide an infput fasta file to  -f|--fa option ...";
  pod2usage(-verbose => 0);
}

unless (defined ($infile_bed)){
  warn "Please provide an infput bed file to  -b|--bed option ...";
  pod2usage(-verbose => 0);
}

unless (-f $infile_fa){
  warn "Cannot find input fasta privided via -f|--fa option ...";
  pod2usage(-verbose => 0);
}
unless (-f $infile_bed){
  warn "Cannot find input bed privided via -b|--bed option ...";
  pod2usage(-verbose => 0);
}
my $cwd = getcwd();
unless ($infile_bed =~ /^\// || $infile_bed =~ /^\.\//){$infile_bed = file($cwd,$infile_bed);}
unless ($infile_fa  =~ /^\// || $infile_fa  =~ /^\.\//){$infile_fa  = file($cwd,$infile_fa);}

#print " ++ outdir: $outdir++\n";
unless(defined $outdir) { # not given via command line
  print " outdir is not defined vi cmdl\n";
  $od = $infile_bed->dir;
}
else { # given via commandline
  my $dir = dir($outdir);
  $od = $dir->absolute;
  unless (-d $od){mkdircheck($od);}
}

# obtain UUID
my $uido = Data::UUID->new();
my $uuid = $uido->create_str();
$bn = basename($infile_bed, ".bed");
$outfile = file($od,$bn.".".$uuid.".plfoldz.out");
open($outfh, ">", $outfile) or die "Cannot open outfile $outfile";
$yamlfile = file($od,$bn.".".$uuid.".plfoldz.yaml");

if (defined $wantlog){
  $logfile = file($od,$bn.".".$uuid.".plfoldz.log");
  #print "logfile $logfile ";
  open($logfh, ">", $logfile) or die "Cannot open logfile $logfile";
  #print "DONE\n";
}

if (defined $dinuc){ $klet = $dinuc };

%params = ( n => $shuffle, # store parameters for YAML dump
	    d => $klet,
	    T => $T,
	    winlen => $winlength,
	    plfoldW =>$plfoldW,
	    plfoldL => $plfoldL,
	    uuid => $uuid,
	    outdir => "$od");

$fastaO    = Bio::ViennaNGS::Fasta->new(fasta=>"$infile_fa");
$intervals = Bio::ViennaNGS::FeatureIO->new(
					    file => "$infile_bed",
					    filetype => 'Bed6',
					    instanceOf => 'Feature',
					   );

foreach my $f (@{$intervals->data}){
  my ($seq,$seq5,$seq3,$motiflength,$winlength5,$winlength3,$z);
  my ($shufseq5,$shufseq3,$ustart,$uend,$start,$dend);
  my %dat_orig = (); # YAML entry for the orig sequence
  my @shuf_oe = ();
  confess "ERROR: ID ".$f->chromosome." not found in Fasta file $infile_fa"
    unless ($fastaO->has_sequid($f->chromosome));
  my $fc = Bio::ViennaNGS::FeatureChain->new(type => "interval",
					     chain => [$f],
					    );
  my $bed6array = $fc->as_bed6_array();
  my $ps = %{$fastaO->primaryseqH}{$f->chromosome};
  my $motifstartpos = $f->start+1;
  if($f->end > $$ps{stop}){
    my $msg = "Motif ".$motifstartpos." -- ".$f->end." is outside of Fasta record ".$f->chromosome." (len= $$ps{stop})";
    warn $msg;
    next;
  }
  $seq = $fastaO->stranded_subsequence($f->chromosome,
				       $f->start+1,
				       $f->end,
				       $f->strand);
  $motiflength = length($seq);

  $ustart = $f->start+1-$winlength;
  $uend = $f->start;
  if($ustart < 1){
    $ustart = 1;
    if ($uend < $ustart){$seq5 = ""}
    else{
      $seq5 = $fastaO->stranded_subsequence($f->chromosome,
					    $ustart,
					    $uend,
					    $f->strand);
    }
  }
  else{
    $seq5 = $fastaO->stranded_subsequence($f->chromosome,
					  $ustart,
					  $uend,
					  $f->strand);
  }
  $winlength5 = length($seq5);
  if ($f->end == $$ps{stop}){$seq3 = ""}
  else{
    $seq3 = $fastaO->stranded_subsequence($f->chromosome,
					  $f->end+1,
					  $f->end+1+$winlength-1,
					  $f->strand);
  }
  $winlength3 = length($seq3);

  if (defined $wantlog){
    print $logfh join "\t", $f->chromosome,$f->name,$f->start+1,$f->end,"\n";
  }

  # compute opening energies for original sequence
  my $seqori = $seq5.$seq.$seq3;
  my $motifendpos = $motiflength+$winlength5;
  my $motifoe = compute_opening_energy($seqori,$motiflength,$plfoldW,$plfoldL,$motifendpos,$logfh);
  %dat_orig = ( seq => $seqori,
		motif => $seq,
		oe => $motifoe,
		endpos => $motifendpos,
		start => $f->start,
		end => $f->end,
		sum => sha1_base64($seqori),
		chr => $f->chromosome,
		strand => $f->strand,
		name => $f->name,
		score => $f->score);

  if (defined $wantlog){
    print $logfh "0 $seqori * (pos $motifendpos): $motifoe\n";
  }
  my @seq_up = split('',$seq5); # sequence upstream of motif
  my @seq_do = split('',$seq3); # sequence sownstream of motif

  if (defined $dinuc){ # initialize output sequence
    for(1..length($seq5)){$shufseq5 .= "a";}
    for(1..length($seq3)){$shufseq3 .= "a";}
  }
  for (my $i=0;$i<$shuffle;$i++){ # loop over n shuffled sequences
    unless (defined $dinuc){ # mononucleotide Fisher-Yates shuffling
      my @shuf5 = @seq_up;
      my @shuf3 = @seq_do;
      if($#seq_up > 0){@shuf5 = fisher_yates_shuffle(\@shuf5);}
      if($#seq_do > 0){@shuf3 = fisher_yates_shuffle(\@shuf3);}
      $shufseq5 = join ("", @shuf5);
      $shufseq3 = join ("", @shuf3);
      # print "M " if (defined $wantlog);
    }
    else { # dinucleotide / k-let shuffling
      $shufseq5 = $seq5;
      $shufseq3 = $seq3;
      if($#seq_up > 0){Ushuffle::shuffle($seq5,$shufseq5,$winlength5,$dinuc);}
      if($#seq_do > 0){Ushuffle::shuffle($seq3,$shufseq3,$winlength3,$dinuc);}
      # print "D "  if (defined $wantlog);
    }
    my $seqshuf = $shufseq5.$seq.$shufseq3;

    # computation of unpair probs / opening energies via ViennaRNA library call:
    my $oe = compute_opening_energy($seqshuf,$motiflength,$plfoldW,$plfoldL,$motifendpos,$logfh);
    if (defined $wantlog) {
      my $j = $i+1;
      if ($verbose){ print $logfh "$j $seqshuf S (pos $motifendpos shuffle $i)\t$oe\n"}
      else{print $logfh "$j\t$oe\n"}
    }
    push @shuf_oe, $oe;
  } # end for
  #print Dumper(\@shuf_oe);

  unless (defined $nozscore){
    if (Statistics::Lite::stddev(@shuf_oe) != 0.){
      $z = ($motifoe - Statistics::Lite::mean(@shuf_oe))/Statistics::Lite::stddev(@shuf_oe);
    }
    else {$z = 'NAN'}

    my $out = $$bed6array[0]."\t$motifoe\t$z";
    print $outfh "$out\n";
    if (defined $wantlog) { print $logfh "z-score: $z\n"}
  }

  my @shuf_clone = @shuf_oe;
  my %yaml = ( origin => \%dat_orig,
	       shufoe => \@shuf_clone,
	     );
  push @data, \%yaml;
  undef @shuf_oe; # memory cleanup
} # end foreach

close($outfh);
if (defined $wantlog){close($logfh)}
my $dmp_yaml = YAML::XS::DumpFile($yamlfile, { params => \%params, data => \@data } );

sub compute_opening_energy {
  my ($seq, $ulength, $window_size,$max_bp_span,$where,$lfh) = @_;
  my @oe = ();
  for (my $i=0;$i<length($seq);$i++){$oe[$i]=0.}
  # compute unpaired probabilities for original sequence
  my $up = RNA::pfl_fold_up($seq, $ulength, $window_size, $max_bp_span);
  # print Dumper($up);
  for (my $i = $ulength;$i <= length($seq);$i++){
    my $up = $$up[$i][$ulength];
    $oe[$i] = -$kT*log($up);
    #print $lfh "|$i|$ulength| up:$up\t oe:$oe[$i]\n";
  }
  return $oe[$where];
}

sub fisher_yates_shuffle {
  my $array = shift;
  my $i;
  for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
  }
  return @$array;
}

__END__

=head1 NAME

plfoldz.pl - Compute z-score for opening energies of intervals

=head1 SYNOPSIS

plfoldz.pl [-f|--fa I<FILE>] [-b|--bed I<FILE>] [-n I<INT>] [-d]
[-o|--outdir I<DIR>] [-w|--window I<INT>] [--log] [--plfoldW I<INT>]
[--plfoldL I<INT>] [-N|--nozscore] [-v|--verbose] options

=head1 DESCRIPTION

This script computes z-scores of opening energies for an RNA sequence
in genomic context related to the same sequence in a shuffled context
(albeit preserving the relative nucleotide frequencies in surrounding
windows, optionally applying dinucleotide shuffling).

The idea is as follows. We are interested in the (mean) opening energy
of interval I in its genomic context (let XXX represent I here)

 ----|augaagaugaXXXaugaagauga|---
 relative to the opening energy in a shuffled context
 ----|uaaggaaaguXXXgauugaaaga|----

The flanking regions upstream and downstream of the motif (interval)
of interest I are shuffled n times, while the motif itself is left
unchanged. Mean opening energies of interval I are computed in a
genomic context (native sequence) as well as for each of the n suffled
sequences by (a library call to) RNAplfold. A z-score is computed for
the mean opening energies of interval I (this can be turned off by the
B<-N|--nozscore> flag, which is essentially useful when the
calculations are done in a distributed cluster environment).

Input RNA sequences are provided as BED6 intervals, together with a
genomic Fasta file (evidently, Fasta and BED IDs should match). The
opening energy of the original interval is reported together with the
z-score of n shufflings as an extended BED6 output file. The output
file is the basename of the input BED6 file, preprended with
".plfoldz.out". An optional log file lists each shuffled sequence
together with the opening energy of the interval of interest.

=head1 OPTIONS

=over

=item B<-f|--fa>

Input Fasta file which should contain the genomic region referenced my
the input BED6 file.

=item B<-b|--bed>

Input BED6 file. Must contain all RNA intervals for which opening
energies and z-scores should be computed.

=item B<-o|--outdir>

Output directory. Output and log files are written there. The name of
the output file is composed on the basename of the input BED file,
with a 'plfoldz.out; suffix.

=item B<--log>

Optional flag, which triggers creation of a log file. The log file has
the same name as the output file, but with a .log extenion.

=item B<-n>

Number of shuffling events (default: 1000).

=item B<-d>

Enable dinucleotide shuffling using the uShuffle library. Default is
Fisher-Yates mononucleotide shuffling.

=item B<-N|--nozscore>

Turn off the default calculation of opening energy z-scores, just
report the mean opening energies for each interval. This option is
typically used in cluster environment, when shuffling is distributed
over several compute nodes. Computation of z-scores is then performed
with F<plfoldz_gather_zscores.pl>.

=item B<-w|--window>

Window size of upstream/downstream nucleotides that will be used for
shuffling.

=item B<--plfoldW>

Size of sliding window for computation of local base pairing
probabilitied . This is essialliy the same as the -W parameter to
RNAplfold.

=item B<--plfoldL>

Maximum base pair span. Corresponds to -L paramater of RNAplfold.

=item B<-v|--verbose>

Be verbose in the log file.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
