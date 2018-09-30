#!/usr/bin/perl -w
use strict;
use Data::Dump::Color qw<dd>;
use DBI;
use Getopt::Long::Descriptive;
use List::Util qw<min max>;
use Sort::Naturally qw<nsort>;
use version 0.77;
our $VERSION = version->declare(1.0.0); ## GIT_VERSION ## our $VERSION = version->declare(<%>);
our $BRANCH = 'master'; ## GIT_BRANCH ## our $BRANCH = '<%>';

=head1 NAME - cluster.pl

=head1 SYNOPSIS

perl cluster.pl -i my_regions.tsv > clustered.tsv

=cut

# get the command line params
my ( $opt, $usage ) = describe_options(
	'cluster %o <some-arg>',
	["
	input file and characteristics:"
	],
	['input|i=s',   "input file with coordinates",            {default => 'demo_sample_male.txt'}],
	['sep|p=s',     "character that separates columns",       {default => "\t"}],
	['exclude|e=s', 'exclude regions that match (_CONTROL_)', {default => '_CONTROL_'}],
	["
	columns in input file (0-based counts):"
	],
	['chrcol=i', 'chromosome (0)', {default => 0}],
	['poscol=i', 'position   (1)', {default => 1}],
	['dircol=i', 'direction  (2)', {default => 2}],
	['bascol=i', 'base       (3)', {default => 3}],
	['namcol=i', 'name       (4)', {default => 4}],
	['cntcol=i', 'counts     (5)', {default => 5}],
	["
	clustering parameters:"
	],
	['min|m=i',      "m) minimal absolute regions per cluster (20)",          {default => 20}],
	['relative|r=i', "r) minimal relative regions to split (0.25)",           {default => 0.25}],
	['logic|l=s',    "minimal regions <-m,-r> logic (OR|AND*)",               {default => 'AND'}],
	['above|a=i',    "a) split regions with gap larger than this (100,000)",  {default => 100_000}],
	['span|s=i',     "s) minimal relative portian a gap spans (0.50)",        {default => 0.50}],
	['gaplogic|g=s', "gap size logic <-a,-s> (OR*|AND)",                      {default => 'OR'}],
	['force|f=i',    "force splitting cluster with this many regions (200)",  {default => 200}],
	['below|b=i',    "never split at gaps smaller than this (1,000)",         {default => 1_000}],
	['distal|q=i',   "separate at gaps > this (regardless of <-m>; 250,000)", {default => 250_000}],
	["
	cluster naming:"
	],
	['names|n=s',    "use this file as source for names",                                       {default => 'refGene.txt'}],
	['type|t=s',     "define the filetype (refgene*|ideo|custom)",                              {default => 'refgene'}],
	['match|z=s',    "only use line that match this (NM_ for refgene)",                         {default => ''}],
	['deadband|d=i', "ignore that many bases of overlap (1000)",                                {default => 1000}],
	['suffix|x=s',   "use this to separate numbered occurence ('.' for refgene; '/' for ideo)", {default => ''}],
	["
	columns (zero-based) in 'custom' name file (presets for refgene/ideo):"
	],
	['nmchr=i',   'chromosome   (2/0)', {default => 2}],
	['nmbegin=i', 'region begin (6/1)', {default => 6}],
	['nmend=i',   'region end   (7/2)', {default => 7}],
	['annot=i',   'annotation  (12/3)', {default => 12}],
	["
	internals (no need to change):"
	],
	["dbh=s",    "handle for SQLite (:memory:)", {default => ":memory:"}],
	["export=s", "column names for exporting",   {default => [qw<chr pos dir base name count orgname>]}],
	["dump",     "dump all options and exit"],
	["
	general:"
	],
	['verbose|v', "print extra stuff"],
	['version|V', 'show the version and exit'],
	['help|h',    "print usage message and exit"],
);
print( $usage->text ), exit if $opt->help;
if ( $opt->version ) {
	my $name = $0;
	$name =~ s!^\.\/!!;
	printf "This is %s %s %s\nUse flags -h or --help for help.\n", $name, $VERSION->stringify, $BRANCH ? qq~($BRANCH)~ : '';
	exit;
}
if ( $opt->{dump} ) {
	dd $opt;
	exit;
}
if ( not -e $opt->{input} ) {
	warn "Input file not specified or not readable!";
	print "\n", $usage->text;
	exit;
}

=head1 DESCIPTION

This script will perform a clustering on a list of coordinates provided via an input file.
The approach will first create one cluster for every molcule (chromosome) found in the file 
and subsequently try to split larger clusters into smaller ones.

Various criteria can be set to control this process, the most important of which are (defaults 
are in round brackets):

=over 3

=item C<--min | -m> - m) minimal absolute regions per cluster (20)

=item C<--relative | -r> - r) minimal relative regions to split (0.25)

=item C<--logic | -l> - minimal regions <m,r> logic (OR | -AND*)

=item C<--above | -a> - a) split regions with gap larger than this (100,000)

=item C<--span | -s> - s) minimal relative portian a gap spans (0.50)

=item C<--gaplogic | -g> - gap size logic <a,s> (OR*|AND)

=item C<--force | -f> - force splitting cluster with this many regions (200)

=item C<--below | -b> - never split at gaps smaller than this (1,000)

=back

=head1 ARGUMENTS

Some advanced settings are available.

=head2 INPUT FILE PARSING

To read the relevant information, the script can be told where to get these from. 
Parameter to extract columns in input file (0-based counts):

=over 3

=item C<--chrcol> - chromosome (0)

=item C<--poscol> - position (1)

=item C<--dircol> - direction (2)

=item C<--bascol> - base (3)

=item C<--namcol> - name (4)

=item C<--cntcol> - counts (5)

=back

Additionally, a string can be specified to exclude lines containing that string, 
e.g. control oligos that should not be reported as missing regions.

=over 3

=item C<--exclude | -e> - disregard lines matching this (_CONTROL_)

=head2 CLUSTER NAMING

To annotate the newly formed clusters, a file with names linked to coordinates can be
used.

=over 3

=item C<--names | -n> - use this file as source for names (e.g. 'refGene.txt')

=item C<--type | -t> - define the filetype <refgene|ideo|custom> ('refgene')

=item C<--match | -z> - only use line that match this (NM_ for refgene)

=item C<--deadband | -d> - ignore that many bases of overlap (1,000)

=item C<--suffix | -x> - use this to separate numbered occurence ('.' for refgene; '/' for ideo)

=back

Where to find the relevant information is adjustable. Columns (zero-based) in 'custom' 
name file (presets for refgene/ideo) can be set. Defaults are given for the two presets 
(C<--type>, see above) refgene and ideo:
 
=over 3
 
=item C<nmchr> - chromosome (2/0)

=item C<nmbegin> - region begin (6/1)

=item C<nmend> - region end (7/2)

=item C<annot> - annotation (12/3)

=back

=cut

# create the sqlite database
my $dbh = DBI->connect( join( ":", 'dbi', 'SQLite', 'dbname=' . $opt->{dbh} ), "", "" );
if ( -e $opt->{names} ) {
	read_namefile($dbh);    # read annotations

	#		dd $dbh->selectall_arrayref( qq~SELECT * FROM `names` LIMIT 10;~, { Slice => {} } );
	#	exit;
}
read_regions($dbh);    # read regions

# idea: merge all region per chromosome and then split into parts if reasonable
# get all chromosome names first
my $Chr = $dbh->selectcol_arrayref(qq~SELECT DISTINCT(chr) FROM primers;~);
$Chr = [nsort @{$Chr}];
my $prefixed = 0;
if ( @{$Chr} and ( $Chr->[0] =~ /^chr/ or $Chr->[-1] =~ /^chr/ ) ) {

	# to create a common order, sort without 'chr' prefix and add that after sorting
	$Chr = [nsort map { substr( $_, 3 ) } @{$Chr}];
	$Chr = [map { 'chr' . $_ } @{$Chr}];
}
my $cluster = 1;

# loop all chromosomes
CHROMOSOME: foreach my $chr ( @{$Chr} ) {

	# set the current cluster count to all regions from chromosome
	$dbh->do(qq~UPDATE primers SET cluster = $cluster WHERE chr = "$chr";~);
	$cluster++;
}
my $max_cluster = $cluster;
my $next        = $max_cluster + 1;
$cluster = 0;
my $gap_require = $opt->{gaplogic} eq 'OR' ? 1 : 2;

# loop all current clusters
CLUSTER: while ( $cluster <= $next ) {
	my $split = 1;
	$cluster++;

	# while the last actions was a split of a cluster, try again
	SPLIT: while ($split) {
		$split = 0;

		# search for the largest gap that has at least $min regions on both sides
		my $All
		= $dbh->selectall_arrayref( qq~SELECT id, pos FROM primers WHERE cluster = $cluster ORDER BY pos;~, {Slice => {}} );
		next CLUSTER unless scalar @{$All} > 2;    # >= $opt->{min};
		my $length = $All->[-1]{'pos'} - $All->[0]{'pos'};
		my $total  = scalar @{$All};

		#next CLUSTER if $length < $opt->{below};
		my %Targets;
		my $first = max( $opt->{min} - 1, int( 0.5 + ( $total * $opt->{relative} ) ) );
		my $last = min( $total - $opt->{min} - 1, $total - int( 0.5 + ( $total * $opt->{relative} ) ) );

		#next CLUSTER if $last < $first;
		# scan through the list and collect distance in the same cluster
		GAP: for my $i ( 1 .. $total - 1 ) {    #$first .. $last ) {
			my $distance = $All->[$i]{pos} - $All->[$i - 1]{pos};
			next GAP unless $distance >= $opt->{below};

			# shortcut to include very large gaps in the list, even if not
			# enough data points are present on both sides of the gap
			if ( $opt->{distal} > 0 and $distance >= $opt->{distal} ) {

				# print STDERR "FORCE!\n";
				$Targets{$i} = $distance;
				next GAP;
			}
			next GAP if $i < $first or $i > $last;
			my $rel_gap = $distance / $length;
			my $valid   = 0;
			$valid++ if $distance >= $opt->{above};
			$valid++ if $rel_gap >= $opt->{span};
			$valid += 2 if $total >= $opt->{force};
			$Targets{$i} = $distance if $valid >= $gap_require;
		} ## end GAP: for my $i ( 1 .. $total - 1 )

		#		select STDERR;
		#		dd \%Targets;
		#		select STDOUT;
		# next cluster if no big gap present
		next CLUSTER if ( not keys %Targets );

		# now get the largest distance
		my @Distance = sort { $Targets{$a} <=> $Targets{$b} } keys %Targets;

		#	dd $Distance[-1];
		$split = 1;

		# now assign the next cluster to all regions RIGHT of the element with biggest distance
		my $cutoff = $All->[$Distance[-1]]{pos};
		$dbh->do(
			qq~UPDATE primers SET cluster = $next 
			WHERE cluster = $cluster AND pos >= $cutoff;~
		);
		$next++;
	} ## end SPLIT: while ($split)
} ## end CLUSTER: while ( $cluster <= $next )

# EXPORTING
my $counter = 0;
my %NAMES;
my %STAT;

# first pass 'dry' to check for ambiguous names
RUN: for my $run ( 'dry', 'wet' ) {

	# loop all chromosomes
	CHR: foreach my $chr ( @{$Chr} ) {
		my $chrcounter = 0;
		my $clusters   = $dbh->selectall_arrayref(
			qq~ 
		SELECT cluster, chr, min(pos) AS begin, max(pos) AS end FROM primers 
		WHERE chr = "$chr" 
		GROUP BY cluster
		ORDER BY pos;~, {Slice => {}}
		);

		# and loop all clusters
		CLUSTER: foreach my $cluster ( @{$clusters} ) {
			$counter++;
			$chrcounter++;
			my $num = $cluster->{cluster};

			#my $name = sprintf "%s:%d-%d", $cluster->{'chr'}, 	$cluster->{'begin'}, $cluster->{'end'};
			my $name = join( ".", $cluster->{'chr'}, $chrcounter );
			if ( $opt->{hasnames} ) {
				my $genes = get_names( $dbh, $cluster->{'chr'}, $cluster->{'begin'}, $cluster->{'end'} );
				if ($genes) {
					$name = $genes;
				}
			}
			$NAMES{$run}{$name}++;    # increase finding for this name - dry / wet
			next CLUSTER if $run eq 'dry';    # done for dry

			# if the sum in the dry run was > 1, use current (wet) count as index
			if ( $NAMES{dry}{$name} > 1 ) {
				$name .= $opt->{suffix} . $NAMES{wet}{$name};
			}
			my $regions = $dbh->selectall_arrayref(
				qq~
		SELECT *, $counter AS counter, "$name" AS name 
		FROM primers WHERE cluster = $num 
		ORDER BY pos;
		;~, {Slice => {}}
			);
			foreach my $line ( @{$regions} ) {
				print join( $opt->{sep}, map { $line->{$_} } @{$opt->{export}} ), "\n";
				$STAT{OLD}{$line->{orgname}}++;
				$STAT{NEW}{$name}++;
				if ( $line->{orgname} ne $name ) {
					$STAT{DIFF}{$line->{orgname}}{$name}++;
				}
			} ## end foreach my $line ( @{$regions} )
		} ## end CLUSTER: foreach my $cluster ( @{$clusters...})
	} ## end CHR: foreach my $chr ( @{$Chr} )
} ## end RUN: for my $run ( 'dry', 'wet' )
if ( $opt->{verbose} ) {
	foreach my $KND (qw<OLD NEW>) {
		open my $statfile, '>', $KND . '_stats.tsv';
		select $statfile;

		#print $KND, "\n", '=' x length($KND), "\n";
		print join( "\t", $KND, 'COUNT' ), "\n";
		my $callable = 0;
		foreach my $annotation ( sort keys %{$STAT{$KND}} ) {
			print join( "\t", $annotation, $STAT{$KND}{$annotation} ), "\n";
			if ( $STAT{$KND}{$annotation} >= 20 ) {
				$callable++;
			}
			else {
				print STDERR "$KND not callable: $annotation (", $STAT{$KND}{$annotation}, ")\n";
			}
		} ## end foreach my $annotation ( sort keys...)
		select STDOUT;
		close $statfile;
		print STDERR "$KND callable = $callable\n";
	} ## end foreach my $KND (qw<OLD NEW>)
	select STDERR;
	dd $STAT{DIFF};
	printf "Original clusters: %5d\n", scalar keys %{$STAT{OLD}};
	printf "Created clusters:  %5d\n", scalar keys %{$STAT{NEW}};
	select STDOUT;
} ## end if ( $opt->{verbose} )
$dbh->disconnect();
exit;

# sub read_regions - this will the internal database with regions to be clustered
sub read_regions {
	my $dbh = shift;
	$dbh->do(
		<<SQL
CREATE TABLE primers (
	id 		INTEGER PRIMARY KEY AUTOINCREMENT,
	chr 	VARCHAR(45),
	pos 	INTEGER,
	dir 	INTEGER,
	base 	CHAR(1),
	orgname VARCHAR(45),
	count 	INTEGER,
	cluster INTEGER DEFAUL 'NULL'
	);
SQL
	);
	my $insert = $dbh->prepare(
		<<'SQL'
INSERT INTO primers	(chr, pos, dir, base, orgname, count) 
VALUES          	(?,   ?,   ?,   ?,    ?,       ?)
SQL
	);
	open my $reader, '<', $opt->{input};
	my $sep     = $opt->{sep};
	my $exc     = $opt->{exclude} ne '' ? $opt->{exclude} : 0;
	my @Columns = ( map { $opt->{$_} } qw<chrcol poscol dircol bascol namcol cntcol> );

	while ( my $line = <$reader> ) {
		next if $line =~ /^#/o;
		if ($exc) {
			next if $line =~ /$exc/;
		}
		chomp($line);
		my @Data = split( /$sep/, $line );
		$insert->execute( map { $Data[$_] } @Columns );
	} ## end while ( my $line = <$reader> )
	$insert->finish();
} ## end sub read_regions

# sub read_namefile - this will fill the internal database with annotations
sub read_namefile {
	my $dbh = shift;
	open my $reader, '<', $opt->{names};
	return unless ref $reader;
	my %Columns = name_columns();
	my @Columns = sort { $a <=> $b } values %Columns;
	$dbh->do(
		<<SQL
CREATE TABLE `names` (
  `chr`   varchar(255) NOT NULL,
  `begin` integer unsigned NOT NULL,
  `end`   integer unsigned NOT NULL,
  `name`  varchar(255) NOT NULL
)	
SQL
	);
	my $insert = $dbh->prepare(
		<<SQL
INSERT INTO `names` (chr, begin, end, name) 
	VALUES          (?,   ?,     ?,   ?)	
SQL
	);
	my $match = '';

	if ( $opt->{match} ) {
		$match = $opt->{match};
		$match = qr/$match/;
	}
	while ( my $line = <$reader> ) {
		next if $line =~ /^#/o;
		if ($match) {
			next unless $line =~ $match;
		}
		chomp($line);
		my @Data = split( "\t", $line );
		$insert->execute( map { $Data[$_] } @Columns );
	} ## end while ( my $line = <$reader> )
	$insert->finish();
	$opt->{hasnames} = 1;
} ## end sub read_namefile

# sub get_names  - this returns a string with annotations for provided coordinates
sub get_names {
	my $dbh   = shift;
	my $chr   = shift;
	my $left  = shift;
	my $right = shift;

	# move margins according to deadband setting
	$left += $opt->{deadband};
	$right -= $opt->{deadband};

	# ensure left is smaller than right
	if ( $left >= $right ) {
		( $left, $right ) = ( $right, $left );
	}

	#dd $chr, $left, $right;
	my $genes = $dbh->selectcol_arrayref(
		qq~
	SELECT DISTINCT(name) FROM names
		WHERE  	chr 	= "$chr" 
		AND 	end 	>= $left 
		AND 	begin 	<= $right
	GROUP BY name
	ORDER BY begin ASC;
	~
	);
	my $count = scalar @{$genes};
	return '' unless $count;
	my $result = $genes->[0];

	if ( $count == 2 ) {
		$result = join( "+", @{$genes} );
	}
	elsif ( $count > 2 ) {
		$result = join( "/", $genes->[0], $genes->[-1] );
	}
	if ( defined $opt->{prepend} and $opt->{prepend} eq 'chr' ) {
		$result = $chr . $result;
		$result =~ s/^chr//;
	}
	return $result;
} ## end sub get_names

# sub name_columns - internal routines to fill the hash with preset values
sub name_columns {

	# preset for refGene files
	if ( $opt->{type} =~ /refg/ ) {
		$opt->{nmchr}   = 2;
		$opt->{nmbegin} = 6;
		$opt->{nmend}   = 7;
		$opt->{annot}   = 12;
		$opt->{match}   = 'NM_';
		$opt->{suffix}  = '.';
	} ## end if ( $opt->{type} =~ /refg/)

	# preset for ideographic data files
	elsif ( $opt->{type} =~ /ideo/ ) {
		$opt->{nmchr}   = 0;
		$opt->{nmbegin} = 1;
		$opt->{nmend}   = 2;
		$opt->{annot}   = 3;
		$opt->{suffix}  = ';';
		$opt->{prepend} = 'chr';
	} ## end elsif ( $opt->{type} =~ /ideo/)
	return (
		'chr'   => $opt->{nmchr},
		'begin' => $opt->{nmbegin},
		'end'   => $opt->{nmend},
		'name'  => $opt->{annot}
	);
} ## end sub name_columns


