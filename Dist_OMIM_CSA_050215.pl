#!/usr/bin/perl -w 

use strict;
use warnings; 
use File::Basename qw/ basename /; 
use IO::Handle;
use Data::Dumper;
use List::Util qw/ any min max /;
use List::MoreUtils qw/ minmax /;
use Math::Complex;
use Carp qw/ confess /; # like die, tells more info that calls error line. route through program. 
use Millie::Distance qw/ calculate_distance /; 

=head1 NAME

COG distances between OMIM and CSA residues.

=cut

##########ACCEPTS UNPARSED DATA OF CSA AND PDB ENTRIES########################

#usuage code?why?

my $PROGNAME = basename($0); #how do i do this for 2 files? WHY DO THIS???
my $DEBUG_FLAG = 1;
my $USAGE =<<"_USAGE"; 

usage: $PROGNAME <CSA_LIT_DATA_121113.txt> <OMIM_PDB_JOIN_sorted.txt>

parses CSA file

_USAGE

################ OPEN DATA FILES #########################
 
my $CSA_file = $ARGV[0] or die "$USAGE\n\nError: expected CSA file argument";
my $OMIM_UNI_file = $ARGV[1] or die "$USAGE\n\nError: expected OMIM file argument";


# INPUT data

############### PARSE OMIM AND CSA DATA ###################

D("Parsing OMIM data...");
my $omim_rows = parse_omim_file( $OMIM_UNI_file );
D("   ...found " . scalar(@$omim_rows). " records");

# $omim_rows = [
#   { pdb_id => '', ... },
#   { pdb_id => '', ... },
# ]

D("Parsing CSA data...");
my $csa_rows = parse_csa_file( $CSA_file );
D("   ...found " . scalar(@$csa_rows). " records");

# $csa_rows = [
#   { pdb_id => '', ... },
#   { pdb_id => '', ... },
# ]

# PROCESS data

################## REORGANISE CSA ROWS BASED ON PDB ################

my %csa_rows_by_pdb_id;
for my $csa_row (@$csa_rows) 
	{

	# the PDB id of the CSA row that we're looking at
	my $pdb_id = $csa_row->{pdb_id};

	# if we haven't seen any entries for this PDB, then 
	# initialise this key in the hash to be an empty array
	$csa_rows_by_pdb_id{ $pdb_id } ||= [];

	# now we push the CSA row onto the array of entries for this PDB. push entries on it.	
	push @{ $csa_rows_by_pdb_id{$pdb_id }}, $csa_row;
	
	}

# just to check: print " CSA rows by PDB_id: ";
# just to checkprint Dumper(\%csa_rows_by_pdb_id); 

# %csa_rows_by_pdb_id = (
# 	"1cuk" => [
#		{ pdb_id => "1cuk", site => "", },
#		{ pdb_id => "1cuk", site => "", },
#	],
# 	"1jpw" => [
#		{ pdb_id => "1jpw", site => "", },
#		{ pdb_id => "1jpw", site => "", },
#	]
# );


################## REORGANISE OMIM ROWS BASED ON PDB ##################

my %omim_rows_by_pdb_id;
for my $omim_row (@$omim_rows)
	{

	my $pdb_id = $omim_row->{pdb_id};

	#print "OMIM:". Dumper($omim_row);

	$omim_rows_by_pdb_id{$pdb_id} ||= [];

	push @{$omim_rows_by_pdb_id{$pdb_id}},$omim_row;
	#die
	}


################## MATCH CSA and OMIM pdb_ids ########################.

my @omim_pdb_id = keys (%omim_rows_by_pdb_id);
my @csa_pdb_id = keys (%csa_rows_by_pdb_id);

#print "OMIM pdb_ids: @omim_pdb_id \n OMIM ";
#print "CSA pdb_ids: @omim_pdb_id \n CSA ";

my @common_omim_csa_pdb = ();

foreach ( keys %omim_rows_by_pdb_id)
{	
	push (@common_omim_csa_pdb, $_) if exists $csa_rows_by_pdb_id{$_};
}

#print " Common PDBs in OMIM and CSA:@common_omim_csa_pdb \n";

#print "  \n The number of OMIM pdbs are: ". scalar(@omim_pdb_id). "\n\n";
#print "  \n The number of CSA pdbs are: ". scalar(@csa_pdb_id). "\n\n";
#print "  \n The number of these common OMIM and CSA PDBs are: ". scalar(@common_omim_csa_pdb)."\n\n";
#print "  \n finished matching CSA PDBs to OMIM PDBs data \n\n"; 


#foreach my $pdb (@common_omim_csa_pdb)
#{
	#print "$pdb \n";
#}

#die;


############### GET RESIDUE ATOMS FOR THE MUTANT AND CSA SITE FROM PDB FILE ##############
	
## IMPORTANT: UNIPROT sequence to PDB structure mapping ##
		# check ($omim_row->{residue_position} = $pdb_row->{res_num})
		# only selected entries where valid = "t" in andrew martin were included.
		# UNIPROT PDB mapping has been done by SIFTS method in EMBL website to get PDB mappings.
	
# get rid of the insertion code ( so no letters) from the $omim_row->{residue_position}. Therefore, substitute the letters for nothing.
		# foreach $omim_row->{pdb_res_num}, s/$[A-Z]//;
		# this can then be entered as a number argument and refer to the amino acid position.

#foreach my $omim_row (@$omim_rows) #get rid of insertion number after pdb_res_num
#{ 
	
	#printf "BEFORE pdbe_res_num: %s \n \n ", $omim_row->{pdb_res_num};
	#$omim_row->{pdb_res_num} =~ s/[A-Z]//; #get rid of the letter on the end to have the res_num
	#printf "AFTER pdbe_res_num: %s \n \n ", $omim_row->{pdb_res_num};
#}

#print " NEXT STEP: Getting residue atoms....\n\n";


#print " matching " . scalar(@$omim_rows). " omim rows to " . scalar(@common_omim_csa_pdb). " csa_omim_pdbs......... \n ";


my @omim_rows_with_csa;
my $ct = 0;
for my $omim_row (@$omim_rows) # create omim rows with matched_common pdbs
{
	#print "  ...OMIM row $ct\n" if $ct++ % 1000;

	#my $found_pdb_id = grep {
	#		my $are_they_equal = $omim_row->{pdb_id} eq $_;
	#		printf "Comparing \$omim_row->{pdb_id}:'%s' to \$_:'%s' ? %s \n", 
	#		$omim_row->{pdb_id}, $_, $are_they_equal ? "YES" : "NO ";
	#		$are_they_equal; #eval as true or false
	#	}
	#	@common_omim_csa_pdb;

	# List::Util provides 'any' which will return true if any of the items match pdb_id

	my $found_pdb_id = any { $omim_row->{pdb_id} eq $_ } @common_omim_csa_pdb;

	if ($found_pdb_id) #true or false/
	{
		push @omim_rows_with_csa, $omim_row;
	}
}

#print " matched " . scalar(@omim_rows_with_csa). " omim rows to common \n ";

# to create a dummie for pdb_id: 3ffg. already done for CSA, do OMIM>


# THIS LAST STEP TAKES AGES ASWELL. 


# @omim_rows_with_csa should contain unique PDB_ids as made from hash: %omim_rows_by_pdb_id. 


###################### PARSE PDB ######### GET RESIDUE ATOMS ###### CALCULATE DISTANCE #######

my $omim_res_atoms;
my $csa_res_atoms;
my $cogdistance;


foreach my $pdb_id(@common_omim_csa_pdb)
{


	my $pdb_file = "/cath/data/current/pdb/$pdb_id"; #...get pdb file name from $pdb_id;
	my $pdb_atoms = parse_pdb_file( $pdb_file ); #as a pdb row is an atom entry, parse it to obtain atoms. created a hash reference
	D("PARSING PDB file: $pdb_file... \n\n\n");
	#D("get_pdb_res_atoms.parsing PDB file: ...found " . scalar (@$pdb_atoms) . " ATOM records");

	my $omim_rows_p = $omim_rows_by_pdb_id{$pdb_id};

	OMIM: foreach my $omim_row (@$omim_rows_p)

			{

			$omim_res_atoms = get_pdb_res_atoms($pdb_atoms,$omim_row->{chain_id}, $omim_row->{pdb_res_num} );
			

			if (not $omim_res_atoms) # not defined
					{
					next OMIM;
					}

			my ($omcogx, $omcogy, $omcogz) = get_centre_of_gravity($omim_res_atoms);
			
			
			#printf "OMIM resisue COG co-ordinates for PDB_ID: $pdb_id and RESIDUE:%s, $omcogx, $omcogy, $omcogz \n\n",$omim_row->{pdb_res_num};
			

			my $csa_rows_p = $csa_rows_by_pdb_id{$pdb_id};

			CSA: foreach my $csa_row (@$csa_rows_p) #get csa_res_atoms
			
				{
				
				$csa_res_atoms = get_pdb_res_atoms( $pdb_atoms,$csa_row->{chain_id}, $csa_row->{res_num} );

				if (not $csa_res_atoms) # not defined
					{
					next CSA;
					}

				my ($cscogx, $cscogy, $cscogz) = get_centre_of_gravity($csa_res_atoms);

				#printf "CSA resisue COG co-ordinates for PDB_ID: $pdb_id and RESIDUE:%s, $cscogx, $cscogy, $cscogz \n\n",$csa_row->{res_num};
				

				#print "NEXT STEP: DISTANCE CALCULATION..........\n"; # need to get shortest distance between the collection of atoms between OMIM and CSA residues.
				$cogdistance = calculate_distance_cog($omcogx, $omcogy, $omcogz,$cscogx, $cscogy, $cscogz);

				my @aadistance = Millie::Distance::calculate_distance($omim_res_atoms,$csa_res_atoms);


				printf " \n PDB: $pdb_id \n ATOM-ATOM Distance for OMIM residue %s %s WT: %s to mutant: %s and CSA residue %s %s %s:\n ", 
				$omim_row->{pdb_res_num}, 
				$omim_row->{chain_id},
				$omim_row->{native_aa},
				$omim_row->{mutant_aa},
				$csa_row->{res_num},
				$csa_row->{res_type},
				$csa_row->{chain_id};

				#foreach my $t (@aadistance)
				#{
					#printf " Atom_Atom Distance:  $t \n "; 
	
				#}

				my $min_aa_dist  = min @aadistance; 
				my $max_aa_dist = max @aadistance;
				print " \n Min atom-atom distance:  $min_aa_dist \n ";
				print " \n Max atom-atom distance:  $max_aa_dist \n ";
				print " \n COG residue distance:  $cogdistance \n ";
				print " \n \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\- \n";
				#printf " %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t $cogdistance \t \n",
				#$omim_row->{pdb_id},
				#$omim_row->{native_aa},
				#$omim_row->{pdb_res_num},
				#$omim_row->{mutant_aa},
				#$omim_row->{chain_id},
				#$csa_row->{res_type},
				#$csa_row->{res_num},
				#$csa_row->{chain_id},	
				#$csa_row->{chem_funct},			
				#$omim_row->{description};    #can i create my own hash called Distance? and refer to it in other parts of my code, as to combine other bits of data?

				}


			}

		
}


#------ DOMAIN IDENTIFIERS---

#my $domain_file = "./DOMAINS/Domain_Common_pdb_omcsal.txt"; #contains info on pdb and domain id.

#my $domain_file = $ARGV[0] or die " can't find argument file"; 

#my $domain_rows = parse_domain_pdb_file($domain_file);




############################### SUBROUTINES############################################

sub D {
	my $str = shift;
	if ($DEBUG_FLAG) {
		print "$str\n";
	}
}

=head2 get_pdb_res_atoms( $pdb_id, $chain_id, $res_label )

Return residue information (including atom x, y, z) in the form of an array reference.from either OMIM and/or CSA residue

	{
		res_aa    => 'G',
		res_label => '134A',
		atoms => [
			{ type => 'CA', x => '', y => '', z => '' },
			{ type => 'CB', x => '', y => '', z => '' },
		]
	}

=cut

sub get_pdb_res_atoms {

	my ($pdb_atoms, $chain_id, $res_num) = @_;
	
	#D("get_pdb_res_atoms($chain_id, $res_num) \n\n");

	my @res_atoms;
	
	RESATOM: for my $atom ( @$pdb_atoms ) 
	{
		if ( $atom->{chain_id} eq $chain_id && $atom->{res_label} eq $res_num ) 
			{
			push @res_atoms, $atom;		# as there are many atom entries for each residue in a chain.
			#print $atom; 
			}
			else 	
			{
			next;
	     		}	
	} 



	if (scalar(@res_atoms) == 0)
	{
	#warn " couldn't find any atom records macthing chain: $chain_id and residue: $res_num ";
	return undef; #returning nothing.
	next RESATOM;
	}
	#D("get_pdb_res_atoms.matched " . scalar(@res_atoms) . " atoms");

	return \@res_atoms;
	
}



=head2 calculate_distance( $res1, $res2 )

Take OMIM residue atom corodinates and calculates distances to CSA residue atoms co-ordinates.

=cut

sub calculate_distance {
	
my ($omim_res_atoms,$csa_res_atoms)= @_;  

my $omim_X_coord; 
my $omim_Y_coord; 
my $omim_Z_coord;
my $csa_X_coord; 
my $csa_Y_coord; 
my $csa_Z_coord;
my @distance_omim_csa;
my $distance_omim_csa;
foreach my $omim_atom (@{$omim_res_atoms})
{
	$omim_X_coord = $omim_atom ->{x_coord};
	#printf " OMIM X-corordinates for residue %s : $omim_X_coord \n", $atom->{pdb_res_num};
        $omim_Y_coord = $omim_atom ->{y_coord};
	$omim_Z_coord = $omim_atom ->{z_coord};
	#printf " OMIM atom: %s \n", $omim_atom->{atom_name};

	foreach my $csa_atom (@{$csa_res_atoms})
	{
		$csa_X_coord = $csa_atom ->{x_coord};
		$csa_Y_coord = $csa_atom ->{y_coord};
		$csa_Z_coord = $csa_atom ->{z_coord};
		#printf " CSA atom: %s \n", $csa_atom->{atom_name};

		#printf " Calculating distance between the OMIM: %s and CSA atom: %s co-ordinates \n\n",$omim_atom->{atom_name},$csa_atom->{atom_name};

		my $X_dist = ( $csa_X_coord - $omim_X_coord)**2;
		#print " X-dist OMIM/CSA:$X_dist ";
		my $Y_dist = ( $csa_Y_coord - $omim_Y_coord)**2; 
		#print " Y-dist OMIM/CSA:$Y_dist ";
		my $Z_dist = ( $csa_Z_coord - $omim_Z_coord)**2;
		#print " Z-dist OMIM/CSA:$Z_dist ";

		my $total = $X_dist + $Y_dist + $Z_dist;

		$distance_omim_csa = sqrt($total);
		#print " distance calculated:$distance_omim_csa \n ";
		push @distance_omim_csa, $distance_omim_csa;

	}

	
	
}

return @distance_omim_csa;


} 
=head2 calculate_distance_cog( $res1, $res2 )

Take OMIM and CSA COG coordinates and works out distances between them. 

=cut

sub calculate_distance_cog {
	
#my $omcogx = $_[0];
#my $omcogy = $_[1]; 
#my $omcogz = $_[2];
#my $cscogx = $_[3];
#my $cscogy = $_[4]; 
#my $cscogz = $_[5]; 

my ($omcogx, $omcogy, $omcogz,$cscogx, $cscogy, $cscogz) = @_; 

#print " 6 arguments for COG-distance: \n"; # the coordinates are the same for OMIM and CSA cog.....not good as not same residue in pdb. 
#print "$omcogx, $omcogy, $omcogz, $cscogx, $cscogy, $cscogz \n";

my $Xcdist = ($omcogx - $cscogx)**2;
my $Ycdist = ($omcogy - $cscogy)**2;
my $Zcdist = ($omcogz - $cscogz)**2;

#print " ...calculating distance: $Xcdist, $Ycdist , $Zcdist...\n"; #NOT FINE keeps giving me 0. 

my $totalc = ($Xcdist + $Ycdist + $Zcdist + 0);

my $cogdistance = sqrt($totalc);

#print "$cogdistance \n\n";
return $cogdistance;

} 


=head1

 COG: Centre of gravity of the residue atoms such as: OMIM_RES_ATOMS and CSA_RES_ATOMS.  

=cut

sub get_centre_of_gravity{

my $res_atom_lines = shift;

confess "error: expected 1 or more atom line." unless scalar(@$res_atom_lines) > 0; 
    
#define x,y,z centre of gravities
my ( $x_sum, $y_sum, $z_sum );

my $num_of_atom_lines = 0;
    
foreach my $res_line ( @$res_atom_lines ) 
{
	#print Dumper($res_line);
	
        #extract x, y and z coordinates (-1 as substr starts from 0)
        my $x = $res_line->{'x_coord'};
        my $y = $res_line->{'y_coord'};
        my $z = $res_line->{'z_coord'};
       
        #print "x: $x, y: $y, z: $z\n";        
       
        #sum the coordinate values for each atom
        $x_sum += $x;
        $y_sum += $y;
        $z_sum += $z;
        
        $num_of_atom_lines++; #increment
}

#calculate centre of gravity, i.e. average x/y/z value
my $x_cog = ( $x_sum / $num_of_atom_lines );
my $y_cog = ( $y_sum / $num_of_atom_lines );
my $z_cog = ( $z_sum / $num_of_atom_lines );
    
#round to 3s.f.
#my $x_cog_3sf = round_to_3sf( $x_cog );
#my $y_cog_3sf = round_to_3sf( $y_cog );
#my $z_cog_3sf = round_to_3sf( $z_cog );
    
return ( $x_cog, $y_cog, $z_cog);

}


=head2 parse_csa_file( $filename )

Parse the CSA data file and return an array of hashrefs

=cut


#--------------------- CSA PARSE AND CREATE DATA STRUCTURE--------------------

sub parse_csa_file 
{
	my $csa_filename = shift;

	open (my $CSA_FH, '<', $CSA_file) or die " failed to open $csa_filename: $!";

	my @CSA;
	my $first_line = $CSA_FH->getline;
	while (my $CSA_line = $CSA_FH->getline)
	{
		chomp $CSA_line;
	
		# for each line split based on the "," character to create columns
		my @cols = split /,/, $CSA_line;
		my %CSA_entry = (	
			'pdb_id'	=> $cols[0],
			'site_num'	=> coerce_str2num($cols[1]),
			'res_type'	=> $cols[2],
			'chain_id'	=> $cols[3],
			'res_num'	=> coerce_str2num($cols[4]),
			'chem_funct'	=> $cols[5],
			'evid_type'	=> $cols[6], 
			'lit_entry'	=> $cols[7]
		);


## GET RID OF HOM ENTRIES ###

		push @CSA, \%CSA_entry;	#put a reference to this hash to our entries in an array- list of pointers to data.

	}
	return \@CSA;
}

#----------------------- OMIM_UNIPROT FILE PARSE AND CREATE DATA STRUCTURE-------------

=head2 parse_omim_file( $omim_filename )

Parse the OMIM data file and return an array of hashrefs

=cut

sub parse_omim_file {
	my $omim_filename = shift;

	open (my $OMIM_UNI_FH, '<', $omim_filename) or die " failed to open $omim_filename: $!";


	my @OMIM_UNI;

	while (my $OMIM_UNI_line = $OMIM_UNI_FH->getline)
	{
		chomp $OMIM_UNI_line;
		
		my @cols = split /\t/, $OMIM_UNI_line;
		my %OMIM_UNI_entry = (		
			'omim_id'		=> $cols[0],
			'sub_code'		=> $cols[1],
			'uniprot_acc'		=> $cols[2], #
			'residue_position'	=> $cols[3], ###
			'native_aa'		=> $cols[4],
			'mutant_aa'		=> $cols[5],
			'valid'			=> $cols[6],
			'description'		=> $cols[7],
			'native_aa_short'	=> $cols[8],
			'pdb_id'		=> $cols[9],
			'chain_id'		=> $cols[10],
			'pdb_res_num'		=> $cols[11],#I WANT! concat of pdb sequence and insertion code (e.g 100A) can't coerce unfortunately.
			'accession'		=> $cols[12],#
			'uniprot_res_num'	=> $cols[13],###
			'pdbe_res_num'		=> $cols[14],
			'aa'			=> $cols[15],###
			'domain_id'		=> $cols[16]
		);

		push @OMIM_UNI, \%OMIM_UNI_entry;
	
	}

	return \@OMIM_UNI;
}


#----------------------- PDB parse--------------------

=head2 parse_pdb_file($pdb_id)

parses pdb file and stores data as an array hash reference.
=cut


sub parse_pdb_file {

	my $PDB_file = $_[0] or die " $USAGE ";

	open (my $PDB_FH, '<', $PDB_file) or die " Error: failed to open file $PDB_file: $!";

	my @PDB; 

	while (my $line = $PDB_FH->getline) { 

		chomp $line;

		next unless $line =~ m/^ATOM/; # next line if dont' contain..
		my %PDB_entry = (
	
			'atom'		=> substr($line,0,6), # position before start, grab how many characters.
			'atom_number'	=> coerce_str2num(substr($line,6,6)),
			'atom_name'	=> substr($line,12,3),
			'alt_loc_id'	=> substr($line,16,1), # alternate location indicator
			'res_name'	=> substr($line,17,3),
			'chain_id'	=> substr($line,21,1), 
			'res_number'	=> coerce_str2num(substr($line,22,4)),
			'insert_code'	=> substr($line,26,1), # some PDBs have it
			'x_coord'	=> coerce_str2num(substr($line,30,8)),
			'y_coord'	=> coerce_str2num(substr($line,38,8)),
			'z_coord'	=> coerce_str2num(substr($line,46,8)),
			'occupancy'	=> coerce_str2num(substr($line,54,6)),
			'b_factor'	=> coerce_str2num(substr($line,60,6)),
			'charge_atom'	=> substr($line,78,2),
		);

	# add residue label which consists of res_number and insert_code
		# res_number => '11', insert_code => 'C'        res_label = '11C'
		# res_number => '11', insert_code => ' '        res_label = '11'
		$PDB_entry{res_label} = $PDB_entry{res_number} . ( $PDB_entry{insert_code} ne ' ' ? $PDB_entry{insert_code} : '' );

	# a new hash for each line, new set of recoreds for each line, no overwrite. atomlines = hash number
	# do checks to see if a chain id is a number.....if not then say!

		#D("Parsing PDB line...");
		#D("$line");
		#D("as:" . Dumper(\%PDB_entry));
		#exit;

	#put a reference to this hash to our entries in an array- list of pointers to data.

		push @PDB, \%PDB_entry;

	}
	#print " Parsed PDB file.....";
	return \@PDB; 
}



#sub parse_domain_pdb_file {

	#my $domain_filename = shift;

	#open (my $DOMAIN_FH, '<', $domain_filename) or die " failed to open $domain_filename: $!";


	#my @DOMAIN_ID;

	#while (my $DOMAIN_line = $DOMAIN_FH->getline)
	#{
		#chomp $DOMAIN_line;
		
		#my @cols = split /\t/, $DOMAIN_line; 
		#my %DOMAIN_entry = (		

			#'pdb_id'		=> $cols[0],
			#'domain_id'		=> $cols[1],
			#'superfam_id'		=> $cols[2],
			#'funfam_id'		=> $cols[3]
			
		#);

		##$DOMAIN_entry{chain_id} = $DOMAIN_entry{domain_id}[4];

		#push @DOMAIN_ID, \%DOMAIN_entry;
	
	#}
	
	#return \@DOMAIN_ID;

#}





sub coerce_str2num{
	my $str = shift;
	$str =~ s/^\s+//; # get rid of whitespace at begining and end of string
	$str =~ s/\s+$//;

	if ( $str !~ m/^-?[0-9]+(\.[0-9]+)?/ ) {
		die " ERROR: failed to coerce string '$str' to a number"
	}
	return $str + 0; #perl knows to return to return a number. math operation. 
}



### LATER PYMOL STUFF ######

#my @csa_pymol_site = get_csa_pymol_site($input_pdb);


#=head2 get_csa_pymol_site ($pdb_id)

#return an array contining the functional site: number, res_type and res_num.

#=cut

#sub get_csa_pymol_site {

#my $pdb = $_[0];
#my @csa_site;

#foreach my $csa_row ( @$csa_rows )
#{
#	if  ($csa_row eq $pdb)
#	{
#	@csa_site = ($csa_rows->{site_num},$csa_rows->{res_type}, $csa_rows->{res_num});
#	return @csa_site;
#	}
#	printf "CSA site is: %s \n ", @csa_site; 
#	printf "In the PDB: %s \n ", $csa_rows->{pdb_id};
#}

#}
#}

#printf "OMIM: %s UNIPROT:%s PDB:%s %s->%s '%s'\n", 
		#$omim_row->{omim_id},
		#$omim_row->{uniprot_acc},
		#$omim_row->{pdb_id},
		#$omim_row->{native_aa},
		#$omim_row->{mutant_aa},
		#$omim_row->{description};



#print <<"PYMOL"; 

#load /cath/data/current/pdb/#PDB_ID

#show cartoon
#hide lines

#select atom1, (residue_position and (name ca))
#select atom2, (@csa_pymol_site)

#color red, atom1
#color blue, atom2

#show sticks, residprintf "CSA site is: %s \n ", @csa_site; 
	
#show sticks, CSA_site

#distance d, atom1, atom2 

#PYMOL








