#!/usr/bin/perl
# 
# This script is to insert CMAP energy term into AMBER topology file
# version 2.01
# Please cite: doi 10.1021/acs.jcim.7b00135
use strict;
use Getopt::Std;
use vars qw ( $perllibdir );
#$perllibdir="$ENV{PDBPTRAJDIR}" if (defined $ENV{PDBPTRAJDIR});
#($perllibdir=$0) =~ s/(\/*)[^\/]+$// if (!defined $perllibdir);
#use lib $perllibdir;
print STDERR "------------------------------------------------------------
Transplanting CMAP into AMBER PRMTOP
Please cite:  Song, et al., Journal of Chemical Information and Modeling 2017 (doi 10.1021/acs.jcim.7b00135).
Programmed by Wei Wang & Wei Ye,
and updated by Song Dong & Ji Dingjue.
------------------------------------------------------------
";
our $version = "Version 2.01
10 JAN 2017\n";
(my $pwd = $0) =~ s/(\/*)[^\/]+$//;
$pwd = "." if $pwd eq "";
my $AMBER_CMAP_parameter = "$pwd\/CMAP_parameter";
our $help = "Usage:
\$ .\/Trans_FF_IDPs [\033[1m\033[34m-p\033[0m AMBER_prmtop_file \033[1m\033[34m-o\033[0m output_file ]\\
\t\033[1m\033[34m-c\033[0m CMAP-parameter-file \033[1m\033[34m-s\033[0m\\
\t\033[1m\033[34m-v\033[0m \033[1m\033[34m-h\033[0m
Explanation:
\033[1m\033[34m-p\033[0m\tRequired.
\tInput AMBER PRMTOP file generated with tleap or xleap.
\033[1m\033[34m-o\033[0m\tRequired.
\tOutput AMBER PRMTOP file with CMAP energy term.
\033[1m\033[34m-c\033[0m\tOptional.
\tSpecify user edited CMAP parameter file.
\tDefault CMAP file is $AMBER_CMAP_parameter.
\033[1m\033[34m-s\033[0m\tSilent Mode, optional.
\tRun the script silently without any interacting UI.
\tNote: In silent mode, all the proteins in the PRMTOP would
\tbe considered as IDP.

\033[1m\033[34m-h\033[0m\tShow this help information.
\033[1m\033[34m-v\033[0m\tShow version information.

";
my %options;
my $exit = 0;
getopts("c:p:o:hvs", \%options);
if (defined $options{h}){
	print $help;
	exit;
}
if (defined $options{v}){
	print STDERR $version;
	exit;
}
unless ($options{p}){
	print $help;
	print STDERR "\033[1m\033[31m!!!\033[0mPlease specify the \033[1m\033[31minput AMBER PRMTOP file\033[0m with \033[1m\033[34m-p\033[0m flag.\n";
	$exit = 1;
}

unless ($options{o}){
	print $help if $exit == 0;
	print STDERR "\033[1m\033[31m!!!\033[0mPlease specify the \033[1m\033[31moutput AMBER PRMTOP file\033[0m with \033[1m\033[34m-o\033[0m flag.\n";
	$exit = 1;
}
if ($exit == 1){
	exit;
}
our $SILENT = 0;
if (defined $options{s}){
	$SILENT = 1;
}
if (defined $options{c}){
	$AMBER_CMAP_parameter = $options{c};
	print "Using user-defined CMAP parameter file. $AMBER_CMAP_parameter\n";	
}else{
	print "Using default CMAP parameter file. $AMBER_CMAP_parameter\n";
}
my $AMBER_TOP_file = $options{p};
my $AMBER_CMAP_TOP_file = $options{o};

#KNOWN RESIDUES & ABBREVIATIONS
our @known_residues = qw /ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HIS/;
our @known_residues_abbr = qw /A R N D C C Q E G H H H I L K M F P S T W Y V H/;
my $n = 0;
our %known_residues_abbr;
our %known_abbr;
foreach (@known_residues){
	$known_residues_abbr{$_} = $known_residues_abbr[$n];
	$known_abbr{$known_residues_abbr[$n]} = $_;
	$n ++;
}
my @Supported_residues = qw /ALA ARG ASN ASP CYS GLN GLU GLY ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HIS HIE HID HIP/;
our %Supported_residues;
foreach (@Supported_residues){
	$Supported_residues{$_} = "yes";
	$Supported_residues{$known_residues_abbr{$_}} = "yes";
}


#####################################################################
#############This part is to load the CMAP parameter################
#####################################################################

my @name_CMAP;
my @content_CMAP;
my %number_CMAP;
(open PAR, $AMBER_CMAP_parameter) || die "\033[1m\033[31m!!!\033[0mCMAP parameter file $AMBER_CMAP_parameter not found!\n";
my $n = 0;
while (<PAR>){
	if (/FLAG\s([\w]+)_MAP/){
		print "Reading CMAP parameter of $1...\n";
		$n ++;
		$name_CMAP[$n] = $known_residues_abbr{$1};
		$number_CMAP{$known_residues_abbr{$1}} = $n;
	}else{
		$content_CMAP[$n] .= $_;
	}
}
close PAR;
if ($n < 1){
	print "\033[1m\033[31m!!!\033[0mCMAP parameter FORMAT error!\n";
	exit;
}
#####################################################################
#############This part is to input the top file######################
#####################################################################
my $temp;
my @tmp;
my @Atom_name;
my @Residues_name_origin;
my @Residues_pointer;
my %LENNARD_JONES_COEF;
my $tt;
my %FLAG_read;
print STDERR "Reading PRMTOP file $AMBER_TOP_file ...";
open TOP, $AMBER_TOP_file;
while(<TOP>){
	if (/FLAG ATOM_NAME/){
		<TOP>;
		$FLAG_read{Atom_Name} = "yes";
		do{
			chomp ($temp = <TOP>);
			last if $temp =~ /^\%/;
			for($n = 0; $n <= 76; $n += 4){
				$tt = substr($temp, $n, 4);
				$tt =~ s/\s+//g;
				last if $tt eq "";
				push @Atom_name, $tt;
			}
		}until(1 == 0);
	}
}
close TOP;
unless (exists $FLAG_read{Atom_Name}){
	print STDERR "\033[1m\033[31m!!!\033[0m%FLAG ATOM_NAME not found!
Please have a check.
Make sure $AMBER_TOP_file is generated with tleap or xleap.
Usually, you don\'t need to modify this file manually.\n";
	exit;
}

open TOP, $AMBER_TOP_file;
while(<TOP>){
	if(/FLAG RESIDUE_LABEL/){
		<TOP>;
		$FLAG_read{Res_label} = "yes";
		do{
			chomp ($temp = <TOP>);
			last if $temp =~ /^\%/;
			for($n = 0; $n <= 76; $n += 4){
				$tt = substr($temp, $n, 4);
				$tt =~ s/\s+//g;
				last if $tt eq "";
				push @Residues_name_origin, $tt;
			}
		}until(1 == 0);
	}
}
close TOP;
unless (exists $FLAG_read{Res_label}){
	print STDERR "\033[1m\033[31m!!!\033[0m%FLAG RESIDUE_LABEL not found!
Please have a check.
Make sure $AMBER_TOP_file is generated with tleap or xleap.
Usually, you don\'t need to modify this file manually.\n";
	exit;
}
open TOP, $AMBER_TOP_file;
while(<TOP>){
	if(/FLAG RESIDUE_POINTER/){
		<TOP>;
		$FLAG_read{Res_pointer} = "yes";
		do{
			chomp ($temp = <TOP>);
			last if $temp =~ /^\%/;
			for($n = 0; $n <= 72; $n += 8){
				$tt = substr($temp, $n, 8);
				$tt =~ s/\s+//g;
				last if $tt eq "";
				push @Residues_pointer, $tt;
			}
		}until(1 == 0);
	}
}
close TOP;
unless (exists $FLAG_read{Res_pointer}){
	print STDERR "\033[1m\033[31m!!!\033[0m%FLAG RESIDUE_POINTER not found!
Please have a check.
Make sure $AMBER_TOP_file is generated with tleap or xleap.
Usually, you don\'t need to modify this file manually.\n";
	exit;
}
open TOP, $AMBER_TOP_file;
while(<TOP>){
	if(/FLAG LENNARD_JONES_ACOEF/){
		<TOP>;
		do{
			$temp = <TOP>;
			last if $temp =~ /^\%/;
			$LENNARD_JONES_COEF{A} .= $temp;
		}until(1 == 0);
	}
}
close TOP;
unless (exists $LENNARD_JONES_COEF{A}){
	print STDERR "\033[1m\033[31m!!!\033[0m%FLAG LENNARD_JONES_ACOEF not found!
Please have a check.
Make sure $AMBER_TOP_file is generated with tleap or xleap.
Usually, you don\'t need to modify this file manually.\n";
	exit;
}
open TOP, $AMBER_TOP_file;
while(<TOP>){
	if (/FLAG LENNARD_JONES_BCOEF/){
		<TOP>;
		do{
			$temp = <TOP>;
			last if $temp =~ /^\%/;
			$LENNARD_JONES_COEF{B} .= $temp;
		}until(1 == 0);
	}
}
close TOP;
unless (exists $LENNARD_JONES_COEF{B}){
	print STDERR "\033[1m\033[31m!!!\033[0m%FLAG LENNARD_JONES_BCOEF not found!
Please have a check.
Make sure $AMBER_TOP_file is generated with tleap or xleap.
Usually, you don\'t need to modify this file manually.\n";
	exit;
}
print "Done.\n";

#
#Specify the intrinsically disordered regions which use CMAP
#

my @Residues_name;
my @AminoAcid_index;
our @known_ligand;
our %known_ligand_template;
my $name;
my $count = 0;
foreach $name(@Residues_name_origin){
	$count ++;
	my $recognized = 0;
	if ($known_residues_abbr{$name}){
		push @Residues_name, $name;
		push @AminoAcid_index, $count;
		next;
	}elsif ($name =~ /[\+\-]/){
		push @Residues_name, "Ion";
		next;
	}elsif($name eq "WAT"){
		push @Residues_name, "Solvent";
		next;
	}elsif($name =~ /^(D|R)?[AGCTU][35]?$/){
		push @Residues_name, "NucAcid";
		next;
	}elsif($SILENT == 1){
		push @Residues_name, "Others";
		print STDERR "\033[1m\033[31m!!!\033[0mUnknown residue $name found, omitted!\n";
		next;
	}elsif(exists $known_ligand_template{$name}){
		push @Residues_name, "m".$known_ligand_template{$name};
		push @AminoAcid_index, $count;
		next;
	}
	my %other_residue_detect = ();
	foreach my $pointer($Residues_pointer[$count-1] .. $Residues_pointer[$count]-1){
		$other_residue_detect{$Atom_name[$pointer]} = "yes";
	}
	unless ($other_residue_detect{N} && $other_residue_detect{CA} && $other_residue_detect{C} && $other_residue_detect{O}){
		push @Residues_name, "Others";
		print "\033[1m\033[31m!!!\033[0mUnknown residue $name found, omitted!\n";
		next;
	}
	printf "\033[1m\033[31m!!!\033[0mAmino acid like unknown residue %s found.\nIs that a modified AA residue?[Y\/N]", $name;
	chomp (my $select = <STDIN>);
	if ($select =~ /^Y(es)?$/i){
LOOP1:	
		printf "Please tell me %s is modification of which AA.\n", $name;
		print "Single and triple character abbreviations are both supported:";
		chomp (my $select = <STDIN>);
		my $tmp = "";
		if (length $select >= 3){
			$tmp =  uc(substr($select, 0, 3));
			if (exists $known_residues_abbr{$tmp}){
				print "Taking $name as modified $tmp...\n";
				push @Residues_name, "m".$tmp;
				push @AminoAcid_index, $count;
				$known_ligand_template{$name} = $tmp;
				sleep (1);
			}else{
				print STDERR "\033[1m\033[31m!!!\033[0mI don't know what is $select\n";
				sleep (1);
				goto LOOP1;
			}
		}elsif (length $select == 1){
			$tmp = uc($select);
			if (exists $known_abbr{$tmp}){
				print "Taking $name as modified $known_abbr{$tmp}...\n";
				push @Residues_name, "m".$known_abbr{$tmp};
				push @AminoAcid_index, $count;
				$known_ligand_template{$name} = $known_abbr{$tmp};
				$recognized = 1;
				sleep (1);
			}else{
				print STDERR "\033[1m\033[31m!!!\033[0mI don't know what is $select\n";
				sleep (1);
				goto LOOP1;
			}
		}else{
			print STDERR "\033[1m\033[31m!!!\033[0mI don't know what is $select\n";
			sleep (1);
			goto LOOP1;
		}
		
	}else{
		printf "Taking %s as other ligand...", $name;
		push @Residues_name, "Others";
		sleep (1);
		print "\n";
	}
}

#Display residue information
print "\nResidues in PRMTOP:
Supported residues are labelled in \033[1m\033[34mbold blue\033[0m.\n";
print "Modified residues are labelled in \033[43m\033[37myellow background\033[0m.\n" if (keys %known_ligand_template) >0;
print "Termini are labelled as \033[1m\033[31mTER\033[0m.\n\n";
my %display_stats = (
	"AA" => 0,
	"Ion" => 0,
	"Solvent" => 0,
	"NucAcid" => 0,
	"Others" => 0,
);
my $display_string = "";
my $display_position = "";
my $position = 0;
my $length_seg = 0;
my $once_position = 0;
foreach $name(@Residues_name){
	chomp $name;
	$position ++;
	if (exists $display_stats{$name}){
		if ($display_stats{$name} == 1){
			next;
		}elsif ($display_stats{AA} == 1){
			print $display_string,"\n";
			print $display_position,"\n\n";
			$display_string = "";
			$display_position = "";
			$length_seg = 0;
			$once_position = 0;
		}
		$display_stats{AA} = 0;
		$display_stats{Ion} = 0;
		$display_stats{Solvent} = 0;
		$display_stats{NucAcid} = 0;
		$display_stats{Others} = 0;
		
		$display_stats{$name} = 1;
		print "===$name===\n";
	}elsif (exists $known_residues_abbr{$name}){
		if ($display_stats{AA} == 0){
			$display_stats{AA} = 1;
		}
		
		$display_stats{Ion} = 0;
		$display_stats{Solvent} = 0;
		$display_stats{NucAcid} = 0;
		$display_stats{Others} = 0;
		
		if (exists $Supported_residues{$name}){
			my $tmp = $known_residues_abbr{$name};
			$display_string .= "\033[1m\033[34m$tmp\033[0m";
		}else{
			$display_string .= $known_residues_abbr{$name};
		}
		$length_seg ++;
		$once_position ++;
		if (($once_position % 5) != 1){
			$display_position .= " ";
		}else{
			$display_position .= sprintf "%s",$position;
			$display_string .= " " x (length ($position)-1);
		}
		foreach my $i($Residues_pointer[$position-1]-1 .. $Residues_pointer[$position]-2){
			if ($Atom_name[$i] eq "OXT"){
				$display_string .= " \033[1m\033[31mTER\033[0m ";
				$display_position .= "     ";
				print $display_string,"\n";
				print $display_position,"\n\n";
				$display_string = "";
				$display_position = "";
				$once_position = 0;
				$length_seg = 0;
				$display_stats{AA} = 0;
				last;
			}
		}
	}else{
	#for the instance of ligand and modified amino acid residues
		if ($display_stats{AA} == 0){
			$display_stats{AA} = 1;
		}
		
		$display_stats{Ion} = 0;
		$display_stats{Solvent} = 0;
		$display_stats{NucAcid} = 0;
		$display_stats{Others} = 0;
		
		$once_position ++;
		$display_string .= sprintf "\033[43m\033[37m%-5s\033[0m",$Residues_name_origin[$position-1];
		$length_seg ++;
		$display_position .= sprintf "%-5s",$position;
		foreach my $i($Residues_pointer[$position-1]-1 .. $Residues_pointer[$position]-2){
			if ($Atom_name[$i] eq "OXT"){
				$display_string .= " \033[1m\033[31mTER\033[0m ";
				$display_position .= "     ";
				print $display_string,"\n";
				print $display_position,"\n\n";
				$display_string = "";
				$display_position = "";
				$once_position = 0;
				$length_seg = 0;
				$display_stats{AA} = 0;
				last;
			}
		}
	}
	if ($length_seg == 50){
		print $display_string,"\n";
		print $display_position,"\n\n";
		$display_string = "";
		$display_position = "";
		$length_seg = 0;
		$once_position = 0;
	}
}

my %exist_CAMP = ();
my @selected_residue_sequence;
my @real_selected_residue;
my @real_selected_supported_residues;
my $number_CMAP_residues = 0;
if ($SILENT == 1){
	foreach (@AminoAcid_index){
		push @selected_residue_sequence, $_-1;
		$number_CMAP_residues ++;
	}
	print "All residues have been taken as IDPs region\n";
}else{
	START:
	@selected_residue_sequence = ();
	$number_CMAP_residues = 0;
	print "\nPlease select amino acids for the IDPs region.
\"-\" and \",\" are suppoted. E.g: 1-20,35,39,40-80
If all the residues should be considered, type \"all\"\n";
	chomp (my $inputInformation = <STDIN>);
	if ($inputInformation =~ /^all$/i){
		foreach (@AminoAcid_index){
			push @selected_residue_sequence, $_-1;
			$number_CMAP_residues ++;
		}
		print "All residues have been taken as IDPs region\n";
	}else{
		my @str = split ",", $inputInformation;
		
		my @tmp = ();
		foreach (@str){
			print $_,"\n";
			unless (@tmp = split "-", $_) {
				print STDERR "Syntax error.\n";
				goto START;
			}
			foreach ($tmp[0] .. $tmp[-1]){
				$Residues_name[$_-1] =~ /m?(\w{3})/;
				unless (exists $known_residues_abbr{$1}){
					print STDERR "Selection Error on $_ .\n";
					goto START;
				}
				push @selected_residue_sequence, $_-1;
				$number_CMAP_residues ++;
			}
			
		}
	}
}
my $i;
my $j;
foreach $i(0..$number_CMAP_residues -1){
	foreach $j(0 .. $number_CMAP_residues - $i -2){
		if ($selected_residue_sequence[$j] > $selected_residue_sequence[$j+1]){
			($selected_residue_sequence[$j],$selected_residue_sequence[$j+1]) =($selected_residue_sequence[$j+1],$selected_residue_sequence[$j]);
		}
	}
}
foreach $i(0..$number_CMAP_residues -1){
	my $strip = 0;
	$j = $selected_residue_sequence[$i];
	foreach ($Residues_pointer[$j]-1 .. $Residues_pointer[$j+1]-2){
		if ($Atom_name[$_] eq "H3" or $Atom_name[$_] eq "OXT"){
			$strip = 1;
		}
	}
	next if $strip == 1;
	push @real_selected_residue, $i;
	#留给第九个参数的接口
	$Residues_name[$j] =~ /m?(\w{3})/;
	push @real_selected_supported_residues, $j if exists $Supported_residues{$1};
}
$number_CMAP_residues = @real_selected_supported_residues;
print "CMAP parameters will add to following residues:\n";
$count = 0;
foreach (@real_selected_supported_residues){
	$count ++;
	$Residues_name[$_] =~ /(m?)(\w{3})/;
	printf "%1s%3s%-4d ",$1,$2,$_+1;
	$exist_CAMP{$known_residues_abbr{$2}} = 1;
	print "\n" unless ($count % 6);
}
if ($SILENT == 0){
	print "\nThe new PRMTOP file with CMAP parameters will be output.
If you want to output a ff14IDPSFF PRMTOP file, please make sure the input PRMTOP file is generated with ff14SB.
Are you sure to continue?[Y\/N]";
	chomp (my $inputInformation = <STDIN>);
	unless ($inputInformation =~ /^y(es)?$/i){
		exit;
	}
}

print "\nChanges will be output to $AMBER_CMAP_TOP_file\n";

my @exist_CAMP = keys %exist_CAMP;
my $number_CMAP_type = @exist_CAMP;
my $AMBER_CMAP_Title = "%FLAG FORCE_FIELD_TYPE
%FORMAT(i2,a78)
 1 CHARMM  31       *>>>>>>>>CHARMM22 All-Hydrogen Topology File for Proteins <<
%COMMENT Transplanted from CHARMM CMAP to AMBER PRMTOP with Trans_FF_IDPs
%COMMENT for Amber ff14IDPSFF\n";
my $AMBER_CMAP_add = "%FLAG CHARMM_UREY_BRADLEY_COUNT
%COMMENT  V(ub) = K_ub(r_ik - R_ub)**2
%COMMENT  Number of Urey Bradley terms and types

%FORMAT(2i8)               
     3      1
%FLAG CHARMM_UREY_BRADLEY
%COMMENT  List of the two atoms and its parameter index
%COMMENT  in each UB term: i,k,index
%FORMAT(10i8)              
       2       5       1       3       5       1       4       5       1       
%FLAG CHARMM_UREY_BRADLEY_FORCE_CONSTANT
%COMMENT  K_ub: kcal/mole/A**2
%FORMAT(5e16.8)            
  0.00000000E+02
%FLAG CHARMM_UREY_BRADLEY_EQUIL_VALUE
%COMMENT  r_ub: A 
%FORMAT(5e16.8)            
  0.00000000E+01\n";
$AMBER_CMAP_add .= "%FLAG CHARMM_NUM_IMPROPERS
%COMMENT  Number of terms contributing to the
%COMMENT  quadratic four atom improper energy term:
%COMMENT  V(improper) = K_psi(psi - psi_0)**2
%FORMAT(10i8)              
     1
%FLAG CHARMM_IMPROPERS
%COMMENT  List of the four atoms in each improper term
%COMMENT  i,j,k,l,index  i,j,k,l,index
%COMMENT  where index is into the following two lists:
%COMMENT  CHARMM_IMPROPER_{FORCE_CONSTANT,IMPROPER_PHASE}
%FORMAT(10i8)              
      15       5      17      16       1      
%FLAG CHARMM_NUM_IMPR_TYPES
%COMMENT  Number of unique parameters contributing to the
%COMMENT  quadratic four atom improper energy term
%FORMAT(i8)                
      1
%FLAG CHARMM_IMPROPER_FORCE_CONSTANT
%COMMENT  K_psi: kcal/mole/rad**2 
%FORMAT(5e16.8)            
  0.00000000E+03  
%FLAG CHARMM_IMPROPER_PHASE
%COMMENT  psi: degrees
%FORMAT(5e16.8)            
  0.00000000E+00\n";
$AMBER_CMAP_add .= "%FLAG LENNARD_JONES_14_ACOEF\n%FORMAT(5E16.8)\n";
$AMBER_CMAP_add .= $LENNARD_JONES_COEF{A};
$AMBER_CMAP_add .= "%FLAG LENNARD_JONES_14_BCOEF\n%FORMAT(5E16.8)\n";
$AMBER_CMAP_add .= $LENNARD_JONES_COEF{B};
$AMBER_CMAP_add .= "%FLAG CHARMM_CMAP_COUNT
%COMMENT  Number of CMAP terms, number of unique CMAP parameters
%FORMAT(2I8)\n";
$AMBER_CMAP_add .= sprintf "%8d%8d\n",$number_CMAP_residues, $number_CMAP_type;
$AMBER_CMAP_add .= "%FLAG CHARMM_CMAP_RESOLUTION
%COMMENT  Number of steps along each phi/psi CMAP axis
%COMMENT  for each CMAP_PARAMETER grid
%FORMAT(20I4)\n";
foreach (1..$number_CMAP_type){
	$AMBER_CMAP_add .= sprintf "%4d",24;
}
$AMBER_CMAP_add .= "\n";

$n = 1;
my %order;
foreach (keys %exist_CAMP){
	$order{$_} = $n;
	$AMBER_CMAP_add .= "%";
	$AMBER_CMAP_add .= sprintf "FLAG CHARMM_CMAP_PARAMETER_%02d\n",$n;
	my $triple = $known_abbr{$_};
	$AMBER_CMAP_add .= "\%COMMENT    $triple   C    N    CA   C    N    CA   C    N    24\n";
	$AMBER_CMAP_add .= $content_CMAP[$number_CMAP{$_}];
	$n++;
}
$AMBER_CMAP_add .= "%FLAG CHARMM_CMAP_INDEX
%COMMENT  Atom index i,j,k,l,m of the cross term
%COMMENT  and then pointer to CHARMM_CMAP_PARAMETER_n
%FORMAT(6I8)\n";
#i,j,k,l,m
#C,N,CA,C,N
#0,1,2,3,4
my @dihedral_pointer = ();
foreach (@real_selected_supported_residues){
	my $i;
	$Residues_name[$_] =~ /m?(\w{3})/;
	my $res = $known_residues_abbr{$1};
	foreach $i($Residues_pointer[$_-1]-1 .. $Residues_pointer[$_]-2){
		#print $i,"||",$Atom_name[$i],"\n";
		$dihedral_pointer[0] = $i+1 if $Atom_name[$i] eq "C";
	}
	foreach $i($Residues_pointer[$_]-1 .. $Residues_pointer[$_+1]-2){
		$dihedral_pointer[1] = $i+1 if $Atom_name[$i] eq "N";
		$dihedral_pointer[2] = $i+1 if $Atom_name[$i] eq "CA";
		$dihedral_pointer[3] = $i+1 if $Atom_name[$i] eq "C";
	}
	foreach $i($Residues_pointer[$_+1]-1 .. $Residues_pointer[$_+2]-2){
		$dihedral_pointer[4] = $i+1 if $Atom_name[$i] eq "N";
	}
	$AMBER_CMAP_add .= sprintf "%8d%8d%8d%8d%8d%8d\n",$dihedral_pointer[0],$dihedral_pointer[1],$dihedral_pointer[2],$dihedral_pointer[3],$dihedral_pointer[4],$order{$res};
}

open TOP,$AMBER_TOP_file;
open OUT,">$AMBER_CMAP_TOP_file";
my $a = <TOP>;
print OUT $a;
print OUT $AMBER_CMAP_Title;
while(<TOP>){
	print OUT $_;
}
print OUT $AMBER_CMAP_add;
print OUT "\n";
close TOP;
close OUT;
print "$AMBER_CMAP_TOP_file written.\n";
