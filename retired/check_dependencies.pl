#!/usr/bin/env perl

# documentation at end of file

use strict;
use Cwd;
use Getopt::Long;
use Pod::Usage;
require CPAN;

my $VERSION = '1.13';

# Initialize values
my (
	$help,
	$print_version,
); # command line variables
my @infiles; 


# Command line options
GetOptions( 
	'help'        => \$help, # print the help
	'version'     => \$print_version, # print the version
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script check_dependencies.pl, version $VERSION\n\n";
	exit;
}


# get current working directory
my $current = Cwd::cwd();


# Check for modules
my @install_list;
foreach ( get_list_of_modules() ) {
	
	my ($module_name, $module_use) = @{$_};
	
	# expand it
	my $mod = CPAN::Shell->expand("Module", $module_name);
	
	# check for it
	printf " checking for %-30s", $module_name; 
	if ($mod->uptodate) {
		# true, it is up to date
		print " up to date\n";
	}
	
	# not up to date or not installed
	else {
		# check if we're just out of date
		if ($mod->inst_version) {
			print " " . $mod->inst_version . " but CPAN has " . 
				$mod->cpan_version . "\n";
			
			# remember for upgrading
			push @install_list, $mod;
		}
		
		# or not present at all
		else {
			print " not installed!\n";
			
			# print explanation for this module
			print "$module_use\n";
			
			# remember for installing
			push @install_list, $mod;
		}
	}
}




# check for Bio::Ensembl
# this is not available through CPAN, why!!!!??????
print " checking for Bio::EnsEMBL modules          ";
my $bio_ensembl = 0;
eval {
	require Bio::EnsEMBL::ApiVersion;
	Bio::EnsEMBL::ApiVersion->import;
	$bio_ensembl = software_version();
};
if ($bio_ensembl) {
	print " core API version $bio_ensembl\n";
}
else {
	print " not installed\n";
	print "  Optional, required for the 'get_ensembl_annotation.pl' script\n";
	print "  It is not available through CPAN, but you can obtain it from\n" .
		"  http://www.ensembl.org/info/docs/api/api_installation.html\n";
}



# offer to install
print "\n\n";
if (@install_list) {
	print "Would you like to install missing dependencies? y or n  "; 
	my $answer = <STDIN>;
	if ($answer =~ /^n/i) {
		print "OK. Here is a list of the missing dependencies I found\n";
		foreach my $mod (@install_list) {
			my $file = $mod->cpan_file;
			print " $file\n";
		}
		exit;
	}
}
else {
	print "Everything appears up to date.\n";
	exit;
}

# installation
my $not_installed = 0;
foreach my $mod (@install_list) {
	my $file = $mod->cpan_file;
	
	# ask user about upgrading
	print "  Install $file? y/n   ";
	my $answer = <STDIN>;
	if ($answer =~ /^y/i) {
		CPAN::Shell->install($file);
		print "\n#####################\n\n"
	}
	else {
		$not_installed++;
	}
	
	# change back to current directory in case CPAN moved us
	chdir $current;
}

# finished
if ($not_installed) {
	print " $not_installed recommended dependencies were not installed\n";
	print " Some biotoolbox programs may not work properly\n";
}
else {
	print " Some dependencies were installed. You may be ready.\n";
	print " Would you like to recheck?  y/n    \n";
	my $answer = <STDIN>;
	if ($answer =~ /^y/i) {
		exec $0;
	}
}





sub get_list_of_modules {
	my @list = (
		['Config::Simple', 
qq(  Required for biotoolbox configuration file
)
		],
		['Statistics::Lite', 
qq(  Required for all data collection and analysis programs
)
		],
		['Statistics::Descriptive', 
qq(  Required for analysis and graphing scripts
)
		],
		['Statistics::LineFit', 
qq(  Required for some conversion scripts
)
		],
		['Bio::Root::Version', 
qq(  Required for all BioPerl database interactions and GFF3 processing
)
		],
		['Parallel::ForkManager', 
qq(  Optional, recommended for multi-threaded execution of some scripts 
)
		],
		['GD', 
qq(  Optional, required for graphing scripts. 
)
		],
		['GD::Graph', 
qq(  Optional, required for graphing scripts
)
		],
		['GD::Graph::smoothlines', 
qq(  Optional, required for graphing smooth bezier curve graphs
)
		],
		['Archive::Zip', 
qq(  Optional, required only for converting Bar files
)
		],
		['Algorithm::Cluster', 
qq(  Optional, required only for cluster analysis
)
		],
		['DBI', 
qq(  Required for BioPerl SQL database functions
)
		],
		['DBD::mysql', 
qq(  Optional, required for interacting with MySQL databases
)
		],
		['DBD::SQLite', 
qq(  Optional, required for interacting with SQLite databases
)
		],
		['Bio::DB::Sam', 
qq(  Optional, required for working with BAM files.
)
		],
		['Bio::DB::BigFile', 
qq(  Optional, required for working with bigWig and bigBed files.
)
		],
		['Bio::DB::USeq', 
qq(  Optional, required for working with useq files.
)
		],
	);
	
	return @list;
}




__END__

=head1 NAME

check_dependencies.pl

A script to check for BioToolBox prerequisites.

=head1 SYNOPSIS

[sudo] ./check_dependencies.pl
  
  Options:
  --version
  --help
  
=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --version

Print the program version number.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program will check for Perl module dependencies used by the biotoolbox 
scripts. It will check for the installed version and the current version 
in CPAN. For missing or out of date modules, it will offer to install them 
through CPAN. Note that any dependencies may not be handled well, and you 
may wish to decline and install manually through CPAN.

If your Perl modules are located in system-owned directories, you may need to 
execute this script with root privilages. Or, if you prefer, check what is 
missing with this program and then install manually with root privilages.

This will require that CPAN is properly configured (proxies, mirrors, etc.).

More information about setting up your computer may be found on the web at
http://code.google.com/p/biotoolbox/wiki/BioToolBoxSetUp.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  

