#!/usr/bin/perl

use strict;
use Cwd;
require CPAN;
my $VERSION = '1.6.4';

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
	require Bio::EnsEMBL::Registry;
	Bio::EnsEMBL::Registry->import;
	$bio_ensembl = 1;
};
if ($bio_ensembl) {
	print " installed\n";
}
else {
	print " not installed\n";
	print "  Optional, only required for the script 'get_ensembl_annotation.pl'\n";
	print "  It is not available through CPAN, but you can obtain if from\n" .
		"  http://www.ensembl.org/info/docs/api/api_installation.html\n";
}



# offer to install
print "\n\n";
if (@install_list) {
	print "Would you like me to help install missing dependencies? y or n  "; 
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
	print "Everything is up to date. You are good to go!\n";
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
	print " You should be good to go!\n";
}





sub get_list_of_modules {
	my @list = (
		['Config::Simple', 
qq(  Required for biotoolbox configuration file
)
		],
		['IO::Zlib', 
qq(  Required for reading and writing gzipped files
)
		],
		['Statistics::Lite', 
qq(  Required for numerous analysis programs
)
		],
		['Statistics::Descriptive', 
qq(  Required for some analysis and graphing scripts
)
		],
		['Bio::Root::Version', 
qq(  Required for BioPerl database functions, data analysis, and GFF3 conversions
)
		],
		['GD', 
qq(  Optional, required for graphing scripts. Requires GD2 libraries to be installed
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
qq(  Optional, required only for script bar2wig.pl
)
		],
		['Algorithm::Cluster', 
qq(  Optional, required only for the script run_cluster.pl
)
		],
		['DBI', 
qq(  Required for BioPerl SQL database functions
)
		],
		['DBD::mysql', 
qq(  Optional, but required for interacting with MySQL databases
)
		],
		['DBD::SQLite', 
qq(  Optional, but required for interacting with SQLite databases
)
		],
		['Bio::DB::Sam', 
qq(  Optional, but required for working with Next Generation Sequencing BAM 
  data files. Requires building the SamTools C libraries.
)
		],
		['Bio::DB::BigFile', 
qq(  Optional, but required and highly recommended for working with bigWig and 
  bigBed files. Requires building Jim Kent's source tree and executables.
)
		],
	);
	
	return @list;
}




__END__

=head1 NAME

check_dependencies.pl

=head1 SYNOPSIS

check_dependencies.pl
  
=head1 OPTIONS

None

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



