#!/usr/bin/perl

use strict;
use Cwd;
require CPAN;

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
	print "I can now help you install missing dependencies if you wish\n\n" 
}
else {
	print "Everything is up to date. You are good to go!\n";
	exit;
}

# installation
foreach my $mod (@install_list) {
	my $file = $mod->cpan_file;
	
	# ask user about upgrading
	print "  Install $file? y/n   ";
	my $answer = <STDIN>;
	if ($answer =~ /^y/i) {
		$mod->install;
		print "\n#####################\n\n"
	}
	
	# change back to current directory in case CPAN moved us
	chdir $current;
}
print " You are now good to go!\n";





sub get_list_of_modules {
	my @list = (
		['Config::Simple', 
qq(  Required for biotoolbox configuration file
)
		],
		['Archive::Zip', 
qq(  Required for script bar2wig.pl
)
		],
		['Statistics::Lite', 
qq(  Required for numerous analysis programs
)
		],
		['Statistics::Descriptive', 
qq( Required for some analysis and graphing scripts
)
		],
		['GD::Graph', 
qq(  Required for graphing scripts
)
		],
		['GD::Graph::smoothlines', 
qq(  Optional, for graphing smooth bezier curve graphs
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
		['Bio::Root::Version', 
qq(  Required for BioPerl database functions, data analysis, and GFF3 conversions
)
		],
		['Algorithm::Cluster', 
qq(  Optional, required only for the script run_cluster.pl
)
		],
		['Bio::Graphics', 
qq(  Optional, but required for working with wig files and GBrowse
)
		],
		['Bio::DB::BigFile', 
qq(  Optional, but required and highly recommended for working with bigWig and 
  bigBed files. Requires building Jim Kent's source tree and executables.
)
		],
		['Bio::DB::Sam', 
qq(  Optional, but required for working with Next Generation Sequencing BAM 
  data files. Requires building the SamTools C libraries.
)
		],
		['Bio::Graphics::Browser2', 
qq(  Optional, but recommended for genome browsing
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

This program will check for module dependencies used by the biotoolbox 
scripts. It will check for the installed version and the current version 
in CPAN. For missing or out of date modules, it will offer to install them 
through CPAN.


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



