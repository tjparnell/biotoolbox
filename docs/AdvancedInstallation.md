# Bio::ToolBox - Installation

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## Advanced Installation

This is an advanced installation guide for getting a complete installation. 

## TLDR Brief guide

This is a no-nonsense, quick guide for those who already know what they're doing on 
an established Linux system with a modern Perl installation, and know how to adjust 
accordingly for their system. If that doesn't describe you, skip ahead to the 
[Detailed guide](#detailed-guide).

- Install external libraries

    These are external libraries that must be compiled and installed prior to installing 
    the Perl adapters. See [External libraries](#external-libraries) below for more 
    information. 

    Install [HTSlib](https://github.com/samtools/htslib) for Bam file support. You may 
    already have it; version 1.9 is officially recommended, but later versions appear
    to work just fine and are preferred, as it includes Cram support. It defaults to 
    installing in `/usr/local`. 

        curl -o htslib-1.19.tar.bz2 -L https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 
        tar -xf htslib-1.19.tar.bz2
        cd htslib-1.19
        make && make install
        cd ..

    Install [libBigWig](https://github.com/dpryan79/libBigWig) for bigWig and bigBed 
    support. It defaults to installing in `/usr/local`.
    
        curl -o libBigWig-0.4.7.tar.gz -L https://github.com/dpryan79/libBigWig/archive/0.4.7.tar.gz 
        tar -xf libBigWig-0.4.7.tar.gz
        cd libBigWig-0.4.7
        make && make install
        cd ..

- Install Perl modules

    These are the Perl modules that should be explicitly installed. A custom 
    [minimal version](https://github.com/tjparnell/bioperl-live/tree/minimal-tjparnell)
    of BioPerl is preferentially used. See [Perl Modules](#perl-modules) below for more 
    information. This assumes [CPAN Minus](https://metacpan.org/pod/App::cpanminus) is 
    installed. 

        curl -O -L https://github.com/tjparnell/bioperl-live/releases/download/minimal-v1.7.8/Minimal-BioPerl-1.7.8.tar.gz
        curl -o bio-db-big-master.tar.gz -L https://github.com/Ensembl/Bio-DB-Big/tarball/master
        cpanm Minimal-BioPerl-1.7.8.tar.gz
        cpanm Template::Tiny bio-db-big-master.tar.gz
        cpanm Bio::DB::HTS
        cpanm Parallel::ForkManager Set::IntervalTree Set::IntSpan::Fast Bio::ToolBox

- External applications

    These are external helper applications for converting to and from bigWig and bigBed 
    formats. Assumes installation to `/usr/local/bin`.

        curl -O -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
        curl -O -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
        curl -O -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
        chmod +x wigToBigWig bigWigToWig bedToBigBed
        mv wigToBigWig bigWigToWig bedToBigBed /usr/local/bin/

## Detailed guide

This assumes installation on a Linux work station with available standard compilation 
tools. Installation on MacOS (x86_64) is also possible with Xcode Command Line Tools 
installed; see see [MacOS Notes](MacOSNotes.md) for additional guidance.

## Perl installations and locations

As a Perl package, BioToolBox needs to be installed under a Perl installation. It is 
not dependent on a specific Perl release version, although later releases (5.16 
or newer) are preferred. Nearly every unix-like OS (Linux, MacOS) includes a system Perl 
installation. If not, one can be installed, either from the OS package manager or 
from source.

Follow one of these options. 

### Home library

When you want to use the system-installed Perl (often `/usr/bin/perl`), but do not have 
write permissions to the system, you can install packages in your home directory. To do 
this, you should first install [local::lib](https://metacpan.org/pod/local::lib), which 
sets up a `perl5` directory for local (home) module installations. The path is set 
appropriately by adding a statement to your home `.profile` or other equivalent file as 
described in the documentation. This can also be used for targeted, standalone 
installations; adjust accordingly. For example, the following command will install 
`local::lib` and the CPAN Minus application

    curl -L https://cpanmin.us | perl - -l $HOME/perl5 local::lib App::cpanminus \
    && echo 'eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.profile \
    && . ~/.profile

### Custom installation

When the system Perl is old (because many vendor OS Perl installations are sadly out of 
date), or you want or need to install a newer, modern Perl, but cannot or do not want to 
overwrite the system Perl, then you can and should install a new Perl version. This can 
be installed anywhere you have read/write access, including your home directory or 
wherever. While a new Perl version can be manually downloaded and installed from the 
main [Perl](https://www.perl.org) site, there are easier ways.

An alternate package manager may be used to install a Perl version in a generally
available location. For example, MacOS users can easily install a modern Perl using
[Homebrew](https://brew.sh). Similarly, Linux (and evidently Microsoft Windows Subsystem
for Linux) users can use [Linuxbrew](http://linuxbrew.sh). These typically install the
latest production release with a single command.

To install a Perl in your home directory (or other location) with a simple, but powerful,
tool, use the excellent [PerlBrew](https://perlbrew.pl). This tool can painlessly compile,
install, and manage one or more Perl release versions side-by-side, allowing you to easily
switch between releases with a simple command. It also manages multiple `local::lib`
installations, in case you want to isolate packages. 

BioToolBox does not utilize threading (it uses forks for parallel execution), so if you 
have a choice, compile a non-threaded Perl for a slight performance gain. 

### System installation

For privileged installations (requiring `root` access or `sudo` privilege) you probably
already know what to do. You can use the `--sudo` or `-S` option to `cpanm`. Note that
installing lots of packages in the OS vendor system perl is generally not recommended, as
it could interfere with other vital OS functions or programs that expect certain versions
or modules to be present. It's best to use one of the other two methods.

## External libraries

There are two external C libraries that are required for reading Bam and BigWig files. 
These are commonly used bioinformatics tools maintained by separate organizations, and 
the Perl modules only provide the XS bindings to these libraries. As such, it's best 
to install these up front separately before attempting the Perl module installation. 
Note that both Perl modules [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS) and 
[Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big) include `INSTALL.pl` scripts within 
their bundles that can compile these external libraries for you in a semi-automated 
fashion. Proceed here if you wish to have more control over what and where these are 
installed.

- [HTSlib](https://github.com/samtools/htslib)

    Follow the directions within for installation. Version 1.9 is officially
    recommended by the `Bio::DB::HTS` authors; however, later versions appear to work
    just fine and should probably be preferred. 
    [Version 1.19](https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2)
    have been used successfully by this author. By default, it installs into
    `/usr/local`, or it may be set to another directory (`$HOME` for example) by
    adding `--prefix=$HOME` option to the `configure` step. This may also be
    available via OS or other package managers. 
    
    Note that the `htslib` package is also included with the `samtools` and `bcftools`
    packages and is compiled therein, but is not installed by default (although it
    could be). In this case, use the commands `make all all-htslib` and
    `make install all all-htslib` to install both.
	
	See above for example code.

- [libBigWig](https://github.com/dpryan79/libBigWig)

    Follow the directions within for installation. By default, it installs into 
    `/usr/local`. To change to a different location, manually edit the `Makefile`
    to change `prefix` to your desired location, and run `make && make install`.
    
    See above for example code.
    
    Note that the Perl module
    [Alien::LibBigWig](https://metacpan.org/pod/Alien::LibBigWig) is available
    to install this library dependency for you automatically for `Bio::DB::Big`.
    although by default it installs an older version and brings along a lot of
    extra `Alien` Perl modules, if that concerns you.
    

## Perl modules

Using a simple CPAN package installer such as [CPAN
Minus](https://metacpan.org/pod/App::cpanminus), i.e. `cpanm`, is highly recommended
for ease and simplicity in installing modules from [CPAN](https://metacpan.org). It
can install directly from CPAN or take a URL or downloaded archive file. Follow the
link to find out how to install `cpanm`. Other CPAN package managers are available
too, if that's your preference. 

The following Perl packages should be explicitly installed. Most of these will 
bring along a number of dependencies (which in turn bring along more dependencies). In 
the end you will have installed dozens of packages. 

- [BioPerl](https://metacpan.org/pod/BioPerl)

    The BioPerl package is a large bundle consisting of hundreds of modules and 
    bundled scripts, the vast majority of which is not needed by Bio::ToolBox.
    However, it is required by Bio::DB::HTS and all of the legacy adapters (see below). 
    
	**NOTE:** For users with no other explicit need for BioPerl, an unofficial, custom, 
	[minimal version](https://github.com/tjparnell/bioperl-live/tree/minimal-tjparnell) 
	has been generated, which is considerably faster to install and has a footprint 
	one-third the size of the full distribution; use this 
	[file](https://github.com/tjparnell/bioperl-live/releases/download/minimal-v1.7.8/Minimal-BioPerl-1.7.8.tar.gz)
	to install.

- [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS)

    This provides a perl interface to the [HTSlib](https://github.com/samtools/htslib) 
    library for working with Bam files. It should be able to identify HTSLIB in standard 
    library locations, such as `/usr/local` for example, on its own. For non-standard 
    locations, specify the location of the HTSlib path to `Build.PL` using the 
    `--htslib` option. 
    
- [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big)

    This provides a perl interface to the UCSC-style bigWig and bigBed formats. 
    As with HTSlib, this should be able to identify the libBigWig library in standard 
    locations, but with non-standard locations, you may specify the path with the 
    `--libbigwig` option to `Build.PL`. 
    
    **NOTE**: The distribution from CPAN will install dozens of unnecessary modules
    for remote URL testing. You may be better off installing directly from
    [source](https://codeload.github.com/Ensembl/Bio-DB-Big/tarball/master).
    
    **NOTE**: This requires Template::Tiny to be installed, which appears to be
    missing from the manifest. Install it first before Bio::DB::Big.
    
- [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

    This is highly recommended to get multi-cpu support for some of the data collection 
    scripts, which can otherwise get slow with a single thread.

- [Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)

    This is necessary for optional functionality (quick intersection of genomic 
    intervals, such as black lists) for a few scripts, namely
    [bam2wig.pl](apps/bam2wig.md).

- [Set::IntSpan::Fast](https://metacpan.org/pod/Set::IntSpan::Fast)

    This is necessary for optional functionality in the
    [bam2wig.pl](apps/bam2wig.md) script. For a slight speed boost, install the
    [XS version](https://metacpan.org/pod/Set::IntSpan::Fast::XS) if possible.
    Note that the XS version requires
    [Data::Swap](https://metacpan.org/pod/Data::Swap); however, this appears to
    [fail to install](http://matrix.cpantesters.org/?dist=Data-Swap+0.08) on
    virtually all recent Perl versions, so unless you're living in the past,
    we're basically stuck with the Perl-only version until someone adopts and
    fixes it.

An example of installing these Perl modules with `cpanm` is below. This assumes that 
you have `local::lib` or a writable Perl installation in your `$PATH`. Adjust accordingly.

    curl -O -L https://github.com/tjparnell/bioperl-live/releases/download/minimal-v1.7.8/Minimal-BioPerl-1.7.8.tar.gz
    cpanm Minimal-BioPerl-1.7.8.tar.gz
    cpanm --configure-args="--htslib /usr/local" Bio::DB::HTS
	cpanm Template::Tiny
	curl -o bio-db-big-master.tar.gz -L https://github.com/Ensembl/Bio-DB-Big/tarball/master
	cpanm --configure-args="--libbigwig /usr/local" bio-db-big-master.tar.gz
    cpanm Parallel::ForkManager Set::IntervalTree Set::IntSpan::Fast Bio::ToolBox

## External applications

Some programs, for example [bam2wig.pl](apps/bam2wig.md) and
[data2wig](apps/data2wig.md), requires external utilities for converting text
formats to binary formats, for example wig files to bigWig. External utilities
are preferred because they're more efficient and spread the load on modern
multi-CPU environments. You may download these from the UCSC Genome Browser
utilities section for either
[Linux](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) or
[macOS](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/). Copy them to
your `bin` directory in your `PATH`, for example `$HOME/bin`, `$HOME/perl5/bin`,
or `/usr/local/bin`. Be sure to make them executable by running `chmod +x` on
each file.

- wigToBigWig
- bedGraphToBigWig
- bigWigToWig
- bedToBigBed

An example for downloading on Linux:

    for name in wigToBigWig bedGraphToBigWig bigWigToWig bigWigToBedGraph bedToBigBed bigBedToBed; \
    do curl -o $HOME/bin/$name http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$name \
    && chmod +x $HOME/bin/$name; done;

**NOTE** Current versions of these utilities do not support directly piping data into
the utility using the `stdin` file name. You will need to either find an older binary
version, compile your own from older source code (see below), or update BioToolBox;
version 2.02 now supports a work-around.

## Legacy Perl modules

These are additional legacy Perl modules that are supported (for example, if you still 
have a [GBrowse](http://gmod.org/wiki/GBrowse) installation), but are either not required 
or have been superseded by other modules. 

- [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store)
- [Bio::DB::Sam](https://metacpan.org/pod/Bio::DB::Sam)
- [Bio::DB::BigWig](https://metacpan.org/pod/Bio::DB::BigWig)
- [Bio::DB::BigBed](https://metacpan.org/pod/Bio::DB::BigBed)
- [Bio::DB::USeq](https://metacpan.org/pod/Bio::DB::USeq)
- [Bio::Graphics::Wiggle](https://metacpan.org/pod/Bio::Graphics::Wiggle)

Some notes are below for anyone who may need to install these. 

### Database support

The Bio::DB::SeqFeature::Store is a convenient SeqFeature annotation database
backed by a SQL engine. It used to be part of the BioPerl distribution prior to
release 1.7.3, but is now split into its own distribution. If you wish to use
annotation databases, you will need a SQL driver, such as
[DBD::SQLite](https://metacpan.org/pod/DBD::SQLite) (recommended for
individuals) or [DBD::mysql](https://metacpan.org/pod/DBD::mysql) (for fancy
multi-user installations). 

### Sam library

The Bio::DB::Sam library _only_ works with the legacy Samtools version, which
included both the C libraries, headers, and executables; use version
[0.1.19](https://github.com/samtools/samtools/archive/0.1.19.tar.gz) for best
results. You will need to compile the Samtools code, but you do not have to install
it (the library is not linked). Before compiling, edit the Makefile to include the
cflags `-fPIC` and (most likely) `-m64` for 64 bit OS. Export the `SAMTOOLS`
environment variable to the path of the Samtools build directory, and then you can
proceed to build the Perl module; it should find the necessary files using the
`SAMTOOLS` environment variable. You may obtain the latest source from
[Github](https://github.com/GMOD/GBrowse-Adaptors/tree/master/Bio-SamTools) or by
downloading a [tarball](https://github.com/GMOD/GBrowse-Adaptors/tarball/master).
**Note** that this project and file contains multiple Perl adapters and cannot be
used directly with `cpanm`, for example.

### UCSC BigFile library

The Bio::DB::BigWig and Bio::DB::BigBed modules are part of the same distribution, 
[Bio-BigFile](https://github.com/GMOD/GBrowse-Adaptors/tree/master/Bio-BigFile). Only 
use the code from the GitHub repository, as it should be compatible with recent UCSC 
libraries, whereas the distribution on CPAN is out of date. 

**NOTE** The UCSC library, when it encounters an error, will immediately terminate
the Perl process, with no chance of trapping the error. The newer `libBigWig` C
library used with Bio::DB::Big (detailed above) does not exhibit this behavior, plus
it's considerably easier to install. Encountering errors rarely happens, however,
because all bioinformatic data is always perfectly formatted and well behaved, right?

You will need the UCSC source code; the
[userApps](http://hgdownload.soe.ucsc.edu/admin/exe/) source code is sufficient, rather
than the entire browser code. Versions 375 and 398, at the time of this writing, works
successfully, but more recent versions appear to have increasing problems with
successful compilation – YMMV. 

**NOTE** If you are compiling the command line utilities, such as `wigToBigWig`, be
aware that in version 439 and later, these utilities no longer accept `stdin` as a
file input. The Bio::ToolBox::big_helper module uses this feature for convenience in
applications such [bam2wig](apps/bam2wig.md). You can compile your own following
these steps, but you do not need to install Bio::DB::BigWig.

This requires at least `OpenSSL` and `libpng` libraries to compile the required
library. For the command line utilities, if desired, you will also need MySQL
libraries; MariaDB, for now, seems adequate as far as I can tell.

On Linux, this is mostly not a problem as these libraries and development files are
readily available through the package manager. If you're building on macOS, see the
notes in the [macOS notes page](MacOSNotes.md).

For purposes of installing the Perl adapter, only the library needs to be compiled.
It does not need to be installed, as nothing is linked. So, you do not need to run
the full `make` command. In the `userApps` folder, run

	cd path/to/userApps
	make installEnvironment

This will generate `kent/src/inc/localEnvironment.mk` for your local machine. Edit this
file to add at the end

	CFLAGS = -fPIC

If you have odd or non-standard locations for some libraries, for example in a
computing cluster where the development files are brought in using environment
[modules](https://modules.readthedocs.io/), you may be able to set additional paths
in this `localEnvironment.mk` file or by directly hacking `kent/src/inc/common.mk`.
**NOTE** Be careful about setting a generic path to the libraries, particularly if
you also have `htslib` installed, since the UCSC userApp provides its own (presumably
modified) `htslib` library, which will conflict with a system available library.

To proceed with library compilation, follow the following steps. Note exporting
environment variables which aid in building the Perl adapter.

	cd path/to/kent/src
	export KENT_SRC=$(PWD)
	cd htslib
	make
	cd ../lib
	make
	cd

To install the Perl adapter, download the tarball from Github (same source as
Bio::DB::Sam above) and follow the steps below.

	curl -o GBrowse-Adaptors.tar.gz -L https://github.com/GMOD/GBrowse-Adaptors/tarball/master
	tar -xf GBrowse-Adaptors.tar.gz
	cd GBrowse-Adaptors-master-85c29de/Bio-BigFile
	export MACHTYPE=local
	perl Build.PL
	./Build
	./Build test
	./Build install

If the environment variables have been set correctly and the library compiled
successfully, then the Perl build process should proceed smoothly. There may be
various warnings emitted during the build process, which can usually be ignored.

Once you have compiled the main `kent/src/lib` library, you can proceed with
compiling the command line utilities, if you desire. You can build just the ones you
want by going into each utility subdirectory within `kent/src/utils/` and issuing
`make`. For example:

	cd /path/to/userApps
	mkdir bin
	cd kent/src/utils/wigToBigWig
	make

The compiled executable should be copied into `userApps/bin`, and you can move it from
there to wherever.



