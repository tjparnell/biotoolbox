# MacOS Notes

While Macs have a Unix-compatible command-line environment (Darwin), there are a few 
issues and solutions that I have encountered that may be useful to someone. Some of 
these pertain to older OS X releases (and left here for posterity), but some, such 
as [rpath](#rpath_errors), are still relevant to current macOS releases.

## Install XCode command line tools

You don't need the full blown XCode installed, just the command line tools. 
Running the following command in terminal will prompt to install them. You will 
need administrator permissions.

	$ xcode-select --install

## Installing your own Perl

Apple is generally a little slow in updating their Perl compared to the latest available
versions, and it is compiled for multiple architectures with backwards compatibility with
32-bit `i386` processors (!), at least as of High Sierra (10.13). Some of the errors below
will go away if you compile your own Perl, but your success may vary.

On Mojave (10.14), the system Perl appears to be ~~broken~~ fixed as of 10.14.3 with 
regards to compiling XS modules. 

For installing your own Perl, the easiest route is to use either [Homebrew](https::/brew.sh) 
or [PerlBrew](https://perlbrew.pl). 

Be sure to use a recent version of PerlBrew. Older versions sometimes fails to set the 
minimum MacOS version. In such cases, I have found the recommendations described 
[here](https://karl.kornel.us/2015/12/perl-osx-1011-warnings/) to be helpful. 

## Linking errors

When linking Perl modules with XS code (compiled C extensions), especially when using 
the system Perl, you may see the following errors.

	ld: warning: object file was built for newer OSX version (10.13) than being linked (10.4)

This is due to Apple compiling their system Perl with far-reaching backwards 
compatibility; the Perl binary was compiled for both `i386` and `x86_64`, but in 
all likelihood your XS was compiled only for `x86_64`. In some cases, this is a 
harmless error; in other cases, it's a deal breaker. The best solution is to 
install your own Perl.

## rpath errors

This is especially notable with the [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big) 
installation as described above, where the `./Build test` fails dramatically because a 
shared library can not be loaded by the bundle, usually with an error message including 
this:

	Reason: unsafe use of relative rpath libBigWig.so in blib/arch/auto/Bio/DB/Big/Big.bundle with restricted binary

The solution is to manually re-link the bundle to the shared library file with the 
following command. See this 
[link](https://stackoverflow.com/questions/33275605/el-capitan-perl-dbd-unsafe-use-of-relative-path) 
for the source of the  solution.

	$ install_name_tool -change libBigWig.so /path/to/lib/libBigWig.so blib/arch/auto/Bio/DB/Big/Big.bundle

## DB_File errors

There are [reports of issues](https://github.com/bioperl/bioperl-live/issues/267) 
regarding certain BioPerl modules that rely on the Berkley database module 
[DB_File](https://metacpan.org/pod/DB_File). This appears to stem from an issue with 
the Apple-supplied library in High Sierra (10.13) as described 
[here](https://discussions.apple.com/thread/8125401). The best solution is to 
install your own `berkley-db` library. 

This issue is not evident on the newer Mojave (10.14) and subsequent releases.

For BioToolBox users, the biggest effect appears to be exceptionally long times 
during `Build` tests, specifically file `04.DB.t` that uses the in-memory database 
adapater (maybe 20-30 seconds instead of 1), and excruciatingly long 
[Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) 
database builds (possibly days or weeks, I give up). 


## Set::IntervalTree failures

Installing [Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree) may lead
failures. This isn't necessarily an issue with per se, but rather one of its
dependencies, [ExtUtils::CppGuess](https://metacpan.org/pod/ExtUtils::CppGuess),
which fails to install on recent versions of macOS. In fact, most tests of it on
Darwin [fail](http://matrix.cpantesters.org/?dist=ExtUtils-CppGuess+0.27).
Fortunately, the failing tests don't appear to be essential for installing
Set::IntervalTree, so force install ExtUtils::CppGuess and try re-installing
Set::IntervalTree again – it will probably work.


## libBigWig

When manually installing libBigWig on recent versions of macOS (observed with Sonoma,
14.x and libBigWig v0.4.7), the compilation may fail at first. To check for libCurl
dependencies, it attempts to compile a small test program and runs the command
`mktemp --suffix=.c`. While that `--suffix` option is available to versions on Linux
platforms, it is not available to the version on macOS, thus breaking the detection
of libCurl. To work around this, we just have to tell it that, yes, we have libCurl.
Comment out the five lines after `# Create a simple test-program...` and add a new line

	HAVE_CURL=YES

Then re-run `make` and it should compile ok with curl support.

However, it's not over yet. If you run `make test`, the `test/testRemote` test seems to
fail. This appears to be an innocuous platform-specific error. If you proceed to
compile the Bio::DB::Big Perl module, it appears to work ok with remote files during
its testing, albeit through a fake remote test with Test::Fake::HTTPD (if it's
installed). However, empirical testing with real remote data (via https) seems to work ok.






