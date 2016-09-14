#!/usr/bin/env perl

# Copyright 2009-2016 Christopher Benner <cbenner@ucsd.edu>
# 
# This file is part of HOMER
#
# HOMER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HOMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

use warnings;
use POSIX;

#use Cwd;

$0 =~ /^(.*)\//;
my $homeDir = $1;
if (!defined($homeDir) || $homeDir eq '') {
	$homeDir = ".";
}
$homeDir .= "/";
my $pwd = '';
unless ($homeDir =~ /^\//) {
	`pwd > .ls`;
	open IN, ".ls";
	while (<IN>) {
		chomp;
		$pwd = $_;
		last;
	}
	close IN;
	$homeDir = $pwd . "/" . $homeDir;
}


print STDERR "\n\tCurrent base directory for HOMER is $homeDir\n\n";


my $baseURL = "http://homer.salk.edu/homer/";
my $makeProgram = 'make';  # gmake for SunOS
my $tarProgram = 'tar';  # gtar for SunOS
my %realCats = ();
$realCats{'ORGANISMS'} = 'o';
$realCats{'PROMOTERS'} = 'p';
$realCats{'GENOMES'} = 'g';
$realCats{'SOFTWARE'} = 's';
my @groups = ('SOFTWARE','ORGANISMS','PROMOTERS','GENOMES');

my %duplicate = ();
my %duplicateLocal = ();
my %dupName = ();

my $listFlag = 0;
my %install = ();
my $installFlag = 0;
my $compileFlag = 0;
my $localFlag = 0;
my %remove = ();
my $removeFlag = 0;
my $checkFlag = 0;
my $updateFlag = 0;
my $reinstallFlag = 0;
my $allFlag = 0;
my $updateScriptFlag = 1;
my $makeFlag = 0;
my $getFactsFlag=0;
my %newSettings = ();

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-list') {
		$listFlag = 1;
	} elsif ($ARGV[$i] eq '-check') {
		$checkFlag = 1;
	} elsif ($ARGV[$i] eq '-make') {
		$makeFlag = 1;
		$localFlag = 1;
	} elsif ($ARGV[$i] eq '-local') {
		$localFlag = 1;
	} elsif ($ARGV[$i] eq '-reinstall') {
		$reinstallFlag = 1;
	} elsif ($ARGV[$i] eq '-all') {
		$allFlag = 1;
	} elsif ($ARGV[$i] eq '-getFacts') {
		$getFactsFlag = 1;
	} elsif ($ARGV[$i] eq '-keepScript') {
		$updateScriptFlag = 0;
	} elsif ($ARGV[$i] eq '-bigWigUrl') {
		$newSettings{'bigWigUrl'} = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bigWigDir') {
		$newSettings{'bigWigDir'} = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-hubsUrl') {
		$newSettings{'hubsUrl'} = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-hubsDir') {
		$newSettings{'hubsDir'} = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-url') {
		$baseURL = $ARGV[++$i];
		print STDERR "\tWill download Homer files from  $baseURL\n";
	} elsif ($ARGV[$i] eq '-sun') {
		$makeProgram = "gmake";
		$tarProgram = "gtar";
		print STDERR "\tUsing gmake and gtar\n";
	} elsif ($ARGV[$i] eq '-update') {
		$updateFlag = 1;
	} elsif ($ARGV[$i] eq '-install') {
		$installFlag = 1;
		my $bail = 0;
		last if ($i == @ARGV-1);
		while ($ARGV[++$i] !~ /^\-/) {
			$install{$ARGV[$i]} = 1;
			print STDERR "\tWill install $ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail=1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-remove') {
		$removeFlag = 1;
		my $bail = 0;
		last if ($i == @ARGV-1);
		while ($ARGV[++$i] !~ /^\-/) {
			$remove{$ARGV[$i]} = 1;
			print STDERR "\tWill remove $ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail=1;
				last;
			}
		}
		last if ($bail == 1);
	} else {
		printCMD();
	}
}

if (@ARGV < 1) {
	printCMD();
}
if ($checkFlag) {
	checkForPrograms();
}
if ($localFlag) {
	compileSoftware();
	exit;
}
if ($getFactsFlag == 1) {
	print STDERR "\tInstalling: Chuck Facts\n";
	print STDERR "\t\tDownloading...\n";
	my $url = $baseURL . "homer.misc.v1.0.zip";
	`wget -O homer.package.zip "$url"`;
	print STDERR "\t\tUnzipping...\n";
	`unzip -o -d "$homeDir" homer.package.zip`;
	`rm homer.package.zip`;
}
my $configFile = $homeDir . "/" . "config.txt";

if ($listFlag || $updateFlag || $removeFlag || $installFlag) { 

	my $updateURL = $baseURL . "update.txt";
	my $updateFile = $homeDir . "/" . "update.txt";
	`wget -O "$updateFile" "$updateURL"`;
	my $config = loadConfigFile($configFile);
	updateSettings($config,\%newSettings);
	my $update = loadConfigFile($updateFile);
	if ($listFlag) {
		printConfigFile($config,"",$update);
	}

	my %installList = ();
	my %removeList = ();
	foreach(@groups) {
		my %a = ();
		my %b = ();
		$installList{$_} = \%a;
		$removeList{$_} = \%a;
	}
	if (scalar(keys %install) < 1 && $installFlag) {
		$installList{'SOFTWARE'}->{'homer'}=1;
	}
	foreach(keys %install) {
		my $pkg = $_;
		if (exists($duplicate{$pkg})) {
			print STDERR "!!! Warning: multiple packages found for \"$pkg\" !!!\n";
			print STDERR "!!! Please select rerun with one of the following:\n";
			foreach(values %{$duplicate{$pkg}}) {
				print STDERR "!!!   $_\n";
			}
			exit;
		}
		if (exists($dupName{$pkg})) {
			my $p = $dupName{$pkg}->{'pkg'};
			my $g = $dupName{$pkg}->{'g'};
			$installList{$g}->{$p}=1;
		} elsif (exists($update->{'packages'}->{$pkg})) {
			my $g = $update->{'packages'}->{$pkg};
			$installList{$g}->{$pkg}=1;
		} else {
			print STDERR "!!! Error: could not find package: \"$pkg\" !!!\n";
			#my @a = keys %{$update->{'packages'}};
			#print STDERR "@a\n";
			exit;
		}
	}
	foreach(keys %remove) {
		my $pkg = $_;
		if (exists($duplicate{$pkg})) {
			print STDERR "!!! Warning: multiple packages found for \"$pkg\" !!!\n";
			print STDERR "!!! Please select rerun with one of the following:\n";
			foreach(values %{$duplicate{$pkg}}) {
				print STDERR "!!!   $_\n";
			}
			exit;
		}
		if (exists($dupName{$pkg})) {
			my $p = $dupName{$pkg}->{'pkg'};
			my $g = $dupName{$pkg}->{'g'};
			$removeList{$g}->{$p}=1;
		} elsif (exists($config->{'packages'}->{$pkg})) {
			my $g = $config->{'packages'}->{$pkg};
			$removeList{$g}->{$pkg}=1;
		} else {
			print STDERR "!!! Error: could not find package: \"$pkg\" !!!\n";
			#my @a = keys %{$update->{'packages'}};
			#print STDERR "@a\n";
			exit;
		}
	}
	if ($removeFlag) {
		my @paks = keys %remove;
		if (@paks < 1) {
			print STDERR "\tNothing needs to be removed\n";
		} else {
			print STDERR "\tThe following packages will be removed:\n";
			foreach(@groups) {
				my $g = $_;
				print STDERR "\t\t$g\n";
				foreach(keys %{$removeList{$g}}) {
					my $pkg = $_;
					my $loc = $homeDir . $config->{$g}->{$pkg}->{'location'};
					my $cmd = "";
					if ($g eq 'ORGANISMS') {
						$cmd = "rm \"$loc/$_\"* \"$loc/../GO/$_.\"*";
					} elsif ($g eq 'PROMOTERS') {
						$cmd = "rm \"$loc/$_.\"*\n";
					} elsif ($g eq 'GENOMES') {
						$cmd = "rm -r \"$loc\"\n";
					} else {
						print STDERR "\t!! $g packages must be removed manually\n";
					}
					print STDERR "\t\t\t$_\tcmd = $cmd\n";
				}
			}
			print STDERR "\n";
			print STDERR "\tWaiting 10s to make sure this is what you want to do (hit ctrl+C to cancel)\n\t\t";
			for (my $i=10;$i>0;$i--) {
				print STDERR "$i ";
				`sleep 1`;
			}
			print STDERR "\n";
			foreach(@groups) {
				my $g = $_;
				foreach(keys %{$removeList{$g}}) {
					my $pkg = $_;
					my $loc = $homeDir . $config->{$g}->{$pkg}->{'location'};
					my $cmd = "";
					if ($g eq 'ORGANISMS') {
						$cmd = "rm \"$loc/$_\"* \"$loc/../GO/$_.\"*";
					} elsif ($g eq 'PROMOTERS') {
						$cmd = "rm \"$loc/$_.\"*\n";
					} elsif ($g eq 'GENOMES') {
						$cmd = "rm -r \"$loc\"\n";
					} else {
						print STDERR "\t!! $g packages must be removed manually\n";
					}

					print STDERR "\tRemoving $pkg ($g)\n\t\tcmd = $cmd\n";
					`$cmd`;
					delete $config->{$g}->{$pkg}
				}
			}
		}
		printConfigFile($config, $homeDir . "/config.txt");
	}


	if ($updateFlag == 1 || $reinstallFlag == 1) {
		foreach(keys %$config) {
			my $cat = $_;
			next if (!exists($realCats{$cat}));
			foreach(keys %{$config->{$cat}}) {
				my $package = $_;
				my $curVer = $config->{$cat}->{$package}->{'version'};
				my $latest = "";
				if (exists($update->{$cat})) {
					if (exists($update->{$cat}->{$package})) {
						$latest = $update->{$cat}->{$package}->{'version'};
					}
				}
				next if ($latest eq '');
				if ($curVer =~ /^v/) {
				} else {
					print STDERR "\tNote: Package $package has non-standard version name, skipping update...\n";
					next;
				}
				if ($reinstallFlag == 1) {
					$install{$package} = 1;
					$installList{$cat}->{$package} = $latest;
				}
				if ($latest ne $curVer) {
					$install{$package} = 1;
					$installList{$cat}->{$package} = $latest;
				}
			}
		}

		my @paks = keys %install;
		if (@paks < 1) {
			print STDERR "\tNothing needs to be updated\n";
		} else {
			print STDERR "\tThe following packages will be updated:\n";
			foreach(@groups) {
				my $g = $_;
				print STDERR "\t\t$g\n";
				foreach(keys %{$installList{$g}}) {
					print STDERR "\t\t\t$_\t$installList{$g}->{$_}\n";
					$installFlag=1;
				}
			}
		}
		print STDERR "\n";
	}
	if ($allFlag) {
		print STDERR "\tIf updating ALL, the following needs to be updated:\n";
		foreach(@groups) {
			my $g = $_;
			print STDERR "\t\t$g\n";
			foreach(keys %{$update->{$g}}) {
				my $pkg = $_;
				my $v = $update->{$g}->{$pkg}->{'version'};
				$installList{$g}->{$pkg} = $v;
				print STDERR "\t\t\t$pkg\t$v\n";
				$install{$_}=1;
				$installFlag = 1;
			}
		}
		print STDERR "\n";
	}

	if ($installFlag) {
		updateConfigureScript($baseURL . "configureHomer.pl", $0 , $updateScriptFlag, @ARGV);

		my @depCheck = ("PROMOTERS","GENOMES");
		foreach(@depCheck) {
			my $g = $_;
			foreach(keys %{$installList{$g}}) {
				my $pkg = $_;
				my $org = $update->{$g}->{$pkg}->{'params'}->[0];
				$installList{'ORGANISMS'}->{$org} = 1;
				print STDERR "\t\tInstalling/Updating Organism $org for package $pkg\n";
			}
		}

		if (!exists($config->{'SOFTWARE'}->{'homer'})) {
			$install{'homer'} = 1;
		}
		print STDERR "\tPackages to Install...\n";
		foreach(@groups) {
			my $g = $_;
			foreach(keys %{$installList{$g}}) {
				my $package = $_;
				next if ($package eq '');
				print STDERR "\t\t$_ -> ";
				my $version = "not available for installation";
				if (exists($update->{$g}->{$package})) {
					$version = $update->{$g}->{$package}->{'version'};
				} else {
					#print STDERR "!!! Could not find $package ($g)\n";
					delete $installList{$g}->{$package};
					next;
				}
				#print STDERR "Checking $package ($g) $version\n";
				print STDERR "$version\n";
				if ($package eq 'homer') {
					checkForPrograms();
				}
			}
		}
		foreach(@groups) {
			my $mode = $_;
			foreach(keys %{$installList{$mode}}) {
				my $package = $_;
				my $url = $update->{$mode}->{$package}->{'url'};
				print STDERR "\tInstalling: $package\n";
				print STDERR "\t\tDownloading...\n";
				`wget -O homer.package.zip "$url"`;
				print STDERR "\t\tUnzipping...\n";
				`unzip -o -d "$homeDir" homer.package.zip`;
				`rm homer.package.zip`;
				$config->{$mode}->{$package} = $update->{$mode}->{$package};
				print STDERR "\t\tFinished Installing $package\n\n";
			}
		}
		if (exists($install{'homer'}) || $installList{'SOFTWARE'}->{'homer'}) {
			compileSoftware();
		}
		printConfigFile($config, $homeDir . "/config.txt");
	}
} elsif (scalar(keys %newSettings) > 0) {

	my $config = loadConfigFile($configFile);
	updateSettings($config,\%newSettings);
	printConfigFile($config, $homeDir . "/config.txt");

}
exit;


sub printCMD {
	print STDERR "\n\tUsage: configureHomer.pl [options]\n";
	print STDERR "\tThis program will install HOMER in the directory containing configureHomer.pl\n";
	print STDERR "\t\ti.e. save this in a directory named homer/ or something like that\n";
	print STDERR "\tTo install the program from scratch, run the following command:\n";
	print STDERR "\t\tperl path-to-homer/configureHomer.pl -install\n";
	print STDERR "\tIf upgrading, make sure HOMER is not running (may prevent replacement of executables)\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-list (prints available and installed packages to screen)\n";
	print STDERR "\t\t-install (install basic Homer software from scratch)\n";
	print STDERR "\t\t-install <package name1> [package name2] ... (installs packages)\n";
	print STDERR "\t\t-remove <package name1> [package name2] ... (safely removes packages)\n";
	print STDERR "\t\t-update (updates all packages to the newest version)\n";
	print STDERR "\t\t-reinstall (forces the reinstallation of all currently installed packages)\n";
	print STDERR "\t\t-all (install everything - will take a while...)\n";
	print STDERR "\t\t-getFacts (add humor to HOMER - to remove delete contents of data/misc/)\n";
	print STDERR "\t\t-check (checks for required 3rd party software)\n";
	print STDERR "\t\t-make (reconfigure files, recompile exectables)\n";
	print STDERR "\t\t-sun (Use with SunOS - i.e. gmake and gtar instead of make and tar)\n";
	print STDERR "\t\t-keepScript (do NOT replace configureHomer.pl script if changes are detected with website)\n";
	print STDERR "\t\t-url <installation url> (For use if the Homer website changes)\n";
	print STDERR "\t\t\tdefault is: $baseURL\n";
	print STDERR "\t\tHubs & BigWig settings (with read existing settings from config.txt if upgrading):\n";
	print STDERR "\t\t\t-bigWigUrl <base urls for bigWigs> (Setting for makeBigWigs.pl)\n";
	print STDERR "\t\t\t-bigWigDir <base directory for bigWigs> (Setting for makeBigWigs.pl)\n";
	print STDERR "\t\t\t-hubsUrl <base urls for hubs> (Setting for makeMultiWigHub.pl)\n";
	print STDERR "\t\t\t-hubsDir <base directory for hubs> (Setting for makeMultiWigHub.pl)\n";
	print STDERR "\tConfiguration files: downloads update.txt from the homer website,\n";
	print STDERR "\t\tand corrects config.txt with the appropriate options\n";
	print STDERR "\n";
	exit;
}

sub updateConfigureScript {

	my $scriptURL = $_[0];
	my $scriptFile = $_[1];
	my $updateScriptFlag = $_[2];
	my $str = "";
	for (my $i=3;$i<@_;$i++) {
		$str .= " " . $_[$i];
	}

	my $tmpFile = rand() . ".tmp";
	print STDERR "`wget -O $tmpFile $scriptURL`;\n";
	`wget -O $tmpFile "$scriptURL"`;
	$delta = `diff $tmpFile $scriptFile | wc -l`;
	if ($delta > 0 && $updateScriptFlag == 1) {
		print STDERR "\tconfigureHomer.pl script differs from one located at:\n\t\t\t$scriptURL\n";
		print STDERR "\tHOMER will now replace the old one - to cancel hit control-C NOW!!!\n";
		print STDERR "\t\t(if you don't want to do this, or if it keeps updating for no reason, use -keepScript)\n";
		for (my $i=10;$i>0;$i--) {
			print STDERR "\t\t$i\n";
			`sleep 1`;
		}
		`mv $tmpFile "$scriptFile"`;
		print STDERR "`perl $scriptFile $str`\n";
		`perl "$scriptFile" $str`;
		exit;
	} else {
		print STDERR "\tconfigureHomer.pl script is up-to-date\n";
	}
	`rm $tmpFile`;
}

sub printConfigFile {
	my ($config, $outputFile,$update) = @_;

	my $file = '';
	if (!defined($outputFile) || $outputFile eq 'stdout' || $outputFile eq '') {
		$file = \*STDOUT;
	} else {
		open CONFIGOUT, ">$outputFile";
		$file = \*CONFIGOUT;
	}
	my %packages = ();
	foreach(@groups) {
		my @packs = (keys %{$config->{$_}});
		$packages{$_} = \@packs;
	}

	if (defined($update)) {
		print $file "Packages with name conflicts have a trailing -o, -p, or -g\n";
		print $file "Version Installed\tPackage\tVersion\tDescription\n";
		foreach(@groups) {
			push(@{$packages{$_}}, keys %{$update->{$_}});
		}
	} else {
		print $file "# Homer Configuration File (automatically generated)\n";
		print $file "#\n";
		print $file "# This file is updated from the Homer website and contains information about data available for\n";
		print $file "# use with the program.\n";
		print $file "#\n";
		print $file "# Each section has the same format, which is <tab> separated values specifying:\n";
		print $file "# package name <tab> version <tab> description <tab> url <tab> optional parameters (, separated)\n";
		print $file "#\n";
	}

	my %done = ();
	foreach(@groups) {
		my $group = $_;	
		print $file "$group\n";
		foreach(@{$packages{$group}}) {
			my $package = $_;
			next if (exists($done{$package . "-" . $group}));
			$done{$package . "-" . $group} = 1;
			if (defined($update)) {
				my $printName = $package;
				if (exists($duplicate{$package})) {
					$printName = $duplicate{$package}->{$group};
				}
				if (exists($config->{$group}->{$package})) {
					if (exists($update->{$group}->{$package})) {
						my $curVer = $config->{$group}->{$package}->{'version'};
						my $upVer = $update->{$group}->{$package}->{'version'};
						if ($curVer ne $upVer) {
							print $file "$curVer";
						} else {
							print $file "+";
						}
					} else {
						print "custom";
					}
				} else {
					print $file "-";
				}
				print $file "\t$printName";
				if (exists($update->{$group}->{$package})) {
					print $file "\t$update->{$group}->{$package}->{'version'}";
					print $file "\t$update->{$group}->{$package}->{'description'}";
				} else {
					print $file "\tNA\tEither Your Custom Addition or No Longer Supported";
				}
				print $file "\n";
			} else {
				next if (!exists($config->{$group}->{$package}));
				print $file "$package";
				print $file "\t$config->{$group}->{$package}->{'version'}";
				print $file "\t$config->{$group}->{$package}->{'description'}";
				print $file "\t$config->{$group}->{$package}->{'url'}";
				print $file "\t$config->{$group}->{$package}->{'location'}";
				print $file "\t";
				my $z = 0;
				foreach(@{$config->{$group}->{$package}->{'params'}}) {
					print $file "," if ($z>0);	
					$z++;	
					print $file "$_";
				}
				print $file "\n";
	
			}
		}
	}
	print $file "SETTINGS\n";
	if (exists($config->{'SETTINGS'})) {
		foreach(keys %{$config->{'SETTINGS'}}) {
			my $var = $_;
			my $val = $config->{'SETTINGS'}->{$var};
			print $file "$var=$val\n";
		}
	}

	if (!defined($outputFile) || $outputFile eq 'stdout' || $outputFile eq '') {
	} else {
		close $file;
	}
}

sub loadConfigFile {
	my ($file) = @_;
	
	my %a = ();
	my %b = ();
	my %c = ();
	my %d = ();
	my %e = ();
	my %f = ();
	my %g = ();
	my $config = {status=>'not installed',
			GENOMES=>\%a,
			PROMOTERS=>\%b,
			SOFTWARE=>\%c,
			ORGANISMS=>\%d,
			packages=>\%e,
			SETTINGS=>\%f
			};
	
	open IN, $file or return $config;
	while (<IN>) {
		chomp;
		s/\r//g;
		s/#.*$//;
		s/^\s+//;
		next if ($_ eq '');

		if (/^SOFTWARE/) {
			$mode = 'SOFTWARE';
			next;
		}
		if (/^ORGANISMS/) {
			$mode = 'ORGANISMS';
			next;
		}
		if (/^PROMOTERS/) {
			$mode = 'PROMOTERS';
			next;
		}
		if (/^GENOMES/) {
			$mode = 'GENOMES';
			next;
		}
		if (/^SETTINGS/) {
			$mode = 'SETTINGS';
			next;
		}
		my @line = split /\t/;
		if (@line > 4) {
			next if ($mode eq '');
			$config->{'status'} = 'installed' if ($mode eq 'SOFTWARE');
			my $package = $line[0];
			my $version = $line[1];
			my $description = $line[2];
			my $url = $line[3];
			my $location = $line[4];
			my @params = ();
			if (@line > 5) {
				@params = split /\,/, $line[5];
			}
			my $p = {package=>$package,version=>$version,description=>$description,url=>$url,
							location=>$location,params=>\@params};
			$config->{$mode}->{$package} = $p;
			if (exists($config->{'packages'}->{$package}) && $config->{'packages'}->{$package} ne $mode) {
				if (!exists($duplicate{$package})) {
					my %a = ();
					$duplicate{$package} = \%a;
					my $oldMode = $config->{'packages'}->{$package};
					$duplicate{$package}->{$oldMode} = $package . "-$realCats{$oldMode}";
					$dupName{$package . "-" . $realCats{$oldMode}} = {g=>$oldMode,pkg=>$package};
				}
				$duplicate{$package}->{$mode} = $package . "-$realCats{$mode}";
				$dupName{$package . "-" . $realCats{$mode}} = {g=>$mode,pkg=>$package};
			}
			$config->{'packages'}->{$package}=$mode;
		} elsif ($mode eq 'SETTINGS') {
			my @val = split /\=/,$line[0];
			next if (@val < 2);
			$config->{$mode}->{$val[0]}=$val[1];
		}
	}
	close IN;
	return $config;
}

sub updateSettings {
	my ($config,$newSettings) = @_;
	print STDERR "\tUpdating Settings...\n";

	foreach(keys %$newSettings) {
		my $var = $_;
		my $val = $newSettings->{$var};
		$config->{'SETTINGS'}->{$var}=$val;
	}
}


sub checkForPrograms {
	my $bad = 0;
	my $warning = 0;
	my %bad = ();
	print STDERR "\n\tChecking for standard utilities and 3rd party software:\n\n";
	my @programs = ("wget","cut","gcc","g++","zip","unzip",$makeProgram,$tarProgram,
										"gunzip","gzip","gs","seqlogo","blat");
	foreach(@programs) {
		my $p = $_;
		print STDERR "\tChecking for $p... ";
		my $result = `which $p`;
		chomp $result;
		if ($result eq "") {
			if ($p eq 'gs' || $p eq 'seqlogo') {
				print STDERR "\t\tThe program $p was not found but is required for making motif logos\n";
				my $warning = 1;
			} elsif ($p eq 'blat') {
				print STDERR "\t\tThe program $p was not found but is required for removing redundant sequences during motif finding\n";
				$warning = 1;
			} else {
				print STDERR "\t\tThe program $p was not found but is strictly required\n";
				$bad{$p} =1 ;
				$bad =1 ;
			}
		} else {
			print STDERR "$result\n";
		}
	}
	print STDERR "\n";

	if ($warning == 1 || $bad == 1) {
		print STDERR "\tInstallation will halt since one or more required programs where not found.\n";
		print STDERR "\t\t(or they are not in your executable path)\n";
		print STDERR "\tMost of these programs are standard Unix Utilities that should be installed by default.\n";
		print STDERR "\tIf you are running Mac OS X, make sure to install developer tools\n";
		if (exists($bad{'wget'})) {
			print STDERR "\n\twget: http://www.gnu.org/software/wget/ or\n";
			print STDERR "\t\tHomer requires a working version of the program \"wget\" in order to automatically download\n";
			print STDERR "\t\tupdates and data.  wget was not detected on your system.  wget is fairly standard on most\n";
			print STDERR "\t\tlinux and unix systems, and should be easily installed through update managers such as \"yum\"\n";
			print STDERR "\t\tIf running Mac OS X, you can either download it from http://www.statusq.org/archives/2008/07/30/1954/\n";
	#		print STDERR "\t\tor create an alias for wget from the command \"curl -O\"\n";
			print STDERR "\n\tIf you prefer to install it from source go to http://www.gnu.org/software/wget/\n"; 
			print STDERR "\n\tAlso, make sure the wget binary is in your executable path\n";
			print STDERR "\n";
		} 
		if (exists($bad{'gs --version'})) {
			print STDERR "\n\tgs (GhostScript): http://pages.cs.wisc.edu/~ghost/\n";
		}
		if (exists($bad{'seqlogo'})) {
			print STDERR "\n\tseqlogo: http://weblogo.berkeley.edu/\n";
		} 
		if (exists($bad{'blat'})) {
			print STDERR "\n\tblat: http://genome-test.cse.ucsc.edu/~kent/exe/ (blatSuite.zip for your OS)\n";
		}


		print STDERR "\n\tAlso, if you think these programs have already been installed, make sure you modify your\n";
		print STDERR "\tconfiguration files to include each progam in the executable path (i.e. ~/.bash_profile)\n";
		print STDERR "\tRefer to the Homer documentation if you need more help: http://homer.salk.edu/homer/\n\n";
		if ($bad == 1) {
			print STDERR "\t--- Cannot Continue until required software is available ---\n";
			exit;
		} else {
			print STDERR "\t--- It is highly recommended that you install the missing software ---\n";
			print STDERR "\t--- Will continue with the installation in 10s ...\n";
			`sleep 10`;
		}
	} else {
		print STDERR "\tAll auxilary programs found.\n\n";
	}

}

sub compileSoftware {

	my @files = (
		"findMotifs.pl", 
		"findMotifsGenome.pl", 
		"compareMotifs.pl", 
		"findKnownMotifs.pl",
		"convertIDs.pl", 
		"findGOtxt.pl", 
		"findGO.pl", 
		"annotatePeaks.pl", 
		"makeBigWig.pl", 
		"makeMultiWigHub.pl",
		"annotateInteractions.pl", 
		"getConservedRegions.pl",
		"getDistalPeaks.pl",
		"annotateTranscripts.pl",
		"checkTagBias.pl", 
		"getGWASoverlap.pl",
		"loadPromoters.pl", 
		"loadGenome.pl", 
		"HomerConfig.pm", 
		"scanMotifGenomeWide.pl",
		"analyzeChIP-Seq.pl", 
		"preparseGenome.pl",
		"GenomeOntology.pl",
		"analyzeRNA.pl",
		"analyzeRepeats.pl",
		"removeOutOfBoundsReads.pl",
		"addGeneAnnotation.pl",
		"prepForR.pl",
		"assignTSStoGene.pl",
		"findHiCInteractionsByChr.pl",
		"getHiCcorrDiff.pl",
		"runHiCpca.pl",
		"findHiCDomains.pl",
		"SIMA.pl",
		"convertOrganismID.pl",
		"../update/updateUCSCGenomeAnnotations.pl",
		"../update/updateGeneIdentifiers.pl",
		"../update/updatePromoters.pl"

	);

	foreach(@files) {
		my $file = $_;
		open IN, "$homeDir/bin/" . $file;
		open OUT, ">.new.pl";
		my $count = 0;
		while (<IN>) {
			$count++;
			if ($count == 1) {
				print OUT "#!/usr/bin/env perl\n";
			} elsif ($count == 2) {
				print OUT "use warnings;\n";
			} elsif ($count == 3) {
				print OUT "use lib \"$homeDir/bin\";\n";
			} elsif ($count == 4) {
				print OUT 'my $homeDir = "' . $homeDir . "\";\n";
			} else {
				print OUT $_;
			}
		
		}
		close IN;
		close OUT;
		`mv .new.pl "$homeDir/bin/$file"`;
		`chmod 755 "$homeDir/bin/$file"`;
	}

	open IN, "$homeDir/cpp/SeqTag.cpp";
	open OUT, ">.new.cpp";
	print OUT "#include \"SeqTag.h\"\n";
	print OUT "//Hard code homer install directory...\n";
	print OUT "const char* HomerConfig::homeDirectory = \"$homeDir\";\n";
	print OUT "\n";
	my $count = 0;
	while (<IN>) {
		$count++;
		next if ($count < 5);
		print OUT $_;
	}
	close IN;
	close OUT;
	`mv .new.cpp "$homeDir/cpp/SeqTag.cpp"`;

	`$makeProgram -C "$homeDir/cpp/" clean`;
	`$makeProgram -C "$homeDir/cpp/"`;

	print STDERR "\n\tSoftware Installed.  If not done so already, add the homer programs to your executable path.\n";
	print STDERR "\n\tAdd this line to your .bash_profile or .bashrc file (or other depending on your shell):\n";
	print STDERR "\t\tPATH=" . '$PATH:' . $homeDir . "/bin/\n\n";
	print STDERR "\n\tSimply typing \"findMotifs.pl\" should work before running Homer.\n\n";
	
}

