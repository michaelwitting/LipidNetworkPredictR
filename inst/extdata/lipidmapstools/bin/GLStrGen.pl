#!/usr/bin/perl -w
#
# File: GLStrGen.pl
# Author: Manish Sud
# Contributor(s): Eoin Fahy
#
# Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the LIPID MAPS (TM) consortium nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# LIPID MAPS (TM) Tools software package is developed by Manish Sud and Eoin Fahy
# at San Diego Supercomputer Center at the University of California, San Diego as
#  a part of LIPID Metabolites and Pathways Strategy consortium efforts.
#
# The LIPID MAPS (TM) consortium is a multi-institutional multi-year effort
# involved in identification and characterization of existing and novel lipids,
# quantifications of changes in their metabolites, and development of biochemical
# pathways and interaction network maps involving lipids. Please acknowledge use
# of the LIPID MAPS (TM) Tools in your work by citing appropriate literature
# references. For further details, please visit www.lipidmaps.org.
#
# Support for this effort is provided by  by National Institutes of Health (NIH)
# and National Institute of General Medical Sciences (NIGMS) Glue Grant
# NIH/NIGMS Grant 1 U54 GM69338.
#

use strict;
use FindBin; use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use LMAPSStr;
use GLStr;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

print "Processing options...\n";
my($SDFileName, $SpecifiedCmpdAbbrevsRef, $WriteSDFile, $AllowArbitraryChainAbbrev);
ProcessOptions();

print "Processing abbreviation(s)...\n";
GLStr::ProcessGLCmpdAbbrevs($SpecifiedCmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileName);

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Process command line options...
sub ProcessOptions {
  # Setup SD file name
  $SDFileName = LMAPSStr::SetupSDFileName('GL', \%Options);
  if (!$Options{overwrite}) {
      if (-e $SDFileName) {
	die "Error: The file $SDFileName already exists\n";
      }
  }
  $SpecifiedCmpdAbbrevsRef = LMAPSStr::SetupCmpdAbbrevs(\%Options);

  $WriteSDFile = ($Options{processmode} =~ /^WriteSDFile$/i) ? 1 : 0;
  $AllowArbitraryChainAbbrev = ($Options{chainabbrevmode} =~ /^Arbitrary$/i) ? 1 : 0;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{chainabbrevmode} = 'MostLikely';
  $Options{mode} = 'Abbrev';
  $Options{processmode} = 'WriteSDFile';

  if (!GetOptions(\%Options, "chainabbrevmode|c=s", "help|h", "mode|m=s", "processmode|p=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{chainabbrevmode} !~ /^(MostLikely|Arbitrary)$/i) {
    die "Error: The value specified, $Options{chainabbrevmode}, for option \"-c, --ChainAbbrevMode\" is not valid. Allowed values: MostLikely or Arbitrary\n";
  }
  if ($Options{mode} !~ /^(Abbrev|AbbrevFileName)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: Abbrev or AbbrevFileName\n";
  }
  if ($Options{processmode} !~ /^(WriteSDFile|CountOnly)$/i) {
    die "Error: The value specified, $Options{processmode}, for option \"-p, --ProcessMode\" is not valid. Allowed values: WriteSDFile or CountOnly\n";
  }
}


__END__

=head1 NAME

GLStrGen.pl - Generate structures for Glycerolipids (GL)

=head1 SYNOPSIS

GLStrGen.pl  GLAbbrev|GLAbbrevFileName ...

GLStrGen.pl [B<-c, --ChainAbbrevMode> I<MostLikely | Arbitrary>]
[B<-h, --help>] [B<-m, --mode> I<Abbrev | AbbrevFileName>]
[B<-p, --ProcessMode> I<WriteSDFile | CountOnly>] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-w, --workingdir> dirname] <arguments>...

=head1 DESCRIPTION

Generate Glycerolipids (GL) structures using compound abbreviations specified on a command line
or in a CSV/TSV Text file. All the command line arguments represent either compound
abbreviations or file name containing abbreviations. Use I<-m, --mode> option to control the type
of command line arguments.

A SD file, containing structures for all GL abbreviations along with ontological information, is
generated as an output.

=head1 SUPPORTED ABBREVIATIONS

Current support for GL structure generation include these main classes and sub classes:

o Monoradylglycerols

    . Monoacylglycerols
    . Monoalkylglycerols
    . Mono-(1Z-alkenyl)-glycerol

o Diradylglycerols

    . Diacylglycerols
    . Alkyl, acylglycerols
    . Dialkylglycerols
    . 1Z-alkenyl, acylglycerols

o Triradylglycerols

    . Triacylglycerols
    . Alkyl, diacylglycerols
    . Dialkyl, monoacylglycerols
    . 1Z-alkenyl, diacylglycerols

=head1 OPTIONS

=over 4

=item B<-c, --ChainAbbrevMode> I<MostLikely|Arbitrary>

Specify what types of acyl chain abbreviations are allowed during processing of complete
abbreviations: allow most likely chain abbreviations containing specific double bond geometry
specifications; allow any acyl chain abbreviation with valid chain length and double bond
geometry specificatios. Possible values: I<MostLikely or Arbitrary>. Default value: I<MostLikely>.

I<Arbitrary> value of B<-c, --ChainAbbrevMode> option is not allowed during processing of
abbreviations containing wild cards.

During I<MostLikely> value of B<-c, --ChainAbbrevMode> option, only the most likely acyl chain
abbreviations specified in ChainAbbrev.pm module are allowed. However, during I<Arbitrary> value
of B<-c, --ChainAbbrevMode> option, any acyl chain abbreviations with valid chain length and
double bond geometry can be specified. The current release of lipidmapstools support chain
lengths from 2 to 50 as specified in ChainAbbev.pm module.

In addition to double bond geometry specifications, valid substituents can be specified for in the acyl
chain abbreviations.

=item B<-h, --help>

Print this help message

=item B<-m, --mode> I<Abbrev|AbbrevFileName>

Controls interpretation of command line arguments. Two different methods are provided:
specify compound abbreviations or a file name containing compound abbreviations. Possible
values: I<Abbrev or AbbrevFileName>. Default: I<Abbrev>

In I<AbbrevFileName> mode, a single line in CSV/TSV files can contain multiple compound
abbreviations. The file extension determines delimiter used to process data lines: comma for
CSV and tab for TSV. For files with TXT extension, only one compound abbreviation per line
is allowed.

Wild card character, *, is also supported in compound abbreviations.

Examples:

    Specific structures: MG(16:0/0:0/0:0) DG(18:1(11E)/16:0/0:0)
                         TG(16:0/16:0/18:1(9Z))
    Specific structures: MG(O-16:0/0:0/0:0) DG(P-16:0/16:0/0:0)
                         TG(O-20:0/16:0/18:1(9Z))
    Specific possibilities: DG(18:*/16:0/0:0) DG(18:1(*)/16:0/0:0)
                            DG(*:*(9Z)/16:0/0:0)
                            DG(*:*(9Z)/*:*(11E)/0:0)
    All TG possibilities: *(*:*/*:*/*:*) or *(*/*/*)
    All MG, DG and TG possibilities: "MG(*:*/0:0/0:0)" "DG(*:*/*:*/0:0)"
                                     "DG(*.*/0:0/*:*)" "TG(*:*/*:*/*:*)"

Along with wild card character, +/- can also be used for chain lengths to indicate even and odd lengths at
sn1/sn2/sn3 positions; additionally > and < qualifiers are also allowed to specify length requirements.
Examples:

    Odd and even number chains at sn1 and sn2: TG(*+:*/*-:*/*:*)
    Odd and even number chains at sn1 and sn2 with length longer than
    10 and 20: TG(*+>10:*/*->20:*/*:*)

Default sn2 stereochemistry is R. However, abbreviation format also supports these additional stereochemistry
specifications for sn2 position: S; U - unknown; rac - racemic mixture. Examples:

    MG(16:0/0:0/0:0)[rac] - racemic mixture
    DG(18:1(11E)/16:0/0:0)[S] - sn2 stereochemistry is S instead of default R
    TG(16:0/16:0/18:1(9Z))[U] - sn2 stereochemistry is unknown

To generate all isomers for specific chains in DG and TG, use of iso designation is also supported.
Stereochemistry specification support is not available with isomeric structure generation.
Examples:

    DG(18:1(11E)/16:0/0:0)[iso2] - Two isomeric structures
    TG(16:0/16:0/18:1(9Z))[iso3] - Three isomeric structures
    TG(16:0/18:0/18:1(9Z))[iso6] - Six isomeric structures

Additionally, all isomeric structures can also be generated by explicit specification of chains at
different positions:

    DG(18:1(11E)/16:0/0:0) DG(16:0/18:1(11E)/0:0)
    TG(16:0/16:0/18:1(9Z)) TG(16:0/18:1(9Z)/16:0)
    TG(18:1(9Z)/16:0/16:0/)

Wild card chain abbreviations are supported with sn2 stereochemistry but not with
isomer abbreviation.

=item B<-p, --ProcessMode> I<WriteSDFile|CountOnly>

Specify how abbreviations are processed: generate structures for specified abbreviations along
with generating a SD file or just count the number of structures corresponding to specified
abbreviations without generating any SD file. Possible values: I<WriteSDFile or CountOnly>.
Default: I<WriteSDFile>.

It can take substantial amount of time for generating all the structures and writing out a SD file
for abbreviations containing wild cards. I<CountOnly> value of B<--ProcessMode> option can
be used to get a quick count of number of structures to be generated without writing out any
SD file.

=item B<-o, --overwrite>

Overwrite existing files

=item B<-r, --root> I<rootname>

New file name is generated using the root: <Root>.sdf. Default for new file names: GLAbbrev.sdf,
<AbbrevFilenName>.sdf, or <FirstAbbrevFileName>1To<Count>.sdf.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory

=back

=head1 EXAMPLES

On some systems, command line scripts may need to be invoked using
I<perl -s GLStrGen.pl>; however, all the examples assume direct invocation
of command line script works.

To generate a GLStructures.sdf file containing a structure specified
by a command line GL abbreviation, type:

    % GLStrGen.pl -r GLStructures -o "MG(16:0/0:0/0:0)"

To generate a GLStructures.sdf file containing structures specified
by a command line GL abbreviations, type:

    % GLStrGen.pl -r GLStructures -o "MG(16:0/0:0/0:0)"
    "DG(18:1(11E)/16:0/0:0)" "TG(16:0/16:0/18:1(9Z))"

To generate a GLStructures.sdf file containing structures specified
by a command line GL abbreviations with specific stereochemistry, type:

    % GLStrGen.pl -r GLStructures -o "MG(16:0/0:0/0:0)[rac]"
    "DG(18:1(11E)/16:0/0:0)[S]" "TG(16:0/16:0/18:1(9Z))[U]"

To generate a GLStructures.sdf file containing all isomeric structures specified
by a command line GL abbreviations, type:

    % GLStrGen.pl -r GLStructures -o "DG(18:1(11E)/16:0/0:0)[iso2]"
    "TG(16:0/16:0/18:1(9Z))[iso3]" "TG(16:0/17:0/18:1(9Z))[iso6]"

To enumerate all possible GL structures and generate a GLStructures.sdf
file, type:

    % GLStrGen.pl -r GLStructures -o "*(*/*/*)"

or

    % GLStrGen.pl -r GLStructures -o "*(*:*/*:*/*:*)"

or

    % GLStrGen.pl -r GLStructures -o "*(*:*(*)/*:*(*)/*:*(*))"

To enumerate all possible Monoradylglycerols structures and generate a
MonoGLStructures.sdf file, type:

    % GLStrGen.pl -r MonoGLStructures -o "MG(*/0:0/0:0)"

To enumerate all possible Diradylglycerols structures and generate a
DiGLStructures.sdf file, type:

    % GLStrGen.pl -r DiGLStructures -o "DG(*/*/0:0)"

To enumerate all possible Monoradylglycerols structures with one double
bond on acyl chain and generate a GLStructures.sdf file, type:

    % GLStrGen.pl -r MonoGLStructures -o "MG(*:1/0:0/0:0)"

To enumerate all possible Monoradylglycerols structures with even chain
lengths and generate a GLStructures.sdf file, type:

    % GLStrGen.pl -r MonoGLStructures -o "MG(*+:*/0:0/0:0)"

To enumerate all possible Diradylglycerols structures with odd chains longer
than 10 at sn1 and even chains longer than 18 at sn2, and generate a
DiGLStructures.sdf file, type:

    % GLStrGen.pl -r DiGLStructures -o "DG(*->10:*/*+>18:*/0:0)"

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

CLStrGen.pl, FAStrGen.pl, GPStrGen.pl, SPStrGen.pl, STStrGen.pl

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
