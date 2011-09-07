#!/usr/bin/perl -w

use strict;
my $libdir = './lib';
opendir (DIR, "$libdir") || die "Could not open $libdir\n";
my @list = readdir(DIR);

my @fl = ();
foreach my $fn(@list){
   if ("$fn" =~ m/\.jar$/){
      push(@fl, "$libdir\/$fn");
   }
}
push(@fl, "./dist/sTASSEL.jar");
my $CP = join(":", @fl);
print $CP;
print "\n";
my @args = @ARGV;

# You may change the the memory allocation (two examples
# are shown at the bottom).
# If you would like to allocate more than 4GB, you need to
# use -d64 to invoke 64-bit version of java

system "java -classpath '$CP' -Xmx512m net.maizegenetics.pipeline.TasselPipeline @args";

# system "java -d64 -classpath $CP -Xms2048m -Xmx8000m net.maizegenetics.pipeline.TasselPipeline";
# qx/java -d64 -classpath $CP -Xms2048m -Xmx8000m net.maizegenetics.pipeline.TasselPipeline;
