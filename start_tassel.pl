#! /usr/bin/perl -w

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

system "java -classpath '$CP' -Xms512m -Xmx1536m net.maizegenetics.tassel.TASSELMainApp @args";
