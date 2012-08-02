#!/usr/bin/perl -w

use strict;
use File::Basename;

my $top = dirname(__FILE__);
$top //= '.';

my $libdir = "$top/lib";
opendir (DIR, "$libdir") || die "Could not open $libdir\n";
my @list = readdir(DIR);

my @fl = ();
foreach my $fn(@list){
   if ("$fn" =~ m/\.jar$/){
      push(@fl, "$libdir\/$fn");
   }
}
push(@fl, "$top/dist/sTASSEL.jar");
my $CP = join(":", @fl);
print $CP . "\n";

# Scan @ARGV for Java memory arguments, and put rest in @args
my $java_mem_default = "-Xms512m -Xmx1536m";
my $java_mem = "";
my @args;
for (my $i=0; $i<=$#ARGV; $i++){
   if ($ARGV[$i] =~ m/Xm/) {
      $java_mem .= "$ARGV[$i] ";
   }
   else{
      push(@args, $ARGV[$i]);
   }
}

if ($java_mem eq "") { $java_mem = $java_mem_default; }

print "Memory Settings: " . $java_mem . "\n";
if (@args != 0){
   print "Tassel Pipeline Arguments: " . "@args\n";
}

system "java -classpath '$CP' $java_mem net.maizegenetics.tassel.TASSELMainApp @args";
