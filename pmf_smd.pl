#!/usr/bin/perl
# interface
#   **************
#   * pmf_smd.pl *
#   **************
#---------------------------------------------------------------------
# This soft computes the potential of mean force (PMF)  from several 
# steered molecular dynamics (SMD) simulations. The PMF is computed
# through the second order cumulant method. Reference of this method
# could be found in:
# "Free energy calculation from steered molecular dynamics simulations
# using Jarzynski's equality". S. Park and F. Khalili-Araghi. 
# J Chem Phys 2003 (119) 3559-66
#
# Usage:
# % pmf_smd.pl ^ basenamefile ^ number
# 
# For example, if you made 20 smd simulations and have files named:
# smd.txt.1 ; smd.txt.2 ... smd.txt.20 
#
# Then:
#
# pmf_smd.pl smd.txt 20
#
# Several information will be asked to the user:
# - The starting value of the X-axis (reaction coordinate)
# - Ending value of the X-axis (reaction coordinate)
# - The number of PMF points
#
# Beware: starting value < ending value
#
# Several files will be created:
# + PMF_curve.dat
# + PMF_gnuplot.gnu
# + PMF_number.png
# 
# PMF_curve.dat is the most important file. This one contains 5 columns
# which represent, respectively, the reaction coordinate, the free energy
# with the first-order cumulant method, the error bar of this last value,
# the free energy with the second-order cumulant method, the error bar of
# this last value and the last column is the Herman error (J Phys Chem 
# 1991 (95) 9029-34. PMF_gnuplot.gnu is a gnuplot script file which is used
# to provide the PMF_number.png.
#
# Enjoy, Florent Barbault
$|;
use strict;
use Math::Trig;
#---------------------------------------------------------------------
# Variable declaration
#
my $basename=$ARGV[0];
my $num_fiche=$ARGV[1];
my $B=1.676;  # value of 1/KbT
my $x0=0;
my $x=0;
my $X=0;
my $pts=0;
my $hist=0;
my $inf=0;
my $sup=0;
my $indice_point=0;
my $choix="y";
my $finput="";
my $fsortie="";
my $fichier="";
my $image="";
my $i=0;
my $serie=0;
my $val_x=0;
my @tab=();
my @tableau=();
my @W=();
my @W2=();
my @W3=();
my @WE=();
my @W2=();
my @WJAR=();
my $we=0.0;
my $ind=0;
my $w=0.0;
my $w2=0.0;
my $w3=0.0;
my $wm=0.0;
my $w2m=0.0;
my $w3m=0.0;
my $dg1=0.0;
my $dg2=0.0;
my $ecg1=0.0;
my $ecg2=0.0;
my $wjar=0.0;
my $herman=0.0;
#---------------------------------------------------------------------
# ask for values of x-axis and number of points
#
print "\n Starting value of x-axis (reaction coordinate)?  ";
$x0=<STDIN>;
chop($x0);
print " Ending value of x-axis (reaction coordinate)?  ";
$x=<STDIN>;
chop($x);
print " How many point for your PMF?  ";
$pts=<STDIN>;
chop($pts);
$hist=($x-$x0)/($pts);
$inf=$x0;
$sup=$x0+$hist;
$X=$x0+$hist/2;
#---------------------------------------------------------------------
#  data collections
#
for ($i=0;$i<$pts;$i++)
{
   $indice_point=$indice_point+1;
   print "\n\t Point number\t$i";
   $fsortie="point_num".$indice_point;
   open (FSOR, "> $fsortie");
   for ($serie=1;$serie<($num_fiche+1);$serie++)
    {
      $fichier=$basename.".".$serie;
      open (FILE, "< $fichier");
      while (<FILE>)
       {
        @tab = split(' ',$_);
        $val_x=$tab[2];
        if ($val_x =~ /^\-?[0-9]+\.?[0-9]*$/)
         {
           if ($val_x >= $inf)
            {
             if ($val_x <= $sup)
              {
                print FSOR "$val_x       $tab[4]\n";
              }
            }
         }
       }
      close(FILE);
    }
   $inf=$inf+$hist;
   $sup=$sup+$hist;
   $X=$X+$hist;
   close(FSOR);
}
#---------------------------------------------------------------------
# calculation
#
open (FOUT, "> PMF_curve.dat");
$X=$x0+$hist/2;
$indice_point=0;
for ($i=0;$i<$pts;$i++)
{
  $ind=0;
  $w=0.0;
  $w2=0.0;
  $w3=0.0;
  $indice_point=$indice_point+1;
  $finput="point_num".$indice_point;
  open (FIN, "< $finput");
  while (<FIN>)
    { 
     @tableau = split(' ',$_);
     $W[$ind]=$tableau[1];
     $W2[$ind]= $tableau[1]*$tableau[1];
     $W3[$ind]= $tableau[1]*$tableau[1]*$tableau[1]; 
     $w=$w+$W[$ind];
     $w2=$w2+$W2[$ind];
     $w3=$w3+$W3[$ind];
     $ind=$ind+1; 
    }
  close (FIN);
  $wm =$w/$ind;
  $w2m=$w2/$ind;
  $w3m=$w3/$ind;
  # computing PMF
  $dg1=$wm;
  $dg2=$wm-($B/2.0)*($w2m-$wm*$wm);
  # Accuracy and errorbars
  $ind=0;
  $we=0.0;
  $wjar=0.0;
  open (FIN, "< $finput");
  while (<FIN>)
    { 
     @tableau = split(' ',$_);
     $WE[$ind]=sqrt(($tableau[1]-$dg1)*($tableau[1]-$dg1));
     $W2[$ind]= $tableau[1]*$tableau[1];
     $WJAR[$ind]=sqrt(($wm-($B/2.0)*($W2[$ind]-$wm*$wm)-$dg2)*($wm-($B/2.0)*($W2[$ind]-$wm*$wm)-$dg2));
     $we=$we+$WE[$ind];
     $wjar=$wjar+$WJAR[$ind];
     $ind=$ind+1; 
   }
  $ecg1=sqrt($we/$ind);
  $ecg2=sqrt($wjar)/$ind;
  $herman=$dg1 - ($B/2.0)*$ecg1*$ecg1;
  close (FIN);
  printf FOUT ("\n %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f",$X,$dg1,$ecg1,$dg2,$ecg2,$herman);
  $X=$X+$hist;
}
close (FOUT);
system ("rm point_num*");
print "\n";
#---------------------------------------------------------------------
# gnuplot PMF
#
open (GNOUT, "> PMF_gnuplot.gnu");
$image="PMF_".$num_fiche.".png";
print GNOUT "set terminal png";
print GNOUT "\nset output \"$image\"";
print GNOUT "\nunset key";
print GNOUT "\nset xlabel \"Reaction Coordinate\"";
print GNOUT "\nset ylabel \"PMF\"";
print GNOUT "\nplot \"PMF_curve.dat\" using 1:4:6 w errorlines lw 2 lc 'red',\\";
for ($serie=1;$serie<($num_fiche);$serie++)
{
  $fichier=$basename.".".$serie;
  print GNOUT "\n \"$fichier\" using 3:5 w li lw 1 lt 0 lc 'black',\\";
}
$fichier=$basename.".".$serie;
print GNOUT "\n \"$fichier\" using 3:5 w li lw 1 lt 0 lc 'black'\n"; #without the last comma
print GNOUT "quit\n";
close(GNOUT);
system("gnuplot < PMF_gnuplot.gnu");




  

