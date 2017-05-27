#!/usr/bin/perl -w

#===============================================================================
#          ___                            _   _      
#         / _ \                          | | (_)     
#     ___/ /_\ \_ __ ___  _ __ ___   __ _| |_ _  ___ 
#    / _ \  _  | '__/ _ \| '_ ` _ \ / _` | __| |/ __|
#   |  __/ | | | | | (_) | | | | | | (_| | |_| | (__ 
#    \___\_| |_/_|  \___/|_| |_| |_|\__,_|\__|_|\___|
#                                                  
#   eAromatic - analysis of aromatic interactions in protein-ligand complexes
#
#   This software is distributed WITHOUT ANY WARRANTY (but with best wishes)
#
#   Report bugs and issues to michal@brylinski.org
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#===============================================================================

 use strict;
 use File::Slurp;
 use Chemistry::Mol;
 use Chemistry::File::MDLMol;
 use Chemistry::Ring::Find ':all';
 use Math::Vector::Real;
 use Math::Trig;
 
 local $| = 1;
 
 print "---------------------------------------------------------------\n";
 print "                           eAromatic\n";
 print "                          version 1.0\n";
 print " analysis of aromatic interactions in protein-ligand complexes\n\n";
 print "        report bugs and issues to michal\@brylinski.org\n";
 print "---------------------------------------------------------------\n\n";
 
 if ($#ARGV < 3)
 {
  print "eAromatic.pl -p <protein in PDB format>\n";
  print "             -l <ligand in MOL format>\n\n";
  print "Optional arguments:\n";
  print "             -c <output from LPC>\n";
  print "             -t <contact distance in A, default 4.5>\n";
  die "\n";
 }
 
 my $fpro1 = '';
 my $flig1 = '';
 my $flpc1 = '';
 
 my $fcut1 = 4.5;
 
 for ( my $i = 0; $i <= $#ARGV; $i++ )
 {
  $fpro1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-p' );
  $flig1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-l' );
  $flpc1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-c' );
  $fcut1 = $ARGV[$i+1] if ( $ARGV[$i] eq '-t' );
 }
 
 die "Provide protein filename\n" if ( !( -e $fpro1 ) or !length($fpro1) );
 die "Provide ligand filename\n" if ( !( -e $flig1 ) or !length($flig1) );
 
 die "Cannot find output from LPC: $flpc1\n" if ( !( -e $flpc1 ) and length($flpc1) );
 
 my @pdb1 = read_file($fpro1); chomp(@pdb1);
 
 my @pdb2 = grep(/PHE|TRP|TYR|HIS/, @pdb1);
 
 push(@pdb2, 'ATOM      1  N  ');
 
 my %pdb3 = ();
 
 my @pdb4 = ();
 
 foreach my $wpdb2 (@pdb2)
 {
  if ( substr($wpdb2, 12, 4) eq ' N  ' )
  {
   if ( @pdb4 )
   {
    if ( substr($pdb4[0], 17, 3) eq 'TRP' )
    {
     my @xyz1 = ( 0.0, 0.0, 0.0 );
     my @xyz2 = ( 0.0, 0.0, 0.0 );
     
     my @pla1 = ();
     my @pla2 = ();
     
     my $n1 = 0;
     my $n2 = 0;
     
     foreach my $wpdb4 (@pdb4)
     {
      my $at1 = substr($wpdb4, 12, 4);
      
      if ( $at1 eq ' CG ' or $at1 eq ' CD1' or $at1 eq ' CD2' or $at1 eq ' NE1' or $at1 eq ' CE2' )
      {
       $xyz1[0] += substr($wpdb4, 30, 8) * 1.0;
       $xyz1[1] += substr($wpdb4, 38, 8) * 1.0;
       $xyz1[2] += substr($wpdb4, 46, 8) * 1.0;
       
       $pla1[$n1][0] = substr($wpdb4, 30, 8) * 1.0;
       $pla1[$n1][1] = substr($wpdb4, 38, 8) * 1.0;
       $pla1[$n1][2] = substr($wpdb4, 46, 8) * 1.0;
       
       $n1++;
      }
      
      if ( $at1 eq ' CD2' or $at1 eq ' CE2' or $at1 eq ' CE3' or $at1 eq ' CZ2' or $at1 eq ' CZ3' or $at1 eq ' CH2' )
      {
       $xyz2[0] += substr($wpdb4, 30, 8) * 1.0;
       $xyz2[1] += substr($wpdb4, 38, 8) * 1.0;
       $xyz2[2] += substr($wpdb4, 46, 8) * 1.0;
       
       $pla2[$n2][0] = substr($wpdb4, 30, 8) * 1.0;
       $pla2[$n2][1] = substr($wpdb4, 38, 8) * 1.0;
       $pla2[$n2][2] = substr($wpdb4, 46, 8) * 1.0;
       
       $n2++;
      }
     }
     
     $xyz1[0] /= $n1;
     $xyz1[1] /= $n1;
     $xyz1[2] /= $n1;
     
     $xyz2[0] /= $n2;
     $xyz2[1] /= $n2;
     $xyz2[2] /= $n2;
     
     my $at2 = substr($pdb4[0], 22, 4) * 1;
     my $at3 = substr($pdb4[0], 17, 3);
     
     $pdb3{"$at2:$at3:1"}[0] = $xyz1[0];
     $pdb3{"$at2:$at3:1"}[1] = $xyz1[1];
     $pdb3{"$at2:$at3:1"}[2] = $xyz1[2];
     
     $pdb3{"$at2:$at3:2"}[0] = $xyz2[0];
     $pdb3{"$at2:$at3:2"}[1] = $xyz2[1];
     $pdb3{"$at2:$at3:2"}[2] = $xyz2[2];
     
     my %dst1 = ();
     my %dst2 = ();
     
     for ( my $xa = 0; $xa < $n1; $xa++ )
     {
      $dst1{"$pla1[$xa][0]:$pla1[$xa][1]:$pla1[$xa][2]"} = sqrt( ( $pla1[$xa][0] - $xyz1[0] )**2 + ( $pla1[$xa][1] - $xyz1[1] )**2 + ( $pla1[$xa][2] - $xyz1[2] )**2 );
     }
     
     for ( my $xa = 0; $xa < $n2; $xa++ )
     {
      $dst2{"$pla2[$xa][0]:$pla2[$xa][1]:$pla2[$xa][2]"} = sqrt( ( $pla2[$xa][0] - $xyz2[0] )**2 + ( $pla2[$xa][1] - $xyz2[1] )**2 + ( $pla2[$xa][2] - $xyz2[2] )**2 );
     }
     
     my @pnt1 = ();
     my @pnt2 = ();
     
     my $pn1 = 0;
     my $pn2 = 0;
     
     foreach my $wdst1 ( sort { $dst1{$b} <=> $dst1{$a} } keys %dst1 )
     {
      my @tt1 = split(/\:/, $wdst1);
      
      $pnt1[$pn1][0] = $tt1[0];
      $pnt1[$pn1][1] = $tt1[1];
      $pnt1[$pn1][2] = $tt1[2];
      
      $pn1++;
     }
     
     foreach my $wdst2 ( sort { $dst2{$b} <=> $dst2{$a} } keys %dst2 )
     {
      my @tt2 = split(/\:/, $wdst2);
      
      $pnt2[$pn2][0] = $tt2[0];
      $pnt2[$pn2][1] = $tt2[1];
      $pnt2[$pn2][2] = $tt2[2];
      
      $pn2++;
     }
     
     my $va1 = V( $pnt1[1][0] - $pnt1[0][0], $pnt1[1][1] - $pnt1[0][1], $pnt1[1][2] - $pnt1[0][2] );
     my $vb1 = V( $pnt1[2][0] - $pnt1[0][0], $pnt1[2][1] - $pnt1[0][1], $pnt1[2][2] - $pnt1[0][2] );
     
     my $vc1 = $va1 x $vb1;
     
     $vc1 /= $vc1->norm();
     
     my $va2 = V( $pnt2[1][0] - $pnt2[0][0], $pnt2[1][1] - $pnt2[0][1], $pnt2[1][2] - $pnt2[0][2] );
     my $vb2 = V( $pnt2[2][0] - $pnt2[0][0], $pnt2[2][1] - $pnt2[0][1], $pnt2[2][2] - $pnt2[0][2] );
     
     my $vc2 = $va2 x $vb2;
     
     $vc2 /= $vc2->norm();
     
     $pdb3{"$at2:$at3:1"}[3] = $vc1;
     
     $pdb3{"$at2:$at3:2"}[3] = $vc2;
    }
    else
    {
     my @xyz1 = ( 0.0, 0.0, 0.0 );
     
     my @pla1 = ();
     
     my $n1 = 0;
     
     foreach my $wpdb4 (@pdb4)
     {
      my $at1 = substr($wpdb4, 12, 4);
      
      if ( $at1 ne ' N  ' and $at1 ne ' CA ' and $at1 ne ' C  ' and $at1 ne ' O  ' and $at1 ne ' CB ' and $at1 ne ' OH ' )
      {
       $xyz1[0] += substr($wpdb4, 30, 8) * 1.0;
       $xyz1[1] += substr($wpdb4, 38, 8) * 1.0;
       $xyz1[2] += substr($wpdb4, 46, 8) * 1.0;
       
       $pla1[$n1][0] = substr($wpdb4, 30, 8) * 1.0;
       $pla1[$n1][1] = substr($wpdb4, 38, 8) * 1.0;
       $pla1[$n1][2] = substr($wpdb4, 46, 8) * 1.0;
       
       $n1++;
      }
     }
     
     $xyz1[0] /= $n1;
     $xyz1[1] /= $n1;
     $xyz1[2] /= $n1;
     
     my $at2 = substr($pdb4[0], 22, 4) * 1;
     my $at3 = substr($pdb4[0], 17, 3);
     
     $pdb3{"$at2:$at3"}[0] = $xyz1[0];
     $pdb3{"$at2:$at3"}[1] = $xyz1[1];
     $pdb3{"$at2:$at3"}[2] = $xyz1[2];
     
     my %dst1 = ();
     
     for ( my $xa = 0; $xa < $n1; $xa++ )
     {
      $dst1{"$pla1[$xa][0]:$pla1[$xa][1]:$pla1[$xa][2]"} = sqrt( ( $pla1[$xa][0] - $xyz1[0] )**2 + ( $pla1[$xa][1] - $xyz1[1] )**2 + ( $pla1[$xa][2] - $xyz1[2] )**2 );
     }
     
     my @pnt1 = ();
     
     my $pn1 = 0;
     
     foreach my $wdst1 ( sort { $dst1{$b} <=> $dst1{$a} } keys %dst1 )
     {
      my @tt1 = split(/\:/, $wdst1);
      
      $pnt1[$pn1][0] = $tt1[0];
      $pnt1[$pn1][1] = $tt1[1];
      $pnt1[$pn1][2] = $tt1[2];
      
      $pn1++;
     }
     
     my $va1 = V( $pnt1[1][0] - $pnt1[0][0], $pnt1[1][1] - $pnt1[0][1], $pnt1[1][2] - $pnt1[0][2] );
     my $vb1 = V( $pnt1[2][0] - $pnt1[0][0], $pnt1[2][1] - $pnt1[0][1], $pnt1[2][2] - $pnt1[0][2] );
     
     my $vc1 = $va1 x $vb1;
     
     $vc1 /= $vc1->norm();
     
     $pdb3{"$at2:$at3"}[3] = $vc1;
    }
   }
   
   @pdb4 = ();
  }
  
  push(@pdb4, $wpdb2);
 }
 
 my $mol1 = Chemistry::Mol->read($flig1);
 
 Chemistry::Ring::aromatize_mol($mol1);
 
 my @rings = find_rings($mol1, sssr => 1);
 
 my $nrings = @rings;
 
 my @aro1 = ();
 
 my $naro1 = 0;
 
 for ( my $xa = 0; $xa < $nrings; $xa++ )
 {
  if ( $rings[$xa]->is_aromatic )
  {
   my @atoms = $rings[$xa]->atoms;
   
   my $natoms = @atoms;
   
   if ( $natoms == 5 or $natoms == 6 )
   {
    my $at5 = '';
    
    my @pla1 = ();
    
    my $npla1 = 0;
    
    foreach my $watoms (@atoms)
    {
     my $at4 = $watoms;
     
     while ( $at4 =~ /a/ ) { $at4 =~ s/a//g; }
     
     if ( length($at5) )
     {
      $at5 .= ':'.$at4;
     }
     else
     {
      $at5 = $at4;
     }
     
     $pla1[$npla1][0] = $watoms->x3();
     $pla1[$npla1][1] = $watoms->y3();
     $pla1[$npla1][2] = $watoms->z3();
     
     $npla1++;
    }
    
    my $cen1 = $rings[$xa]->centroid;
    
    $aro1[$naro1][0] = $cen1->x();
    $aro1[$naro1][1] = $cen1->y();
    $aro1[$naro1][2] = $cen1->z();
    
    $aro1[$naro1][4] = $at5;
    
    my %dst1 = ();
    
    for ( my $xb = 0; $xb < $npla1; $xb++ )
    {
     $dst1{"$pla1[$xb][0]:$pla1[$xb][1]:$pla1[$xb][2]"} = sqrt( ( $pla1[$xb][0] - $cen1->x() )**2 + ( $pla1[$xb][1] - $cen1->y() )**2 + ( $pla1[$xb][2] - $cen1->z() )**2 );
    }
     
    my @pnt1 = ();
    
    my $pn1 = 0;
    
    foreach my $wdst1 ( sort { $dst1{$b} <=> $dst1{$a} } keys %dst1 )
    {
     my @tt1 = split(/\:/, $wdst1);
     
     $pnt1[$pn1][0] = $tt1[0];
     $pnt1[$pn1][1] = $tt1[1];
     $pnt1[$pn1][2] = $tt1[2];
     
     $pn1++;
    }
    
    my $va1 = V( $pnt1[1][0] - $pnt1[0][0], $pnt1[1][1] - $pnt1[0][1], $pnt1[1][2] - $pnt1[0][2] );
    my $vb1 = V( $pnt1[2][0] - $pnt1[0][0], $pnt1[2][1] - $pnt1[0][1], $pnt1[2][2] - $pnt1[0][2] );
    
    my $vc1 = $va1 x $vb1;
    
    $vc1 /= $vc1->norm();
    
    $aro1[$naro1][3] = $vc1;
    
    $naro1++;
   }
  }
 }
 
 if ( $naro1 )
 {
  my @con1 = ();
  
  if ( -e $flpc1 and length($flpc1) )
  {
   my @lpc1 = read_file($flpc1); chomp(@lpc1);
   
   my $w1 = 0;
   
   foreach my $wlpc1 (@lpc1)
   {
    $w1 = 0 if ( $wlpc1 =~ /Complementarity values for the ligand/ );
    
    if ( $w1 == 2 and length($wlpc1) and !( $wlpc1 =~ /----------/ ) )
    {
     my $lig1 = substr($wlpc1, 0, 3) * 1;
     
     my $pro1 = substr($wlpc1, 21, 3);
     my $pro2 = substr($wlpc1, 25, 4) * 1;
     my $pro3 = substr($wlpc1, 34, 5);
     
     while ( $pro3 =~ /\ / ) { $pro3 =~ s/\ //g; }
     
     push(@con1, "$lig1:$pro1:$pro2:$pro3") if ( $pro1 eq 'PHE' or $pro1 eq 'TRP' or $pro1 eq 'TYR' or $pro1 eq 'HIS' );
    }
    
    $w1++ if ( $wlpc1 =~ /  N   Name   Class    Residue       Name   Class/ );
   }
  }
  
  else
  {
   foreach my $wpdb2 (@pdb2)
   {
    if ( length($wpdb2) > 53 )
    {
     my $x1 = substr($wpdb2, 30, 8) * 1.0;
     my $y1 = substr($wpdb2, 38, 8) * 1.0;
     my $z1 = substr($wpdb2, 46, 8) * 1.0;
     
     my $pro1 = substr($wpdb2, 17, 3);
     my $pro2 = substr($wpdb2, 22, 4) * 1;
     my $pro3 = substr($wpdb2, 12, 4);
     
     while ( $pro3 =~ /\ / ) { $pro3 =~ s/\ //g; }
     
     my @mol2 = read_file($flig1); chomp(@mol2);
     
     my $m1 = substr($mol2[3], 0, 3) * 1;
     
     for ( my $xa = 4; $xa < $m1 + 4; $xa++ )
     {
      my $lig1 = $xa - 3;
      
      my $x2 = substr($mol2[$xa],  0, 10) * 1.0;
      my $y2 = substr($mol2[$xa], 10, 10) * 1.0;
      my $z2 = substr($mol2[$xa], 20, 10) * 1.0;
      
      my $r1 = sqrt( ( $x1 - $x2 )**2 + ( $y1 - $y2 )**2 + ( $z1 - $z2 )**2 );
      
      if ( $r1 <= $fcut1 )
      {
       push(@con1, "$lig1:$pro1:$pro2:$pro3") if ( $pro1 eq 'PHE' or $pro1 eq 'TRP' or $pro1 eq 'TYR' or $pro1 eq 'HIS' );
      }
     }
    }
   }
  }
  
  my %aro2 = ();
  
  my $naro2 = 0;
  
  foreach my $wcon1 ( sort @con1 )
  {
   my @con2 = split(/\:/, $wcon1);
   
   my $lig1 = $con2[0];
   my $pro1 = $con2[1];
   my $pro2 = $con2[2];
   my $pro3 = $con2[3];
   
   if ( $pro3 ne 'N' and $pro3 ne 'CA' and $pro3 ne 'C' and $pro3 ne 'O' and $pro3 ne 'CB' and $pro3 ne 'OH' )
   {
    my $at7 = 0;
    my $at8 = 0;
    
    if ( $pro1 eq 'TRP' )
    {
     if ( $pro3 eq 'CG' or $pro3 eq 'CD1' or $pro3 eq 'CD2' or $pro3 eq 'NE1' or $pro3 eq 'CE2' )
     {
      $at7 = 1;
     }
     
     if ( $pro3 eq 'CD2' or $pro3 eq 'CE2' or $pro3 eq 'CE3' or $pro3 eq 'CZ2' or $pro3 eq 'CZ3' or $pro3 eq 'CH2' )
     {
      $at8 = 1;
     }
    }
    
    for ( my $xa = 0; $xa < $naro1; $xa++ )
    {
     my @at6 = split(/\:/, $aro1[$xa][4]);
     
     my $w2 = 0;
     
     foreach my $wat6 (@at6)
     {
      $w2 = 1 if ( $lig1 == $wat6 );
     }
     
     if ( $w2 )
     {
      if ( $pro1 eq 'TRP' )
      {
       if ( $at7 )
       {
        $aro2{"$xa?$pro2:$pro1:1"} = 1;
        
        $naro2++;
       }
       if ( $at8 )
       {
        $aro2{"$xa?$pro2:$pro1:2"} = 1;
        
        $naro2++;
       }
      }
      else
      {
       $aro2{"$xa?$pro2:$pro1"} = 1;
       
       $naro2++;
      }
     }
    }
   }
  }
  
  if ( $naro2 )
  {
   foreach my $waro2 ( keys %aro2 )
   {
    my @aro3 = split(/\?/, $waro2);
    
    my $x1 = $aro1[$aro3[0]][0];
    my $y1 = $aro1[$aro3[0]][1];
    my $z1 = $aro1[$aro3[0]][2];
    
    my $x2 = $pdb3{$aro3[1]}[0];
    my $y2 = $pdb3{$aro3[1]}[1];
    my $z2 = $pdb3{$aro3[1]}[2];
    
    my @aro4 = split(/\:/, $aro3[1]);
    
    my @aro5 = split(/\:/, $aro1[$aro3[0]][4]);
    
    my $naro5 = @aro5;
    
    my $aa1 = $aro4[1];
    
    my $aa2 = '';
    
    if ( $aa1 eq 'TRP' )
    {
     if ( $aro4[2] == 1 )
     {
      $aa2 = 'TR5';
     }
     else
     {
      $aa2 = 'TR6';
     }
    }
    else
    {
     $aa2 = $aa1;
    }
    
    my $r1 = sqrt( ( $x1 - $x2 )**2 + ( $y1 - $y2 )**2 + ( $z1 - $z2 )**2 );
    
    my $a1 = $aro1[$aro3[0]][3] * $pdb3{$aro3[1]}[3];
    
    my $a2 = rad2deg(acos($a1));
    
    $a2 = 180.0 - $a2 if ( $a2 > 90.0 );
    
    # aa (TR5: TRP 1st ring, TR6: TRP 2nd ring) - residue index - ligand ring number - distance ( ring centers ) - angle ( ring planes, normal vector product ) - n-member ring - ligand ring atoms
    
    printf("AROMATIC %s %4d %2d %7.4f %7.2f %2d %s\n", $aa2, $aro4[0], $aro3[0], $r1, $a2, $naro5, $aro1[$aro3[0]][4]);
   }
  }
 }
 
 exit(0);
