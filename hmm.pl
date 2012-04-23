# function calculate_emission_probabilities
# calculates emission probabilities
# from input text file
# the useful text is marked as <u>...</u>
# number of hidden states = 2

# This script is doing following
# calsulates emission probabilities
# read text message for decoding
# applies viterbi algorithm for
# decoding text message
#!/usr/bin/perl -w

calculate_emission_probabilities();



for($i=0; $i<2; $i++)
{
  for($j=0; $j<@words; $j++)
  {
     print  $ep[$i][$j];
     print " ";
     print "\n";
  }
  print "\n";
}



 # Read new text message for decoding

 $dr="C:\\Users\\Owner\\Desktop\\A\\HMM";
 $file="115.txt";
 
 $file="$dr\\$file";
 $count_word=0;

       open (FH, $file) || die "$!" ;
       
  
       while($ln = <FH>) 
       {
          $ln =~ s/!|	|\n|chr(10)/ /g;            
          $ln =~ s/( )*/ /;
        

          @words=split(" ", $ln);
          
         

          foreach $w (@words)
          {

              $y[$count_word]=$w_index{$w};
              $yw[$count_word]=$w;
              $count_word=$count_word+1;
          }
       }
   close (FH);



forward_viterbi();



sub calculate_emission_probabilities ()
{
 $word_count=0;

$MAX_words=1000;
$number_of_states=2;

for($i=0; $i<$number_of_states; $i++)
{
  for($j=0; $j<$MAX_words; $j++)
  {
     $emission_probability[$i][$j]=0;
  }
}

 $dr="C:\\Users\\Owner\\Desktop\\A\\HMM";
 $file="111.txt";
 
 $file="$dr\\$file";

       open (FH, $file) || die "$!" ;
       
  
       while($ln = <FH>) 
       {

          


          $ln =~ s/!|	|\n|chr(10)/ /g;            
          $ln =~ s/( )*/ /;
        

          @words=split(" ", $ln);
          
          $current_hidden_state_index=0;

          foreach $w (@words)
          {
               
             if ($w =~ m/<u>/)
             {
               $w =~ s/<u>//;
               $current_hidden_state_index=1;
             } 

                      $w1=$w;
                      $w1 =~ m/<\/u>/;


                      if (!exists $w_index{$w1} )
                      {
                             $w_index{$w1}=$word_count;
                             $word_count=$word_count+1;
                      }
                     
                            
      $emission_probability[$current_hidden_state_index][$w_index{$w1}]=$emission_probability[$current_hidden_state_index][$w_index{$w1}]+1;
                       
 

             if ($w =~ m/<\/u>/)
             {
               $w =~ s/<\/u>//;
               $current_hidden_state_index=0;
             } 
         }
       }

close(FH);

$temp[0]=0;
$temp[1]=0;

$temp_total=0;
for($i=0; $i<$number_of_states; $i++)
{
  for($j=0; $j<$word_count-1; $j++)
  {
     
      $temp[$i]=$temp[$i]+$emission_probability[$i][$j];

      $temp_total=$temp_total + $emission_probability[$i][$j];

  }
}


for($i=0; $i<$number_of_states; $i++)
{
  for($j=0; $j<$word_count-1; $j++)
  {
    $ep[$i][$j]=$emission_probability[$i][$j] / $temp[$i];
  }
}


}




sub forward_viterbi ()
{


#  INPUT 
#  y - the sequence of observation
#  X - the set of hidden states
#  sp - start probability
#  tp - transition probabilities
# 
#  states: Bull, Bear, Even
#  ep - emission probabilities
#  
#  



@X=("Useless", "Useful");
@sp=(0.5, 0.5);

@tp[0]=[0.5, 0.5];
@tp[1]=[0.5,  0.5];


for ($i=0; $i<@X; $i++)
{
   @T[$i]=[$sp[$i], $X[$i], $sp[$i]];
}


for ($i=0; $i<@y; $i++)
{
   for ($next_state=0; $next_state < @X; $next_state++)
   {
       $total=0;
       $argmax=None;
       $valmax=0;
       for ($state=0; $state<@X; $state++)
       {
           $prob=$T[$state][0];
           $v_path=$T[$state][1];
           $v_prob=$T[$state][2];


       print $y[$i];
       print " ";
       print $ep[$state][$y[$i]];
       print " ";

           $p=$ep[$state][$y[$i]] * $tp[$state][$next_state];
           $prob=$prob * $p;
           $v_prob=$v_prob * $p;
           $total=$total + $prob;
           if ($v_prob > $valmax)
               {
                  $argmax = $v_path."_".$X[$next_state];
                  $valmax=$v_prob;
               }
       }
       @U[$next_state]= [$total, $argmax, $valmax];
   }
   @T=@U;
}

$total=0;
$argmax="";
$valmax=0;
for ($i=0; $i<@X; $i++)
{


   $prob=$T[$i][0];
   $v_path=$T[$i][1];
   $v_prob=$T[$i][2];


   $total=$total + $prob ;
   if ($v_prob > $valmax)
      {
          $argmax=$v_path;
          $valmax=$v_prob
      }
}

print "Total        = $total\n";
print "Viterbi prob.= $valmax\n";
print "Viterbi path = $argmax\n";


$p="";
@w=split ("_", $argmax);



for ($i=0; $i < @w; $i++)
{
   if ($w[$i] eq "Useful")
   {
       $p=$p.$yw[$i]." ";
   }
}

print "\n\n$p";

return ($total, $argmax, $valmax);
}       
            