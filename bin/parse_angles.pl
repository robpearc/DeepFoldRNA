chomp(my @pred = `cat $ARGV[0]`);
chomp(my $datadir = "$ARGV[1]");

open(FILE1,">$datadir/eta.txt");
open(FILE2,">$datadir/theta.txt");
open(FILE3,">$datadir/chi.txt");
foreach my $line (@pred){
  chomp(my @splitLine = split(/\s+/,$line));
  next if($splitLine[0] eq "No." || $splitLine[0] eq "");
  print FILE1 "$splitLine[0] $splitLine[9]\n";
  print FILE2 "$splitLine[0] $splitLine[10]\n";
  print FILE3 "$splitLine[0] $splitLine[8]\n";
}
close(FILE1);
close(FILE2);
close(FILE3);

