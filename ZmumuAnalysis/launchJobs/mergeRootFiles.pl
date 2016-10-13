#!/usr/bin/perl
$step = $ARGV[0];
$dataset = $ARGV[1];
print "-------------------------------\n";
print "step: $step - dataset: $dataset \n";
print "-------------------------------\n";


chdir $dataset;
#system("cd $dataset");
@FILES= <*>; 
$tmp = 1;
@toMerge = ();
#maximum number of tmp
$tmp_limit = int(scalar(@FILES)/$step)+1;
print $tmp_limit;

for ($i = 0; $i < scalar(@FILES); $i++) {
	push(@toMerge, "@FILES[$i]");
	if (scalar(@toMerge) == $step) {
	    $in = ($tmp-1)*$step;
	    $fn = $tmp*$step;
	    print "Merging from $in to $fn...\n";
	    system("hadd /home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs/DYJets_$tmp.root @toMerge");
	    @toMerge = ();
	    $tmp++;
	}
        if (($tmp_limit-$tmp+1) == 1 && $i == scalar(@FILES)-1){
	    $in = ($tmp-1)*$step;
	    $fn = $tmp*$step;
	    print "Merging last bunch of files...\n";
	    system("hadd /home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs/DYJets_$tmp.root @toMerge");
	    @toMerge = ();
        }
}

@TMP=<TMP*>;
@toMergeTMP = ();
for ($i = 0; $i < scalar(@TMP); $i++) {
   push(@toMergeTMP, "@TMP[$i]");
}

#system("rm -f avoidText");
