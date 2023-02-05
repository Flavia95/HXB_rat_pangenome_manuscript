while (<>){
	$line=$_; # one line at a time
	chomp($line); # remove the line break
	if ($line =~/>/){ # when the line is the fasta name 
		$line=~/>.+-(.+)_\w\d\d.ab1/; # extract run name from middle part of the chars 
		$id=$1;
		if (length($id)==8) { 	# these are reruns	
			$rerun="_rerun";
			$id=~s/_R$//; # remove the  "_R" rerun indicator 
		} else {
			$rerun="";
		}
		print "\n", $id, "\t",$line, $rerun, "\t";
	} else  { # when the line is fasta sequence
		print $line;
	}
}

