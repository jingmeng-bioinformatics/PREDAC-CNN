#! usr/bin/perl -w
use strict;
use lib "./";
use naiveBayesian2;



##############################################################################
#This is to predict the antigenic relationship based on the model built
#三个参数，如下：
#modelFile: 模型文件所在路径，如model
#seqFile: 用于预测病毒间关系的病毒序列，如H3N2_seq_forTest
#antigenicGeneticFile: 输出病毒间抗原关系的文件，如antigenicRelation_forTest
##############################################################################
my ($modelFile,$seqFile,$antigenicGeneticFile)=@ARGV;
if(scalar @ARGV < 3){
	print "Please input the modelfile, the sequence file andd the outputFile\n";die;
}

##Model
print "Step 1 : Start reading model...\n";
my $infile=$modelFile;
my ($groups,$site_groups,$chars,$character,$receptor,$threshold,$prior,$distribution)=Model($infile);

#Read the sequences in standard format;
my $seqRef=readSeq($seqFile);

##	predict the Glycosylation Sites
print "Step 2 : Start reading glycosylation sites...\n";
my $glyc={}; #$glyc->{sequence}=\@glycosylation sites
foreach my$id(keys %{$seqRef}){
	my $oneseq=join("",@{$seqRef->{$id}});
	my @pos=ngly($oneseq);
	push(@{$glyc->{$id}},@pos);
}


##	Prediction For each File
print "Step 3 : Predciton ...\n";
open(OUT,">$antigenicGeneticFile") or die "$!";
my @id=keys %{$seqRef};
my @indexes;
my $dataset;
my @aaDiff;
for(my $num1=0;$num1<scalar @id;$num1++){
	my @seq1=@{$seqRef->{$id[$num1]}};
	for(my $num2=$num1+1;$num2<scalar @id;$num2++){
		my @seq2=@{$seqRef->{$id[$num2]}};
		my $diff=0;
		my @compare;
		for(my $i=0;$i<scalar @seq1;$i++){
			if($seq1[$i]=~/[^ACDEFGHIKLMNPQRSTVWY]/ or $seq2[$i]=~/[^ACDEFGHIKLMNPQRSTVWY]/){
				next;
			}
			if($seq1[$i] ne $seq2[$i]){
				$diff++;
				push(@compare,$seq1[$i].($i+1).$seq2[$i]);
			}
		}
		if($diff <=2){                                 #小于等于2个aa差异，直接认为抗原相似
			print OUT $id[$num1],"\t",$id[$num2],"\tMin\t";
			if($diff ==0){
				print OUT "NO\n";
			}else{
				print OUT join("\t",@compare),"\n";
			}
			next;
		}
		if($diff > 9){                                 #大于9个aa差异，直接认为抗原差异
			print OUT $id[$num1],"\t",$id[$num2],"\tMax\t",join("\t",@compare),"\n";
			next;
		}
		push(@aaDiff,join("\t",@compare));
		push(@indexes,$id[$num1]."\t".$id[$num2]);
		push(@{$dataset->[$#indexes]},@compare);
		if(scalar @indexes == 100000){
			my @prediction=AV($dataset,\@indexes,$groups,$site_groups,$chars,$character,$receptor,$threshold,$prior,$distribution,$glyc);
			for(my $i=0;$i<scalar @aaDiff;$i++){
				print OUT $prediction[$i],"\t",$aaDiff[$i],"\n";
			}
			@indexes=();
			$dataset=[];
			@aaDiff=();
		}
	}
}

if(scalar @indexes > 0){
	my @prediction=AV($dataset,\@indexes,$groups,$site_groups,$chars,$character,$receptor,$threshold,$prior,$distribution,$glyc);
	for(my $i=0;$i<scalar @aaDiff;$i++){
		print OUT $prediction[$i],"\t",$aaDiff[$i],"\n";
	}
}	
close OUT or die"$!";

