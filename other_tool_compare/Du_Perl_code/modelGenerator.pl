#! usr/bin/perl -w
use strict;
use lib "./";
use naiveBayesian2;


#####################################################################
#This is to generate the naive bayes model as that in Du's paper
#输入参数为5个文件路径以及1个输出的模型文件路径，如下：
#epitopeFile: 抗原表位文件所在路径，如epitope_H3N2
#aaIndexFile:aaIndex所在路径，如aaIndex
#receptorFile: 受体文件所在路径，如receptor_H3N2
#seqFile: 训练数据涉及的病毒序列文件，如H3N2_seq
#seqCompareFile:训练数据中病毒间aa差异文件，如seqCompare_H3N2
#modelFile: 输出的模型文件，如model
####################################################################
my ($epitopeFile,$aaIndexFile,$receptorFile,$seqFile,$seqCompareFile,$modelFile)=@ARGV;
if(scalar @ARGV < 6){
	print "Please input the epitopeFile, aaIndexFile,receptorFile,sequence File, sequence compare file and model file\n";die;
}


##	Feature : Site Groups
my @groups;
my %groups;
my $group_sites;

#my $infile="epitope_H3N2";
my $infile=$epitopeFile;
open(IN,"<$infile")||die "$!";
while(<IN>){
	chomp;
	next if(/^$/);
	next if(/^\#/);
	if(/^[A-Z]/i){
		@groups=split(/\t/,$_);
	}elsif(/^[0-9]/){
		my @line=split(/\t/);
		push(@{$group_sites->{$line[1]}},$line[0]);
		$groups{$line[0]}=$line[1];
	}
}
close(IN)||die "$!";
@groups=sort @groups;


##	AA characters
my @chars=qw(FASG890101 GRAR740103 ZIMJ680104 JANJ780101 CHAM820101); 
my $character; #$character->{index}->{aa}=value

#$infile="aaIndex";
$infile=$aaIndexFile;
open(IN,"<$infile")||die "$!";
$_=<IN>;
chomp;
my @aa=split(/\t/);
shift @aa;
pop @aa;
while(<IN>){
	chomp;
	next if(/^$/);
	my @lineArray = split(/\t/);
	my $char=shift @lineArray;
	pop @lineArray;
	die "Bad line formate in file properties : $_\n" if($#lineArray != $#aa);
	for(my $i=0;$i<=$#aa;$i++){
		$character->{$char}->{$aa[$i]}=$lineArray[$i];
	}
}
close(IN)||die "$!";

##	Receptor Influence
my %receptor;
#$infile="receptor_H3N2";
$infile=$receptorFile;
open(IN,"<$infile")||die "$!";
while(<IN>){
	chomp;
	next if(/^$/);
	my @lineArray = split(/\t/);
	$receptor{ $lineArray[0] } = $lineArray[1];
}
close(IN)||die "$!";

#Read the glycosylated sites;
my $glyc={}; #$glyc->{sequence}=\@glycosylation sites
#my $seqRef=readSeq("H3N2_seq");
my $seqRef=readSeq($seqFile);
foreach my$id(keys %{$seqRef}){
	my $oneseq=join("",@{$seqRef->{$id}});
	my @pos=ngly($oneseq);
	push(@{$glyc->{$id}},@pos);
}

##	Dataset
print "Start storing dataset ...\n";
my @antigenicity; #$train_antigenicity->{$substitution number}->[case index]=$value
my @indexes; #$train_indexes->{$substitution number}->[case index]=$pair

my $dataset=[]; #$dataset->[case index]=\@changes
#$infile="H3N2_seqCompare";
$infile=$seqCompareFile;
open(IN,"<$infile")||die "$!";
while(<IN>){
	chomp;
	next if(/^$/);
	my @lineArray = split(/\t/);
	die "Bad line formate in sequence comparison file : $_\n" if($#lineArray < 2);
	my $seqA=shift @lineArray;
	my $seqB=shift @lineArray;
	my $mark=shift @lineArray;

	my @new;  #
	foreach my$one(@lineArray){  #只考虑329以内的编号
		my $two=$one;
		$two=~s/[A-Z]//gi;
		my $temp=$one;
		$temp=~s/[0-9]//g;
		next if($temp=~/[^ACDEFGHIKLMNPQRSTVWY]/i);
		if($two <= 329){
			push(@new,$one);
		}
	}
	next if(scalar @new ==0);
	next if($mark ==0 && (scalar @new) >10);  #只考虑小于等于9个aa差异的病毒对
	next if($mark ==1 && (scalar @new) < 3);  #只考虑等于3个aa差异的病毒对

	push(@antigenicity,$mark);
	push(@indexes,$seqA."\t".$seqB );
	@{ $dataset->[$#antigenicity] }=@new;
}
close(IN)||die "$!";

my $change=0;
foreach(@antigenicity){
	$change+=$_;
}
my $prior=(1+$change)/(1+scalar(@antigenicity)-$change);  

##	Matrix
print "Start extract matrix...\n";
my $data;   #$data->[特征]->[所有病毒对在该特征上的改变]
my @attribute;

my $groups; #$groups->[group index]->[case index]=changed sites number。group即表位
$groups=Groups($dataset,\@groups,\%groups);
for(my $i=0;$i<=$#{$groups};$i++){
	push(@attribute,$groups[$i]);
	@{$data->[$#attribute]}=@{$groups->[$i]};
}


my $chars; #$chars->[char index]->[case index]=changed aa character value
$chars=Chars($dataset,\@chars,$character);
for(my $i=0;$i<=$#{$chars};$i++){
	push(@attribute,$chars[$i]);
	@{$data->[$#attribute]}=@{$chars->[$i]};
}

my @receptor;
@receptor=Receptor($dataset,\%receptor);
push(@attribute,"receptor");
@{$data->[$#attribute]}=@receptor;

my @glyc;
@glyc=Glyc(\@indexes,$glyc);
push(@attribute,"glycosylation");
@{$data->[$#attribute]}=@glyc;
print scalar @glyc,"\n";


##	Discretion
print "Start discretion...\n";
my @threshold = Discretion($data,@antigenicity);

##	Learning
print "Start parameter learning...\n";
my $distribution = NaiveBayesian($data,\@antigenicity,\@threshold);

##	Output Model
print "Start output model...\n";
#my $outfile="model";
my $outfile=$modelFile;
open(OUT,">$outfile")||die "$!";

#Output	Attribute:Group
print OUT "<Attribute:Group>\n";
foreach(@groups){
	print OUT "$_\t".join("\t",@{$group_sites->{$_}})."\n";
}
print OUT "//\n";

#Output	Attribute:Char
print OUT "<Attribute:Char>\n";
print OUT "AA\t".join("\t",@aa)."\n";
foreach my $one(@chars){
	print OUT "$one";
	foreach my $two(@aa){
		print OUT "\t$character->{$one}->{$two}";
	}
	print OUT "\n";
}
print OUT "//\n";

#Output	Attribute:Receptor
print OUT "<Attribute:Receptor>\n";
foreach(sort {$a <=> $b} keys %receptor){
	print OUT "$_\t$receptor{$_}\n";
}
print OUT "//\n";

#Output	Attribute:Glycosylation
print OUT "<Attribute:Glycosylation>\n//\n";

#Output	Parameter
print OUT "<Parameter>\n";
for(my $i=0;$i<=$#threshold;$i++){
	print OUT "$i\t$threshold[$i]\n";
}
print OUT "//\n";

#Output	Model
print OUT "<Model>\n";
print OUT "Prior\t$prior\n";
for(my $i=0;$i<=$#threshold;$i++){
	print OUT "$i\t$distribution->[$i]->[0]->[0]\t$distribution->[$i]->[0]->[1]\t$distribution->[$i]->[1]->[0]\t$distribution->[$i]->[1]->[1]\n";
}
print OUT "//\n";

close(OUT)||die "$!";
print "Done!\n";
