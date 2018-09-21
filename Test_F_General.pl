#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use Getopt::Long;
use Data::Compare;
use JSON;
use Getopt::Long;
use File::Compare;


use lib "$Bin/";
use GBuffer;
use GenBoPatient;
use Test::More;
plan tests => 38;


print "\n\n";
print "### EXECUTE GENERAL FUNCTIONAL TESTS ###\n";
print "\n\n";

my $project_name = 'GLOBAL_TESTS_TMP';
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
$project->version( 'HG19' );
$project->{year} = '2016';


my $project_dir = prepare_datas();
setPatients();

print "### LAUNCH TESTS\n";

# TESTS sur les principaux fichiers / dossiers du projet
test_getSomaticFile();
test_getPedigreeFile();
test_getAlignmentDir();
test_getVariationsDir();
test_getGvcfDir();
test_getDejaVuProjectDir();
test_getRawCallingDir();
test_getIndelsDir();
test_getLargeIndelsDir();

# TESTS sur les principaux fichiers / dossiers pour chaque patient
test_getBamFiles();
test_getVariationsFiles();
test_getIndelsFiles();
test_getLargeIndelsFiles();
test_getCoverageFile();


print "\n\n";
supress_tmp_datas($project_dir);



#### TESTS_F #####


sub test_getCoverageFile {
	my $test_name = 'test_getCoverageFile';
	foreach my $patient (@{$project->getPatients()}) {
		my $pat_name = $patient->name();
		my $obs = $patient->getCoverageFile();
		my $exp = $project_dir.'/HG19/align/coverage/'.$pat_name.'.cov.gz';
		ok(comparePaths($obs, $exp), "[$pat_name] ".$test_name);
		ok(isExistsPath($obs, 'f'), "[$pat_name] ".$test_name.' and exists');
	}
}

sub test_getLargeIndelsFiles {
	my $test_name = 'test_getLargeIndelsFiles';
	foreach my $patient (@{$project->getPatients()}) {
		my @lFiles = @{$patient->getLargeIndelsFiles()};
		my $pat_name = $patient->name();
		my $obs = $lFiles[0];
		my $exp = $project_dir.'/HG19/large_indels/lifinder/'.$pat_name.'.bed.gz';
		ok(comparePaths($obs, $exp), "[$pat_name] ".$test_name);
		ok(isExistsPath($obs, 'f'), "[$pat_name] ".$test_name.' and exists');
	}
}

sub test_getIndelsFiles {
	my $test_name = 'test_getIndelsFiles';
	foreach my $patient (@{$project->getPatients()}) {
		my @lFiles = @{$patient->getIndelsFiles()};
		my $pat_name = $patient->name();
		my $obs = $lFiles[0];
		my $exp = $project_dir.'/HG19/indels/unifiedgenotyper/'.$pat_name.'.vcf.gz';
		ok(comparePaths($obs, $exp), "[$pat_name] ".$test_name);
		ok(isExistsPath($obs, 'f'), "[$pat_name] ".$test_name.' and exists');
	}
}

sub test_getVariationsFiles {
	my $test_name = 'test_getVariationsFiles';
	foreach my $patient (@{$project->getPatients()}) {
		my @lFiles = @{$patient->getVariationsFiles()};
		my $pat_name = $patient->name();
		my $obs = $lFiles[0];
		my $exp = $project_dir.'/HG19/variations/unifiedgenotyper/'.$pat_name.'.vcf.gz';
		ok(comparePaths($obs, $exp), "[$pat_name] ".$test_name);
		ok(isExistsPath($obs, 'f'), "[$pat_name] ".$test_name.' and exists');
	}
}

sub test_getBamFiles {
	my $test_name = 'test_getBamFiles';
	foreach my $patient (@{$project->getPatients()}) {
		my @lFiles = @{$patient->getBamFiles()};
		my $pat_name = $patient->name();
		my $obs = $lFiles[0];
		my $exp = $project_dir.'/HG19/align/bwa/'.$pat_name.'.bam';
		ok(comparePaths($obs, $exp), "[$pat_name] ".$test_name);
		ok(isExistsPath($obs, 'f'), "[$pat_name] ".$test_name.' and exists');
	}
}

sub test_getLargeIndelsDir {
	my $test_name = 'test_getLargeIndelsDir';
	my $exp = $project_dir.'/HG19/large_indels/lifinder/';
	my $obs = $project->getLargeIndelsDir( 'lifinder' );
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'f'), $test_name.' and exists');
}

sub test_getIndelsDir {
	my $test_name = 'test_getIndelsDir';
	my $exp = $project_dir.'/HG19/indels/unifiedgenotyper/';
	my $obs = $project->getIndelsDir( 'unifiedgenotyper' );
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'f'), $test_name.' and exists');
}

sub test_getRawCallingDir {
	my $test_name = 'test_getRawCallingDir';
	my $exp = $project_dir.'/HG19/raw_calling/';
	my $obs = $project->getRawCallingDir();
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'd'), $test_name.' and exists');
}

sub test_getSomaticFile {
	my $test_name = 'test_getSomaticFile';
	my $exp = $project_dir.'/'.$project_name.'.somatic';
	my $obs = $project->getSomaticFile();
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'f'), $test_name.' and exists');
}

sub test_getPedigreeFile {
	my $test_name = 'test_getPedigreeFile';
	my $exp = $project_dir.'/'.$project_name.'.ped';
	my $obs = $project->getPedigreeFile();
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'f'), $test_name.' and exists');
}

sub test_getAlignmentDir {
	my $test_name = 'test_getAlignmentDir';
	my $exp = $project_dir.'/HG19/align/bwa/';
	my $obs = $project->getAlignmentDir( 'bwa' );
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'd'), $test_name.' and exists');
}

sub test_getVariationsDir {
	my $test_name = 'test_getVariationsDir';
	my $exp = $project_dir.'/HG19/variations/unifiedgenotyper/';
	my $obs = $project->getVariationsDir( 'unifiedgenotyper' );
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'd'), $test_name.' and exists');
}

sub test_getGvcfDir {
	my $test_name = 'test_getGvcfDir';
	my $exp = $project_dir.'/HG19/gvcf/haplotypecaller/';
	my $obs = $project->getGvcfDir( 'haplotypecaller' );
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'd'), $test_name.' and exists');
}

sub test_getDejaVuProjectDir {
	my $test_name = 'test_getDejaVuProjectDir';
	my $exp = $project_dir.'/HG19/deja_vu/';
	my $obs = $project->getDejaVuProjectDir();
	ok(comparePaths($obs, $exp), $test_name);
	ok(isExistsPath($obs, 'd'), $test_name.' and exists');
}



##### METHODS #####



sub isExistsPath {
	my ($path, $type) = @_;
	if ($type eq 'd') {
		return 1 if (-d $path);
	}
	elsif ($type eq 'f') {
		return 1 if (-e $path);
	}
	return;
}

sub comparePaths {
	my ($path1, $path2) = @_;
	my $c = Data::Compare->new(esthetic_correct_path($path1), esthetic_correct_path($path2));
	return $c->Cmp;
}

sub setPatients {
	print "### CREATE NEW FICTIFS PATIENTS IN PROJECT\n";
	my %names;
	my $id = 0;
	my @lNames;
	push(@lNames, 'PATIENT_A');
	push(@lNames, 'PATIENT_B');
	foreach my $pat_name (@lNames) {
		$id++;
		my $h = {
	            'status' => '1',
	            'flowcell' => 'A',
	            'genbo_id' => '0',
	            'patient_id' => $id,
	            'capture_id' => '100',
	            'project_id' => '1000',
	            'name' => $pat_name,
	            'origin' => $pat_name,
	            'mother' => '',
	            'sex' => '1',
	            'description' => '',
	            'run_id' => '891',
	            'bar_code' => 'B03',
	            'creation_date' => '1986-05-26 00:15:00',
	            'father' => '',
	            'project_id_dest' => '0',
	            'family' => 'FAMILY'
		};
		$h->{id} = $h->{patient_id};
		$names{$h->{id}} = undef;
		$h->{project} = $project;
		my $obj = new GenBoPatient($h);
		$obj->{alignmentMethods} = ['bwa'];
		$obj->{callingMethods} = ['unifiedgenotyper', 'lifinder'];
		$project->{objects}->{patients}->{$h->{id}} = $obj;
		print "# CREATED $pat_name: ".ref($obj)." object\n";
	}
	$project->{patients_object} = \%names;
	print "\n\n"
}

sub touch_file {
	my $file = shift;
	eval {
		`touch $file`;
	};
	if ($@) {
		warn "\n\nERROR: can't create this file: $file. Die.\n\n";
		warn Dumper $@;
		warn "\n\n";
		die;
	}
	unless (-e $file) {
		warn "\n\nERROR: $file not found... Die !\n\n";
		die;
	}
}

sub create_path {
	my $path = shift;
	eval {
		mkdir $path;
	};
	if ($@) {
		warn "\n\nERROR: can't create this path: $path. Die.\n\n";
		warn Dumper $@;
		warn "\n\n";
		die;
	}
	unless (-d $path) {
		warn "\n\nERROR: $path not found... Die !\n\n";
		die;
	}
}

sub esthetic_correct_path {
	my $path = shift;
	$path =~ s/\/\//\//g;
	return $path;
}

sub supress_tmp_datas {
	my $project_dir = shift;
	print "### SUPRESS DATAS\n";
	unless (-d $project_dir) {
		warn "\n\nERROR: $project_dir not found... Very strnage ???? Die.\n\n";
		die;
	}
	print '# SUPRESS PROJECT DIR: '.$project_dir."\n";
	`rm -r $project_dir`;
	if (-d $project_dir) {
		warn "\n\nERROR: can't supress $project_dir... Very strnage ???? Die.\n\n";
		die;
	}
	print "\n\n";
}

sub prepare_datas {
	print "### PREPARE DATAS WITH FICTIFS DIR / FILES\n";
	my $root_dir = esthetic_correct_path( $buffer->config->{project_data}->{root}.'/' );
	print '# CHECK ROOT DIR: '.$root_dir."\n";
	unless (-d $root_dir) {
		warn "\n\nERROR: $root_dir not found... Die !\n\n";
		die;
	}
	
	my $ngs_dir = esthetic_correct_path( $root_dir.'/'.$buffer->config->{project_data}->{ngs}.'/' );
	print '# CHECK NGS DIR: '.$ngs_dir."\n";
	unless (-d $ngs_dir) {
		warn "\n\nERROR: $ngs_dir not found... Die !\n\n";
		die;
	}
	
	my $project_dir = esthetic_correct_path( $ngs_dir.'/'.$project_name.'/' );
	print '# CREATE PROJECT DIR: '.$project_dir."\n";
	create_path($project_dir);
	
	my $ref_dir = esthetic_correct_path( $project_dir.'/'.$project->getVersion().'/' );
	print '# CREATE '.$project->getVersion().' DIR: '.$ref_dir."\n";
	create_path($ref_dir);
	
	my $seq_dir = esthetic_correct_path( $project_dir.'/sequences/' );
	print '# CREATE SEQUENCES DIR: '.$seq_dir."\n";
	create_path($seq_dir);
	
	my $seq_machine_dir = esthetic_correct_path( $seq_dir.'/miseq/' );
	print '# CREATE SEQUENCES MACHINE DIR: '.$seq_machine_dir."\n";
	create_path($seq_machine_dir);
	
	my $dejavu_dir = esthetic_correct_path( $ref_dir.'/deja_vu/' );
	print '# CREATE DEJA_VU DIR: '.$dejavu_dir."\n";
	create_path($dejavu_dir);
	
	my $align_dir = esthetic_correct_path( $ref_dir.'/align/' );
	print '# CREATE ALIGN DIR: '.$align_dir."\n";
	create_path($align_dir);
	
	my $align_bwa_dir = esthetic_correct_path( $align_dir.'/bwa/' );
	print '# CREATE ALIGN BWA DIR: '.$align_bwa_dir."\n";
	create_path($align_bwa_dir);
	
	my $align_cov_dir = esthetic_correct_path( $align_dir.'/coverage/' );
	print '# CREATE ALIGN COVERAGE DIR: '.$align_cov_dir."\n";
	create_path($align_cov_dir);
	
	my $var_dir = esthetic_correct_path( $ref_dir.'/variations/' );
	print '# CREATE VARIATIONS DIR: '.$var_dir."\n";
	create_path($var_dir);
	
	my $var_uni_dir = esthetic_correct_path( $var_dir.'/unifiedgenotyper/' );
	print '# CREATE VARIATIONS UNIFIEDGENOTYPER DIR: '.$var_uni_dir."\n";
	create_path($var_uni_dir);
	
	my $indels_dir = esthetic_correct_path( $ref_dir.'/indels/' );
	print '# CREATE INDELS DIR: '.$indels_dir."\n";
	create_path($indels_dir);
	
	my $indels_uni_dir = esthetic_correct_path( $indels_dir.'/unifiedgenotyper/' );
	print '# CREATE INDELS UNIFIEDGENOTYPER DIR: '.$indels_uni_dir."\n";
	create_path($indels_uni_dir);
	
	my $large_indels_dir = esthetic_correct_path( $ref_dir.'/large_indels/' );
	print '# CREATE LARGE_INDELS DIR: '.$large_indels_dir."\n";
	create_path($large_indels_dir);
	
	my $large_indels_lifinder_dir = esthetic_correct_path( $large_indels_dir.'/lifinder/' );
	print '# CREATE LARGE_INDELS UNIFIEDGENOTYPER DIR: '.$large_indels_lifinder_dir."\n";
	create_path($large_indels_lifinder_dir);
	
	my $ped_file = esthetic_correct_path( $project_dir.'/'.$project_name.'.ped' );
	print '# ADD PED FILE: '.$ped_file."\n";
	touch_file($ped_file);
	
	my $somatic_file = esthetic_correct_path( $project_dir.'/'.$project_name.'.somatic' );
	print '# ADD SOMATIC FILE: '.$somatic_file."\n";
	touch_file($somatic_file);
	
	my $pat_A_fastq_1 = esthetic_correct_path( $seq_machine_dir.'/PATIENT_A_S30_L001_R1_001.fastq.gz' );
	print '# ADD PATIENT A FASTQ FILE 1: '.$pat_A_fastq_1."\n";
	touch_file($pat_A_fastq_1);
	
	my $pat_A_fastq_2 = esthetic_correct_path( $seq_machine_dir.'/PATIENT_A_S30_L001_R2_001.fastq.gz' );
	print '# ADD PATIENT A FASTQ FILE 2: '.$pat_A_fastq_2."\n";
	touch_file($pat_A_fastq_2);
	
	my $pat_B_fastq_1 = esthetic_correct_path( $seq_machine_dir.'/PATIENT_B_S32_L001_R1_001.fastq.gz' );
	print '# ADD PATIENT B FASTQ FILE 1: '.$pat_B_fastq_1."\n";
	touch_file($pat_B_fastq_1);
	
	my $pat_B_fastq_2 = esthetic_correct_path( $seq_machine_dir.'/PATIENT_B_S32_L001_R2_001.fastq.gz' );
	print '# ADD PATIENT B FASTQ FILE 2: '.$pat_B_fastq_2."\n";
	touch_file($pat_B_fastq_2);
	
	my $pat_A_bam = esthetic_correct_path( $align_bwa_dir.'/PATIENT_A.bam' );
	print '# ADD PATIENT A BAM FILE: '.$pat_A_bam."\n";
	touch_file($pat_A_bam);
	
	my $pat_A_bai = esthetic_correct_path( $align_bwa_dir.'/PATIENT_A.bam.bai' );
	print '# ADD PATIENT A BAM.BAI FILE: '.$pat_A_bai."\n";
	touch_file($pat_A_bai);
	
	my $pat_B_bam = esthetic_correct_path( $align_bwa_dir.'/PATIENT_B.bam' );
	print '# ADD PATIENT B BAM FILE: '.$pat_B_bam."\n";
	touch_file($pat_B_bam);
	
	my $pat_B_bai = esthetic_correct_path( $align_bwa_dir.'/PATIENT_B.bam.bai' );
	print '# ADD PATIENT B BAM.BAI FILE: '.$pat_B_bai."\n";
	touch_file($pat_B_bai);
	
	my $pat_A_cov = esthetic_correct_path( $align_cov_dir.'/PATIENT_A.cov.gz' );
	print '# ADD PATIENT A COV FILE: '.$pat_A_cov."\n";
	touch_file($pat_A_cov);
	
	my $pat_A_cov_tbi = esthetic_correct_path( $align_cov_dir.'/PATIENT_A.cov.tbi' );
	print '# ADD PATIENT A COV.TBI FILE: '.$pat_A_cov_tbi."\n";
	touch_file($pat_A_cov_tbi);
	
	my $pat_B_cov = esthetic_correct_path( $align_cov_dir.'/PATIENT_B.cov.gz' );
	print '# ADD PATIENT B COV FILE: '.$pat_B_cov."\n";
	touch_file($pat_B_cov);
	
	my $pat_B_cov_tbi = esthetic_correct_path( $align_cov_dir.'/PATIENT_B.cov.tbi' );
	print '# ADD PATIENT B COV.TBI FILE: '.$pat_B_cov_tbi."\n";
	touch_file($pat_B_cov_tbi);
	
	my $pat_A_var_vcf = esthetic_correct_path( $var_uni_dir.'/PATIENT_A.vcf.gz' );
	print '# ADD PATIENT A VARIATIONS VCF FILE: '.$pat_A_var_vcf."\n";
	touch_file($pat_A_var_vcf);
	
	my $pat_A_var_vcf_tbi = esthetic_correct_path( $var_uni_dir.'/PATIENT_A.vcf.gz.tbi' );
	print '# ADD PATIENT A VARIATIONS VCF TBI FILE: '.$pat_A_var_vcf_tbi."\n";
	touch_file($pat_A_var_vcf_tbi);
	
	my $pat_B_var_vcf = esthetic_correct_path( $var_uni_dir.'/PATIENT_B.vcf.gz' );
	print '# ADD PATIENT B VARIATIONS VCF FILE: '.$pat_B_var_vcf."\n";
	touch_file($pat_B_var_vcf);
	
	my $pat_B_var_vcf_tbi = esthetic_correct_path( $var_uni_dir.'/PATIENT_B.vcf.gz.tbi' );
	print '# ADD PATIENT B VARIATIONS VCF TBI FILE: '.$pat_B_var_vcf_tbi."\n";
	touch_file($pat_B_var_vcf_tbi);
	
	my $pat_A_indels_vcf = esthetic_correct_path( $indels_uni_dir.'/PATIENT_A.vcf.gz' );
	print '# ADD PATIENT A INDELS VCF FILE: '.$pat_A_indels_vcf."\n";
	touch_file($pat_A_indels_vcf);
	
	my $pat_A_indels_vcf_tbi = esthetic_correct_path( $indels_uni_dir.'/PATIENT_A.vcf.gz.tbi' );
	print '# ADD PATIENT A INDELS VCF TBI FILE: '.$pat_A_indels_vcf_tbi."\n";
	touch_file($pat_A_indels_vcf_tbi);
	
	my $pat_B_indels_vcf = esthetic_correct_path( $indels_uni_dir.'/PATIENT_B.vcf.gz' );
	print '# ADD PATIENT B INDELS VCF FILE: '.$pat_B_indels_vcf."\n";
	touch_file($pat_B_indels_vcf);
	
	my $pat_B_indels_vcf_tbi = esthetic_correct_path( $indels_uni_dir.'/PATIENT_B.vcf.gz.tbi' );
	print '# ADD PATIENT B INDELS VCF TBI FILE: '.$pat_B_indels_vcf_tbi."\n";
	touch_file($pat_B_indels_vcf_tbi);
	
	my $pat_A_large_indels_vcf = esthetic_correct_path( $large_indels_lifinder_dir.'/PATIENT_A.bed.gz' );
	print '# ADD PATIENT A LARGE_INDELS VCF FILE: '.$pat_A_large_indels_vcf."\n";
	touch_file($pat_A_large_indels_vcf);
	
	my $pat_A_large_indels_vcf_tbi = esthetic_correct_path( $large_indels_lifinder_dir.'/PATIENT_A.bed.gz.tbi' );
	print '# ADD PATIENT A LARGE_INDELS VCF TBI FILE: '.$pat_A_large_indels_vcf_tbi."\n";
	touch_file($pat_A_large_indels_vcf_tbi);
	
	my $pat_B_large_indels_vcf = esthetic_correct_path( $large_indels_lifinder_dir.'/PATIENT_B.bed.gz' );
	print '# ADD PATIENT B LARGE_INDELS VCF FILE: '.$pat_B_large_indels_vcf."\n";
	touch_file($pat_B_large_indels_vcf);
	
	my $pat_B_large_indels_vcf_tbi = esthetic_correct_path( $large_indels_lifinder_dir.'/PATIENT_B.bed.gz.tbi' );
	print '# ADD PATIENT B LARGE_INDELS VCF TBI FILE: '.$pat_B_large_indels_vcf_tbi."\n";
	touch_file($pat_B_large_indels_vcf_tbi);
	
	print "\n\n";
	return $project_dir;
}