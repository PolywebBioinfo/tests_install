#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Getopt::Long;
use Data::Compare;
use JSON;
use Getopt::Long;
use File::Compare;
use DateTime;
use GBuffer;

my ($update_data, $with_db, $single_chr, $multi_chr, $create_json, $my_project_name, $no_test, $not_exact_test, $only_test, $print_cmd, $print_url);
GetOptions(
	'single_chr!' => \$single_chr,
	'multi_chr!' => \$multi_chr,
	'update_data!' => \$update_data,
	'with_db!' => \$with_db,
	'print_cmd!' => \$print_cmd,
	'print_url!' => \$print_url,
	'create_json=s' => \$create_json,
	'project_name=s' => \$my_project_name,
	'no_test!' => \$no_test,
	'not_exact_test!' => \$not_exact_test,
	'only_test=s' => \$only_test,
);

if ($not_exact_test) {
	warn "\n\n### TYPE TEST: compare values and common genes only ###\n";
}

if ($create_json) {
	unless (-d $create_json) {
		die("\n\nERROR: with -create_json option,you need defined a valid directory to write JSON results. Die...\n\n");
	}
}

use lib "$Bin/";
my $project_name = 'TESTS_F';
$project_name = $my_project_name if ($my_project_name);
my $cmd_interface = "$Bin/../../../../../polymorphism-cgi/json_output_nodb/polyquery.pl ";
if ($with_db) {
	$cmd_interface.= " project=$project_name";
	$cmd_interface.= " test_with_db=1" unless ($my_project_name);
}
else {
	$cmd_interface.= " project=$project_name";
	$cmd_interface.= " test=1" unless ($my_project_name);
}

#warn $cmd_interface; die;

use Test::More;
plan tests => 242;

# pour monitorer le script: perl -d:NYTProf test.pl puis nytprofhtml.
# Ex: perl -d:NYTProf /var/www/cgi-bin/mbras/polymorphism-cgi/json_output_nodb/interface_json.pl project=$project_name stat=all filter_type_variation=silent+coding+intergenic+intronic+ncrna+maturemirna+pseudogene+frameshift+ mode=ind level_ind=variation level_fam=variation
# Ex: export PATH=$PATH:/bip-d/activePerl/site/bin/
# Ex: /bip-d/activePerl/site/bin/nytprofhtml -o nytprof -f nytprof.out



my @lFiltersChr = ('1', '1,2,3,Y');
@lFiltersChr = ('1') if ($single_chr);
@lFiltersChr = ('1,2,3,Y') if ($multi_chr);
foreach my $filter_chromosome (@lFiltersChr) {
	print "\n\n";
	print "### EXECUTE POLYQUERY FUNCTIONAL TESTS ";
	if ($with_db) { print "[USING DATABASE] [CHR $filter_chromosome] ### "; }
	else { print "[NOT USING DATABASE] [CHR $filter_chromosome] ###"; }
	print "\n\n";
	
#	test_ind_compound($filter_chromosome);
#	test_fam_compound($filter_chromosome);
#	test_ind_compound_by_gene($filter_chromosome);
#	test_model_recessif_compound($filter_chromosome);
#	die;
	
	
	#test_compare_store_ids($my_project_name);
	test_default($filter_chromosome);
	test_only_intergenic($filter_chromosome);
	test_all_databases($filter_chromosome);
	test_ind_utr_splicing_esssplicing_stop_startstop_noframeshift($filter_chromosome);
#	test_all_large_deletion($filter_chromosome);
	test_ind_in_the_attic($filter_chromosome); 
	test_ind_exclude($filter_chromosome);
	test_ind_intersect($filter_chromosome);
	test_ind_filter_he($filter_chromosome);
	test_ind_filter_ho($filter_chromosome);
	test_ind_in_the_attic_exclude_intersect_he_ho($filter_chromosome);
	test_ind_in_the_attic_by_gene($filter_chromosome);
	test_ind_exclude_by_gene($filter_chromosome);
	test_ind_intersect_by_gene($filter_chromosome);
	test_ind_filter_he_by_gene($filter_chromosome);
	test_ind_filter_ho_by_gene($filter_chromosome);
	test_ind_in_the_attic_exclude_intersect_he_ho_by_gene($filter_chromosome);
	test_patients_fam_exclude($filter_chromosome);
	test_patients_fam_intersect($filter_chromosome);
	test_fam_intersect($filter_chromosome);
	test_fam_exclude($filter_chromosome);
	test_fam_intersect_familial($filter_chromosome);
	test_fam_exclude_familial($filter_chromosome);
	test_fam_intersect_familial_by_gene($filter_chromosome);
	test_fam_exclude_familial_by_gene($filter_chromosome);
	test_fam_filter_he($filter_chromosome);
	test_fam_filter_ho($filter_chromosome);
	test_ind_recessif($filter_chromosome);
	test_ind_compound($filter_chromosome);
	test_fam_dominant($filter_chromosome);
	test_fam_denovo($filter_chromosome);
	test_fam_recessif($filter_chromosome);
	test_fam_compound($filter_chromosome);
#	test_ind_disease($filter_chromosome);
	test_ind_atleast($filter_chromosome);
	test_fam_atleast($filter_chromosome);
	test_ind_atleast_byGenes($filter_chromosome);
	test_fam_atleast_byGenes($filter_chromosome);
	test_filter_region_and_stats_all($filter_chromosome);
	test_filter_region_and_stats_region($filter_chromosome);
	test_filter_region_exclude_stats_all($filter_chromosome);
	test_filter_region_exclude_stats_region($filter_chromosome);
	test_filter_default_text_search($filter_chromosome);
	test_ind_recessif_by_gene($filter_chromosome);
	test_ind_compound_by_gene($filter_chromosome);
	test_fam_recessif_by_gene($filter_chromosome);
	test_fam_dominant_by_gene($filter_chromosome);
	test_fam_denovo_by_gene($filter_chromosome);
	test_fam_strict_denovo($filter_chromosome);
	test_dejavu_nb15($filter_chromosome);
	test_dejavu_nb1($filter_chromosome);
	test_dejavu_nb0($filter_chromosome);
	test_dejavu_ho_nb15($filter_chromosome);
	test_dejavu_ho_nb1($filter_chromosome);
	test_dejavu_ho_nb0($filter_chromosome);
	###test_export_xls_gene($filter_chromosome);
	###test_export_xls_variants($filter_chromosome);
	test_ind_var_var_nomodel_intersect($filter_chromosome);
	test_ind_var_var_recessif_intersect($filter_chromosome);
	test_ind_var_var_nomodel_exclude($filter_chromosome);
	test_ind_var_var_recessif_exclude($filter_chromosome);
	test_ind_var_var_nomodel_intheattic($filter_chromosome);
	test_ind_var_var_recessif_intheattic($filter_chromosome);
	test_ind_var_var_nomodel_he($filter_chromosome);
	test_ind_var_var_recessif_he($filter_chromosome);
	test_ind_var_var_nomodel_ho($filter_chromosome);
	test_ind_var_var_recessif_ho($filter_chromosome);
	test_ind_gene_var_nomodel_intersect($filter_chromosome);
	test_ind_gene_var_recessif_intersect($filter_chromosome); # OK JE PENSE
	test_ind_gene_var_nomodel_exclude($filter_chromosome);
	test_ind_gene_var_recessif_exclude($filter_chromosome);
	test_ind_gene_var_nomodel_intheattic($filter_chromosome);
	test_ind_gene_var_recessif_intheattic($filter_chromosome);
	test_fam_var_var_nomodel_intersect($filter_chromosome);
	test_fam_var_var_recessif_intersect($filter_chromosome);
	test_fam_var_var_nomodel_exclude($filter_chromosome);
	test_fam_var_var_recessif_exclude($filter_chromosome);
	test_fam_var_var_nomodel_intheattic($filter_chromosome);
	test_fam_var_var_recessif_intheattic($filter_chromosome);
	test_fam_var_var_nomodel_he($filter_chromosome);
	test_fam_var_var_recessif_he($filter_chromosome);
	test_fam_var_var_nomodel_ho($filter_chromosome);
	test_fam_var_var_recessif_ho($filter_chromosome);
	test_fam_gene_var_nomodel_intersect($filter_chromosome);
	test_fam_gene_var_recessif_intersect($filter_chromosome);#to check
	test_fam_gene_var_nomodel_exclude($filter_chromosome);
	test_fam_gene_var_recessif_exclude($filter_chromosome);
	test_fam_gene_var_nomodel_intheattic($filter_chromosome);
	test_fam_gene_var_recessif_intheattic($filter_chromosome);
	test_fam_var_var_nomodel_intersectfam($filter_chromosome);
	test_fam_var_var_recessif_intersectfam($filter_chromosome);
	test_fam_var_gene_nomodel_intersectfam($filter_chromosome);
	test_fam_var_gene_recessif_intersectfam($filter_chromosome);
	test_model_recessif_compound($filter_chromosome);
	test_divers_1($filter_chromosome);
	test_divers_2($filter_chromosome);
	test_divers_3($filter_chromosome);
	test_divers_4($filter_chromosome);
	test_divers_5($filter_chromosome);
	test_divers_6($filter_chromosome);
	test_divers_7($filter_chromosome);
	test_divers_8($filter_chromosome);
	test_divers_9($filter_chromosome);
	test_divers_10($filter_chromosome);
	test_divers_11($filter_chromosome);
	test_divers_12($filter_chromosome);
	test_divers_13($filter_chromosome);
	test_divers_14($filter_chromosome);
	test_divers_15($filter_chromosome);
	test_divers_16($filter_chromosome);
	test_divers_17($filter_chromosome);
	test_divers_18($filter_chromosome);
	test_divers_19($filter_chromosome);
#	test_divers_20($filter_chromosome);
	test_only_genes($filter_chromosome);
	test_only_genes_recessif($filter_chromosome);
	test_fam_mosaique_1($filter_chromosome);
	test_fam_mosaique_2($filter_chromosome);
	test_ind_region_ho_ind_25($filter_chromosome);
	test_ind_region_ho_ind_50($filter_chromosome);
	test_ind_region_ho_ind_75($filter_chromosome);
	test_ind_region_ho_ind_100($filter_chromosome);
	test_ind_region_rec_fam_25($filter_chromosome);
	test_ind_region_rec_fam_50($filter_chromosome);
	test_ind_region_rec_fam_75($filter_chromosome);
	test_ind_region_rec_fam_100($filter_chromosome);
}



##### METHODS #####



sub test_ind_region_rec_fam_25 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_rec_fam_25';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=25";
	$cmd .= " mode=fam";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_rec_fam_50 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_rec_fam_50';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=50";
	$cmd .= " mode=fam";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_rec_fam_75 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_rec_fam_75';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=75";
	$cmd .= " mode=fam";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_rec_fam_100 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_rec_fam_100';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=100";
	$cmd .= " mode=fam";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_ho_ind_25 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_ho_ind_25';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=25";
	$cmd .= " mode=ind";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_ho_ind_50 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_ho_ind_50';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=50";
	$cmd .= " mode=ind";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_ho_ind_75 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_ho_ind_75';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=75";
	$cmd .= " mode=ind";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_region_ho_ind_100 {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_region_ho_ind_100';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=+";
	$cmd .= " filter_nbvar_regionho=25";
	$cmd .= " mode=ind";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_mosaique_1 {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_mosaique_1';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=mosaic";
	$cmd .= " mode=fam";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_mosaique_2 {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_mosaique_2';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=intergenic+";
	$cmd .= " model=mosaic";
	$cmd .= " mode=fam";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_only_genes_recessif {
	my $filter_chromosome = shift;
	my $test_name = 'test_only_genes_recessif';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " view_all_genes=1";
	$cmd .= " only_genes=RNASEL,FMO3,PCSK9";
	$cmd .= " view_all_genes=1";	
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_only_genes {
	my $filter_chromosome = shift;
	my $test_name = 'test_only_genes';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " view_all_genes=1";
	$cmd .= " only_genes=RHCE,MTF1,SDC3,CASZ1,ATP13A2";
	$cmd .= " view_all_genes=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_model_recessif_compound {
	my $filter_chromosome = shift;
	my $test_name = 'test_model_recessif_compound';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif_compound";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}


sub test_fam_var_var_nomodel_intersectfam {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_nomodel_intersectfam';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_recessif_intersectfam {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_recessif_intersectfam';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_gene_nomodel_intersectfam {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_gene_nomodel_intersectfam';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=gene";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_gene_recessif_intersectfam {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_gene_recessif_intersectfam';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=gene";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_gene_var_nomodel_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_gene_var_nomodel_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_gene_var_recessif_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_gene_var_recessif_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_gene_var_nomodel_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_gene_var_nomodel_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_gene_var_recessif_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_gene_var_recessif_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_gene_var_nomodel_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_gene_var_nomodel_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_gene_var_recessif_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_gene_var_recessif_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_nomodel_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_nomodel_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_recessif_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_recessif_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_nomodel_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_nomodel_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_recessif_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_recessif_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_nomodel_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_nomodel_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_recessif_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_recessif_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_nomodel_he {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_nomodel_he';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_recessif_he {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_recessif_he';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_nomodel_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_nomodel_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_var_var_recessif_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_var_var_recessif_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_gene_var_nomodel_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_gene_var_nomodel_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_gene_var_recessif_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_gene_var_recessif_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_gene_var_nomodel_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_gene_var_nomodel_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_gene_var_recessif_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_gene_var_recessif_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_gene_var_nomodel_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_gene_var_nomodel_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_gene_var_recessif_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_gene_var_recessif_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_nomodel_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_nomodel_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_recessif_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_recessif_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_nomodel_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_nomodel_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_recessif_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_recessif_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_nomodel_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_nomodel_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_recessif_intheattic {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_recessif_intheattic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_nomodel_he {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_nomodel_he';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_recessif_he {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_recessif_he';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_nomodel_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_nomodel_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_var_var_recessif_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_var_var_recessif_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_1 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_1';
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	if ($only_test) { return unless ($only_test eq $test_name); }
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+splicing+coding+intergenic+intronic+ncrna+pseudogene+non-frameshift+freq_001+freq_01+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=10";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_10'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_10'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_2 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_2';
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	if ($only_test) { return unless ($only_test eq $test_name); }
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=large_deletion+deletion+insertion+silent+utr+splicing+coding+intergenic+intronic+ncrna+pseudogene+non-frameshift+freq_01+freq_05+freq_1+prediction_0+prediction_1";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=20";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_20'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_20'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_3 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_3';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=intergenic+freq_1+prediction_0+prediction_1+prediction_2";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=1";
	$cmd .= " dejavu_ho=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_1ho'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_1ho'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_4 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_4';
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=intergenic+freq_1+prediction_0+prediction_1+prediction_2";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_1'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_1'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_5 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_5';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=large_deletion+deletion+insertion+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=20";
	$cmd .= " dejavu_ho=1";
	$cmd .= " filter_patient=LEF_MEL+";
	$cmd .= " filter_not_patient=GUI_REF+";
	$cmd .= " filter_attic=BOU_CHR+GUI_NAT+GUI_YAE+GUI_YOS+LEF_JER+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_20ho'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_20ho'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_6 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_6';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=large_deletion+deletion+insertion+";
	$cmd .= " model=recessif";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=20";
	$cmd .= " filter_attic=GUI_NAT+GUI_REF+GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_20'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_20'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_7 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_7';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=large_deletion+silent+utr+splicing+intergenic+intronic+ncrna+pseudogene+freq_05+freq_1+prediction_0";
	$cmd .= " model=denovo";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=10";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_10'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_10'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_8 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_8';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=large_deletion+intergenic+";
	$cmd .= " model=recessif";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=20";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_20'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_20'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_9 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_9';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_10 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_10';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_not=GUI+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_11 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_11';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_NAT+GUI_REF+GUI_YAE+GUI_YOS+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_12 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_12';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_13 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_13';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_14 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_14';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+splicing+essential_splicing+coding+intergenic+intronic+ncrna+maturemirna+pseudogene+frameshift+non-frameshift+prediction_0";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_REF+LEF_MEL+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_15 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_15';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+splicing+essential_splicing+coding+intergenic+intronic+ncrna+maturemirna+pseudogene+frameshift+non-frameshift+prediction_0";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_REF+LEF_MEL+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_16 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_16';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+splicing+essential_splicing+coding+intergenic+intronic+ncrna+maturemirna+pseudogene+frameshift+non-frameshift+prediction_0";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_17 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_17';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+splicing+essential_splicing+coding+intergenic+intronic+ncrna+maturemirna+pseudogene+frameshift+non-frameshift+prediction_0";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=gene";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_18 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_18';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_REF+LEF_MEL+";
	$cmd .= " filter_ho=BOU_CHR+GUI_NAT+GUI_YAE+GUI_YOS+LEF_JER+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_19 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_19';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=prediction_0+prediction_1";
	$cmd .= " filter_region=1:100:10000000:0+2:100:10000000:0+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=30";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_30'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_30'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_divers_20 {
	my $filter_chromosome = shift;
	my $test_name = 'test_divers_20';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=freq_05+freq_1+prediction_0";
	$cmd .= " filter_region=1:100:10000000:-1+2:100:10000000:-1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=30";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.'_with_dejavu_30'.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.'_with_dejavu_30'.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_export_xls_variants {
	my $filter_chromosome = shift;
	my $test_name = 'test_export_xls_variants';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " type=genes";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " resume_vector=2:all:13551:28,53,182,295-297,322,422,506,533,833,837,888,944,945,967,969,979,980,982,991,995,1001,1007,1015,1044,1059,1068,1088,1118,1239,1257,1273,1277,1281,1296,1297,1362,1390,1395,1418,1465,1469,1480,1495,1516,1575,1656,1895,2113,2237,2365,2405,2471,2505,2535,2553,2672,2677,2686,2687,2715,2756,2776,2777,2814,2949,2980,3058,3094,3150,3353,3415,3592,3644,3657,3813,3863,3864,3874,3925,3926,3965,3975,4160,4164,4205,4228,4235,4251,4252,4279,4327,4427,4472-4475,4501,4505,4537,4572,4590,4732,4735,4737,4742,4744,4745,4747,4748,4750,4754,4757,4759,4761-4765,4768,4795,4797,4867,4868,4890-4893,4940,4956,5016,5064-5066,5105,5116,5118,5119,5126,5127,5138,5147-5151,5153,5154,5247,5291,5297,5329,5333,5349,5359,5482,5516,5625,5683,5816,5818,5862,5887,5959,5969,5980,5982-5984,5998,6051,6154,6160,6321,6443,6593,6621,6648,6662,6679,6681,6691,6692,6706,6743,6754,6759,6765,6773,6775,6783,6789,6793,6794,6809,6826-6828,6833,6874,6875,6950,6954,6955,7040,7058,7059,7115,7178,7199,7278,7298,7328,7446,7473,7528,7647,7649,7650,7653,7679,7736,7765,7816,7833,7862-7864,7913,7923,7929,7935,7974,8031,8040,8097,8107,8134,8160,8167,8227,8256,8267,8715,8716,8757,8796,8838,8839,8920,8921,8942,8943,8950,8952,8959,8970,8999,9022,9027,9028,9035,9037,9042,9057,9220,9222,9265,9282,9353,9356,9393,9407,9515,9527,9536,9573,9684,9722,9802,10010,10041,10070,10108,10142,10149,10174,10213,10256,10356,10520,10522,10602,10703,10748,10851,10860,10943,10990,10991,11111,11157,11179,11184,11185,11218,11349,11371,11426,11471,11540,11575,11717,11741,11806,11841-11843,11860,11899,11934,11982,12037,12049,12057,12061,12064,12231,12312,12477,12566,12635,12645,12709,12724,12770,12787,12809,12810,12889,12951,12979,12981,12983,12985,12997,13075,13087,13137,13138,13256,13272,13273,13349,13359,13409,13443,13491,13495";
	$cmd .= " xls_by_variants=1	";
	if ($with_db) { $test_name = '[WITH DB] '.$test_name; }
	else { $test_name = '[NO DB] '.$test_name; }
	ok(launchCmd($cmd, $test_name), $test_name);
}

sub test_export_xls_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_export_xls_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " type=genes";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " resume_vector=2:all:13551:28,53,182,295-297,322,422,506,533,833,837,888,944,945,967,969,979,980,982,991,995,1001,1007,1015,1044,1059,1068,1088,1118,1239,1257,1273,1277,1281,1296,1297,1362,1390,1395,1418,1465,1469,1480,1495,1516,1575,1656,1895,2113,2237,2365,2405,2471,2505,2535,2553,2672,2677,2686,2687,2715,2756,2776,2777,2814,2949,2980,3058,3094,3150,3353,3415,3592,3644,3657,3813,3863,3864,3874,3925,3926,3965,3975,4160,4164,4205,4228,4235,4251,4252,4279,4327,4427,4472-4475,4501,4505,4537,4572,4590,4732,4735,4737,4742,4744,4745,4747,4748,4750,4754,4757,4759,4761-4765,4768,4795,4797,4867,4868,4890-4893,4940,4956,5016,5064-5066,5105,5116,5118,5119,5126,5127,5138,5147-5151,5153,5154,5247,5291,5297,5329,5333,5349,5359,5482,5516,5625,5683,5816,5818,5862,5887,5959,5969,5980,5982-5984,5998,6051,6154,6160,6321,6443,6593,6621,6648,6662,6679,6681,6691,6692,6706,6743,6754,6759,6765,6773,6775,6783,6789,6793,6794,6809,6826-6828,6833,6874,6875,6950,6954,6955,7040,7058,7059,7115,7178,7199,7278,7298,7328,7446,7473,7528,7647,7649,7650,7653,7679,7736,7765,7816,7833,7862-7864,7913,7923,7929,7935,7974,8031,8040,8097,8107,8134,8160,8167,8227,8256,8267,8715,8716,8757,8796,8838,8839,8920,8921,8942,8943,8950,8952,8959,8970,8999,9022,9027,9028,9035,9037,9042,9057,9220,9222,9265,9282,9353,9356,9393,9407,9515,9527,9536,9573,9684,9722,9802,10010,10041,10070,10108,10142,10149,10174,10213,10256,10356,10520,10522,10602,10703,10748,10851,10860,10943,10990,10991,11111,11157,11179,11184,11185,11218,11349,11371,11426,11471,11540,11575,11717,11741,11806,11841-11843,11860,11899,11934,11982,12037,12049,12057,12061,12064,12231,12312,12477,12566,12635,12645,12709,12724,12770,12787,12809,12810,12889,12951,12979,12981,12983,12985,12997,13075,13087,13137,13138,13256,13272,13273,13349,13359,13409,13443,13491,13495";
	$cmd .= " xls_by_genes=1";
	if ($with_db) { $test_name = '[WITH DB] '.$test_name; }
	else { $test_name = '[NO DB] '.$test_name; }
	ok(launchCmd($cmd, $test_name), $test_name);
}

sub test_dejavu_ho_nb0 {
	my $filter_chromosome = shift;
	my $test_name = 'test_dejavu_ho_nb0';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_none+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=0";
	$cmd .= " dejavu_ho=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_dejavu_ho_nb1 {
	my $filter_chromosome = shift;
	my $test_name = 'test_dejavu_ho_nb1';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_none+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=1";
	$cmd .= " dejavu_ho=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_dejavu_ho_nb15 {
	my $filter_chromosome = shift;
	my $test_name = 'test_dejavu_ho_nb15';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_none+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=15";
	$cmd .= " dejavu_ho=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_dejavu_nb0 {
	my $filter_chromosome = shift;
	my $test_name = 'test_dejavu_nb0';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_none+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=0";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_dejavu_nb1 {
	my $filter_chromosome = shift;
	my $test_name = 'test_dejavu_nb1';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_none+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=1";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_dejavu_nb15 {
	my $filter_chromosome = shift;
	my $test_name = 'test_dejavu_nb15';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_none+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " dejavu=15";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_all_large_deletion {
	my $filter_chromosome = shift;
	my $test_name = 'test_all_large_deletion';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=deletion+insertion+substitution";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_compound {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_compound';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " model=compound";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_denovo_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_denovo_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " model=denovo";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_dominant_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_dominant_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " model=dominant";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_recessif_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_recessif_by_gene';
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	if ($only_test) { return unless ($only_test eq $test_name); }
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " model=recessif";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_compound_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_compound_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " model=compound";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_recessif_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_recessif_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " model=recessif";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_filter_default_text_search {
	my $filter_chromosome = shift;
	my $test_name = 'test_filter_default_text_search';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_text='not CDK + nuclease'";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_filter_region_exclude_stats_region {
	my $filter_chromosome = shift;
	my $test_name = 'test_filter_region_exclude_stats_region';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=region";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_region=1:100:1000000:-1+1:2000000:3000000:-1+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_filter_region_exclude_stats_all {
	my $filter_chromosome = shift;
	my $test_name = 'test_filter_region_exclude_stats_all';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_region=1:100:1000000:-1+1:2000000:3000000:-1+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_filter_region_and_stats_region {
	my $filter_chromosome = shift;
	my $test_name = 'test_filter_region_and_stats_region';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=region";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_region=1:100:1000000:0+1:2000000:3000000:0+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_filter_region_and_stats_all {
	my $filter_chromosome = shift;
	my $test_name = 'test_filter_region_and_stats_all';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_region=1:100:1000000:0+1:2000000:3000000:0+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_atleast_byGenes {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_atleast_byGenes';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=gene";
	$cmd .= " nb=2";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_atleast_byGenes {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_atleast_byGenes';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " nb=2";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_atleast {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_atleast';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " nb=2";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_atleast {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_atleast';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " nb=2";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_disease {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_disease';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_diseases=2054";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_recessif {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_recessif';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " model=recessif";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_strict_denovo {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_strict_denovo';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=intergenic+intronic+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " model=strict-denovo";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_denovo {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_denovo';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " model=denovo";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_dominant {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_dominant';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " model=dominant";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_compound {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_compound';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " model=compound";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_recessif {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_recessif';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+intergenic+intronic+pseudogene+freq_05+freq_1 ";
	$cmd .= " mode=ind";
	$cmd .= " level_fam=variation";
	$cmd .= " level_ind=variation";
	$cmd .= " model=recessif";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_filter_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_filter_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_filter_he {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_filter_he';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_patients_fam_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_patients_fam_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_patients_fam_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_patients_fam_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_not=GUI+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_exclude_familial_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_exclude_familial_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=gene";
	$cmd .= " fam_not=LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_intersect_familial_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_intersect_familial_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=gene";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_exclude_familial {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_exclude_familial';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_not=LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_fam_intersect_familial {
	my $filter_chromosome = shift;
	my $test_name = 'test_fam_intersect_familial';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=fam";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " fam_and=GUI+LEF+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_in_the_attic_exclude_intersect_he_ho_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_in_the_attic_exclude_intersect_he_ho_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_NAT+GUI_YAE+GUI_REF+";
	$cmd .= " filter_not_patient=BOU_CHR+LEF_JER+LEF_MEL+";
	$cmd .= " filter_attic=GUI_YOS+";
	$cmd .= " filter_he=GUI_REF+";
	$cmd .= " filter_ho=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_filter_ho_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_filter_ho_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_filter_he_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_filter_he_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_intersect_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_intersect_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_exclude_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_exclude_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_in_the_attic_by_gene {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_in_the_attic_by_gene';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=gene";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_in_the_attic_exclude_intersect_he_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_in_the_attic_exclude_intersect_he_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_NAT+GUI_YAE+GUI_REF+";
	$cmd .= " filter_not_patient=BOU_CHR+LEF_JER+LEF_MEL+";
	$cmd .= " filter_attic=GUI_YOS+";
	$cmd .= " filter_he=GUI_REF+";
	$cmd .= " filter_ho=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_filter_ho {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_filter_ho';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_ho=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_filter_he {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_filter_he';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_he=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_intersect {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_intersect';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_patient=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_exclude {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_exclude';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_not_patient=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_in_the_attic {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_in_the_attic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	if ($only_test) { return unless ($only_test eq $test_name); }
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	$cmd .= " filter_attic=GUI_NAT+GUI_YAE+";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_ind_utr_splicing_esssplicing_stop_startstop_noframeshift {
	my $filter_chromosome = shift;
	my $test_name = 'test_ind_utr_splicing_esssplicing_stop_startstop_noframeshift';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+coding+intergenic+intronic+ncrna+maturemirna+pseudogene+frameshift+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_all_databases {
	my $filter_chromosome = shift;
	my $test_name = 'test_all_databases';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_only_intergenic {
	my $filter_chromosome = shift;
	my $test_name = 'test_only_intergenic';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+utr+splicing+essential_splicing+coding+stop+phase+upstream_downstream+intronic+ncrna+maturemirna+pseudogene+frameshift+non-frameshift+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_default {
	my $filter_chromosome = shift;
	my $test_name = 'test_default';
	if ($only_test) { return unless ($only_test eq $test_name); }
	my $cmd = $cmd_interface;
	$cmd .= " filter_chromosome=$filter_chromosome";
	$cmd .= " stat=all";
	$cmd .= " filter_type_variation=silent+intergenic+intronic+pseudogene+";
	$cmd .= " mode=ind";
	$cmd .= " level_ind=variation";
	$cmd .= " level_fam=variation";
	my ($obsJson, $time) = launchCmd($cmd, $test_name, $filter_chromosome);
	my $expJson = getExpJson($test_name, $filter_chromosome);
	if ($with_db) { $test_name = '[WITH DB] '.$test_name.' ['.$time.'s]'; }
	else { $test_name = '[NO DB] '.$test_name.' ['.$time.'s]'; }
	ok(compareTwoJson($obsJson, $expJson), $test_name) unless ($update_data);
}

sub test_compare_store_ids {
	my $project_name = shift;
	my $test_name = 'test_compare_store_ids';
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name );
	my ($hIds_from_vcf, $hIds_from_cache);
	foreach my $var (@{$project->getChromosome('1')->getStructuralVariations()}) {
		$hIds_from_vcf->{$var->id()} = undef;
	}
	my $buffer_2 = new GBuffer;
	my $project_2 = $buffer_2->newProjectCache( -name 			=> $project_name,
											-cache 		=> '1',
											-typeFilters 	=> '', );
	foreach my $var (@{$project_2->getChromosome('1')->getStructuralVariations()}) {
		$hIds_from_cache->{$var->id()} = undef;
	}
	if ($with_db) { $test_name = '[WITH DB] '.$test_name; }
	else { $test_name = '[NO DB] '.$test_name; }
	ok(compareTwoHashes($hIds_from_cache, $hIds_from_vcf), $test_name) unless ($update_data);
	exit(0);
}

sub compareTwoHashes {
	my ($h1, $h2) = @_;
	if ($not_exact_test) {
		my ($h1_new, $h2_new) = transformHashesWithCommonGenes($h1, $h2);
		$h1 = $h1_new;
		$h2 = $h2_new;
	}
	my $c = Data::Compare->new($h1, $h2);
	return $c->Cmp;
}

sub compareTwoJson {
	my ($json1, $json2) = @_;
	$json1 =~ s/.+\n//;
	$json2 =~ s/.+\n//;
	return 1 if ($no_test);
	my $h1 = JSON->new->utf8->decode($json1);
	my $h2 = JSON->new->utf8->decode($json2);
	if ($not_exact_test) {
		my ($h1_new, $h2_new) = transformHashesWithCommonGenes($h1, $h2);
		$h1 = $h1_new;
		$h2 = $h2_new;
	}
	if ($create_json) {
		my $f1 = $create_json.'/JSON_OBS.tab';
		my $f2 = $create_json.'/JSON_EXP.tab';
		open (FILE, ">$f1");
		print FILE Dumper $h1;
		close(FILE);
		open (FILE2, ">$f2");
		print FILE2 Dumper $h2;
		close(FILE2);
		exit(1);
	}
	my $c = Data::Compare->new($h1, $h2);
	return $c->Cmp;
}

sub transformHashesWithCommonGenes {
	my ($h1, $h2) = @_;
	my ($h1_new, $h2_new);
	
	$h1_new = $h1;
	delete $h1_new->{progress};
	my @lGenes;
	my $hGenesStats_1;
	my $hGenesStats_2;
	foreach my $hGene (@{$h1->{items}->{genes}}) {
		delete $hGene->{dejavu_capture_diag};
		delete $hGene->{is_omim};
		$hGenesStats_1->{$hGene->{id}} = $hGene;
	}
	delete $h1_new->{items}->{genes};
	foreach my $id (sort keys %{$hGenesStats_1}) { push(@lGenes, $hGenesStats_1->{$id}); }
	$h1_new->{items}->{genes} = \@lGenes;
	
	$h2_new = $h2;
	delete $h2_new->{progress};
	my @lGenes;
	foreach my $hGene (@{$h2->{items}->{genes}}) {
		delete $hGene->{dejavu_capture_diag};
		delete $hGene->{is_omim};
		$hGenesStats_2->{$hGene->{id}} = $hGene;
	}
	delete $h2_new->{items}->{genes};
	foreach my $id (sort keys %{$hGenesStats_2}) { push(@lGenes, $hGenesStats_2->{$id}); }
	$h2_new->{items}->{genes} = \@lGenes;
	return ($h1_new, $h2_new);
}


sub compareTwoFiles {
	my ($exp_file, $obs_file) = @_;
	if (compare($exp_file, $obs_file) == 0) {
		`rm $obs_file`;
		return 1;
	}
	else {
		warn "\n\nERROR: files are different !\n";
		warn " -> Exp: $exp_file\n";
		warn " -> Obs: $obs_file\n\n";
		return 0; 
	}
}

sub getExpJsonFileName {
	my ($test_name, $filter_chromosome) = @_;
	my ($file, $suffix);
	my $suffix;
	unless ($filter_chromosome eq '1') {
		$suffix = $filter_chromosome;
		$suffix =~ s/,/_/g;
	}
	$file = "$Bin/exp_json/$test_name.json";
	$file = "$Bin/exp_json/$test_name"."_"."$suffix.json" if ($suffix);
	return $file;
	
}

sub getExpJson {
	my ($test_name, $filter_chromosome) = @_;
	return 1 if ($no_test);
	my $file = getExpJsonFileName($test_name, $filter_chromosome);
	my $json;
	open (EXP, $file);
	while (<EXP>) { $json = $_; }
	close (EXP);
	return $json;
}

sub convertCmdToUrl {
	my $cmd = shift;
	my $url = "http://www.polyweb.fr/mbras/polyweb/vector/gene.html?project=$my_project_name";
	my @lArgs = split(' ', $cmd);
	my @lNewArgs;
	my $i = 0;
	for ($i = 1; $i<scalar(@lArgs); $i++) {
		my ($name, $value) = split('=', $lArgs[$i]);
		next if ($name eq 'project');
		next if ($name eq 'stat');
		next if ($name eq 'filter_chromosome');
		push (@lNewArgs, $lArgs[$i]);
	} 
	$url .= '&use_args='.join(';', @lNewArgs);
	return $url; 
}

sub launchCmd {
	my ($cmd, $test_name, $filter_chromosome) = @_;
	if ($update_data) {
		warn "\nWARN: $test_name json expected is changing !!\n";
		my $file = getExpJsonFileName($test_name, $filter_chromosome);
		my $cmd_test = $cmd;
		$cmd_test .= " >$file";
		`$cmd_test`;
		`chmod 777 $file`;
	}
	else {
		if ($print_url) {
			print "\n";
			print 'URL: '.convertCmdToUrl($cmd)."\n";
		}
		if ($print_cmd) {
			print "\n";
			print 'CMD: '.$cmd."\n";
		}
		my $sttime = DateTime->now();
		my $obs = `$cmd`;
		my $entime = DateTime->now;
		my $elapse = $entime - $sttime;
		return ($obs, 'NO TEST '.$elapse->in_units('seconds')) if ($no_test);
		return ($obs, $elapse->in_units('seconds'));
	}
}
