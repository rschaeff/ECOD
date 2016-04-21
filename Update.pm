package ECOD::Update;
require Exporter;

use warnings;
use strict;
use Carp;
$Carp::Verbose = 1;

use lib('/home/rschaeff/lib');
use XML::LibXML;
use XML::Grishin;
use Domains::PDBml;
use Domains::Range;
use Domains::SGE;
use SQL::Util;
use ECOD::Reference;
load_references();

use Text::CSV;

our @ISA = qw(Exporter);

our @EXPORT = (
    "&cluster_pfgroups_by_pfam",,
    "&cluster_pfgroups_by_ecodf",           "&convert_hmmer_txt_to_xml",
    "&hmmer_annotation",                    "&convert_hmmer_to_pfam_annotation",
    "&calculate_divergence",                "&convert_hmmer_to_ecodf_annotation",
    "&convert_manual_range_to_seqid_range", "&convert_mc_manual_range_to_seqid_range",
    "&determine_update_step",               "&generate_blast_library",
    "&assemble_hmmer_profile_library",      "&generate_distributable_files",
    "&assemble_hh_profile_library",         "&assemble_domain_fasta_library",
    "&assemble_chainwise_fasta_library",    "&generate_hmmer_jobs",
    "&generate_domain_seqinput_file",       "&generate_domain_range_cache_file",
    "&generate_domain_fasta_jobs",          "&generate_domain_data_directories",
    "&generate_chainwise_fasta_jobs",       "&generate_hhblits_profile_jobs",
    "&generate_ecod_domain_dump",           "&generate_ecod_fasta_dump",
    "&generate_ecod_pdb_tarball",           "&generate_ecod_rep_list",
    "&generate_distributable_files",        "&generate_ecod_rep_list",
    "&ecodf_acc_to_pfam_acc",               "&get_pdbs_by_latest_week",
    "&generate_ligand_annotation_jobs",     "&generate_pdb_structure_jobs",
    "&generate_image_jobs",                 "&ligand_annotation",
    "&build_chainwise",                     "&find_overlaps",
    "&delete_chainwise_domains",            "&pf_reorder",
    "&rep_check",                           "&rep_check_ecodf",
    "&add_ecod_domain_ids",                 "&fix_manual_range_nodes",
    "&fix_implicit_range_nodes",            "&apply_merge",
    "&apply_merge_v2",                      "&merge_new_ecod_pre_xml_to_old_ecod_xml",
    "&parse_ecod_txt_to_xml",               "&process_run_list_summary_to_ecod",
    "&convert_side_load_domains_to_xml",    "&process_side_load_domains_to_ecod",
    "&isMultiChainDomain",                  "&update_release_statistics",
    "&write_stats_gnuplot_script",          "&basic_ecod_stats",
);

my $DOMAIN_IMAGE_EXE = '/data/ecod/domain_image_new.py';
if ( !-f $DOMAIN_IMAGE_EXE ) { die "ERROR! Could not find domain image exe $DOMAIN_IMAGE_EXE\n"; }
my $PYMOL_EXE = '/home/rschaeff/bin/pymol/pymol';
if ( !-f $PYMOL_EXE ) { die "ERROR! Could not find pymol image exe $PYMOL_EXE\n"; }
my $XML_PDB_EXE = '/data/ecod/database_versions/bin/xml_pdb_process_v3.pl';
if ( !-f $XML_PDB_EXE ) { die "ERROR! Could not find pdb formatter $XML_PDB_EXE\n"; }

my $GENERATE_SEQDB_EXE = '/data/ecod/bin/generate_seqdb_from_text.pl';
if ( !-f $GENERATE_SEQDB_EXE ) { die "ERROR! file not found $GENERATE_SEQDB_EXE\n"; }
my $GENERATE_CHAINWISE_SEQDB_EXE = '/data/ecod/bin/chainwise_seq_generate_para.pl';
if ( !-f $GENERATE_CHAINWISE_SEQDB_EXE ) { die "ERROR! file not found $GENERATE_CHAINWISE_SEQDB_EXE\n"; }

my $HHBLITS_EXE = '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhblits';
if ( !-f $HHBLITS_EXE ) { die "ERROR! HHblits exe not found $HHBLITS_EXE\n" }
my $HHMAKE_EXE = '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhmake';
if ( !-f $HHMAKE_EXE ) { die "ERROR! HHmake exe not found $HHMAKE_EXE\n" }
my $MAKEBLAST_EXE = "/usr1/ncbi-blast-2.2.25+/bin/makeblastdb";
if ( !-f $MAKEBLAST_EXE ) { die "ERROR! makeblastdb not found $MAKEBLAST_EXE\n" }
my $HMMSCAN_EXE = "/usr1/local/bin/hmmscan";
if ( !-f $HMMSCAN_EXE ) { die "ERROR! hmmscan not found $HMMSCAN_EXE\n"; }
my $SINGLE_LIGAND_ANNOT_EXE = '/data/ecod/database_versions/bin/single_annotate_ecod_ligand.pl';
if ( !-f $SINGLE_LIGAND_ANNOT_EXE ) { die "ERROR! File not found: $SINGLE_LIGAND_ANNOT_EXE\n" }
my $HH_LIB = '/home/rschaeff/hhsuite-2.0.15-linux-x86_x64/lib';

#HMM parameters need to be offloaded, this is bad
my $HMM_LABEL = 'ECODf';

my $HMM_VNUM = 'v28';
my $HMM_LIB  = "/data/ecodf/ecodf.$HMM_VNUM.hmm";
if ( !-f $HMM_LIB ) { die "ERROR! HMM lib not found: $HMM_LIB\n" }

my $ECODF_XML_FN = "/data/ecodf/ecodf.$HMM_VNUM.xml";
if ( !-f $ECODF_XML_FN ) { die "Could not find ECODf xml $ECODF_XML_FN\n" }

my $DEBUG = 0;

my $SCOP_REFERENCE  = 'v1.75';
my $GAP_TOL         = 20;
my $INCLUDE_LIGANDS = 1;

my $LIGAND_INCLUSION_THRESHOLD = 4;    #Angstroms
my $USE_OBSOLETE               = 1;
my $FORCE_OVERWRITE            = 0;

my $FORCE_REPLACE = 0;

umask(0002);
my $UID = 1219;                        #rschaeff
my $GID = 902;                         #ecod

my $UNIPROT_LIB =
  '/home/rschaeff_1/side/wtf_hh/uniprot20_2012_10_klust20_dc_2012_12_10/uniprot20_2012_10_klust20_dc_2012_12_10';
my $NR_LIB           = '/home/rschaeff_1/side/wtf_hh/nr20_12Aug11';
my $HH_BLITS_CPU     = 8;
my $PSIPRED_DIR      = '/usr1/psipred/bin';
my $PSIPRED_DATA_DIR = '/usr1/psipred/data';
my $DOMAIN_DATA_DIR  = '/data/ecod/domain_data';
my $CHAIN_DATA_DIR   = '/data/ecod/chain_data';

sub basic_ecod_stats {
    my $sub = 'basic_ecod_stats';
    my ( $ref_xml_fn, $ref_version ) = @_;

    open( my $xml_fh, $ref_xml_fn ) or die "ERROR! $sub: Could not open $ref_xml_fn for reading:$!\n";

    my $ecod_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );
    close $xml_fh;

    my %stats;
    my $arch      = $ecod_xml_doc->findnodes('//architecture')->size();
    my $x_groups  = $ecod_xml_doc->findnodes('//x_group')->size();
    my $h_groups  = $ecod_xml_doc->findnodes('//h_group')->size();
    my $f_groups  = $ecod_xml_doc->findnodes('//f_group')->size();
    my $pf_groups = $ecod_xml_doc->findnodes('//pf_group')->size();
    my $domains   = $ecod_xml_doc->findnodes('//domain')->size();

    my %pdbs;
    my %peptide_chain;
    my %chain;
    foreach my $structure_node ( $ecod_xml_doc->findnodes('//structure') ) {
        my $pdb_id = $structure_node->findvalue('@pdb_id');
        $pdbs{$pdb_id}++;
    }
    my $peptides = $ecod_xml_doc->findnodes('//peptide')->size();
    foreach my $peptide_node ( $ecod_xml_doc->findnodes('//peptide')->get_nodelist() ) {
        my $pdb_id = lc( $peptide_node->findvalue('@pdb_id') );
        $pdbs{$pdb_id}++;
        my $chain_id  = $peptide_node->findvalue('@chain_id');
        my $pdb_chain = $pdb_id . "_" . $chain_id;
        $chain{$pdb_chain}++;
        $peptide_chain{$pdb_chain}++;
    }

    #synthetics
    my $synth = $ecod_xml_doc->findnodes('//synthetic')->size();
    foreach my $synth_node ( $ecod_xml_doc->findnodes('//synth')->get_nodelist() ) {
        my $pdb_id = lc( $synth_node->findvalue('@pdb_id') );
        $pdbs{$pdb_id}++;
        my $chain_id  = $synth_node->findvalue('@chain_id');
        my $pdb_chain = $pdb_id . "_" . $chain_id;
        $chain{$pdb_chain}++;
        $peptide_chain{$pdb_chain}++;
    }

    #pss
    my $pss = $ecod_xml_doc->findnodes('//pss')->size();
    foreach my $pss_node ( $ecod_xml_doc->findnodes('//pss')->get_nodelist() ) {
        my $pdb_id = lc( $pss_node->findvalue('@pdb_id') );
        $pdbs{$pdb_id}++;
        my $chain_id  = $pss_node->findvalue('@chain_id');
        my $pdb_chain = $pdb_id . "_" . $chain_id;
        $chain{$pdb_chain}++;
        $peptide_chain{$pdb_chain}++;
    }

    #nonpep_poly
    my $npp = $ecod_xml_doc->findnodes('//nonpeptide_poly')->size;
    foreach my $npp_node ( $ecod_xml_doc->findnodes('//nonpeptide_poly')->get_nodelist() ) {

        my $pdb_id = lc( $npp_node->findvalue('@pdb') );
        $pdbs{$pdb_id}++;
        my $chain_str = $npp_node->findvalue('@chain');
        my @chains = split( /,/, $chain_str );
        foreach my $chain_id (@chains) {
            my $pdb_chain = $pdb_id . "_" . $chain_id;
            $chain{$pdb_chain}++;
        }
    }

    #Coiled-coil
    my $coil = $ecod_xml_doc->findnodes('//coil')->size();
    foreach my $coil_node ( $ecod_xml_doc->findnodes('//coil')->get_nodelist() ) {
        my $pdb_id = lc( $coil_node->findvalue('@pdb_id') );
        $pdbs{$pdb_id}++;
        my $chain_id  = $coil_node->findvalue('@chain_id');
        my $pdb_chain = $pdb_id . "_" . $chain_id;
        $chain{$pdb_chain}++;
        $peptide_chain{$pdb_chain}++;
    }

    #MCC
    my $mcc = $ecod_xml_doc->findnodes('//mcc')->size();
    foreach my $mcc_node ( $ecod_xml_doc->findnodes('//mcc')->get_nodelist() ) {

        my $pdb_id = lc( $mcc_node->findvalue('@pdb_id') );
        $pdbs{$pdb_id}++;
        my $chain_id  = $mcc_node->findvalue('@chain');
        my $pdb_chain = $pdb_id . "_" . $chain_id;
        $chain{$pdb_chain}++;
        $peptide_chain{$pdb_chain}++;
    }

    #CreatedOn

    my $xml_date;
    my $graph_date;
    $ref_version =~ /develop((\d+)\w?)/;
    my $ref_int        = $2;
    my $ref_subversion = $1;
    my $ref_v          = 'v' . $ref_int;
    print "$ref_version $ref_int $ref_v $ref_subversion\n";

    my %month_translation = (
        "Jan" => "01",
        "Feb" => "02",
        "Mar" => "03",
        "Apr" => "04",
        "May" => "05",
        "Jun" => "06",
        "Jul" => "07",
        "Aug" => "08",
        "Sep" => "09",
        "Oct" => "10",
        "Nov" => "11",
        "Dec" => "12"
    );

    if ( $ecod_xml_doc->exists(qq{//createdOn[\@version="$ref_version"]}) ) {
        my $createdOn = $ecod_xml_doc->findnodes(qq{//createdOn[\@version="$ref_version"]})->get_node(1);
        $xml_date = $createdOn->textContent;
    }
    elsif ( $ecod_xml_doc->exists(qq{//createdOn[\@version="$ref_int"]}) ) {
        my $createdOn = $ecod_xml_doc->findnodes(qq{//createdOn[\@version="$ref_int"]})->get_node(1);
        $xml_date = $createdOn->textContent;
    }
    elsif ( $ecod_xml_doc->exists(qq{//createdOn[\@version="$ref_v"]}) ) {
        my $createdOn = $ecod_xml_doc->findnodes(qq{//createdOn[\@version="$ref_v"]})->get_node(1);
        $xml_date = $createdOn->textContent;
    }
    elsif ( $ecod_xml_doc->findnodes(qq{//createdOn})->size() == 1 ) {
        my $createdOn = $ecod_xml_doc->findnodes(qq{//createdOn})->get_node(1);
        $xml_date = $createdOn->textContent;
    }
    else {
        print "WARNING! No date found for $ref_version\n";
        $xml_date = "NOT_FOUND";
    }

    if ( $xml_date =~ /(\w{3})\s(\w{3})\s+(\d+)\s+(\d+\:\d+\:\d+)\s+(\w+)\s+(\w{4})/ ) {
        my $weekday  = $1;
        my $month    = $month_translation{$2};
        my $day      = sprintf "%02i", $3;
        my $time     = $4;
        my $timezone = $5;
        my $year     = $6;
        $graph_date = $year . $month . $day;
    }
    else {
        $graph_date = "NOT_FOUND";
    }
    print "$xml_date $graph_date\n";

    $stats{arch} = $arch;
    $stats{xgrp} = $x_groups;
    $stats{hgrp} = $h_groups;
    $stats{tgrp} = $f_groups;
    $stats{fgrp} = $pf_groups;
    $stats{domn} = $domains;
    $stats{pdbs} = keys %pdbs;
    $stats{date} = $graph_date;

    return \%stats;
}

sub write_stats_gnuplot_script {
    my $sub = 'write_stats_gnuplot_script';
    my ( $ref_aref, $out_fn ) = @_;

    my $out_gplt_fn = $out_fn;
    $out_gplt_fn =~ s/txt/gplt/;

    open( OUT, ">$out_gplt_fn" );

    print OUT "set term svg size 800,600\n";
    print OUT "set output 'stats_groups.svg'\n";

    print OUT "set ylabel 'Groups'\n";
    print OUT "set xlabel 'Version'\n";

    my @ref_labels;
    for ( my $i = 0 ; $i < scalar(@$ref_aref) ; $i++ ) {
        push( @ref_labels, "\"$$ref_aref[$i]\" $i" );
    }

    my $string = join( ",", @ref_labels );

    print OUT "set xtics rotate ($string)\n";
    print OUT "set key below\n";
    print OUT
"plot '$out_fn' u 3 t 'X-groups', '$out_fn' u 4 t 'H-groups', '$out_fn' u 5 t 'T-groups', '$out_fn' u 6 t 'F-groups'\n";

    print OUT "set output 'stats_domains.svg'\n";
    print OUT "set ylabel 'Domain/Structure Counts'\n";

    print OUT "plot '$out_fn' u 7 t 'Domains', '$out_fn' u 8 t 'PDBs'\n";
}

sub update_release_statistics {
    my $sub = 'update_release_statistics';

    my ( $ecod_master_ref_xml, $stats_xml_fn, $stats_txt_fn ) = @_;

    my %references;
    foreach my $reference_node ( $ecod_master_ref_xml->findnodes('//reference') ) {

        my $version = $reference_node->findvalue('@version');
        if ( $version =~ /repsonly/ ) { next }        #This is a hack, remove those versions from master index
        if ( $reference_node->exists('ref_xml') ) {
            my $ref_xml = $reference_node->findvalue('ref_xml');
            $references{$version} = $ref_xml;
        }

    }

    my @references =
      sort { $a =~ /develop(\d+)/; my $anum = $1; $b =~ /develop(\d+)/; my $bnum = $1; $anum <=> $bnum }
      keys %references;

    open( my $xml_fh, $stats_xml_fn ) or die "ERROR! Could not open $stats_xml_fn for reading:$!\n";
    my $stats_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );

    my $stats_list_node = $stats_xml_doc->findnodes('//release_statistics_list')->get_node(1);

    foreach my $ref (@references) {
        print "#$ref $references{$ref}\n";
        if ( $stats_xml_doc->exists(qq{//release_statistics[\@version="$ref"]}) ) {
            next;
        }
        if ( !-f $references{$ref} ) { print "WARNING! $references{$ref} not found\n"; next }

        my $stats_ref = basic_ecod_stats( $references{$ref}, $ref );

        my $version    = $ref;
        my $arch_count = $$stats_ref{arch};
        my $x_count    = $$stats_ref{xgrp};
        my $h_count    = $$stats_ref{hgrp};
        my $t_count    = $$stats_ref{tgrp};
        my $f_count    = $$stats_ref{fgrp};

        my $dom_count  = $$stats_ref{domn};
        my $pdb_count  = $$stats_ref{pdbs};
        my $week_label = $$stats_ref{date};
        my $stats_node = $stats_xml_doc->createElement('release_statistics');

        $stats_node->setAttribute( 'version',    $version );
        $stats_node->setAttribute( 'arch_count', $arch_count );
        $stats_node->setAttribute( 'x_count',    $x_count );
        $stats_node->setAttribute( 'h_count',    $h_count );
        $stats_node->setAttribute( 't_count',    $t_count );
        $stats_node->setAttribute( 'f_count',    $f_count );
        $stats_node->setAttribute( 'dom_count',  $dom_count );
        $stats_node->setAttribute( 'pdb_count',  $pdb_count );
        $stats_node->setAttribute( 'week_label', $week_label );

        $stats_list_node->appendChild($stats_node);

    }

    my %dates_by_ref;
    foreach my $ref (@references) {
        my $stats_node = $stats_xml_doc->findnodes(qq{//release_statistics[\@version="$ref"]})->get_node(1);
        my $week_label = $stats_node->findvalue('@week_label');
        $dates_by_ref{$ref} = $week_label;
    }

    my $out_fn = $stats_xml_fn;
    open( OUT, ">$out_fn" ) or die "ERROR! Could not open $out_fn for writing:$!\n";
    my $doc_string = $stats_xml_doc->toString(1);
    print OUT $doc_string;
    close OUT;
    my @dates;
    foreach my $ref (@references) {
        push( @dates, $dates_by_ref{$ref} );
    }

    my $stats_fn = $stats_txt_fn;
    open( OUT, ">$stats_fn" ) or die "ERROR! could not open $stats_fn for writing:$!\n";
    print OUT "#version\tarch_count\tx_count\th_count\tt_count\tf_count\tdom_count\tpdb_count\tdate\n";
    foreach my $stats_node ( $stats_xml_doc->findnodes('//release_statistics') ) {

        my $version    = $stats_node->findvalue('@version');
        my $arch_count = $stats_node->findvalue('@arch_count');
        my $x_count    = $stats_node->findvalue('@x_count');
        my $h_count    = $stats_node->findvalue('@h_count');
        my $t_count    = $stats_node->findvalue('@t_count');
        my $f_count    = $stats_node->findvalue('@f_count');
        my $dom_count  = $stats_node->findvalue('@dom_count');
        my $pdb_count  = $stats_node->findvalue('@pdb_count');
        my $week_label = $stats_node->findvalue('@week_label');

        printf OUT "%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\n",
          $version,
          $arch_count,
          $x_count,
          $h_count,
          $t_count,
          $f_count,
          $dom_count,
          $pdb_count,
          $week_label;
    }
    close OUT;

    write_stats_gnuplot_script( \@dates, $stats_fn );

}

sub ecodf_acc_to_pfam_acc {
    my $sub = 'ecodf_acc_to_pfam_acc';

    my $ecodf_xml_fn = "/data/ecodf/ecodf.$HMM_VNUM.xml";
    if ( !-f $ecodf_xml_fn ) {
        die "ERROR! $sub: Could not find ECODf file $ecodf_xml_fn\n";
    }
    my $ecodf_xml_doc = xml_open($ecodf_xml_fn);

    my %ecodf2pfam_accs;
    foreach my $family_node ( $ecodf_xml_doc->findnodes('//family') ) {
        my $ecodf_acc = $family_node->findvalue('@ecodf_acc');
        $ecodf_acc =~ /(EF\d+)\./;
        $ecodf_acc = $1;
        if ( $family_node->exists('@pfam_acc') ) {
            my $pfam_acc = $family_node->findvalue('@pfam_acc');
            $pfam_acc =~ /(PF\d+)\./;
            $pfam_acc = $1;
            $ecodf2pfam_accs{$ecodf_acc} = $pfam_acc;
        }
    }
    return \%ecodf2pfam_accs;
}
my $FGROUP_CLUSTER_METHOD = 'blastclust';

sub generate_distributable_files {
    my $sub = 'generate_distributable_files';

    my ( $version, $master_ecod_xml, $tag ) = @_;
    print "$sub: $version\n";
    my $recent_week = $master_ecod_xml->findvalue('//@mostRecentWeek');

    #current ECOD XML;
    if ( !$master_ecod_xml->exists(qq{//reference[\@version="$version"]}) ) {
        die "ERROR! $sub: No $version found in master ref\n";
    }
    my $ecod_xml_fn = $master_ecod_xml->findvalue(qq{//reference[\@version="$version"]/ref_xml});
    open( my $xml_fh, $ecod_xml_fn ) or die "ERROR! Could not open $ecod_xml_fn for reading:$!\n";
    my $ecod_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );

    my $ecod_fasta_fn = $master_ecod_xml->findvalue(qq{//reference[\@version="$version"]/domain_fa});
    my $top_dir       = $master_ecod_xml->findvalue(qq{//reference[\@version="$version"]/local_dir});

    my $current_ecod_domain_dump = "$top_dir/ecod.$version.domains.txt";
    if ($tag) {
        $current_ecod_domain_dump = "$top_dir/ecod.$version.$tag.domains.txt";
    }
    my %pf_lookup;
    if ( !-f $current_ecod_domain_dump ) {
        generate_ecod_domain_dump( $ecod_xml_doc, $current_ecod_domain_dump, $version, \%pf_lookup, $tag,
            $FGROUP_CLUSTER_METHOD );
        printf "keys %i\n", scalar keys %pf_lookup;
    }

    if ( !-f $ecod_fasta_fn ) { die "ERROR! Could not find ecod fasta $ecod_fasta_fn\n"; }
    my $current_ecod_fasta_dump = "$top_dir/ecod.$version.fasta.txt";
    if ($tag) {
        $current_ecod_fasta_dump = "$top_dir/ecod.$version.$tag.fasta.txt";
    }
    if ( !-f $current_ecod_fasta_dump ) {
        generate_ecod_fasta_dump( $ecod_fasta_fn, $current_ecod_fasta_dump, \%pf_lookup, $tag );
    }

    #current ECOD pdb tarball
    my $current_ecod_pdb_tarball = "$top_dir/ecod.$version.pdb.tar.gz";
    if ($tag) {
        $current_ecod_pdb_tarball = "$top_dir/ecod.$version.$tag.pdb.tar.gz";
    }
    if ( !-f $current_ecod_pdb_tarball ) {
        generate_ecod_pdb_tarball( $ecod_xml_doc, $current_ecod_pdb_tarball, $tag );
    }

    #current ECOD rep txt list
    my $current_ecod_pdb_replist = "$top_dir/reps.$version.txt";
    if ( !-f $current_ecod_pdb_replist ) {
        generate_ecod_rep_list( $ecod_xml_doc, $current_ecod_pdb_replist );
    }

}

sub generate_ecod_rep_list {
    my ( $ecod_xml_doc, $out_fn ) = @_;

    open( OUT, ">$out_fn" ) or die "ERROR! COuld not open $ARGV[0] for writing:$!\n";

    foreach my $rep_domain ( $ecod_xml_doc->findnodes('//domain[@manual_rep="true"]') ) {
        my $uid            = $rep_domain->findvalue('@uid');
        my $ecod_domain_id = $rep_domain->findvalue('@ecod_domain_id');
        my $seqid_range    = $rep_domain->findvalue('seqid_range');
        print OUT "$uid\t$ecod_domain_id\t$seqid_range\n";
    }
    close OUT;
}

sub generate_ecod_domain_dump {
    my $sub = 'generate_ecod_domain_dump';
    my ( $ecod_xml_doc, $output_fn, $current_version, $pf_lookup, $use_pfam_acc, $tag, $method ) = @_;

    #print "o: $output_fn, c: $current_version\n";
    open( OUT, ">$output_fn" ) or die "ERROR! $sub: Could not open $output_fn for writing:$!\n";

    my $ecodf_translate;
    if ($use_pfam_acc) {
        $ecodf_translate = ecodf_acc_to_pfam_acc();
    }
    if ( !$tag ) { $tag = 'none' }
    my %nodes = (
        'none' => \&find_strict_domain_nodes,
        F40    => \&find_strict_domain_reps,
        F70    => \&find_strict_domain_reps,
        F99    => \&find_strict_domain_reps,
    );

    #Header needs to be moved outside of code
    print OUT "#$output_fn\n";
    print OUT "#ECOD version $current_version\n";
    print OUT "#Domain list version 1.4\n";
    print OUT "#Grishin lab (http://prodata.swmed.edu/ecod)\n";
    print OUT
"#uid\tecod_domain_id\tmanual_rep\tf_id\tpdb\tchain\tpdb_range\tseqid_range\tarch_name\tx_name\th_name\tt_name\tf_name\tasm_status\tligand\n";
    foreach my $arch_node ( $ecod_xml_doc->findnodes('//architecture')->get_nodelist() ) {
        my $arch_name = $arch_node->findvalue('@arch_name');
        foreach my $x_node ( $arch_node->findnodes('x_group')->get_nodelist() ) {
            my $x_name = 'NO_X_NAME';
            if ( $x_node->exists('@name') ) {
                $x_name = $x_node->findvalue('@name');
            }
            foreach my $h_node ( $x_node->findnodes('h_group') ) {
                my $h_name = 'NO_H_NAME';
                if ( $h_node->exists('@name') ) {
                    $h_name = $h_node->findvalue('@name');
                }

                foreach my $f_node ( $h_node->findnodes('f_group') ) {
                    my $f_name = 'NO_T_NAME';
                    if ( $f_node->exists('@name') ) {
                        $f_name = $f_node->findvalue('@name');
                    }
                    my $f_id = $f_node->findvalue('@f_id');
                    foreach my $domain_node ( $nodes{$tag}->( $f_node, $tag, $method ) ) {
                        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
                        my $uid            = $domain_node->findvalue('@uid');
                        my $pdb_range      = get_pdb_range($domain_node);
                        my $seqid_range    = get_seqid_range($domain_node);

                        #my $pdb_range           = $domain_node->findvalue('derived_range');
                        my $pdb   = $domain_node->findvalue('structure/@pdb_id');
                        my $chain = $domain_node->findvalue('structure/@chain_id');

                        my $manual_rep = "AUTO_NONREP";
                        if ( $domain_node->findvalue('@manual_rep') eq 'true' ) {
                            $manual_rep = "MANUAL_REP";
                        }

                        my $pf_type         = 'F_UNCLASSIFIED';
                        my $domain_assembly = 'NOT_DOMAIN_ASSEMBLY';
                        my $ligand_str      = 'NO_LIGANDS_4A';

                        if ( $domain_node->findvalue('ligand_str') =~ /\w+/ ) {
                            $ligand_str = $domain_node->findvalue('ligand_str');
                        }
                        $$pf_lookup{$uid} = $f_id;
                        print OUT
"$uid\t$ecod_domain_id\t$manual_rep\t$f_id\t$pdb\t$chain\t$pdb_range\t$seqid_range\t$arch_name\t\"$x_name\"\t\"$h_name\"\t\"$f_name\"\t$pf_type\t$domain_assembly\t$ligand_str\n";
                    }
                    foreach my $domain_assembly_node ( $f_node->findnodes('domain_assembly') ) {
                        foreach my $domain_node ( $nodes{$tag}->( $domain_assembly_node, $tag, $method ) ) {
                            my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
                            my $uid            = $domain_node->findvalue('@uid');

                            #my $pdb_range           = $domain_node->findvalue('derived_range');
                            my $pdb_range   = get_pdb_range($domain_node);
                            my $seqid_range = get_seqid_range($domain_node);
                            my $pdb         = $domain_node->findvalue('structure/@pdb_id');
                            my $chain       = $domain_node->findvalue('structure/@chain_id');

                            my $manual_rep = "AUTO_NONREP";
                            if ( $domain_node->findvalue('@manual_rep') eq 'true' ) {
                                $manual_rep = "MANUAL_REP";
                            }

                            my $pf_type         = 'F_UNCLASSIFIED';
                            my $domain_assembly = 'NOT_DOMAIN_ASSEMBLY';
                            my $ligand_str      = 'NO_LIGANDS_4A';
                            $domain_assembly = $domain_node->parentNode->findvalue('@uid');
                            if ( $domain_node->findvalue('ligand_str') =~ /\w+/ ) {
                                $ligand_str = $domain_node->findvalue('ligand_str');
                            }
                            $$pf_lookup{$uid} = $f_id;
                            print OUT
"$uid\t$ecod_domain_id\t$manual_rep\t$f_id\t$pdb\t$chain\t$pdb_range\t$seqid_range\t$arch_name\t\"$x_name\"\t\"$h_name\"\t\"$f_name\"\t$pf_type\t$domain_assembly\t$ligand_str\n";
                        }
                    }

                    foreach my $pf_node ( $f_node->findnodes('pf_group') ) {
                        my $pf_type = 'UNK_F_TYPE';
                        my $pf_id   = $pf_node->findvalue('@pf_id');
                        if ( $pf_node->findvalue('@pfam_cluster') eq 'true' ) {

                            #$pf_type = 'PFAM_F_GROUP';
                            $pf_type = $pf_node->findvalue('pfam/@pfam_acc');
                        }
                        elsif ( $pf_node->findvalue('@hh_cluster') eq 'true' ) {
                            $pf_type = 'HH_F_GROUP';
                        }
                        elsif ( $pf_node->findvalue('@jb_cluster') eq 'true' ) {
                            $pf_type = 'F_UNCLASSIFIED';
                        }
                        elsif ( $pf_node->findvalue('@ecodf_cluster') eq 'true' ) {

                            #$pf_type = 'ECODF_F_GROUP';
                            $pf_type = $pf_node->findvalue('ecodf/@ecodf_acc');
                            if ($use_pfam_acc) {
                                my @f_names = split( /\,/, $pf_type );
                                my @new_f_names;
                                foreach my $f_name (@f_names) {
                                    if ( exists $$ecodf_translate{$f_name} ) {
                                        push( @new_f_names, $$ecodf_translate{$f_name} );
                                    }
                                    else {
                                        push( @new_f_names, $f_name );
                                    }
                                }
                                $pf_type = join( ",", @new_f_names );
                            }
                        }

                        foreach my $domain_node ( $nodes{$tag}->( $pf_node, $tag, $method ) ) {
                            my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
                            my $pdb_range;
                            if ( $domain_node->exists('range') ) {
                                $pdb_range = $domain_node->findvalue('range');
                            }
                            elsif ( $domain_node->exists('derived_range') ) {
                                $pdb_range = $domain_node->findvalue('derived_range');
                            }
                            else {
                                $pdb_range = 'NO_RANGE';
                            }
                            my $seqid_range = get_seqid_range($domain_node);
                            my $manual_rep  = "AUTO_NONREP";
                            if ( $domain_node->findvalue('@manual_rep') eq 'true' ) {
                                $manual_rep = "MANUAL_REP";
                            }

                            my $ligand_str      = 'NO_LIGANDS_4A';
                            my $pdb             = $domain_node->findvalue('structure/@pdb_id');
                            my $chain           = $domain_node->findvalue('structure/@chain_id');
                            my $domain_assembly = 'NOT_DOMAIN_ASSEMBLY';
                            if ( $domain_node->parentNode->nodeName eq 'domain_assembly' ) {

                                #$domain_assembly  = 'IS_DOMAIN_ASSEMBLY';
                                $domain_assembly = $domain_node->parentNode->findvalue('@uid');
                            }
                            if ( $domain_node->findvalue('ligand_str') =~ /\w+/ ) {
                                $ligand_str = $domain_node->findvalue('ligand_str');
                            }
                            $$pf_lookup{$uid} = $pf_id;

                            #print "$ecod_domain_id\t$f_id\t$pdb\t$chain\t$pdb_range\n";
                            print OUT
"$uid\t$ecod_domain_id\t$manual_rep\t$pf_id\t$pdb\t$chain\t$pdb_range\t$seqid_range\t$arch_name\t\"$x_name\"\t\"$h_name\"\t\"$f_name\"\t$pf_type\t$domain_assembly\t$ligand_str\n";
                        }
                        foreach my $domain_assembly_node ( $pf_node->findnodes('domain_assembly') ) {
                            foreach my $domain_node ( $nodes{$tag}->( $domain_assembly_node, $tag, $method ) ) {
                                my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
                                my $uid            = $domain_node->findvalue('@uid');

                                #my $pdb_range           = $domain_node->findvalue('derived_range');
                                my $pdb_range   = get_pdb_range($domain_node);
                                my $seqid_range = get_seqid_range($domain_node);
                                my $pdb         = $domain_node->findvalue('structure/@pdb_id');
                                my $chain       = $domain_node->findvalue('structure/@chain_id');

                                my $manual_rep = "AUTO_NONREP";
                                if ( $domain_node->findvalue('@manual_rep') eq 'true' ) {
                                    $manual_rep = "MANUAL_REP";
                                }

                                my $domain_assembly = 'IS_DOMAIN_ASSEMBLY';
                                $domain_assembly = $domain_node->parentNode->findvalue('@uid');

                                my $ligand_str = 'NO_LIGANDS_4A';
                                if ( $domain_node->findvalue('ligand_str') =~ /\w+/ ) {
                                    $ligand_str = $domain_node->findvalue('ligand_str');
                                }
                                $$pf_lookup{$uid} = $f_id;
                                print OUT
"$uid\t$ecod_domain_id\t$manual_rep\t$f_id\t$pdb\t$chain\t$pdb_range\t$seqid_range\t$arch_name\t\"$x_name\"\t\"$h_name\"\t\"$f_name\"\t$pf_type\t$domain_assembly\t$ligand_str\n";
                            }
                        }
                    }
                }
            }
        }
    }
}

sub generate_ecod_fasta_dump {
    my $sub = 'generate_ecod_fasta_dump';
    my ( $ecod_fasta_fn, $ecod_fasta_output_fn, $pf_lookup_href ) = @_;

    printf "e: $ecod_fasta_fn eo: $ecod_fasta_output_fn keys: %i\n", scalar keys %$pf_lookup_href;

    open( IN,  $ecod_fasta_fn )           or die "ERROR! $sub: Could not open $ecod_fasta_fn for reading:$!\n";
    open( OUT, ">$ecod_fasta_output_fn" ) or die "ERROR! $sub: Could not open $ecod_fasta_output_fn for writing:$!\n";

    my $start = 0;
    while ( my $ln = <IN> ) {
        if ( $ln =~ /^>/ ) {
            $ln =~ s/^>//;
            my @F = split( /\s+/, $ln );

            my $uid            = $F[2];
            my $range          = $F[1];
            my $ecod_domain_id = $F[0];
            if ( !$$pf_lookup_href{$uid} ) {

                #print "No fid for $uid?\n";
                $start = 0;
            }
            else {
                print OUT ">$uid|$ecod_domain_id|$$pf_lookup_href{$uid}|$range\n";
                $start = 1;
            }
        }
        else {
            if ($start) {
                print OUT $ln;
            }
        }
    }
}

sub generate_ecod_pdb_tarball {
    my $sub = 'generate_ecod_pdb_tarball';
    my ( $ecod_xml_doc, $output_tarball, $tag ) = @_;

    my @pdbs;
    my $missing = 0;
    my $total   = 0;
    my $xpath   = '//domain[@manual_rep="true"]';
    if ($tag) {
        $xpath = "//domain/cluster[\@level='$tag'][\@domain_rep='true']";
    }
    foreach my $domain_node ( $ecod_xml_doc->findnodes($xpath) ) {
        if ($tag) {
            $domain_node = $domain_node->parentNode;
        }

        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
        my $uid            = $domain_node->findvalue('@uid');
        my $short_uid      = substr( $uid, 2, 5 );

        my $pdb_fn = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.pdbnum.pdb";

        #my $pdb_fn = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.pdb";

        if ( !-f $pdb_fn ) {

            #	print "WARNING! No pdb for $pdb_fn $ecod_domain_id/$uid\n";
            $missing++;
        }
        else {
            push( @pdbs, $pdb_fn );
        }
        $total++;
    }
    printf "#%i missing pdbs of %i total\n", $missing, $total;
    open( OUT, ">pdb_lst.tmp.fn" );
    foreach my $pdb_fn (@pdbs) {
        print OUT "$pdb_fn\n";
    }
    close OUT;

    my $tar_cmd = "tar -czf $output_tarball -T pdb_lst.tmp.fn\n";
    system($tar_cmd);

}

sub determine_update_step {
    my ($update_node) = @_;
    my $STEP;
    if (  !has_premerge_ecod_fn($update_node)
        && has_manual_table_fn($update_node)
        && -f get_manual_table_fn($update_node) )
    {
        $STEP = 1;
    }
    elsif ( !has_merge_fn($update_node)
        && has_premerge_ecod_fn($update_node)
        && -f get_premerge_ecod_fn($update_node) )
    {
        $STEP = 2;
    }
    elsif ( !has_ecod_fn($update_node) ) {
        if (   has_previous_ecod_fn($update_node)
            && has_premerge_ecod_fn($update_node)
            && has_merge_fn($update_node)
			&& !has_ecod_tmp_fn($update_node)
            && -f get_previous_ecod_fn($update_node)
            && -f get_premerge_ecod_fn($update_node)
            && -f get_merge_fn($update_node) )
        {
            $STEP = 3;
        }
        elsif ( has_ecod_tmp_fn($update_node)
            && -f get_ecod_tmp_fn($update_node) )
        {
            $STEP = 4;
        }
        else {
            die "ERROR! Confusion case 1\n";
        }
    }
    elsif ( has_ecod_fn($update_node) && -f get_ecod_fn($update_node) ) {
        my $ecod_node = get_ecod_xml_node($update_node);
        if ( !has_chainwise_ecod_fn($update_node) ) {
            $STEP = 5;
        }
        elsif ( !has_derived_libraries($update_node) ) {
            if ( has_chainwise_ecod_fn($update_node) && -f get_chainwise_ecod_fn($update_node) ) {
                if ( overlaps_cleaned($ecod_node)
                    && !ligands_annotated($ecod_node) )
                {
                    $STEP = 6;
                }
                elsif ( ligands_annotated($ecod_node) ) {
                    $STEP = 7;
                }
                else {
                    die "ERROR! Confusion case 1a\n";
                }
            }
            else {
                die "ERROR! Confusion case 2\n";
            }
        }
        elsif ( has_derived_libraries($update_node) ) {
            my $lib_node = get_derived_libraries_node($update_node);
            if (  !pfam_clustered($ecod_node)
                && has_update_hmmer_db($lib_node)
                && -f get_update_hmmer_db($lib_node) )
            {
                $STEP = 8;
            }
            elsif ( pfam_clustered($ecod_node) && !version_registered($update_node) ) {
                $STEP = 9;
            }
            elsif ( version_registered($update_node)
                && !divergence_calculated($update_node) )
            {
                $STEP = 10;
            }
            elsif ( divergence_calculated($update_node)
                && !divergence_repaired($update_node) )
            {
                $STEP = 11;
            }
            elsif ( divergence_repaired($update_node)
                && !representatives_calculated($update_node) )
            {
                $STEP = 12;
            }
            elsif (divergence_repaired($update_node)
                && representatives_calculated($update_node) )
            {
                $STEP = 13;
            }
            elsif ( statistics_updated($update_node) ) {
                $STEP = 14;
            }
            else {
                die "ERROR! Confusion case 3\n";
            }
        }
    }
    else {
        die "ERROR! Could not determine update step, aborting...\n";
    }

    return $STEP;
}

sub update_dali_summ_cache {

    my $sub = 'update_dali_summ_cache';
    my ( $old_rep_cache_fn, $ecod_xml_doc, $version ) = @_;

    if ( !-f $old_rep_cache_fn ) {
        die "ERROR! Previous rep cache not found: $old_rep_cache_fn\n";
    }

    my $rep_cache = retrieve($old_rep_cache_fn);

    printf "#Found %i items in previous dali_summ rep_cache\n", scalar keys %$rep_cache;

    my $man_uids = get_manual_domain_ids($ecod_xml_doc);
    printf "#Found %i man_reps in ecod xml\n", scalar keys %$man_uids;

    my %domains;
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
        $domains{$uid} = $domain_node;
    }

    foreach my $uid ( keys %$man_uids ) {
        my ( $pdb_id, $chain_id ) = get_pdb_chain( $domains{$uid} );
        next if $chain_id eq '.';
        if ( !exists $$rep_cache{$uid} ) {
            add_to_rep_cache( $uid, $rep_cache, \%domains );
        }
    }
    foreach my $uid ( keys %$rep_cache ) {
        if ( !exists $$man_uids{$uid} ) {
            drop_from_rep_cache( $uid, $rep_cache );
        }
    }

    my $cache_fn = "$version.dali_summ.rep_cache";
    if ( !nstore( $rep_cache, $cache_fn ) ) {
        die "ERROR! Could not store rep_stats in $cache_fn\n";
    }

}

sub drop_from_rep_cache {
    my ( $uid, $rep_cache ) = @_;

    delete $$rep_cache{$uid};
}

sub add_to_rep_cache {
    my ( $uid, $rep_cache, $domain_index ) = @_;

    my $rep_domain_node = $$domain_index{$uid};

    my $seqid_range = get_seqid_range($rep_domain_node);

    my ( $ecod_pdb, $ecod_chain ) = get_pdb_chain($rep_domain_node);
    if ( $ecod_chain eq '.' ) { next }    #This is hacky and should be fixed 1/11/2013

    my $ecod_domain_id = get_ecod_domain_id( $$domain_index{$uid} );

    #my $ecod_pdb = $1;
    my $short_uid = substr( $uid, 2, 5 );
    if ( !-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb" ) {
        print "WARNING! No clean domain file for $ecod_domain_id\n";
    }
    else {
        #push (@ecod_reps, $ecod_domain_id);
        #$ecod_rep_range{$ecod_domain_id} = $seqid_range;
        $$rep_cache{$uid}{range}          = $seqid_range;
        $$rep_cache{$uid}{ecod_domain_id} = $ecod_domain_id;
    }
    my ( $ecod_seqid_aref, $ecod_struct_seqid_aref, $ecod_pdbnum_href, $ecod_asym_id ) =
      pdbml_seq_parse( $ecod_pdb, $ecod_chain );
    if ( !$ecod_seqid_aref ) {
        print "WARNING! $ecod_pdb, $ecod_chain obsolete in reps\n";
        next;
    }
    print "$ecod_pdb, $ecod_chain, $ecod_domain_id\n";
    my ($ecod_seq) = pdbml_seq_fetch( $ecod_pdb, $ecod_asym_id, $ecod_chain, $ecod_seqid_aref );
    if ($ecod_seqid_aref) {
        $$rep_cache{$uid}{seqid_aref}        = $ecod_seqid_aref;
        $$rep_cache{$uid}{struct_seqid_aref} = $ecod_struct_seqid_aref;
        $$rep_cache{$uid}{pdbnum_href}       = $ecod_pdbnum_href;
        $$rep_cache{$uid}{asym_id}           = $ecod_asym_id;
        $$rep_cache{$uid}{seq}               = $ecod_seq;
        $$rep_cache{$uid}{pdb}               = $ecod_pdb;
        $$rep_cache{$uid}{chain}             = $ecod_chain;
    }
}

sub generate_dali_summ_cache {
    my $sub = 'generate_dali_summ_cache';
    my ( $dali_summ_cache_fn, $ecod_xml_doc ) = @_;

    my $rep_domains = $ecod_xml_doc->findnodes(qq{//domain[\@manual_rep="true"]});

    my @ecod_reps;
    my %rep_stats;

    foreach my $rep_domain_node ( $rep_domains->get_nodelist() ) {
        my ( $uid, $ecod_domain_id ) = get_ids($rep_domain_node);
        my $short_uid = substr( $uid, 2, 5 );

        print "$uid $ecod_domain_id\n";

        next if isObsolete($rep_domain_node);

        my $seqid_range = get_seqid_range($rep_domain_node);

        my ( $ecod_pdb, $ecod_chain ) = get_pdb_chain($rep_domain_node);
        if ( $ecod_chain eq '.' ) { next }    #This is hacky and should be fixed 1/11/2013

        if ( $ecod_domain_id =~ /e(\w(\w{2})\w)/ ) {

            #my $ecod_pdb = $1;
            my $ecod_two = $2;
            if ( !-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb" ) {
                print "WARNING! No clean domain file for $ecod_domain_id\n";
            }
            else {
                push( @ecod_reps, $ecod_domain_id );

                #$ecod_rep_range{$ecod_domain_id} = $seqid_range;
                $rep_stats{$uid}{range} = $seqid_range;
            }
            my ( $ecod_seqid_aref, $ecod_struct_seqid_aref, $ecod_pdbnum_href, $ecod_asym_id ) =
              pdbml_seq_parse( $ecod_pdb, $ecod_chain );
            if ( !$ecod_seqid_aref ) {
                print "WARNING! $ecod_pdb, $ecod_chain obsolete in reps\n";
                next;
            }
            print "$ecod_pdb, $ecod_chain, $ecod_domain_id\n";
            my ($ecod_seq) = pdbml_seq_fetch( $ecod_pdb, $ecod_asym_id, $ecod_chain, $ecod_seqid_aref );
            if ($ecod_seqid_aref) {
                $rep_stats{$uid}{seqid_aref}        = $ecod_seqid_aref;
                $rep_stats{$uid}{struct_seqid_aref} = $ecod_struct_seqid_aref;
                $rep_stats{$uid}{pdbnum_href}       = $ecod_pdbnum_href;
                $rep_stats{$uid}{asym_id}           = $ecod_asym_id;
                $rep_stats{$uid}{seq}               = $ecod_seq;
                $rep_stats{$uid}{ecod_domain_id}    = $ecod_domain_id;
            }
        }
    }

    if ( !nstore( \%rep_stats, $dali_summ_cache_fn ) ) {
        die "ERROR! Could not store rep_stats in $dali_summ_cache_fn\n";
    }
}

sub calculate_divergence {
    my $sub              = 'calculate_divergence';
    my $current_version  = $LATEST_REFERENCE;
    my $previous_version = $PREVIOUS_VERSIONS{$current_version};

    print "cV: $current_version\n";

    my $current_pdb_doc_fn;
    if ( $PDB_XML_DOC{$current_version} && -f $PDB_XML_DOC{$current_version} ) {
        print "Current pdb doc found\n";
        $current_pdb_doc_fn = $PDB_XML_DOC{$current_version};
    }
    elsif ( $PDB_XML_DOC{$previous_version} && -f $PDB_XML_DOC{$previous_version} ) {
        print "Updating from previous pdb doc\n";
        open( my $xml_fh, $PDB_XML_DOC{$previous_version} )
          or die "ERROR! Could not open $PDB_XML_DOC{$previous_version} for reading:$!\n";
        my $previous_pdb_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );
        update_pdb_doc( $previous_pdb_xml_doc, $current_version );
        $current_pdb_doc_fn = $PDB_XML_DOC{$current_version};
    }
    else {
        print "No previous found, creating new\n";
        create_pdb_doc($current_version);
        $current_pdb_doc_fn = $PDB_XML_DOC{$current_version};
    }

    #Compare

    my $ecod_xml_fn = $REF_XML{$current_version};
    open( my $xml_fh1, "<", $ecod_xml_fn ) or die "ERROR! Could not open $ecod_xml_fn for reading:$!\n";

    #my $ecod_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );
    my $ecod_pull_parse = XML::LibXML::Reader->new( IO => $xml_fh1 );

    open( my $xml_fh2, "<", $current_pdb_doc_fn ) or die "ERROR! Could not open $current_pdb_doc_fn for reading:$!\n";

    #my $pdb_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );
    my $pdb_pull_parse = XML::LibXML::Reader->new( IO => $xml_fh2 );

    my $repair_job_xml = compare_pdb_xml_doc_to_ecod( $pdb_pull_parse, $ecod_pull_parse, $current_version );
    return $repair_job_xml;

}

sub update_pdb_doc {
    my ( $pdb_chain_doc, $current_version ) = @_;

    my %old_pdbs;
    foreach my $pdb_node ( $pdb_chain_doc->findnodes('/pdb_chain_document/pdb_list/pdb')->get_nodelist() ) {
        my $pdb_id = $pdb_node->findvalue('@pdb_id');
        $old_pdbs{$pdb_id}++;
    }
    printf "Found %i old pdbs\n", scalar keys %old_pdbs;

    my $local_dir = $LOCAL_DIR{$current_version};

    my $pdb_top_dir = '/usr2/pdb/data/structures/divided/XML-noatom';
    if ( !-d $pdb_top_dir ) { die "ERROR! $pdb_top_dir not found\n" }
    my $pdb_list_node = $pdb_chain_doc->findnodes('//pdb_list')->get_node(1);
    my @pdbs          = glob("$pdb_top_dir/??/*.xml.gz");

    #my $pdb_large_structures_top_dir = '/usr2/pdb/data/large_structures/XML-noatom';

    printf "Found %i total pdbs\n", scalar(@pdbs);

    foreach my $pdb_ln (@pdbs) {

        $pdb_ln =~ /(\w{4})-noatom.xml/;
        my $pdb = $1;
        if ( $old_pdbs{$pdb} ) { next }
        my ( $chains, $chain_seqid_aref, $chain_struct_seqid_aref, $chain_pdbnum_href ) = pdbml_fetch_chains($pdb);

        my $pdb_node = $pdb_chain_doc->createElement('pdb');
        $pdb_node->setAttribute( 'pdb_id', $pdb );

        my ( $date, $method_aref ) = pdbml_date_method($pdb);
        my $method = join( ",", @$method_aref );

        $pdb_node->setAttribute( 'date',   $date );
        $pdb_node->setAttribute( 'method', $method );

        foreach my $chain (@$chains) {
            my $chain_node = $pdb_chain_doc->createElement('chain');
            $chain_node->setAttribute( 'chain_id', $chain );
            $pdb_node->appendChild($chain_node);

            my $seqid_range_node = $pdb_chain_doc->createElement('seqid_range');
            my $range_node       = $pdb_chain_doc->createElement('range');

            my ( $seqid_range, $range );
            if ( scalar( @{ $$chain_struct_seqid_aref{$chain} } ) == 0 ) {
                $seqid_range_node->setAttribute( 'unstructured', 'true' );
                $range_node->setAttribute( 'unstructured', 'true' );
            }
            else {
                $seqid_range = scopify_range( rangify( @{ $$chain_struct_seqid_aref{$chain} } ), $chain );
                $seqid_range_node->appendTextNode($seqid_range);
                $range = scopify_range( pdb_rangify( $$chain_struct_seqid_aref{$chain}, $$chain_pdbnum_href{$chain} ),
                    $chain );
                $range_node->appendTextNode($range);
            }

            $chain_node->appendChild($seqid_range_node);
            $chain_node->appendChild($range_node);
        }
        $pdb_list_node->appendChild($pdb_node);
    }

    my $doc_string = $pdb_chain_doc->toString(1);
    open( OUT, ">$local_dir/pdb.$current_version.xml" ) or die;
    print OUT $doc_string;
    close OUT;

    register_pdb_doc( "$local_dir/pdb.$current_version.xml", $current_version );
    load_references();

}

sub get_pdbs_by_latest_week {

    my $fn = '/usr2/pdb/data/status/latest/added.sf';
    !-f $fn && return 0;

    my %pdbs;
    open( my $fh, "<", $fn ) or die "ERROR! Could not open $fn for reading:$!\n";
    while ( my $ln = <$fh> ) {
        chomp $ln;
        $pdbs{$ln}++;
    }
    return \%pdbs;

}

sub compare_pdb_xml_doc_to_ecod {
    my ( $pdb_pull_parse, $ecod_pull_parse, $version ) = @_;

    my %chains;
    my %skip;

    my ( $compare_repair_xml_doc, $compare_repair_top_node ) = xml_create('job_set_top');

    my $local_dir = $LOCAL_DIR{$version};

    my $repair_dir = "$local_dir/repair.$version";
    if ( !-d $repair_dir ) {
        mkdir($repair_dir);
    }

    my $dump_dir = "$repair_dir/ecod_dump";
    if ( !-d $dump_dir ) {
        mkdir($dump_dir);
    }

    my $latest_week_pdbs_href = get_pdbs_by_latest_week();

    my $job_list_dir_node = $compare_repair_xml_doc->createElement('job_list_dir');
    $job_list_dir_node->appendTextNode($repair_dir);

    my $job_dump_dir_node = $compare_repair_xml_doc->createElement('job_dump_dir');
    $job_dump_dir_node->appendTextNode($dump_dir);
    my $job_list_node = $compare_repair_xml_doc->createElement('job_list');

    $compare_repair_top_node->appendChild($job_list_dir_node);
    $compare_repair_top_node->appendChild($job_dump_dir_node);
    $compare_repair_top_node->appendChild($job_list_node);

    print "Index ECOD domains...\n";
    while ( $ecod_pull_parse->read() ) {

        my $node_name = $ecod_pull_parse->name;
        my $depth     = $ecod_pull_parse->depth;
        my $nodeType  = $ecod_pull_parse->nodeType;
        my $isempty   = $ecod_pull_parse->isEmptyElement;

        #print "$node_name $depth $nodeType $isempty\n";

        if ( $nodeType == 1 ) {
            if ( $node_name eq 'domain' ) {
                my $domain_node = $ecod_pull_parse->copyCurrentNode(1);

                my $pdb = $domain_node->findvalue('structure/@pdb_id');

                #my $chain 	= $domain_node->findvalue('structure/@chain_id');
                my ( $uid, $ecod_domain_id ) = get_ids($domain_node);

                my $seqid_range = get_seqid_range($domain_node);
                if ( $seqid_range =~ /\d+/ ) {
                    my @segs = split( /,/, $seqid_range );
                    foreach my $seg (@segs) {
                        my ( $seg_str, $chain_str ) = scop_range_split($seg);
                        my $seg_aref = range_expand($seg_str);
                        range_include( $seg_aref, \@{ $chains{$pdb}{$chain_str} } );
                    }
                }
                else {
                    print "WARNING! $ecod_domain_id has non-numeric range\n";
                }
            }
        }
    }

    print "Index ECOD special architectures...\n";
    my $dbh = db_connect('ecod');
    my $sth = $dbh->prepare('SELECT uid, seqid_range, pdb, chain, type FROM special');
    $sth->execute();
    while ( my $row_aref = $sth->fetchrow_arrayref() ) {
        my $uid         = $$row_aref[0];
        my $seqid_range = $$row_aref[1];
        my $pdb         = $$row_aref[2];
        my $chain       = $$row_aref[3];
        my $type        = $$row_aref[4];
        print "DEBUG: $uid $pdb $chain $seqid_range $type\n";
        if ( $seqid_range =~ /\d+\-/ && $type ne 'nonpeptide_poly' ) {
            my @segs = split( /,/, $seqid_range );
            foreach my $seg (@segs) {
                my ( $seg_str, $chain_str ) = scop_range_split($seg);
                my $seg_aref = range_expand($seg_str);
                range_include( $seg_aref, \@{ $chains{$pdb}{$chain_str} } );
            }
        }
        else {
            print "WARNING $uid has non-numeric range\n";
        }
    }

    print "Compare to pdb\n";
    my $job_id = 1;
    while ( $pdb_pull_parse->read() ) {
        my $node_name = $pdb_pull_parse->name();
        if ( $node_name eq 'pdb' ) {

            my $pdb_node = $pdb_pull_parse->copyCurrentNode(1);

            my $pdb_id = $pdb_node->findvalue('@pdb_id');
            next if $$latest_week_pdbs_href{$pdb_id};    #Skip this pdb if it appears in the latest PDB release.

            foreach my $chain_node ( $pdb_node->findnodes('chain')->get_nodelist() ) {
                my $chain_id    = $chain_node->findvalue('@chain_id');
                my $seqid_range = $chain_node->findvalue('seqid_range');
                if ( $skip{$pdb_id}{$chain_id} ) { next }
                if ( !$chains{$pdb_id}{$chain_id} && $seqid_range =~ /\d+\-/ ) {

                    #print "MISSING $pdb_id $chain_id 0 \n";
                    my $job_node = $compare_repair_xml_doc->createElement('job');
                    $job_node->appendTextChild( 'query_pdb',   $pdb_id );
                    $job_node->appendTextChild( 'query_chain', $chain_id );
                    $job_node->appendTextChild( 'reference',   $version );
                    $job_node->setAttribute( 'id', $job_id++ );

                    $job_node->appendTextChild( 'repair_type', 'missing' );

                    my $mode_node = $compare_repair_xml_doc->createElement('mode');
                    $mode_node->setAttribute( 'input', 'struct_seqid' );
                    $mode_node->appendTextNode('seq_iter');
                    $job_node->appendChild($mode_node);
                    $job_list_node->appendChild($job_node);

                    next;
                }

                #print "$pdb_id $chain_id $seqid_range\n";
                if ( $seqid_range =~ /\d+\-/ ) {
                    my $range_aref = range_expand($seqid_range);

                    if ( ref( $chains{$pdb_id}{$chain_id} ) ne 'ARRAY' ) {
                        die "ERROR! $pdb_id $chain_id has non-array ref range\n";
                    }

                    my $c1 = region_coverage( $range_aref, $chains{$pdb_id}{$chain_id} );
                    my $c2 = region_coverage( $chains{$pdb_id}{$chain_id}, $range_aref );

                    my $residue_coverage = residue_coverage( $chains{$pdb_id}{$chain_id}, $range_aref );

                    #my $total_residues = scalar(@{$chains{$pdb_id}{$chain_id}});
                    my $total_residues = scalar(@$range_aref);

                    #print "COVERAGE $pdb_id $chain_id $c1 $c2 $residue_coverage $total_residues\n";

                    if ( $total_residues - $residue_coverage >= 30 ) {
                        my $job_node = $compare_repair_xml_doc->createElement('job');
                        $job_node->appendTextChild( 'query_pdb',   $pdb_id );
                        $job_node->appendTextChild( 'query_chain', $chain_id );
                        $job_node->appendTextChild( 'reference',   $version );
                        $job_node->setAttribute( 'id', $job_id++ );

                        $job_node->appendTextChild( 'repair_type', 'low_coverage' );
                        $job_node->setAttribute( 'missing_res', $total_residues - $residue_coverage );

                        my $mode_node = $compare_repair_xml_doc->createElement('mode');
                        $mode_node->setAttribute( 'input', 'struct_seqid' );
                        $mode_node->appendTextNode('seq_iter');
                        $job_node->appendChild($mode_node);
                        $job_list_node->appendChild($job_node);

                    }
                }
            }
        }
    }

    my $out_fn     = "$repair_dir/repair.$version.job.xml";
    my $doc_string = $compare_repair_xml_doc->toString(1);
    open( OUT, ">", $out_fn ) or die "ERROR! Could not open $out_fn for writing:$!\n";
    print OUT $doc_string;
    close OUT;

    return $out_fn;
}

sub compare_pdb_xml_doc_to_ecod_obsolete {
    my ( $pdb_xml_doc, $ecod_xml_doc, $version ) = @_;

    my %chains;
    my $compare_repair_xml_doc  = XML::LibXML->createDocument();
    my $compare_repair_top_node = $compare_repair_xml_doc->createElement('job_set_top');
    $compare_repair_xml_doc->setDocumentElement($compare_repair_top_node);

    my $local_dir = $LOCAL_DIR{$version};

    my $repair_dir = "$local_dir/repair.$version";
    if ( !-d $repair_dir ) {
        mkdir($repair_dir);
    }
    my $out_fn = "$repair_dir/repair.$version.job.xml";

    my $dump_dir = "$repair_dir/ecod_dump";
    if ( !-d $dump_dir ) {
        mkdir($dump_dir);
    }

    my $job_list_dir_node = $compare_repair_xml_doc->createElement('job_list_dir');
    $job_list_dir_node->appendTextNode($repair_dir);
    my $job_dump_dir_node = $compare_repair_xml_doc->createElement('job_dump_dir');
    $job_dump_dir_node->appendTextNode($dump_dir);
    my $job_list_node = $compare_repair_xml_doc->createElement('job_list');
    $compare_repair_top_node->appendChild($job_list_dir_node);
    $compare_repair_top_node->appendChild($job_dump_dir_node);
    $compare_repair_top_node->appendChild($job_list_node);

    print "Index existing domains...\n";
    foreach my $domain_node (
        $ecod_xml_doc->findnodes('/ecod_document/domain_dictionary/architecture/x_group/h_group/f_group//domain')
        ->get_nodelist() )
    {

        my $pdb = $domain_node->findvalue('structure/@pdb_id');

        #my $chain 	= $domain_node->findvalue('structure/@chain_id');
        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');

        my $seqid_range;
        if ( $domain_node->exists('seqid_range') ) {
            $seqid_range = $domain_node->findvalue('seqid_range');
        }
        elsif ( $domain_node->exists('derived_seqid_range') ) {
            $seqid_range = $domain_node->findvalue('derived_seqid_range');
        }
        else {
            print "WARNING! No range for $ecod_domain_id, skipping...\n";
            next;
        }
        if ( $seqid_range =~ /\d+/ ) {
            my @segs = split( /,/, $seqid_range );
            foreach my $seg (@segs) {
                my ( $seg_str, $chain_str ) = scop_range_split($seg);
                my $seg_aref = range_expand($seg_str);
                range_include( $seg_aref, \@{ $chains{$pdb}{$chain_str} } );
            }
        }
        else {
            print "WARNING! $ecod_domain_id has non-numeric range\n";
        }
    }

    print "Index existing peptides\n";
    foreach my $peptide_node ( $ecod_xml_doc->findnodes('/ecod_document/peptide_list/peptide')->get_nodelist() ) {

        my $pdb   = $peptide_node->findvalue('@pdb_id');
        my $chain = $peptide_node->findvalue('@chain_id');

        if ( $peptide_node->exists('seqid_range') ) {
            my $seqid_range = $peptide_node->findvalue('seqid_range');

            if ( $seqid_range =~ /\d+/ ) {
                my $range_aref = range_expand($seqid_range);
                range_include( $range_aref, \@{ $chains{$pdb}{$chain} } );
            }
            else {
                print "WARNING! Non-numeric peptide range for $pdb $chain\n";
            }
        }
        else {
            print "WARNING! No seqid_range for peptide $pdb $chain\n";
        }
    }

    print "Indexing existing pss\n";
    foreach my $pss_node ( $ecod_xml_doc->findnodes('/ecod_document/poor_sequence_structure_list/pss')->get_nodelist() )
    {

        my $pdb   = $pss_node->findvalue('@pdb_id');
        my $chain = $pss_node->findvalue('@chain_id');

        if ( $pss_node->exists('seqid_range') ) {
            my $seqid_range = $pss_node->findvalue('seqid_range');
            if ( $seqid_range =~ /\d+/ ) {
                my $range_aref = range_expand($seqid_range);
                range_include( $range_aref, \@{ $chains{$pdb}{$chain} } );
            }
        }
        else {
            print "WARNING! No seqid_range for pss $pdb $chain\n";
        }

    }

    print "Indexing existing coiled-coils\n";
    foreach my $coiled_coil_node ( $ecod_xml_doc->findnodes('/ecod_document/coil_list/coil')->get_nodelist() ) {

        my $pdb   = $coiled_coil_node->findvalue('@pdb_id');
        my $chain = $coiled_coil_node->findvalue('@chain_id');

        if ( $coiled_coil_node->exists('seqid_range') ) {
            my $seqid_range = $coiled_coil_node->findvalue('seqid_range');
            if ( $seqid_range =~ /\d+/ ) {
                my $range_aref = range_expand($seqid_range);
                range_include( $range_aref, \@{ $chains{$pdb}{$chain} } );
            }
        }
        else {
            print "WARNING! NO seqid range for coil $pdb $chain\n";
        }
    }

    print "Indexing existing synth\n";
    foreach
      my $synth_node ( $ecod_xml_doc->findnodes('/ecod_document/synthetic_peptide_list/synthetic')->get_nodelist() )
    {

        my $pdb   = $synth_node->findvalue('@pdb_id');
        my $chain = $synth_node->findvalue('@chain_id');

        if ( $synth_node->exists('seqid_range') ) {
            my $seqid_range = $synth_node->findvalue('seqid_range');
            if ( $seqid_range =~ /\d+/ ) {
                my $range_aref = range_expand($seqid_range);
                range_include( $range_aref, \@{ $chains{$pdb}{$chain} } );
            }
        }
        else {
            print "WARNING! No seqid_range for synth $pdb $chain\n";
        }

    }

    #Test

    print "Compare to pdb\n";
    my $job_id = 1;
    foreach my $pdb_node ( $pdb_xml_doc->findnodes('/pdb_chain_document/pdb_list/pdb')->get_nodelist() ) {
        my $pdb_id = $pdb_node->findvalue('@pdb_id');

        foreach my $chain_node ( $pdb_node->findnodes('chain')->get_nodelist() ) {
            my $chain_id    = $chain_node->findvalue('@chain_id');
            my $seqid_range = $chain_node->findvalue('seqid_range');
            if ( !$chains{$pdb_id}{$chain_id} && $seqid_range =~ /\d+/ ) {
                print "MISSING $pdb_id $chain_id 0 \n";
                my $job_node = $compare_repair_xml_doc->createElement('job');
                $job_node->appendTextChild( 'query_pdb',   $pdb_id );
                $job_node->appendTextChild( 'query_chain', $chain_id );
                $job_node->appendTextChild( 'reference',   $version );
                $job_node->setAttribute( 'id', $job_id++ );

                $job_node->appendTextChild( 'repair_type', 'missing' );

                my $mode_node = $compare_repair_xml_doc->createElement('mode');
                $mode_node->setAttribute( 'input', 'struct_seqid' );
                $mode_node->appendTextNode('seq_iter');
                $job_node->appendChild($mode_node);
                $job_list_node->appendChild($job_node);

                next;
            }

            #print "$pdb_id $chain_id $seqid_range\n";
            if ( $seqid_range =~ /\d+/ ) {
                my $range_aref = range_expand($seqid_range);

                if ( ref( $chains{$pdb_id}{$chain_id} ) ne 'ARRAY' ) {
                    die "ERROR! $pdb_id $chain_id has non-array ref range\n";
                }

                my $c1 = region_coverage( $range_aref, $chains{$pdb_id}{$chain_id} );
                my $c2 = region_coverage( $chains{$pdb_id}{$chain_id}, $range_aref );

                my $residue_coverage = residue_coverage( $chains{$pdb_id}{$chain_id}, $range_aref );

                #my $total_residues = scalar(@{$chains{$pdb_id}{$chain_id}}); #what
                my $total_residues = scalar(@$range_aref);
                print "COVERAGE $pdb_id $chain_id $c1 $c2 $residue_coverage $total_residues\n";
                if ( $total_residues - $residue_coverage >= 30 ) {
                    my $job_node = $compare_repair_xml_doc->createElement('job');
                    $job_node->appendTextChild( 'query_pdb',   $pdb_id );
                    $job_node->appendTextChild( 'query_chain', $chain_id );
                    $job_node->appendTextChild( 'reference',   $version );
                    $job_node->setAttribute( 'id', $job_id++ );

                    $job_node->appendTextChild( 'repair_type', 'low_coverage' );
                    $job_node->setAttribute( 'missing_res', $total_residues - $residue_coverage );

                    my $mode_node = $compare_repair_xml_doc->createElement('mode');
                    $mode_node->setAttribute( 'input', 'struct_seqid' );
                    $mode_node->appendTextNode('seq_iter');
                    $job_node->appendChild($mode_node);
                    $job_list_node->appendChild($job_node);

                }
            }

        }
    }

    #my $out_fn = "$repair_dir/repair.$version.job.xml";
    my $doc_string = $compare_repair_xml_doc->toString(1);
    open( OUT, ">", $out_fn ) or die "ERROR! Could not open $out_fn for writing:$!\n";
    print OUT $doc_string;
    close OUT;

    return $out_fn;
}

sub create_pdb_doc {
    my ($current_version) = @_;
    my $pdb_chain_doc     = XML::LibXML->createDocument();
    my $pdb_root_node     = $pdb_chain_doc->createElement('pdb_chain_document');
    $pdb_chain_doc->setDocumentElement($pdb_root_node);

    my $pdb_list_node = $pdb_chain_doc->createElement('pdb_list');
    $pdb_root_node->appendChild($pdb_list_node);

    my $local_dir = $LOCAL_DIR{$current_version};

    my $pdb_top_dir = '/usr2/pdb/data/structures/divided/XML-noatom';
    my @pdbs        = glob("$pdb_top_dir/??/*.xml.gz");

    printf "Found %i pdbs\n", scalar(@pdbs);

    my $date         = `date`;
    my $created_node = $pdb_chain_doc->createElement('createdOn');
    $created_node->appendTextNode($date);
    $pdb_root_node->appendChild($created_node);

    foreach my $pdb_ln (@pdbs) {

        $pdb_ln =~ /(\w{4})-noatom.xml/;
        my $pdb = $1;
        my ( $chains, $chain_seqid_aref, $chain_struct_seqid_aref, $chain_pdbnum_href ) = pdbml_fetch_chains($pdb);

        my $pdb_node = $pdb_chain_doc->createElement('pdb');
        $pdb_node->setAttribute( 'pdb_id', $pdb );
        foreach my $chain (@$chains) {
            my $chain_node = $pdb_chain_doc->createElement('chain');
            $chain_node->setAttribute( 'chain_id', $chain );
            $pdb_node->appendChild($chain_node);
            my $seqid_range = scopify_range( rangify( @{ $$chain_struct_seqid_aref{$chain} } ), $chain );
            my $range =
              scopify_range( pdb_rangify( $$chain_struct_seqid_aref{$chain}, $$chain_pdbnum_href{$chain} ), $chain );

            my $seqid_range_node = $pdb_chain_doc->createElement('seqid_range');
            $seqid_range_node->appendTextNode($seqid_range);
            my $range_node = $pdb_chain_doc->createElement('range');
            $range_node->appendTextNode($range);

            if ( scalar( @{ $$chain_struct_seqid_aref{$chain} } ) == 0 ) {
                $seqid_range_node->setAttribute( 'unstructured', 'true' );
                $range_node->setAttribute( 'unstructured', 'true' );
            }

            $chain_node->appendChild($seqid_range_node);
            $chain_node->appendChild($range_node);
        }
        $pdb_list_node->appendChild($pdb_node);
    }

    my $doc_string = $pdb_chain_doc->toString(1);
    open( OUT, ">$local_dir/pdb.$current_version.xml" ) or die;
    print OUT $doc_string;
    close OUT;

    register_pdb_doc( "$local_dir/pdb.$current_version.xml", $current_version );
    load_references();

}

sub generate_domain_data_directories {
    my $sub = 'generate_domain_data_directories';

    my ($ecod_xml_doc) = @_;

    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
        my $short_uid = substr( $uid, 2, 5 );

        if ( !-d "$DOMAIN_DATA_DIR" ) {
            mkdir "$DOMAIN_DATA_DIR";
        }
        if ( !-d "$DOMAIN_DATA_DIR/$short_uid" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid");
        }
        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/$uid/" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid/$uid");
        }
    }
}

sub generate_ligand_annotation_jobs {
    my $sub = 'generate_ligand_annotation_jobs';

    my ($ecod_xml_doc) = @_;

    my @jobs;
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc))
    {
		
		my ($uid, $ecod_domain_id) = get_ids($domain_node);
        my $short_uid      = substr( $uid, 2, 5 );

		my ($pdb, $chain) = get_pdb_chain($domain_node);

		next if isObsolete($domain_node);
		next unless has_ligand_annotation($domain_node);

        my $domain_seqid_range = get_seqid_range($domain_node);

        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/$uid/" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid/$uid");
        }
        if ( -f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.ligand_contact.txt" ) { next }

        #my $ligand_href = find_ligands_by_range($pdb, $domain_seqid_range);
        my $job = "$SINGLE_LIGAND_ANNOT_EXE $uid $pdb $domain_seqid_range\n";
        job_create( "ligand.$uid.job", $job );
        push( @jobs, "ligand.$uid.job" );

    }
    return \@jobs;
}

sub ligand_annotation {
    my $sub            = 'ligand_annotation';
    my ($ecod_xml_doc) = @_;
    my $count          = 0;
  DOMAIN:
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {

		if (has_ligand_annotation($domain_node) && !$FORCE_OVERWRITE ) { next } 

        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);

        my $short_uid = substr( $uid, 2, 5 );
        my ($pdb, $chain) = get_pdb_chain($domain_node);

        my $domain_path = "$DOMAIN_DATA_DIR/$short_uid/$uid";
        if ( !-d $domain_path ) {
            print "WARNING! domain path $domain_path not found for $uid $ecod_domain_id\n";
            next;
        }

        my $ligand_contact_fn = "$domain_path/$uid.ligand_contact.txt";
        if ( -f $ligand_contact_fn ) {
            open( IN, $ligand_contact_fn ) or die "ERROR! Could not open ligand contact file: $ligand_contact_fn\n";
            my $found = 0;
            while ( my $ln = <IN> ) {
                my @F = split( /\s+/, $ln );
                if ( scalar(@F) == 5 ) {
                    my $file_uid = $F[0];
                    if ( $file_uid != $uid ) {
                        print "WARNING! $file_uid != $uid in $ligand_contact_fn\n";
                        print "rm $ligand_contact_fn\n";
                        next DOMAIN;
                    }
                    my $file_pdb = $F[1];
                    if ( $file_pdb ne $pdb ) {
                        print "WARNING! $file_pdb ne $pdb in $ligand_contact_fn\n";
                        print "rm $ligand_contact_fn\n";
                        next DOMAIN;
                    }

                    my $file_seqid_range       = $F[2];
                    my $file_ligand_str        = $F[3];
                    my $file_ligand_pdbnum_str = $F[4];

                    #print "$uid $file_ligand_str $file_ligand_pdbnum_str\n";
                    my $ligand_str = $ecod_xml_doc->createElement('ligand_str');
                    $ligand_str->appendTextNode($file_ligand_str);
                    $ligand_str->setAttribute( 'inclusion_threshold', $LIGAND_INCLUSION_THRESHOLD );
                    $domain_node->appendChild($ligand_str);

                    my $ligand_pdbnum_str = $ecod_xml_doc->createElement('ligand_pdbnum_str');
                    $ligand_pdbnum_str->appendTextNode($file_ligand_pdbnum_str);
                    $ligand_pdbnum_str->setAttribute( 'inclusion_threshold', $LIGAND_INCLUSION_THRESHOLD );
                    $domain_node->appendChild($ligand_pdbnum_str);
                    $found = 1;

                }
            }
        }
        else {
            my $ligand_str = $ecod_xml_doc->createElement('ligand_str');
            $ligand_str->setAttribute( 'xsi:nil', 'true' );
            $domain_node->appendChild($ligand_str);
            my $ligand_pdbnum_str = $ecod_xml_doc->createElement('ligand_pdbnum_str');
            $ligand_pdbnum_str->setAttribute( 'xsi:nil', 'true' );
            $domain_node->appendChild($ligand_pdbnum_str);
            print "WARNING! No ligand contact file for $uid $ecod_domain_id\n";
        }
        $count++;
    }
    printf "#%i domains ligand_annotated\n", $count;
}

sub generate_image_jobs {
    my $sub = "generate_image_jobs";
    my ($ecod_xml_doc) = @_;
    my @jobs;
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {
		my ($uid, $ecod_domain_id) = get_ids($domain_node);
        my $short_uid      = substr( $uid, 2, 5 );

		my ($pdb, $chain) = get_pdb_chain($domain_node);

        my $pdb_range;
		my $pdb_range = get_pdb_range($domain_node);
	  	my $run_str = "$PYMOL_EXE -qc $DOMAIN_IMAGE_EXE -- $uid $pdb $pdb_range ";

        if ( $domain_node->parentNode->nodeName eq 'domain_assembly' ) {
            my @sib_pdb_ranges;
            foreach my $sibling_domain_node ( $domain_node->parentNode->findnodes('domain') ) {
                my $sib_uid = $sibling_domain_node->findvalue('@uid');
                if ( $sib_uid eq $uid ) { next }
				my $sib_pdb_range = get_pdb_range($sibling_domain_node);
				push( @sib_pdb_ranges, $sib_pdb_range );
            }
            my $sib_range_str = join( ",", @sib_pdb_ranges );
            $run_str .= "$sib_range_str ";
        }

        if ( $domain_node->exists('ligand_pdbnum_str') ) {
            my $ligand_range_str = $domain_node->findvalue('ligand_pdbnum_str');
            $run_str .= "$ligand_range_str ";
        }

        if ( !-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.png" || $FORCE_OVERWRITE ) {
            my $job = "$run_str\n";
            job_create( "image_gen.$uid.job", $run_str );
            push( @jobs, "image_gen.$uid.job" );
        }

    }
    return \@jobs;
}

sub generate_pdb_structure_jobs {
    my $sub = "generate_pdb_structure_jobs";
    my ( $ecod_xml_doc, $mode ) = @_;

    my @jobs;
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc)) ) {

        #if ($domain_node->findvalue('structure/@structure_obsolete') eq 'true') { next }

        my $pdb       = $domain_node->findvalue('structure/@pdb_id');
		my ($uid, $ecod_domain_id) = get_ids($domain_node);
        my $short_uid = substr( $uid, 2, 5 );

        if ( !-d "$DOMAIN_DATA_DIR/$short_uid" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid");
        }
        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/$uid" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid/$uid");
        }

        $ecod_domain_id =~ /e\w(\w{2})/;
        my $two = $1;

		my $range = get_seqid_range($domain_node);

        my $ligand_range;
        if ( $INCLUDE_LIGANDS && $domain_node->exists('ligand_pdbnum_str') ) {
            $ligand_range = $domain_node->findvalue('ligand_pdbnum_str');
        }

        #my $MODE = 'seqres';
        #my $MODE = 'pdbnum';
        my $input_mode;
        if ( $mode eq 'pdbnum' ) {
            $input_mode = "pdb";
        }
        elsif ( $mode eq 'seqres' ) {
            $input_mode = "seq";
        }
        else {
            die "unk mode $mode\n";
        }
        my $cmd;
        my $pdb_fn;
        $pdb_fn = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.$mode.pdb";

        if ( !-f $pdb_fn || $FORCE_OVERWRITE ) {
            my $job;
            if ($ligand_range) {
                $job =
"$XML_PDB_EXE $pdb $range $DOMAIN_DATA_DIR/$short_uid/$uid $uid --ligand_range=$ligand_range --mode=$input_mode\n";
            }
            else {
                $job =
                  "$XML_PDB_EXE $pdb $range $DOMAIN_DATA_DIR/$short_uid/$uid $uid --ligand_range=0 --mode=$input_mode\n"
                  ;    #this is super hacky, plz fix
            }
            job_create( "gen_pdb.$mode.$uid.job", $job );
            push( @jobs, "gen_pdb.$mode.$uid.job" );
        }
    }
    return \@jobs;
}

sub convert_side_load_domains_to_xml {
    my $sub = 'convert_side_load_domains_to_xml';
    my ($side_load_fn) = @_;

    open( IN, $side_load_fn ) or die "ERROR! Could not open $side_load_fn for reading:$!\n";
    my $uid = 1;

    my @domains;
    my $i = 0;
    while ( my $ln = <IN> ) {

        my @F = split( /\t/, $ln );

        my $pdb            = lc( $F[0] );
        my $chain          = $F[1];
        my $range          = $F[2];
        my $hit_domain_uid = $F[3];
        my $hit_domain_id  = $F[4];
        if ( $hit_domain_id !~ /e\w{4}.\d+/ ) { next }
        my $reference = $F[5];
        chomp $reference;

        $domains[$i]{pdb}            = $pdb;
        $domains[$i]{chain}          = $chain;
        $domains[$i]{range}          = $range;
        $domains[$i]{hit_domain_id}  = $hit_domain_id;
        $domains[$i]{hit_domain_uid} = $hit_domain_uid;
        $domains[$i]{reference}      = $reference;

        $i++;

    }

    my $domain_xml_doc  = XML::LibXML->createDocument();
    my $domain_doc_node = $domain_xml_doc->createElement('domain_doc');
    $domain_xml_doc->setDocumentElement($domain_doc_node);

    my $domain_list_node = $domain_xml_doc->createElement('domain_list');
    $domain_doc_node->appendChild($domain_list_node);

    my $ref = $LATEST_REFERENCE;

    my %seen;
    for ( my $i = 0 ; $i < scalar(@domains) ; $i++ ) {

        my $domain_node = $domain_xml_doc->createElement('domain');

        $domain_node->setAttribute( 'yx_uid', sprintf "%09i", $uid++ );    #not a true uid
        $domain_node->setAttribute( 'manual_rep',   'true' );
        $domain_node->setAttribute( 'yx_side_load', 'true' );

        $domain_list_node->appendChild($domain_node);

        $ref = $domains[$i]{reference};
        my $pdb   = $domains[$i]{pdb};
        my $chain = $domains[$i]{chain};

        my $pdb_chain = lc($pdb) . $chain;
        $seen{$pdb_chain}++;
        my $ecod_domain_id = 'e' . lc($pdb) . $chain . $seen{$pdb_chain};

        my $manual_range = scopify_range( $domains[$i]{range}, $chain );
        my $hit_domain   = $domains[$i]{hit_domain_id};
        my $hit_uid      = $domains[$i]{hit_domain_uid};

        my $ecod_domain_node = $domain_xml_doc->createElement('ecod_domain');
        $ecod_domain_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
        $ecod_domain_node->setAttribute( 'reference',      $ref );
        $domain_node->setAttribute( 'manual_range', 'true' );

        my $structure_node = $domain_xml_doc->createElement('structure');
        $structure_node->setAttribute( 'db',                   'PDB' );
        $structure_node->setAttribute( 'pdb_id',               $pdb );
        $structure_node->setAttribute( 'chain_id',             $chain );
        $structure_node->setAttribute( 'chain_case_ambiguous', 'false' );

        my $manual_range_node = $domain_xml_doc->createElement('manual_range');
        $manual_range_node->appendTextNode($manual_range);

        $domain_node->appendChild($ecod_domain_node);
        $domain_node->appendChild($manual_range_node);
        $domain_node->appendChild($structure_node);

        if ( $hit_domain =~ /(e\w{4}.\d+)/ ) {
            my $hit_ecod_domain_id = $1;
            my $hit_domain_node    = $domain_xml_doc->createElement('ecod_representative_domain');
            $hit_domain_node->setAttribute( 'ecod_domain_id', $hit_ecod_domain_id );
            $hit_domain_node->setAttribute( 'uid',            sprintf "%09i", $hit_uid );
            $hit_domain_node->setAttribute( 'reference',      $ref );
            $domain_node->appendChild($hit_domain_node);
        }
        elsif ( $hit_domain =~ /(d\w{4}.(\d+|_))/ ) {
            my $hit_scop_domain_id = $1;
            my $hit_domain_node    = $domain_xml_doc->createElement('scop_representative_domain');
            $hit_domain_node->setAttribute( 'scop_domain_id', $hit_scop_domain_id );
            $hit_domain_node->setAttribute( 'reference',      $ref );
            $domain_node->appendChild($hit_domain_node);
        }
        else {
            die "ERROR! hit domain regexp fail on $hit_domain\n";
        }

    }

    return $domain_xml_doc;
}

sub process_side_load_domains_to_ecod {
    my $sub = 'process_side_load_domains_to_xml';

    my ( $yx_domain_xml_doc, $ecod_xml_doc ) = @_;
    my $domain_XPath = '//domain_list/domain';

    my $maxUID = $ecod_xml_doc->findvalue('//@maxUID');

    my $dict_node = $ecod_xml_doc->findnodes('//domain_dictionary')->get_node(1);

    my %seen;
    my %dict_ranges;
    my %chain_ranges;
    my %yx_side_load;

    foreach my $ecod_domain_node ( find_domain_nodes($ecod_xml_doc)) ) {
        my $ecod_domain_id = $ecod_domain_node->findvalue('@ecod_domain_id');
        my $seqid_range = get_seqid_range($ecod_domain_node);

        $seen{$ecod_domain_id}++;

        if ( $ecod_domain_node->findvalue('@yx_side_load') eq 'true' ) {
            $yx_side_load{$ecod_domain_id}++;
        }

		my ($pdb_id, $chain_id) = get_pdb_chain($ecod_domain_node);

        push( @{ $chain_ranges{$pdb_id}{$chain_id} }, $seqid_range );

    }
    my $size = $yx_domain_xml_doc->findnodes($domain_XPath)->size();

    #printf "DEBUG $sub: size $size\n";
    if ( $size > 0 ) {
        my $uid = $ecod_xml_doc->findvalue('//@maxUID');
      YX_DOMAIN:
        foreach my $yx_domain_node ( find_domain_nodes($yx_domain_xml_doc) ) {

            my $ecod_domain_id = $yx_domain_node->findvalue('ecod_domain/@ecod_domain_id');
            if ( $yx_side_load{$ecod_domain_id} ) { next }

			my ($pdb_id, $chain_id) = get_pdb_chain($domain_node);

            #my $uid = $yx_domain_node->findvalue('@uid');

            #my $seqid_range = $yx_domain_node->findvalue('seqid_range');
            my $seqid_range      = $yx_domain_node->findvalue('manual_range');    #Really?
            my $seqid_range_aref = range_expand($seqid_range);

            foreach my $ref_range_str ( @{ $chain_ranges{$pdb_id}{$chain_id} } ) {
                my ( $ref_range, $ref_chain ) = scop_range_split($ref_range_str);
                my $ref_range_aref = range_expand($ref_range);
                print "DEBUG: 1:$ref_range_str 2:$seqid_range\n";
                if ( residue_coverage( $ref_range_aref, $seqid_range_aref ) > 25 ) {
                    print "\tBAIL OUT\n";
                    next YX_DOMAIN;
                }
            }

            while ( $seen{$ecod_domain_id} ) {
                print "WARNING! ID collision on domain id $ecod_domain_id\n";
                $ecod_domain_id =~ /e(\w{4})(.)(\d+)/;
                my $pdb   = $1;
                my $chain = $2;
                my $dnum  = $3;
                $dnum++;
                $ecod_domain_id = "e" . $pdb . $chain . $dnum;
            }
            $seen{$ecod_domain_id}++;

            if ( !$yx_domain_node->exists('@ecod_domain_id') ) {
                $yx_domain_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
            }

            if ( $yx_domain_node->exists('ecod_representative_domain/@ecod_domain_id') ) {
                my $hit_ecod_domain_id = $yx_domain_node->findvalue('ecod_representative_domain/@ecod_domain_id');
                my $ecod_hit_domain_domain_nodes =
                  $ecod_xml_doc->findnodes(qq{//domain[\@ecod_domain_id="$hit_ecod_domain_id"]});
                printf "Found %i nodes for hit domain $hit_ecod_domain_id\n", $ecod_hit_domain_domain_nodes->size();
                if ( $ecod_hit_domain_domain_nodes->size() == 1 ) {
                    my $f_node;
                    if ( $ecod_hit_domain_domain_nodes->get_node(1)->parentNode->exists('@f_id') ) {
                        $f_node = $ecod_hit_domain_domain_nodes->get_node(1)->parentNode;
                    }
                    elsif ( $ecod_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode->exists('@f_id') ) {
                        $f_node = $ecod_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode;
                    }
                    else {
                        die "FAIL1 $ecod_domain_id $hit_ecod_domain_id\n";
                    }

                    #my $f_node = $ecod_hit_domain_domain_nodes->get_node(1)->parentNode();
                    printf "F:%s %s\n", $f_node->findvalue('@name'), $f_node->findvalue('@f_id');
                    $f_node->appendChild($yx_domain_node);
                    my $set_uid = sprintf "%09i", $uid++;
                    my $sql_uid = sprintf "%09i", sql_uid_gen();

                    #$yx_domain_node->setAttribute('uid', $set_uid);
                    $yx_domain_node->setAttribute( 'uid', $sql_uid );    #10/9/2015 XML->SQL uid switch
                    if ( $uid > $maxUID ) { $maxUID = $uid }
                }
                else {
                    print "WARNING! Could not find F-group for ECOD hit domain $hit_ecod_domain_id\n";
                }
            }
            if ( $yx_domain_node->exists('scop_representative_domain/@scop_domain_id') ) {
                my $hit_scop_domain_id = $yx_domain_node->findvalue('scop_representative_domain/@scop_domain_id');
                my $scop_hit_domain_domain_nodes =
                  $ecod_xml_doc->findnodes(qq{//scop_domain[\@scop_domain_id="$hit_scop_domain_id"]});
                printf "Found %i nodes for hit domain $hit_scop_domain_id\n", $scop_hit_domain_domain_nodes->size();
                if ( $scop_hit_domain_domain_nodes->size() == 1 ) {
                    my $f_node;
                    if ( $scop_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode->exists('@f_id') ) {
                        $f_node = $scop_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode;
                    }
                    elsif (
                        $scop_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode->parentNode->exists('@f_id')
                      )
                    {
                        $f_node = $scop_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode->parentNode;
                    }
                    else {
                        die "FAIL2\n";
                    }

                    #my $f_node = $scop_hit_domain_domain_nodes->get_node(1)->parentNode->parentNode;
                    printf "F:%s %s\n", $f_node->findvalue('@name'), $f_node->findvalue('@f_id');
                    $f_node->appendChild($yx_domain_node);
                    if ( $uid > $maxUID ) { $maxUID = $uid }
                }
                else {
                    print "WARNING! Could not find F-group for SCOP hit domain $hit_scop_domain_id\n";
                }
            }
        }
    }
    else {
        die "ERROR! $sub: No domains in side load xml doc? ($size)\n";
    }

    $dict_node->setAttribute( 'maxUID', $maxUID );

}

sub process_run_list_summary_to_ecod {
    my $sub = 'process_run_list_summary_to_ecod';
    my ( $domain_summ_doc, $ecod_xml_doc, $force_replace ) = @_;

    my $pdb_chain_node_XPath = '//pdb_chain';

    my $doc_node = $ecod_xml_doc->findnodes('//domain_dictionary')->get_node(1);
    my $UID      = $doc_node->findvalue('@maxUID') + 1;

    #my $UID = $ecod_xml_doc->findvalue('//domain_dictionary/@maxUID') + 1;

    my %dict_ecod_ids;
    my %dict_ecod_uids;
    my %new_ecod_ids;
    my %dict_ranges;
    my %chain_ranges;
    foreach my $dict_d_node ( find_domain_nodes($ecod_xml_doc)) ) {
        my $uid = sprintf "%09i", $dict_d_node->findvalue('@uid');
        if ( !$dict_d_node->exists('@ecod_domain_id') ) {
            print "WARNING! Missing ecod domain id for $uid\n";
            next;
        }
        my $ecod_domain_id = $dict_d_node->findvalue('@ecod_domain_id');
        my $range;
        $dict_ecod_ids{$ecod_domain_id}++;
        $dict_ecod_uids{$ecod_domain_id} = $uid;

        if ( $dict_d_node->findvalue('@derived_range') eq 'true' ) {
            $range = $dict_d_node->findvalue('derived_seqid_range');
        }
        else {
            $range = $dict_d_node->findvalue('seqid_range');
        }
        $dict_ranges{$ecod_domain_id} = $range;
        my $pdb_id   = $dict_d_node->findvalue('structure/@pdb_id');
        my $chain_id = $dict_d_node->findvalue('structure/@chain_id');

        push( @{ $chain_ranges{$pdb_id}{$chain_id} }, $range );
    }

    my $restrict_pdbs = 0;

    #my %pdbs;
    my %domains;

    my ( $total_domains, $total_chains, $added_domain, $good_chains, $found_domains );

    my %domain_nodes;
    foreach my $ds_pc_node ( $domain_summ_doc->findnodes($pdb_chain_node_XPath)->get_nodelist() ) {

        my $query_pdb   = lc( $ds_pc_node->findvalue('@pdb') );
        my $query_chain = $ds_pc_node->findvalue('@chain');

        #print "p: $query_pdb c: $query_chain\n";

        #if ($restrict_pdbs && !$pdbs{$F[0]}{$F[1]}) { next }

        if ( !$ds_pc_node->exists('domain_parse[@type="merge"]/optimized_chain_domain_coverage') ) {

            #print "opt fail 1 $query_pdb $query_chain\n";
            next;
        }
        if ( $ds_pc_node->findvalue('domain_parse[@type="merge"]/optimized_chain_domain_coverage/@used_res') == 0 ) {

            #print "opt fail 2 $query_pdb $query_chain\n";
            next;
        }
        my $chain_coverage =
          $ds_pc_node->findvalue('domain_parse[@type="merge"]/optimized_chain_domain_coverage/@unused_res');

        if ( $chain_coverage < 25 ) {
            $good_chains++;

            my $domain_XPath = 'domain_parse[@type="merge"]/job/domain_list/domain';

          DOMAIN:
            foreach my $d_node ( $ds_pc_node->findnodes($domain_XPath)->get_nodelist() ) {

                my $ecod_domain_id = $d_node->findvalue('@domain_id');

                #my $derived_seqid_range	= scopify_range($d_node->findvalue('ungapped_seqid_range'), $query_chain);
                my $derived_seqid_range;
                if ( !$d_node->exists('optimized_seqid_range') ) {

                    #print "WARNING! No optimized seqid range for $query_pdb $query_chain $ecod_domain_id\n";
                    if ( $d_node->exists('ungapped_seqid_range') ) {
                        $derived_seqid_range =
                          scopify_range( $d_node->findvalue('ungapped_seqid_range'), $query_chain );
                    }
                    else {
                        die "ERROR! No range found for $ecod_domain_id\n";
                    }
                }
                else {
                    $derived_seqid_range = scopify_range( $d_node->findvalue('optimized_seqid_range'), $query_chain );
                }

                #my $derived_pdb_range	= scopify_range($d_node->findvalue('ungapped_pdb_range'), $query_chain);
                my $derived_pdb_range;
                if ( !$d_node->exists('optimized_pdb_range') ) {

                    #print "WARNING! No optimized pdb range for $query_pdb $query_chain $ecod_domain_id\n";
                    if ( $d_node->exists('ungapped_pdb_range') ) {
                        $derived_pdb_range = scopify_range( $d_node->findvalue('ungapped_pdb_range'), $query_chain );
                    }
                    else {
                        die "ERROR! No range found for $ecod_domain_id\n";
                    }
                }
                else {
                    $derived_pdb_range = scopify_range( $d_node->findvalue('optimized_pdb_range'), $query_chain );
                }

                my $derived_from;
                if ( ( $d_node->findvalue('derived_range/@derivedFrom') =~ /\d+/ ) ) {
                    $derived_from = sprintf "%09i", $d_node->findvalue('derived_range/@derivedFrom');

                    #print "T1: $ecod_domain_id $derived_from\n";
                }
                my $derived_from_domain_id = $d_node->findvalue('hit_domain/@ecod_domain_id');
                if ( !$derived_from ) {
                    if ( $dict_ecod_uids{$derived_from_domain_id} ) {
                        $derived_from = $dict_ecod_uids{$derived_from_domain_id};

                        #print "T2: $ecod_domain_id $derived_from\n";
                    }
                    else {
                        #print "T3: $ecod_domain_id NONE\n";
                        print "Can't find hit $ecod_domain_id\n";
                        next;
                    }
                }
                my $derived_from_reference = $d_node->findvalue('hit_domain/@reference');

                my $f_id = $d_node->findvalue('hit_domain/classification/f_group/@f_id');

                #print "f1 $f_id\n";
                #my $uid = sprintf "%09i", $UID++;
                my $sql_uid = sprintf "%09i", sql_uid_gen();

                my $ecod_domain_node = $ecod_xml_doc->createElement('domain');
                $ecod_domain_node->setAttribute( 'manual_rep', 'false' );
                if ( $chain_coverage >= 25 ) {
                    $ecod_domain_node->setAttribute( 'cluster_derivation_warning', 'true' );
                    $ecod_domain_node->setAttribute( 'unused_res',                 $chain_coverage );
                }

                if ( $restrict_pdbs && !$domains{$ecod_domain_id} ) { next }

                if ( $dict_ecod_ids{$ecod_domain_id} ) {

                    #Renumber, fix overlap manually, impossible to rule which domain is correct automatically.

                    if ( $dict_ranges{$ecod_domain_id} eq $derived_seqid_range ) {
                        print "WARNING! range equiv on $ecod_domain_id, skipping\n";
                        next;
                    }

                    my $new_ecod_domain_id = $ecod_domain_id;
                    $new_ecod_domain_id =~ /e\w{4}\w{1,4}(\d+)/ or die;
                    my $domain_number = $1;
                    while ( $dict_ecod_ids{$new_ecod_domain_id} || $new_ecod_ids{$new_ecod_domain_id} ) {
                        $domain_number++;
                        $new_ecod_domain_id = "e$query_pdb$query_chain$domain_number";
                    }
                    print "WARNING! ID collision, new domain $ecod_domain_id renumbered to $new_ecod_domain_id\n";
                    $ecod_domain_id = $new_ecod_domain_id;

                    #die "ERROR! $ecod_domain_id already in domain dictionary? ID collision?\n";
                }
                my $total_res_conflict = 0;
                foreach my $range ( @{ $chain_ranges{$query_pdb}{$query_chain} } ) {
                    my $c1 = region_coverage( range_expand($range), range_expand($derived_seqid_range) );
                    my $c2 = region_coverage( range_expand($derived_seqid_range), range_expand($range) );
                    my $r1 = residue_coverage( range_expand($derived_seqid_range), range_expand($range) );
                    $total_res_conflict += $r1;

                    #print "$range $derived_seqid_range $c1 $c2\n";
                    if ( $range eq $derived_seqid_range ) {
                        print "WARNING! identical pc range for $ecod_domain_id,skipping...\n";
                        next DOMAIN;
                    }

                    if ( $c1 > 0.8 && $c2 > 0.8 && !$force_replace ) {
                        print "WARNING: Loose correspondence for $ecod_domain_id, skipping...\n";
                        next DOMAIN;
                    }
                }
                if ( $total_res_conflict > 10 && !$force_replace ) {
                    print "WARNING: Residue conflict > 10 for $ecod_domain_id, skipping...\n";
                    next DOMAIN;
                }
                $dict_ecod_ids{$ecod_domain_id}++;
                $dict_ranges{$ecod_domain_id} = $derived_seqid_range;
                $ecod_domain_node->setAttribute( 'derived_range',  'true' );
                $ecod_domain_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );

                #$ecod_domain_node->setAttribute('uid', $uid);
                $ecod_domain_node->setAttribute( 'uid', $sql_uid );    #10/9/2015 XML->SQL uid switch
                $new_ecod_ids{$ecod_domain_id}++;

                my $derived_seqid_range_node = $ecod_xml_doc->createElement('derived_seqid_range');
                $derived_seqid_range_node->setAttribute( 'derivedFromManual',  'false' );
                $derived_seqid_range_node->setAttribute( 'derivedByPartition', 'true' );
                $derived_seqid_range_node->setAttribute( 'derivedFrom',        $derived_from );
                $derived_seqid_range_node->appendTextNode($derived_seqid_range);
                $ecod_domain_node->appendChild($derived_seqid_range_node);

                my $derived_range_node = $ecod_xml_doc->createElement('derived_range');
                $derived_range_node->setAttribute( 'derivedFromManual',  'false' );
                $derived_range_node->setAttribute( 'derivedByPartition', 'true' );
                $derived_range_node->setAttribute( 'derivedFrom',        $derived_from );
                $derived_range_node->appendTextNode($derived_pdb_range);
                $ecod_domain_node->appendChild($derived_range_node);

                my $structure_node = $ecod_xml_doc->createElement('structure');
                $structure_node->setAttribute( 'chain_case_ambiguous', 'false' );
                $structure_node->setAttribute( 'chain_id',             $query_chain );
                $structure_node->setAttribute( 'pdb_id',               $query_pdb );
                $structure_node->setAttribute( 'db',                   'PDB' );
                $ecod_domain_node->appendChild($structure_node);

                my $ecod_domain_ecod_domain_node = $ecod_xml_doc->createElement('ecod_domain');  #This is a little ridic
                $ecod_domain_ecod_domain_node->setAttribute( 'chain_case_ambiguous', 'false' );
                $ecod_domain_ecod_domain_node->setAttribute( 'ecod_domain_id',       $ecod_domain_id );
                $ecod_domain_node->appendChild($ecod_domain_ecod_domain_node);

                my $ecod_representative_domain_node = $ecod_xml_doc->createElement('ecod_representative_domain');
                $ecod_representative_domain_node->setAttribute( 'ecod_domain_id', $derived_from_domain_id );
                $ecod_representative_domain_node->setAttribute( 'uid',            $derived_from );
                $ecod_representative_domain_node->setAttribute( 'reference',      $derived_from_reference );
                $ecod_domain_node->appendChild($ecod_representative_domain_node);

                if ( $d_node->exists('hit_domain/domain/scop_domain') ) {
                    my $scop_representative_domain_node = $ecod_xml_doc->createElement('scop_representative_domain');

                    my $derived_from_scop_domain_id =
                      $d_node->findvalue('hit_domain/domain/scop_domain/@scop_domain_id');
                    $scop_representative_domain_node->setAttribute( 'scop_domain_id', $derived_from_scop_domain_id );
                    $scop_representative_domain_node->setAttribute( 'uid',            $derived_from );
                    my $derived_from_scop_reference = $d_node->findvalue('hit_domain/domain/scop_domain/@reference');
                    $scop_representative_domain_node->setAttribute( 'reference', $derived_from_scop_reference );
                }

                if ( $d_node->exists('hit_domain/domain/ecod_domain') ) {
                    my $ecod_representative_domain_node = $ecod_xml_doc->createElement('ecod_representative_domain');

                    my $derived_from_ecod_domain_id =
                      $d_node->findvalue('hot_domain/domain/ecod_domain/@ecod_domain_id');
                    $ecod_representative_domain_node->setAttribute( 'scop_domain_id', $derived_from_ecod_domain_id );
                    $ecod_representative_domain_node->setAttribute( 'uid',            $derived_from );
                    if ( $d_node->exists('hit_domain/domain/ecod_domain/@reference') ) {
                        my $derived_from_ecod_reference =
                          $d_node->findvalue('hit_domain/domain/ecod_domain/@reference');
                        $ecod_representative_domain_node->setAttribute( 'reference', $derived_from_ecod_reference );
                    }
                    else {

                        #my $derived_from_ecod_reference		= $DEFAULT_ECOD_REFERENCE;
                        die "ERROR! No reference for $ecod_domain_id\n";
                    }
                }

                #push (@{$domain_nodes{$f_id}}, \$ecod_domain_node);
                push( @{ $domain_nodes{$derived_from} }, \$ecod_domain_node );
                $found_domains++;
                $total_domains++;

            }
        }
        else {

            my $domain_XPath = 'domain_parse[@type="merge"]/job/domain_list/domain';
            $total_domains += $ds_pc_node->findnodes($domain_XPath)->size();
        }
        $total_chains++;
    }
    my $pdb_chains_node_XPath = '//pdb_chains';
    foreach my $ds_pc_node ( $domain_summ_doc->findnodes($pdb_chains_node_XPath)->get_nodelist() ) {

        my $query_pdb = lc( $ds_pc_node->findvalue('@pdb') );

        #my $query_chain	= $ds_pc_node->findvalue('@chain');
        my $query_chains = $ds_pc_node->findvalue('@chains');
        my @query_chains = split( "", $query_chains );

        #print "p: $query_pdb c: $query_chain\n";

        #if ($restrict_pdbs && !$pdbs{$F[0]}{$F[1]}) { next }

        if ( !$ds_pc_node->exists('domain_parse[@type="merge"]/optimized_chain_domain_coverage') ) {

            #print "opt fail asm  1 $query_pdb $query_chains\n";
            next;
        }
        if ( $ds_pc_node->findvalue('domain_parse[@type="merge"]/optimized_chain_domain_coverage/@used_res') == 0 ) {

            #print "opt fail asm 2 $query_pdb $query_chains\n";
            next;
        }
        my $chain_coverage =
          $ds_pc_node->findvalue('domain_parse[@type="merge"]/optimized_chain_domain_coverage/@unused_res');

        if ( $chain_coverage < 250 ) {
            $good_chains++;

            my $domain_XPath = 'domain_parse[@type="merge"]/job/domain_list/domain';

          DOMAIN:
            foreach my $d_node ( $ds_pc_node->findnodes($domain_XPath)->get_nodelist() ) {

                my $ecod_domain_id = $d_node->findvalue('@domain_id');

                my $derived_seqid_range;
                if ( !$d_node->exists('optimized_seqid_range') ) {

                    if ( $d_node->exists('ungapped_seqid_range') ) {

                        $derived_seqid_range = $d_node->findvalue('ungapped_seqid_range');
                    }
                    else {
                        die "ERROR! No range found for $ecod_domain_id\n";
                    }
                }
                else {
                    $derived_seqid_range = $d_node->findvalue('optimized_seqid_range');
                }

                my $derived_pdb_range;
                if ( !$d_node->exists('optimized_pdb_range') ) {

                    if ( $d_node->exists('ungapped_pdb_range') ) {
                        $derived_pdb_range = $d_node->findvalue('ungapped_pdb_range');
                    }
                    else {
                        die "ERROR! No range found for $ecod_domain_id\n";
                    }
                }
                else {
                    $derived_pdb_range = $d_node->findvalue('optimized_pdb_range');
                }

                my $derived_from;
                if ( ( $d_node->findvalue('derived_range/@derivedFrom') =~ /\d+/ ) ) {
                    $derived_from = sprintf "%09i", $d_node->findvalue('derived_range/@derivedFrom');

                    #print "T1: $ecod_domain_id $derived_from\n";
                }
                my $derived_from_domain_id = $d_node->findvalue('hit_domain/@ecod_domain_id');
                if ( !$derived_from ) {
                    if ( $dict_ecod_uids{$derived_from_domain_id} ) {
                        $derived_from = $dict_ecod_uids{$derived_from_domain_id};

                        #print "T2: $ecod_domain_id $derived_from\n";
                    }
                    else {
                        #print "T3: $ecod_domain_id NONE\n";
                        next;
                    }
                }
                my $derived_from_reference = $d_node->findvalue('hit_domain/@reference');

                my $f_id = $d_node->findvalue('hit_domain/classification/f_group/@f_id');

                #print "f1 $f_id\n";
                #my $uid = sprintf "%09i", $UID++;
                my $sql_uid = sprintf "%09i", sql_uid_gen();

                my $ecod_domain_node = $ecod_xml_doc->createElement('domain');
                $ecod_domain_node->setAttribute( 'manual_rep', 'false' );
                if ( $chain_coverage >= 25 ) {
                    $ecod_domain_node->setAttribute( 'cluster_derivation_warning', 'true' );
                    $ecod_domain_node->setAttribute( 'unused_res',                 $chain_coverage );
                }

                if ( $restrict_pdbs && !$domains{$ecod_domain_id} ) { next }

                if ( $dict_ecod_ids{$ecod_domain_id} ) {

                    #Renumber, fix overlap manually, impossible to rule which domain is correct automatically.

                    if ( $dict_ranges{$ecod_domain_id} eq $derived_seqid_range ) {
                        print "WARNING! range equiv on $ecod_domain_id, skipping\n";
                        next;
                    }

                    my $new_ecod_domain_id = $ecod_domain_id;
                    $new_ecod_domain_id =~ /e\w{4}\w{1,4}(\d+)/ or die;
                    my $domain_number = $1;
                    while ( $dict_ecod_ids{$new_ecod_domain_id} || $new_ecod_ids{$new_ecod_domain_id} ) {
                        $domain_number++;
                        $new_ecod_domain_id = "e$query_pdb.$domain_number";
                    }
                    print "WARNING! ID collision, new domain $ecod_domain_id renumbered to $new_ecod_domain_id\n";
                    $ecod_domain_id = $new_ecod_domain_id;

                    #die "ERROR! $ecod_domain_id already in domain dictionary? ID collision?\n";
                }
                my $total_res_conflict = 0;
                foreach my $query_chain (@query_chains) {
                    foreach my $range ( @{ $chain_ranges{$query_pdb}{$query_chain} } ) {
                        my $c1 = region_coverage( range_expand($range), range_expand($derived_seqid_range) );
                        my $c2 = region_coverage( range_expand($derived_seqid_range), range_expand($range) );
                        my $r1 = residue_coverage( range_expand($derived_seqid_range), range_expand($range) );
                        $total_res_conflict += $r1;

                        #print "$range $derived_seqid_range $c1 $c2\n";

                        if ( $c1 > 0.8 && $c2 > 0.8 ) {
                            print "WARNING: Loose correspondence for $ecod_domain_id, skipping...\n";

                            #next DOMAIN;
                            $ecod_domain_node->setAttribute( 'multi_chain_clash', 'true' );
                        }
                    }
                }
                if ( $total_res_conflict > 10 ) {
                    print "WARNING: Residue conflict > 10 for $ecod_domain_id, skipping...\n";
                    $ecod_domain_node->setAttribute( 'multi_chain_clash', 'true' );

                    #next DOMAIN;
                }
                $dict_ecod_ids{$ecod_domain_id}++;
                $dict_ranges{$ecod_domain_id} = $derived_seqid_range;
                $ecod_domain_node->setAttribute( 'derived_range',  'true' );
                $ecod_domain_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );

                #$ecod_domain_node->setAttribute('uid', $uid);
                $ecod_domain_node->setAttribute( 'uid', $sql_uid );    #10/9/2015 XML->SQL uid switch
                $new_ecod_ids{$ecod_domain_id}++;

                my $derived_seqid_range_node = $ecod_xml_doc->createElement('derived_seqid_range');
                $derived_seqid_range_node->setAttribute( 'derivedFromManual',  'false' );
                $derived_seqid_range_node->setAttribute( 'derivedByPartition', 'true' );
                $derived_seqid_range_node->setAttribute( 'derivedFrom',        $derived_from );
                $derived_seqid_range_node->appendTextNode($derived_seqid_range);
                $ecod_domain_node->appendChild($derived_seqid_range_node);

                my $derived_range_node = $ecod_xml_doc->createElement('derived_range');
                $derived_range_node->setAttribute( 'derivedFromManual',  'false' );
                $derived_range_node->setAttribute( 'derivedByPartition', 'true' );
                $derived_range_node->setAttribute( 'derivedFrom',        $derived_from );
                $derived_range_node->appendTextNode($derived_pdb_range);
                $ecod_domain_node->appendChild($derived_range_node);

                my $structure_node = $ecod_xml_doc->createElement('structure');
                $structure_node->setAttribute( 'chain_case_ambiguous', 'false' );
                if ( scalar(@query_chains) == 1 ) {
                    $structure_node->setAttribute( 'chain_id', $query_chains[0] );
                }
                else {
                    $structure_node->setAttribute( 'chain_id', "." );
                }
                $structure_node->setAttribute( 'pdb_id', $query_pdb );
                $structure_node->setAttribute( 'db',     'PDB' );
                $ecod_domain_node->appendChild($structure_node);

                my $ecod_domain_ecod_domain_node = $ecod_xml_doc->createElement('ecod_domain');  #This is a little ridic
                $ecod_domain_ecod_domain_node->setAttribute( 'chain_case_ambiguous', 'false' );
                $ecod_domain_ecod_domain_node->setAttribute( 'ecod_domain_id',       $ecod_domain_id );
                $ecod_domain_node->appendChild($ecod_domain_ecod_domain_node);

                my $ecod_representative_domain_node = $ecod_xml_doc->createElement('ecod_representative_domain');
                $ecod_representative_domain_node->setAttribute( 'ecod_domain_id', $derived_from_domain_id );
                $ecod_representative_domain_node->setAttribute( 'uid',            $derived_from );
                $ecod_representative_domain_node->setAttribute( 'reference',      $derived_from_reference );
                $ecod_domain_node->appendChild($ecod_representative_domain_node);

                if ( $d_node->exists('hit_domain/domain/scop_domain') ) {
                    my $scop_representative_domain_node = $ecod_xml_doc->createElement('scop_representative_domain');

                    my $derived_from_scop_domain_id =
                      $d_node->findvalue('hit_domain/domain/scop_domain/@scop_domain_id');
                    $scop_representative_domain_node->setAttribute( 'scop_domain_id', $derived_from_scop_domain_id );
                    $scop_representative_domain_node->setAttribute( 'uid',            $derived_from );
                    my $derived_from_scop_reference = $d_node->findvalue('hit_domain/domain/scop_domain/@reference');
                    $scop_representative_domain_node->setAttribute( 'reference', $derived_from_scop_reference );
                }

                if ( $d_node->exists('hit_domain/domain/ecod_domain') ) {
                    my $ecod_representative_domain_node = $ecod_xml_doc->createElement('ecod_representative_domain');

                    my $derived_from_ecod_domain_id =
                      $d_node->findvalue('hot_domain/domain/ecod_domain/@ecod_domain_id');
                    $ecod_representative_domain_node->setAttribute( 'scop_domain_id', $derived_from_ecod_domain_id );
                    $ecod_representative_domain_node->setAttribute( 'uid',            $derived_from );
                    if ( $d_node->exists('hit_domain/domain/ecod_domain/@reference') ) {
                        my $derived_from_ecod_reference =
                          $d_node->findvalue('hit_domain/domain/ecod_domain/@reference');
                        $ecod_representative_domain_node->setAttribute( 'reference', $derived_from_ecod_reference );
                    }
                    else {
                        #my $derived_from_ecod_reference		= $DEFAULT_ECOD_REFERENCE;
                        #$ecod_representative_domain_node->setAttribute('reference', $derived_from_ecod_reference);
                        die "ERROR! No reference for $ecod_domain_id\n";
                    }
                }

                #push (@{$domain_nodes{$f_id}}, \$ecod_domain_node);
                push( @{ $domain_nodes{$derived_from} }, \$ecod_domain_node );
                $found_domains++;
                $total_domains++;

            }
        }
        else {

            my $domain_XPath = 'domain_parse[@type="merge"]/job_asm/domain_list/domain';
            $total_domains += $ds_pc_node->findnodes($domain_XPath)->size();
        }
        $total_chains++;
    }

    #$doc_node->setAttribute('maxUID', $UID);

    #my $f_group_XPath = '//f_group';
    #
    #foreach my $f_group_node ($ecod_xml_doc->findnodes($f_group_XPath)->get_nodelist() ) {
    #
    #	my $f_id = $f_group_node->findvalue('@f_id');
    #	if (!$domain_nodes{$f_id}) { next }
    #	#printf "f: $f_id\n", scalar(@{$domain_nodes{$f_id}});
    #	foreach my $domain_node (@{$domain_nodes{$f_id}}) {
    #		$f_group_node->appendChild($$domain_node);
    #	}
    #
    #}

    my $domain_XPath = '//domain';

    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {

        my $domain_uid = get_uid($domain_node);

		my $f_node = get_f_node_from_domain_node($domain_node);
    
        foreach my $derived_domain_node ( @{ $domain_nodes{$domain_uid} } ) {
            $f_node->appendChild($$derived_domain_node);
            $added_domain++;
        }
    }

    printf "#Added %i, Found %i, Total %i, Good Chains %i, Total Chains %i\n", $added_domain, $found_domains,
      $total_domains, $good_chains, $total_chains;

}

sub parse_ecod_txt_to_xml {
    my $sub = 'parse_ecod_txt_to_xml';

    my ( $ecod_manual_txt_fn, $ecod_version ) = @_;

    if ( !-f $ecod_manual_txt_fn ) { die "ERROR! $sub: $ecod_manual_txt_fn file not found\n"; }
    open( IN, "<:encoding(UTF-8)", $ecod_manual_txt_fn )
      or die "ERROR! COuld not open $ecod_manual_txt_fn for reading:$!\n";

    my $arch;
    my $arch_comment;

    my ( $x_group, $x_group_name, $x_group_comment );
    my ( $h_group, $h_group_name, $h_group_comment, $h_group_scop_comment, $h_group_alert_comment );
    my ( $f_group, $f_group_name, $f_group_comment, $f_group_scop_comment, $f_group_alert_comment );

    my ( %x_groups, %h_groups, %f_groups );

    my ( $ecod_xml_doc, $ecod_xml_doc_node ) = xml_create('ecod_document');

    my $date = `date`;
    chomp($date);
    my $doc_creation_node = $ecod_xml_doc->createElement('createdOn');
    $doc_creation_node->appendTextNode($date);
    $doc_creation_node->setAttribute( 'version', $ecod_version );
    $ecod_xml_doc_node->appendChild($doc_creation_node);

    my $dictionary_node = $ecod_xml_doc->createElement('domain_dictionary');
    $ecod_xml_doc_node->appendChild($dictionary_node);

    my ( $arch_node, $x_node, $h_node, $f_node );
    my $arch_id      = 1;
    my $arch_ordinal = 0;
    my $x_ordinal    = 0;
    my $h_ordinal    = 0;
    my $f_ordinal    = 0;

    my $UID = sprintf "%09i", 0;

    while ( my $ln = <IN> ) {
        if ( $ln =~ /^\x{FEFF}/ ) {
            print "WARNING! BOM detected in $ln\n";
        }
        $ln =~ s/^\x{FEFF}//;

        #Arch
        if ( $ln =~ /^#([\w\+\/\-\s]+)/ ) {
            $arch = $1;
            if ( $ln =~ /\{(.*)\}/ ) {
                $arch_comment = $1;
            }
            $arch_comment =~ s/\s+$//;
            $arch =~ s/\s+$//;
            chomp($arch);
            if ($DEBUG) {
                print "DEBUG $sub: arch_ln $arch $arch_comment\n";
            }

            $arch_node = $ecod_xml_doc->createElement('architecture');
            $arch_node->setAttribute( 'arch_name',    $arch );
            $arch_node->setAttribute( 'arch_id',      $arch_id );
            $arch_node->setAttribute( 'arch_ordinal', $arch_ordinal++ );
            $arch_id++;

            #reset ordinals
            $x_ordinal = $h_ordinal = $f_ordinal = 0;

            if ($arch_comment) {
                $arch_node->appendTextChild( 'arch_comment', $arch_comment );
            }
            $dictionary_node->appendChild($arch_node);
            next;
        }

        #X_group
        if (
            $ln =~ /^(\d+)\s+ 			#x_id
			     ([\w\s\-\"\'\/\,\+\.\:]+)  #x_name
			    /x
          )
        {

            $x_group      = $1;
            $x_group_name = $2;
            $x_group_name =~ s/\s+$//;    #Trailing whitespace only,
            if ($DEBUG) {
                print "DEBUG: x_group_ln $x_group $x_group_name $x_group_comment\n";
            }

            my $x_node = $ecod_xml_doc->createElement('x_group');
            $x_node->setAttribute( 'x_id',      $x_group );
            $x_node->setAttribute( 'name',      $x_group_name );
            $x_node->setAttribute( 'x_ordinal', $x_ordinal++ );
            $h_ordinal = $f_ordinal = 0;    #Reset ordinals

            if ( !$arch ) { die "No arch, $x_group\n" }

            if ( $ln =~ /\{(.*)\}/ ) {
                my $x_group_comment = $1;
                chomp $x_group_comment;
                $x_group_comment =~ s/\s+$//;
                $x_node->appendTextChild( 'comment', $x_group_comment );
            }

            $arch_node->appendChild($x_node);
            $x_groups{$x_group} = $x_node;
            next;
        }

        #H_group
        if ( $ln =~ /^((\d+)\.\d+)\s+([\w\s\-\"\'\/\,\+\.\:]+)/ ) {
            $h_group = $1;
            my $h_x_group = $2;
            $h_group_name = $3;
            $h_group_name =~ s/\s+$//;
            chomp($h_group_name);

            my $h_x_node;
            if ( exists $x_groups{$h_x_group} ) {
                $h_x_node = $x_groups{$h_x_group};
            }
            else {
                $h_x_node = $ecod_xml_doc->createElement('x_group');
                $h_x_node->setAttribute( 'x_id',      $h_x_group );
                $h_x_node->setAttribute( 'x_ordinal', $x_ordinal++ );
                $h_ordinal = 0;

                #$h_x_node->setAttribute('architecture', $arch);
            }

            my $h_node = $ecod_xml_doc->createElement('h_group');
            $h_node->setAttribute( 'h_id',      $h_group );
            $h_node->setAttribute( 'h_ordinal', $h_ordinal++ );
            $h_node->setAttribute( 'name',      $h_group_name );

            if ($DEBUG) {
                print "DEBUG: h_group_ln $h_group $h_group_name\n";
            }

            if ( $ln =~ /\{([\w\s\-\:\'\+\,\.]+)\}/ ) {
                $h_group_comment = $1;
                $h_node->appendTextChild( 'comment', $h_group_comment );
            }

            if ( $ln =~ /\@([\w\s\-\:\'\.\+\,\.]+)\@/ ) {
                $h_group_scop_comment = $1;
                $h_node->appendTextChild( 'scop_comment', $h_group_scop_comment );
            }
            if ( $ln =~ /\!([\w\s\-\:\'\+\,\.]+)\!/ ) {
                $h_group_alert_comment = $1;
                $h_node->appendTextChild( 'alert_comment', $h_group_alert_comment );
            }

            $h_groups{$h_group} = $h_node;
            $h_x_node->appendChild($h_node);
            $arch_node->appendChild($h_x_node);
            next;
        }

        #F_group
        if (
            $ln =~ /^(((\d+)\.\d+)\.\d+)\s+ 		#f_id
				([\w\s\-\"\'\(\)\>\/\+\,\.\:]+) #f_group_name
				/x
          )
        {
            $f_group = $1;
            my $f_h_group = $2;
            my $f_x_group = $3;

            my $f_h_node;
            if ( $h_groups{$f_h_group} ) {
                $f_h_node = $h_groups{$f_h_group};
            }
            elsif ( $x_groups{$f_x_group} ) {
                my $f_x_node = $x_groups{$f_x_group};

                $f_h_node = $ecod_xml_doc->createElement('h_group');
                $f_h_node->setAttribute( 'h_id',      $f_h_group );
                $f_h_node->setAttribute( 'h_ordinal', $h_ordinal++ );
                $f_ordinal = 0;

                $f_x_node->appendChild($f_h_node);

                #$dictionary_node->appendChild($f_x_node);
                $arch_node->appendChild($f_x_node);
                $h_groups{$f_h_group} = $f_h_node;
            }
            else {

                my $f_x_node = $ecod_xml_doc->createElement('x_group');
                $f_x_node->setAttribute( 'x_id', $f_x_group );

                #$f_x_node->setAttribute('architecture', $arch);
                $f_x_node->setAttribute( 'x_ordinal', $x_ordinal++ );

                #$dictionary_node->appendChild($f_x_node);
                $arch_node->appendChild($f_x_node);
                $x_groups{$f_x_group} = $f_x_node;
                $h_ordinal = 0;

                $f_h_node = $ecod_xml_doc->createElement('h_group');
                $f_h_node->setAttribute( 'h_id',      $f_h_group );
                $f_h_node->setAttribute( 'h_ordinal', $h_ordinal++ );
                $h_groups{$f_h_group} = $f_h_node;
                $f_ordinal = 0;
                $f_x_node->appendChild($f_h_node);

                #$dictionary_node->appendChild($f_x_node);
                $arch_node->appendChild($f_x_node);

            }

            $f_node = $ecod_xml_doc->createElement('f_group');
            $f_node->setAttribute( 'f_id',      $f_group );
            $f_node->setAttribute( 'f_ordinal', $f_ordinal++ );
            $f_groups{$f_group} = $f_node;

            $f_group_name = $4;
            chomp($f_group_name);
            $f_group_name =~ s/\s+$//;
            $f_node->setAttribute( 'name', $f_group_name );

            if ( $ln =~ /\{([\w\s\-\:\'\+\,\.]+)\}/ ) {
                $f_group_comment = $1;
                $f_node->appendTextChild( 'comment', $f_group_comment );
            }
            if ( $ln =~ /\@([\w\s\-\:\'\.\+\,\.]+)\@/ ) {
                $f_group_scop_comment = $1;
                $f_node->appendTextChild( 'scop_comment', $f_group_scop_comment );
            }
            if ($DEBUG) {
                print "DEBUG: f_group_ln $f_group $f_group_name\n";
            }
            if ( $ln =~ /\!([\w\s\-\:\'\+\,\.]+)\!/ ) {
                $f_group_alert_comment = $1;
                $f_node->appendTextChild( 'alert_comment', $f_group_alert_comment );
            }
            $f_h_node->appendChild($f_node);
            next;
        }

        #Domain
        #my $dom_regexp = qr/(d|e)[0-9]\w{3}.(?:\d+|\_)(\[.*?\])?/; #The good old day of one character chain ids
        my $dom_regexp = qr/(d|e)[0-9]\w{3}\w+(?:\d+|\_)(\[.*?\])?/;
        my $rep_regexp = qr/%$dom_regexp%/;
        if ( $ln =~ /^(($dom_regexp\s?($rep_regexp)?\s?[\,\&]\s?)+\s?$dom_regexp\s?($rep_regexp)?\s?)/ ) {

            #print "$sub: WARNING! swap on $ln\n";
            my $domain_assembly_node = $ecod_xml_doc->createElement('domain_assembly');

            #my $sql_uid = sprintf "%09i", sql_uid_gen();
            $domain_assembly_node->setAttribute( 'uid', $UID++ );

            #This needs to be improved.
            my $domain_chunk = $1;
            my @domains;
            my $domain_assembly_type;
            if ( $domain_chunk =~ /&/ ) {
                @domains = split( /&/, $domain_chunk );
                $domain_assembly_type = 'obligate';
                $domain_assembly_node->setAttribute( 'assembly_type', $domain_assembly_type );
            }
            else {
                @domains = split( /\, /, $domain_chunk );
                $domain_assembly_type = 'display';
                $domain_assembly_node->setAttribute( 'assembly_type', $domain_assembly_type );
            }

            my $prime = 0;
            foreach my $domain_bit (@domains) {
                chomp $domain_bit;
                my $domain_node = $ecod_xml_doc->createElement('domain');

                if ( $domain_bit =~ /((d|e)\w{4}\w+[\d\_]+)\[([\w\d\-\,\:]*?)\]/ ) {
                    my $domain_id    = $1;
                    my $cbit         = $2;
                    my $domain_range = $3;

                    if ( $cbit eq 'd' ) {

                        my $scop_domain_node = $ecod_xml_doc->createElement('scop_domain');
                        $scop_domain_node->setAttribute( 'scop_domain_id', $domain_id );
                        $domain_node->setAttribute( 'uid', $UID );
                        $UID++;

                        $scop_domain_node->setAttribute( 'reference', $SCOP_REFERENCE );

                        my $range_node = $ecod_xml_doc->createElement('manual_range');
                        $range_node->appendTextNode($domain_range);
                        $domain_node->appendChild($range_node);
                        $domain_node->setAttribute( 'manual_range', 'true' );
                        $domain_node->appendChild($scop_domain_node);
                        if ( !$prime ) {
                            $domain_node->setAttribute( 'primary', 'true' );
                            $prime = 1;
                        }
                    }
                    elsif ( $cbit eq 'e' ) {

                        my $ecod_domain_node = $ecod_xml_doc->createElement('ecod_domain');
                        $domain_id =~ /(e\w{4}\w+\d+)/;
                        $domain_id = $1;
                        $ecod_domain_node->setAttribute( 'ecod_domain_id', $domain_id );
                        $domain_node->setAttribute( 'uid', $UID );
                        $UID++;
                        $ecod_domain_node->setAttribute( 'reference', $ecod_version );

                        my $range_node = $ecod_xml_doc->createElement('manual_range');
                        $range_node->appendTextNode($domain_range);
                        $domain_node->appendChild($range_node);
                        $domain_node->setAttribute( 'manual_range', 'true' );
                        $domain_node->appendChild($ecod_domain_node);
                        if ( !$prime ) {
                            $domain_node->setAttribute( 'primary', 'true' );
                            $prime = 1;
                        }
                    }
                    else {
                        die "ERROR! Something terrible has happened in dbit regexp\n";
                    }

                }
                elsif ( $domain_bit =~ /((d|e)\w{4}.[\d\_]+)/ ) {
                    my $domain_id = $1;
                    my $cbit      = $2;

                    if ( $cbit eq 'd' ) {
                        my $scop_domain_node = $ecod_xml_doc->createElement('scop_domain');
                        $scop_domain_node->setAttribute( 'scop_domain_id', $domain_id );
                        $domain_node->setAttribute( 'uid', $UID );
                        $UID++;
                        $scop_domain_node->setAttribute( 'reference', $SCOP_REFERENCE );
                        $domain_node->setAttribute( 'scop_implicit_range', 'true' );
                        $domain_node->appendChild($scop_domain_node);
                        if ( !$prime ) {
                            $domain_node->setAttribute( 'primary', 'true' );
                            $prime = 1;
                        }
                    }
                    elsif ( $cbit eq 'e' ) {    #this should never happen;
                        die "ERROR! ECOD domains with implicit range are not allowed\n";
                        my $ecod_domain_node = $ecod_xml_doc->createElement('ecod_domain');
                        $ecod_domain_node->setAttribute( 'ecod_domain_id', $domain_id );
                        $domain_node->setAttribute( 'uid', $UID );
                        $UID++;
                        $ecod_domain_node->setAttribute( 'reference', $ecod_version );
                        $domain_node->setAttribute( 'ecod_implicit_range', 'true' );
                        print "$sub: WARNING! $domain_id is ECOD and has no range ...\n";
                        $domain_node->appendChild($ecod_domain_node);

                        if ( !$prime ) {
                            $domain_node->setAttribute( 'primary', 'true' );
                            $prime = 1;
                        }
                    }
                }

                if ( $domain_bit =~ /%((e|d)[0-9]\w{3}\w+(_|\d+))\[([\:\w\d\-\,]*?)\]?%/ ) {
                    my $manual_rep = 'false';
                    my $rep_domain = $1;
                    my $cbit       = $2;
                    my $rep_range  = $4;
                    if ( $cbit eq 'd' ) {
                        my $manual_scop_rep_node = $ecod_xml_doc->createElement('scop_representative_domain');
                        if ($rep_range) {
                            $manual_scop_rep_node->setAttribute( 'scop_domain_id', $rep_domain );
                            $manual_scop_rep_node->appendTextChild( 'rep_manual_range', $rep_range );
                        }
                        else {
                            $manual_scop_rep_node->setAttribute( 'scop_domain_id', $rep_domain );
                        }

                        $manual_scop_rep_node->setAttribute( 'reference', $SCOP_REFERENCE );
                        $domain_node->appendChild($manual_scop_rep_node);
                        $domain_node->setAttribute( 'manual_rep', $manual_rep );
                    }
                    elsif ( $cbit eq 'e' ) {
                        my $manual_ecod_rep_node = $ecod_xml_doc->createElement('ecod_representative_domain');
                        $manual_ecod_rep_node->setAttribute( 'ecod_domain_id', $rep_domain );
                        $manual_ecod_rep_node->setAttribute( 'reference',      $ecod_version );
                        $domain_node->appendChild($manual_ecod_rep_node);
                        $domain_node->setAttribute( 'manual_rep', $manual_rep );

                    }

                }
                else {
                    $domain_node->setAttribute( 'manual_rep', 'true' );
                }

                $domain_assembly_node->appendChild($domain_node);
            }

            if ( $ln =~ /\{([\w\s\-\,]+)\}/ ) {
                my $domain_comment = $1;
                $domain_assembly_node->appendTextChild( 'comment', $domain_comment );
            }
            if ( $ln =~ /\!([\w\s\-\:\,\']+)\!/ ) {
                my $domain_alert_comment = $1;
                $domain_assembly_node->appendTextChild( 'alert_comment', $domain_alert_comment );
            }
            $f_node->appendChild($domain_assembly_node);
            next;
        }
        elsif ( $ln =~ /^(d\w{4}(.).)/ ) {    #Scop domain case, maybe implicit range
            my $domain      = $1;
            my $domain_node = $ecod_xml_doc->createElement('domain');
            $domain_node->setAttribute( 'uid', $UID );
            $UID++;

            my $scop_domain_node = $ecod_xml_doc->createElement('scop_domain');
            $scop_domain_node->setAttribute( 'scop_domain_id', $domain );
            $scop_domain_node->setAttribute( 'reference',      $SCOP_REFERENCE );
            $domain_node->appendChild($scop_domain_node);

            if ( $ln =~ /%(.*)%/ )
            { #Doing some dangerous assumptions here, we need to de-reference this into a uid in the context of the current version futher down the pipeline
                my $manual_rep           = 'false';
                my $rep_domain           = $1;
                my $manual_scop_rep_node = $ecod_xml_doc->createElement('scop_representative_domain');
                $manual_scop_rep_node->setAttribute( 'scop_domain_id', $rep_domain );
                $manual_scop_rep_node->setAttribute( 'reference',      $SCOP_REFERENCE );
                $domain_node->appendChild($manual_scop_rep_node);
                $domain_node->setAttribute( 'manual_rep', $manual_rep );
            }
            else {
                $domain_node->setAttribute( 'manual_rep', 'true' );
            }

            if ( $ln =~ /\[([\w\d\,\-\:]*)\]/ ) {
                my $domain_range = $1;
                $domain_node->setAttribute( 'manual_range', 'true' );
                $domain_node->appendTextChild( 'manual_range', $domain_range );
            }
            else {
                $domain_node->setAttribute( 'scop_implicit_range', 'true' );
            }

            if ( $ln =~ /\{([\w\s\-\:]+)\}/ ) {
                my $domain_comment = $1;
                $domain_node->appendTextChild( 'comment', $domain_comment );
            }

            if ( $ln =~ /\!([\w\s\-\:\']+)\!/ ) {
                my $domain_alert_comment = $1;
                $domain_node->appendTextChild( 'alert_comment', $domain_alert_comment );
            }

            $f_node->appendChild($domain_node);
            next;
        }
        elsif ( $ln =~ /^(e\w{4}(\.|\w+)\d+)/ ) {
            my $domain      = $1;
            my $domain_node = $ecod_xml_doc->createElement('domain');
            $domain_node->setAttribute( 'uid', $UID );
            $UID++;

            my $ecod_domain_node = $ecod_xml_doc->createElement('ecod_domain');
            $ecod_domain_node->setAttribute( 'ecod_domain_id', $domain );
            $ecod_domain_node->setAttribute( 'reference',      $ecod_version );

            $domain_node->appendChild($ecod_domain_node);

            if ( $ln =~ /%((e|d)[0-9]\w{3}.(_|\d+))\[([\w\:\d\-\,]*?)\]?%/ ) {
                my $manual_rep = 'false';
                my $rep_domain = $1;
                my $cbit       = $2;
                my $dnum       = $3;
                my $rep_range  = $4;
                if ( $rep_domain =~ /^e/ ) {
                    my $manual_ecod_rep_node = $ecod_xml_doc->createElement('ecod_representative_domain');
                    $manual_ecod_rep_node->setAttribute( 'ecod_domain_id', $rep_domain );
                    $manual_ecod_rep_node->setAttribute( 'reference',      $ecod_version );
                    $domain_node->appendChild($manual_ecod_rep_node);
                    $domain_node->setAttribute( 'manual_rep', $manual_rep );
                }
                elsif ( $rep_domain =~ /^d/ ) {
                    my $manual_scop_rep_node = $ecod_xml_doc->createElement('scop_representative_domain');
                    $manual_scop_rep_node->setAttribute( 'scop_domain_id', $rep_domain );
                    $manual_scop_rep_node->setAttribute( 'reference',      $SCOP_REFERENCE );
                    if ($rep_range) {
                        $manual_scop_rep_node->appendTextChild( 'rep_manual_range', $rep_range );
                    }
                    $domain_node->appendChild($manual_scop_rep_node);
                    $domain_node->setAttribute( 'manual_rep', $manual_rep );
                }
            }
            else {
                $domain_node->setAttribute( 'manual_rep', 'true' );
            }

            if ( $ln =~ /\[([\w\d\-\,\:]+)\]/ ) {
                my $domain_range = $1;
                $domain_node->setAttribute( 'manual_range', 'true' );
                my $domain_range_node = $ecod_xml_doc->createElement('manual_range');
                $domain_range_node->appendTextNode($domain_range);
                $domain_node->appendChild($domain_range_node);

            }
            else {
                $domain_node->setAttribute( 'scop_implicit_range', 'true' );
            }

            if ( $ln =~ /\{([\w\s\-\:]+)\}/ ) {
                my $domain_comment      = $1;
                my $domain_comment_node = $ecod_xml_doc->createElement('comment');
                $domain_comment_node->appendTextNode($domain_comment);
                $domain_node->appendChild($domain_comment_node);
            }
            if ( $ln =~ /\!([\w\s\-\:\']+)\!/ ) {
                my $domain_alert_comment      = $1;
                my $domain_alert_comment_node = $ecod_xml_doc->createElement('alert_comment');
                $domain_alert_comment_node->appendTextNode($domain_alert_comment);
                $domain_node->appendChild($domain_alert_comment_node);
            }
            $f_node->appendChild($domain_node);
            next;
        }
        if ( $ln =~ /^\s+$/ ) {
            next;
        }
        die "ERROR! $sub: UNKNOWN LINE: $ln\n";
    }
    return ($ecod_xml_doc);
}

sub merge_new_ecod_pre_xml_to_old_ecod_xml {
    my $sub = 'merge_new_ecod_pre_xml_to_old_ecod_xml';
    my ( $ecod_pre_xml_fn, $old_ecod_xml_fn, $test ) = @_;

	my $old_ecod_xml_doc = xml_open($old_ecod_xml_fn);
	my $new_ecod_xml_doc = xml_open($ecod_pre_xml_fn);

	my ($merge_xml_doc, $merge_root_node) = xml_create('merge_document');

    print "$sub: OLD summary\n";
    ecod_group_summary($old_ecod_xml_doc);
    print "$sub: NEW summary\n";
    ecod_group_summary($new_ecod_xml_doc);

    my $old_version = $old_ecod_xml_doc->findvalue('//createdOn/@version');
    my $new_version = $new_ecod_xml_doc->findvalue('//createdOn/@version');

    my $obsolete_node = $merge_xml_doc->createElement('obsolete_elements');
    $merge_root_node->appendChild($obsolete_node);
    my $modify_node = $merge_xml_doc->createElement('modify_elements');
    $merge_root_node->appendChild($modify_node);
    my $new_node = $merge_xml_doc->createElement('new_elements');
    $merge_root_node->appendChild($new_node);

    #Arch
    my %old_arch_groups;
    my %old_arch_seen;
    my %old_arch_name;
    my %old_arch_name_seen;
    my %old_arch_comment;
    my %new_arch_groups;
    my %new_arch_seen;
    my %new_arch_name;
    my %new_arch_name_seen;
    my %new_arch_comment;

    my %old_arch_ordinal;
    my %new_arch_ordinal;

    print "$sub: read old arch\n";
    foreach my $arch ( find_architecture_nodes($old_ecod_xml_doc) ) {

        my $arch_id   = get_arch_id($arch);
		my $arch_name = has_arch_name($arch) ? get_arch_name($arch) : 'NA';
		my $arch_commment = has_arch_comment($arch) ? get_arch_comment($arch) : 'NA';
		my $arch_ordinal = has_arch_ordinal($arch) ? get_arch_ordinal($arch) : 'NA';

        $old_arch_seen{$arch_id}++;

        #$old_arch_name_seen{$arch_name}++;
        $old_arch_name_seen{$arch_name} = $arch_id;
        $old_arch_name{$arch_id}        = $arch_name;
        $old_arch_comment{$arch_id}     = $arch_comment;
        $old_arch_ordinal{$arch_id}     = $arch_ordinal;
    }

    print "$sub: read new arch\n";
    foreach my $arch ( find_architecture_nodes($new_ecod_xml_doc) ) {

        my $arch_id   = get_arch_id($arch);
		my $arch_name = has_arch_name($arch) ? get_arch_name($arch) : 'NA';
		my $arch_commment = has_arch_comment($arch) ? get_arch_comment($arch) : 'NA';
		my $arch_ordinal = has_arch_ordinal($arch) ? get_arch_ordinal($arch) : 'NA';

		$new_arch_seen{$arch_id}++;
        $new_arch_name{$arch_id}    = $arch_name;
        $new_arch_comment{$arch_id} = $arch_comment;
        $new_arch_name_seen{$arch_name} = $arch_id;
        $new_arch_ordinal{$arch_id}     = $arch_ordinal;
    }

    print "generate old arch rules (obsolete or modify)\n";
    foreach my $old_arch_id ( sort { $a <=> $b } keys %old_arch_seen ) {

        #Architectures do not have a persistent ID
        if ( !$new_arch_name_seen{ $old_arch_name{$old_arch_id} } ) {
            my $obsolete_arch_node = $merge_xml_doc->createElement('obsolete');
            $obsolete_arch_node->setAttribute( 'type',    'arch' );
            $obsolete_arch_node->setAttribute( 'arch_id', $old_arch_id );
            $obsolete_arch_node->setAttribute( 'name',    $old_arch_name{$old_arch_id} );
            $obsolete_node->appendChild($obsolete_arch_node);
        }

        #Architectures cannot be modified, we will assume name change = obsolete + new

        if ( $old_arch_ordinal{$old_arch_id} ne $new_arch_ordinal{ $new_arch_name_seen{ $old_arch_name{$old_arch_id} } }
          )
        {
            my $modify_arch_node = $merge_xml_doc->createElement('modify');
            $modify_arch_node->setAttribute( 'type',        'arch' );
            $modify_arch_node->setAttribute( 'mod_type',    'ordinal_shift' );
            $modify_arch_node->setAttribute( 'arch_id',     $old_arch_id );
            $modify_arch_node->setAttribute( 'old_ordinal', $old_arch_ordinal{$old_arch_id} );
            $modify_arch_node->setAttribute( 'new_ordinal',
                $new_arch_ordinal{ $new_arch_name_seen{ $old_arch_name{$old_arch_id} } } );
            $modify_node->appendChild($modify_arch_node);

        }
    }

    print "generate new arch rules (create)\n";
    foreach my $new_arch_id ( sort { $a <=> $b } keys %new_arch_seen ) {

        if ( !$old_arch_name_seen{ $new_arch_name{$new_arch_id} } ) {

            my $new_arch_node = $merge_xml_doc->createElement('new');
            $new_arch_node->setAttribute( 'type',    'arch' );
            $new_arch_node->setAttribute( 'arch_id', $new_arch_id );
            $new_arch_node->setAttribute( 'name',    $new_arch_name{$new_arch_id} );

            $new_node->appendChild($new_arch_node);
        }
    }

    #X_group;
    my %old_x_groups;
    my %old_x_seen;
    my %old_x_name;
    my %old_x_arch;
    my %old_x_arch_id;
    my %new_x_arch_id;
    my %new_x_ordinal;
    my %old_x_ordinal;
    my %old_x_arch_comment;
    my %new_x_seen;
    my %new_x_name;
    my %new_x_arch;
    my %new_x_arch_comment;
    print "Read old x\n";

    foreach my $x_group ( find_x_group_nodes($old_ecod_xml_doc)) ) {

        my $x_id   = get_x_id($x_group);
		my $x_name = has_name($x_group) ? get_name($x_group) : 'NA';
        #my $x_arch = $x_group->findvalue('@architecture');

        my $x_arch = has_arch_name($x_group->parentNode) ? get_arch_name($x_group)->parentNode : 'NA';
        my $x_arch_comment = has_arch_comment($x_group) ? get_arch_comment($x_group) : 'NA';
		my $x_ordinal = has_x_ordinal($x_group) ? get_x_ordinal($x_group) : 'NaN';

        $old_x_seen{$x_id}++;
        $old_x_name{$x_id}         = $x_name;
        $old_x_arch{$x_id}         = $x_arch;
        $old_x_ordinal{$x_id}      = $x_ordinal;
        $old_x_arch_comment{$x_id} = $x_arch_comment;

        $old_x_arch_id{$x_id} = $x_arch_id;
    }

    print "Read new x\n";
    foreach my $x_group ( $new_ecod_xml_doc->findnodes('//x_group')->get_nodelist() ) {

        my $x_id   = get_x_id($x_group);
        my $x_name = has_name($x_group) ? get_name($x_group) : 'NA';

        #my $x_arch = $x_group->findvalue('@architecture');
        my $x_arch         = get_arch_name($x_group->parentNode);
        my $x_arch_id      = get_arch_id($x_group->parentNode);
        my $x_arch_comment 	= has_arch_comment($x_group) ? get_arch_comment($x_group) : 'NA';
        my $x_ordinal 		= has_x_ordinal($x_group) ? get_x_ordinal($x_group) 'NaN';

        $new_x_seen{$x_id}++;
        $new_x_name{$x_id}         = $x_name;
        $new_x_arch{$x_id}         = $x_arch;
        $new_x_arch_comment{$x_id} = $x_arch_comment;
        $new_x_ordinal{$x_id}      = $x_ordinal;

        $new_x_arch_id{$x_id} = $x_arch_id;
    }

    print "generate old x rules\n";
    foreach my $old_x_id ( sort { $a <=> $b } keys %old_x_seen ) {

        if ( !$new_x_seen{$old_x_id} ) {
            my $obsolete_x_group_node = $merge_xml_doc->createElement('obsolete');
            $obsolete_x_group_node->setAttribute( 'type', 'x_group' );
            $obsolete_x_group_node->setAttribute( 'x_id', $old_x_id );
            if ( $old_x_name{$old_x_id} ) {
                $obsolete_x_group_node->setAttribute( 'name', $old_x_name{$old_x_id} );
            }
            else {
                print "WARNING! No name for $old_x_id\n";
            }
            $obsolete_node->appendChild($obsolete_x_group_node);
            next;
        }

        if ( $old_x_name{$old_x_id} ne $new_x_name{$old_x_id} ) {

            my $modify_x_group_node = $merge_xml_doc->createElement('modify');
            $modify_x_group_node->setAttribute( 'type',     'x_group' );
            $modify_x_group_node->setAttribute( 'x_id',     $old_x_id );
            $modify_x_group_node->setAttribute( 'mod_type', 'name_change' );

            my $old_name_node = $merge_xml_doc->createElement('name');
            $old_name_node->appendTextNode( $old_x_name{$old_x_id} );
            $old_name_node->setAttribute( 'version', $old_version );
            $modify_x_group_node->appendChild($old_name_node);

            my $new_name_node = $merge_xml_doc->createElement('name');
            $new_name_node->appendTextNode( $new_x_name{$old_x_id} );
            $new_name_node->setAttribute( 'version', $new_version );
            $modify_x_group_node->appendChild($new_name_node);

            $modify_node->appendChild($modify_x_group_node);

        }
        if ( $old_x_ordinal{$old_x_id} ne $new_x_ordinal{$old_x_id} ) {

            my $modify_x_group_node = $merge_xml_doc->createElement('modify');
            $modify_x_group_node->setAttribute( 'type',     'x_group' );
            $modify_x_group_node->setAttribute( 'x_id',     $old_x_id );
            $modify_x_group_node->setAttribute( 'mod_type', 'ordinal_shift' );
            if ( $old_x_ordinal{$old_x_id} ) {
                $modify_x_group_node->setAttribute( 'old_ordinal', $old_x_ordinal{$old_x_id} );
            }
            else {
                print "WARNING! No old x_ordinal for $old_x_id\n";
            }
            if ( $new_x_ordinal{$old_x_id} ) {
                $modify_x_group_node->setAttribute( 'new_ordinal', $new_x_ordinal{$old_x_id} );
            }
            else {
                print "WARNING! No new x_ordinal for $old_x_id\n";
            }

            $modify_node->appendChild($modify_x_group_node);
        }

        if ( $old_x_arch{$old_x_id} ne $new_x_arch{$old_x_id} ) {

            if ($DEBUG) {
                print
"DEBUG $sub: arch_shift $old_x_id $old_x_arch{$old_x_id} $new_x_arch{$old_x_id} $old_x_arch_id{$old_x_id} $new_x_arch_id{$old_x_id}\n";
            }

            my $modify_x_group_node = $merge_xml_doc->createElement('modify');
            $modify_x_group_node->setAttribute( 'type',     'x_group' );
            $modify_x_group_node->setAttribute( 'x_id',     $old_x_id );
            $modify_x_group_node->setAttribute( 'mod_type', 'arch_shift' );

            my $old_arch_name_node = $merge_xml_doc->createElement('arch_name');
            $old_arch_name_node->setAttribute( 'version', $old_version );
            if ( $old_x_arch_id{$old_x_id} ) {
                $old_arch_name_node->setAttribute( 'arch_id', $old_x_arch_id{$old_x_id} );
            }
            else {
                print "WARNING! No arch_id for $old_x_id\n";
            }
            $old_arch_name_node->appendTextNode( $old_x_arch{$old_x_id} );
            $modify_x_group_node->appendChild($old_arch_name_node);

            my $new_arch_name_node = $merge_xml_doc->createElement('arch_name');
            $new_arch_name_node->setAttribute( 'version', $new_version );
            if ( $new_x_arch_id{$old_x_id} ) {
                $new_arch_name_node->setAttribute( 'arch_id', $new_x_arch_id{$old_x_id} );
            }
            else {
                print "WARNING! No arch_id for $old_x_id\n";
            }
            if ( $new_x_arch{$old_x_id} ) {
                $new_arch_name_node->appendTextNode( $new_x_arch{$old_x_id} );
            }
            else {
                print "WARNING! No new arch_name for $old_x_id\n";
            }
            $modify_x_group_node->appendChild($new_arch_name_node);

            $modify_node->appendChild($modify_x_group_node);
        }

    }

    print "generate new x rules\n";
    foreach my $new_x_id ( sort { $a <=> $b } keys %new_x_seen ) {
        if ( !$old_x_seen{$new_x_id} ) {

            my $new_x_group_node = $merge_xml_doc->createElement('new');
            $new_x_group_node->setAttribute( 'type',         'x_group' );
            $new_x_group_node->setAttribute( 'x_id',         $new_x_id );
            $new_x_group_node->setAttribute( 'name',         $new_x_name{$new_x_id} );
            $new_x_group_node->setAttribute( 'architecture', $new_x_arch{$new_x_id} );
            if ( $new_x_arch_comment{$new_x_id} ) {
                $new_x_group_node->setAttribute( 'arch_comment', $new_x_arch_comment{$new_x_id} );
            }
            else {
                print "WARNING! No arch_comment for $new_x_id\n";
            }

            $new_node->appendChild($new_x_group_node);

        }
    }

    #H_group
    my %old_h_groups;
    my %old_h_seen;
    my %old_h_name;
    my %new_h_seen;
    my %new_h_name;

    my %old_h_ordinal;
    my %new_h_ordinal;

    print "Read old h\n";
    foreach my $h_group ( find_h_group_nodes($old_ecod_xml_doc) ) {

        my $h_id   = get_h_id($h_group);
        my $h_name = has_name($h_group) ? get_name($h_group) : 'NA';
        my $h_ordinal = has_h_ordinal($h_group) ? get_h_ordinal($h_group) : 'NaN';

        $old_h_seen{$h_id}++;
        $old_h_name{$h_id}    = $h_name;
        $old_h_ordinal{$h_id} = $h_ordinal;
    }

    print "Read new h\n";
    foreach my $h_group ( $new_ecod_xml_doc->findnodes('//h_group') ) {

        my $h_id   = get_h_id($h_group);
        my $h_name = has_name($h_group) ? get_name($h_group) : 'NA';
        if ( $h_group->exists('@name') ) {
            $h_name = $h_group->findvalue('@name');
        }
        my $h_ordinal = 'NaN';
        if ( $h_group->exists('@h_ordinal') ) {
            $h_ordinal = $h_group->findvalue('@h_ordinal');
        }

        $new_h_seen{$h_id}++;
        $new_h_name{$h_id}    = $h_name;
        $new_h_ordinal{$h_id} = $h_ordinal;
    }

    print "generate old h rules\n";
    foreach my $old_h_id ( sort { $a <=> $b } keys %old_h_seen ) {

        if ( !$new_h_seen{$old_h_id} ) {
            if ($DEBUG) {
                print "OBS h_group? $old_h_id $old_h_name{$old_h_id}\n";
            }
            my $obsolete_h_group_node = $merge_xml_doc->createElement('obsolete');
            $obsolete_h_group_node->setAttribute( 'type', 'h_group' );
            $obsolete_h_group_node->setAttribute( 'h_id', $old_h_id );
            if ( $old_h_name{$old_h_id} ) {
                $obsolete_h_group_node->setAttribute( 'name', $old_h_name{$old_h_id} );
            }
            else {
                print "WARNING! No h_name for $old_h_id\n";
            }
            $obsolete_node->appendChild($obsolete_h_group_node);

            next;
        }

        if ( $old_h_name{$old_h_id} ne $new_h_name{$old_h_id} ) {
            if ($DEBUG) {
                print "NAME hold->hnew $old_h_id $old_h_name{$old_h_id} => $new_h_name{$old_h_id}\n";
            }

            my $modify_h_group_node = $merge_xml_doc->createElement('modify');
            $modify_h_group_node->setAttribute( 'type',     'h_group' );
            $modify_h_group_node->setAttribute( 'h_id',     $old_h_id );
            $modify_h_group_node->setAttribute( 'mod_type', 'name_change' );

            my $old_name_node = $merge_xml_doc->createElement('name');
            $old_name_node->setAttribute( 'version', $old_version );
            $old_name_node->appendTextNode( $old_h_name{$old_h_id} );
            $modify_h_group_node->appendChild($old_name_node);

            my $new_name_node = $merge_xml_doc->createElement('name');
            $new_name_node->setAttribute( 'version', $new_version );
            $new_name_node->appendTextNode( $new_h_name{$old_h_id} );
            $modify_h_group_node->appendChild($new_name_node);

            $modify_node->appendChild($modify_h_group_node);

        }

        if ( $old_h_ordinal{$old_h_id} ne $new_h_ordinal{$old_h_id} ) {

            my $modify_h_group_node = $merge_xml_doc->createElement('modify');
            $modify_h_group_node->setAttribute( 'type',     'h_group' );
            $modify_h_group_node->setAttribute( 'h_id',     $old_h_id );
            $modify_h_group_node->setAttribute( 'mod_type', 'ordinal_shift' );
            if ( $old_h_ordinal{$old_h_id} ) {
                $modify_h_group_node->setAttribute( 'old_ordinal', $old_h_ordinal{$old_h_id} );
            }
            else {
                print "WARNING! No old h_ordinal for $old_h_id\n";
            }
            if ( $new_h_ordinal{$old_h_id} ) {
                $modify_h_group_node->setAttribute( 'new_ordinal', $new_h_ordinal{$old_h_id} );
            }
            else {
                print "WARNING! No new h_ordinal for $old_h_id\n";
            }

            $modify_node->appendChild($modify_h_group_node);
        }

    }

    print "Generate new h rules\n";
    foreach my $new_h_id ( sort { $a <=> $b } keys %new_h_seen ) {
        if ( !$old_h_seen{$new_h_id} ) {
            if ($DEBUG) {
                print "NEW h_group $new_h_id $new_h_name{$new_h_id}\n";
            }

            my $new_h_group_node = $merge_xml_doc->createElement('new');
            $new_h_group_node->setAttribute( 'type', 'h_group' );
            $new_h_group_node->setAttribute( 'h_id', $new_h_id );
            if ( $new_h_name{$new_h_id} ) {
                $new_h_group_node->setAttribute( 'name', $new_h_name{$new_h_id} );
            }
            else {
                print "WARNING! No new h_name for $new_h_id\n";
            }

            $new_node->appendChild($new_h_group_node);
        }
    }

    #F_group
    my %old_f_groups;
    my %old_f_seen;
    my %old_f_name;
    my %new_f_seen;
    my %new_f_name;
    my %new_f_ordinal;
    my %old_f_ordinal;

    print "Read old f\n";
    foreach my $f_group ( $old_ecod_xml_doc->findnodes('//f_group') ) {

        my $f_id   = $f_group->findvalue('@f_id');
        my $f_name = 'NA';
        if ( $f_group->exists('@name') ) {
            $f_name = $f_group->findvalue('@name');
        }

        my $f_ordinal = 'NaN';
        if ( $f_group->exists('@f_group') ) {
            $f_ordinal = $f_group->findvalue('@f_ordinal');
        }

        $old_f_seen{$f_id}++;
        $old_f_name{$f_id}    = $f_name;
        $old_f_ordinal{$f_id} = $f_ordinal;
    }

    print "Read new f\n";
    foreach my $f_group ( $new_ecod_xml_doc->findnodes('//f_group') ) {

        my $f_id   = $f_group->findvalue('@f_id');
        my $f_name = 'NA';
        if ( $f_group->exists('@name') ) {
            $f_name = $f_group->findvalue('@name');
        }
        my $f_ordinal = 'NaN';
        if ( $f_group->exists('@f_ordinal') ) {
            $f_ordinal = $f_group->findvalue('@f_ordinal');
        }

        $new_f_seen{$f_id}++;
        $new_f_name{$f_id}    = $f_name;
        $new_f_ordinal{$f_id} = $f_ordinal;
    }

    print "Generate old f rules\n";
    foreach my $old_f_id ( sort { $a cmp $b } keys %old_f_seen ) {

        if ( !$new_f_seen{$old_f_id} ) {
            if ($DEBUG) {
                print "OBS f_group? $old_f_id $old_f_name{$old_f_id}\n";
            }
            my $obsolete_f_group_node = $merge_xml_doc->createElement('obsolete');
            $obsolete_f_group_node->setAttribute( 'type', 'f_group' );
            $obsolete_f_group_node->setAttribute( 'f_id', $old_f_id );
            if ( $old_f_name{$old_f_id} ) {
                $obsolete_f_group_node->setAttribute( 'name', $old_f_name{$old_f_id} );
            }
            else {
                print "WARNING! No old f_name for $old_f_id\n";
            }
            $obsolete_node->appendChild($obsolete_f_group_node);
            next;
        }

        if ( $old_f_name{$old_f_id} ne $new_f_name{$old_f_id} ) {
            if ($DEBUG) {
                print "NAME fold->fnew $old_f_id $old_f_name{$old_f_id} => $new_f_name{$old_f_id}\n";
            }

            my $modify_f_group_node = $merge_xml_doc->createElement('modify');
            $modify_f_group_node->setAttribute( 'type',     'f_group' );
            $modify_f_group_node->setAttribute( 'f_id',     $old_f_id );
            $modify_f_group_node->setAttribute( 'mod_type', 'name_change' );

            my $old_name_node = $merge_xml_doc->createElement('name');
            $old_name_node->setAttribute( 'version', $old_version );
            $old_name_node->appendTextNode( $old_f_name{$old_f_id} );
            $modify_f_group_node->appendChild($old_name_node);

            my $new_name_node = $merge_xml_doc->createElement('name');
            $new_name_node->setAttribute( 'version', $new_version );
            $new_name_node->appendTextNode( $new_f_name{$old_f_id} );
            $modify_f_group_node->appendChild($new_name_node);

            $modify_node->appendChild($modify_f_group_node);
        }

        if ( $old_f_ordinal{$old_f_id} ne $new_f_ordinal{$old_f_id} ) {

            my $modify_f_group_node = $merge_xml_doc->createElement('modify');
            $modify_f_group_node->setAttribute( 'type',     'f_group' );
            $modify_f_group_node->setAttribute( 'f_id',     $old_f_id );
            $modify_f_group_node->setAttribute( 'mod_type', 'ordinal_shift' );
            if ( $old_f_ordinal{$old_f_id} ) {
                $modify_f_group_node->setAttribute( 'old_ordinal', $old_f_ordinal{$old_f_id} );
            }
            else {
                print "WARNING! No old f_ordinal for $old_f_id\n";
            }
            if ( $new_f_ordinal{$old_f_id} ) {
                $modify_f_group_node->setAttribute( 'new_ordinal', $new_f_ordinal{$old_f_id} );
            }
            else {
                print "WARNING! No new f_ordinal for $old_f_id\n";
            }

            $modify_node->appendChild($modify_f_group_node);
        }
        else {
            if ($DEBUG) {
                print "DEBUG: f_ordinal $old_f_id $old_f_ordinal{$old_f_id} $new_f_ordinal{$old_f_id}\n";
            }
        }

    }

    print "Generate new f rules\n";
    foreach my $new_f_id ( sort { $a cmp $b } keys %new_f_seen ) {
        if ( !$old_f_seen{$new_f_id} ) {
            if ($DEBUG) {
                print "NEW f_group $new_f_id $new_f_name{$new_f_id}\n";
            }

            my $new_f_group_node = $merge_xml_doc->createElement('new');
            $new_f_group_node->setAttribute( 'type', 'f_group' );
            $new_f_group_node->setAttribute( 'f_id', $new_f_id );
            $new_f_group_node->setAttribute( 'name', $new_f_name{$new_f_id} );
            if ( $new_f_ordinal{$new_f_id} ) {
                $new_f_group_node->setAttribute( 'f_ordinal', $new_f_ordinal{$new_f_id} );
            }
            else {
                print "WARNING! No new f_ordinal for $new_f_id\n";
            }

            $new_node->appendChild($new_f_group_node);

        }
    }

    #Domain assemblies

    my (
        %new_domain_assembly_seen, %new_domain_assembly_parent, %new_domain_assembly_range,
        %new_domain_assembly_uid,  %new_domain_assembly_unhash, %new_domain_assembly_parent_uid,
        %new_domain_assembly_fid,  %new_domain_assembly_consist
    );

    print "Read new domain assemblies\n";
    foreach my $domain_assembly ( $new_ecod_xml_doc->findnodes('//domain_assembly')->get_nodelist() ) {

        my $uid = $domain_assembly->findvalue('@uid');

        my $f_id;
        if ( $domain_assembly->parentNode->nodeName eq 'f_group' ) {
            $f_id = $domain_assembly->parentNode->findvalue('@f_id');
        }
        elsif ( $domain_assembly->parentNode->parentNode->nodeName eq 'f_group' ) {
            $f_id = $domain_assembly->parentNode->parentNode->findvalue('@f_id');
        }
        else {
            die "Confused domain assembly $uid\n";
        }

        my $assembly_type = $domain_assembly->findvalue('@assembly_type');

        my $primary_domain_id;
        my $primary_domain;
        if ( $domain_assembly->exists('domain[@primary="true"]') ) {
            $primary_domain = $domain_assembly->findnodes('domain[@primary="true"]')->get_node(1);
            if ( $primary_domain->exists('scop_domain') ) {
                $primary_domain_id = $primary_domain->findvalue('scop_domain/@scop_domain_id');
            }
            elsif ( $primary_domain->exists('ecod_domain') ) {
                $primary_domain_id = $primary_domain->findvalue('ecod_domain/@ecod_domain_id');
            }
            else {
                die "?da $uid\n";
            }
        }
        else {
            print "WARNING! Domain assembly $uid has no primary domainid\n";
        }

        my %domain_consist;
        foreach my $domain ( $domain_assembly->findnodes('domain')->get_nodelist() ) {
            my $domain_id;
            if ( $domain->exists('scop_domain') ) {
                $domain_id = $domain->findvalue('scop_domain/@scop_domain_id');
            }
            elsif ( $domain->exists('ecod_domain') ) {
                $domain_id = $domain->findvalue('ecod_domain/@ecod_domain_id');
            }
            else {
                die "?da $uid\n";
            }
            $domain_consist{$domain_id}++;
        }
        my $range;
        my $type;
        if ( $primary_domain->findvalue('@scop_implicit_range') eq 'true' ) {
            $range = 'imp';
            $type  = 'implicit';
        }
        elsif ( $primary_domain->findvalue('@manual_range') eq 'true' ) {
            $range = $primary_domain->findvalue('manual_range');
            $type  = 'explicit';
        }
        else {
            $range = $primary_domain->findvalue('derived_range');
            $type  = 'derived';
        }

        my $hash = $primary_domain_id . $range;
        $new_domain_assembly_seen{$hash}++;
        $new_domain_assembly_fid{$hash}     = $f_id;
        $new_domain_assembly_range{$hash}   = $range;
        $new_domain_assembly_uid{$hash}     = $uid;
        $new_domain_assembly_unhash{$hash}  = $primary_domain_id;
        $new_domain_assembly_consist{$hash} = \%domain_consist;

        #Number of domains?
    }

    my (
        %old_domain_assembly_seen, %old_domain_assembly_parent,  %old_domain_assembly_range,
        %old_domain_assembly_uid,  %old_domain_assembly_unhash,  %old_domain_assembly_parent_uid,
        %old_domain_assembly_fid,  %old_domain_assembly_uid_fid, %old_domain_assembly_consist
    );
    my ( %old_domain_assembly_type, %old_domain_assembly_derivedFrom, %old_domain_assembly_derivedFrom_uid );
    print "Read old domain assemblies\n";
    foreach my $domain_assembly ( $old_ecod_xml_doc->findnodes('//domain_assembly')->get_nodelist() ) {

        my $uid = $domain_assembly->findvalue('@uid');

        my $f_id;
        if ( $domain_assembly->parentNode->nodeName eq 'f_group' ) {
            $f_id = $domain_assembly->parentNode->findvalue('@f_id');
        }
        elsif ( $domain_assembly->parentNode->parentNode->nodeName eq 'f_group' ) {
            $f_id = $domain_assembly->parentNode->parentNode->findvalue('@f_id');
        }
        else {
            die "Confused domain assembly $uid\n";
        }
        $old_domain_assembly_uid_fid{$uid} = $f_id;
        my $assembly_type = $domain_assembly->findvalue('@assembly_type');

        my $primary_domain_id;
        my $primary_domain;
        if ( $domain_assembly->exists('domain[@primary="true"]') ) {
            $primary_domain = $domain_assembly->findnodes('domain[@primary="true"]')->get_node(1);
            if ( $primary_domain->exists('scop_domain') ) {
                $primary_domain_id = $primary_domain->findvalue('scop_domain/@scop_domain_id');
            }
            elsif ( $primary_domain->exists('ecod_domain') ) {
                $primary_domain_id = $primary_domain->findvalue('ecod_domain/@ecod_domain_id');
            }
            else {
                die "?da $uid\n";
            }
        }
        else {
            print "WARNING! Domain assembly $uid has no primary domainid\n";
        }
        my %domain_consist;
        foreach my $domain ( $domain_assembly->findnodes('domain')->get_nodelist() ) {
            my $domain_id;
            if ( $domain->exists('scop_domain') ) {
                $domain_id = $domain->findvalue('scop_domain/@scop_domain_id');
            }
            elsif ( $domain->exists('ecod_domain') ) {
                $domain_id = $domain->findvalue('ecod_domain/@ecod_domain_id');
            }
            else {
                die "?da $uid\n";
            }
            $domain_consist{$domain_id}++;
        }
        my $range;
        my $type;
        if ( $primary_domain->findvalue('@scop_implicit_range') eq 'true' ) {
            $range = 'imp';
            $type  = 'implicit';
        }
        elsif ( $primary_domain->findvalue('@manual_range') eq 'true' ) {
            $range = $primary_domain->findvalue('manual_range');
            $type  = 'explicit';
        }
        else {
            $range = $primary_domain->findvalue('derived_range');
            $type  = 'derived';
        }

        my $derivedFrom;
        my $assemblyDerived;
        if ( $domain_assembly->findvalue('@domain_assembly_representative') ne 'true' ) {
            $assemblyDerived = 1;
            $derivedFrom     = $domain_assembly->findvalue('@derivedFrom');
        }

        my $hash = $primary_domain_id . $range;
        $old_domain_assembly_seen{$hash}++;
        $old_domain_assembly_fid{$hash}     = $f_id;
        $old_domain_assembly_range{$hash}   = $range;
        $old_domain_assembly_uid{$hash}     = $uid;
        $old_domain_assembly_unhash{$hash}  = $primary_domain_id;
        $old_domain_assembly_consist{$hash} = \%domain_consist;
        if ($assemblyDerived) {
            $old_domain_assembly_type{$hash}           = 'derived';
            $old_domain_assembly_derivedFrom{$hash}    = $derivedFrom;
            $old_domain_assembly_derivedFrom_uid{$uid} = $derivedFrom;
        }
        else {
            $old_domain_assembly_type{$hash} = 'rep';
        }

    }

    print "Generate old domain assembly rules\n";
    foreach my $old_domain_assembly ( sort { $a cmp $b } keys %old_domain_assembly_seen ) {

        my $primary_domain_id = $old_domain_assembly_unhash{$old_domain_assembly};

        #if (!$old_domain_da{$primary_domain_id}) {
        #	die "ERROR! $primary_domain_id?\n";
        #}

        if (   $new_domain_assembly_seen{$old_domain_assembly}
            && $old_domain_assembly_range{$old_domain_assembly} eq $new_domain_assembly_range{$old_domain_assembly} )
        {
#print "ndas: $old_domain_assembly $new_domain_assembly_seen{$old_domain_assembly} $old_domain_assembly_fid{$old_domain_assembly} $new_domain_assembly_fid{$old_domain_assembly}\n";
            if ( $old_domain_assembly_fid{$old_domain_assembly} eq $new_domain_assembly_fid{$old_domain_assembly} ) {
                if ($DEBUG) {
                    print
"DOMAIN_ASSEMBLY_CORRESPONDENCE $old_domain_assembly $new_domain_assembly_uid{$old_domain_assembly} $new_domain_assembly_range{$old_domain_assembly} $old_domain_assembly_range{$old_domain_assembly}\n";
                }

                my $modify_domain_assembly_node = $merge_xml_doc->createElement('modify');
                $modify_domain_assembly_node->setAttribute( 'type',      'domain_assembly' );
                $modify_domain_assembly_node->setAttribute( 'mod_type',  'inherit_annotation' );
                $modify_domain_assembly_node->setAttribute( 'domain_id', $primary_domain_id );
                $modify_domain_assembly_node->setAttribute( 'new_uid', $new_domain_assembly_uid{$old_domain_assembly} );
                $modify_domain_assembly_node->setAttribute( 'old_uid', $old_domain_assembly_uid{$old_domain_assembly} );

                $modify_node->appendChild($modify_domain_assembly_node);
            }
            else {
                if ($DEBUG) {
                    print
"DOMAIN_ASSEMBLY_FSHIFT $old_domain_assembly $new_domain_assembly_uid{$old_domain_assembly} $old_domain_assembly_uid{$old_domain_assembly} $new_domain_assembly_range{$old_domain_assembly} $old_domain_assembly_range{$old_domain_assembly}\n";
                }

                my $modify_domain_assembly_node = $merge_xml_doc->createElement('modify');
                $modify_domain_assembly_node->setAttribute( 'type',      'domain_assembly' );
                $modify_domain_assembly_node->setAttribute( 'mod_type',  'f_shift' );
                $modify_domain_assembly_node->setAttribute( 'domain_id', $primary_domain_id );
                $modify_domain_assembly_node->setAttribute( 'new_uid', $new_domain_assembly_uid{$old_domain_assembly} );
                $modify_domain_assembly_node->setAttribute( 'old_uid', $old_domain_assembly_uid{$old_domain_assembly} );
                $modify_domain_assembly_node->setAttribute( 'new_fid', $new_domain_assembly_fid{$old_domain_assembly} );
                $modify_domain_assembly_node->setAttribute( 'old_fid', $old_domain_assembly_fid{$old_domain_assembly} );

                $modify_node->appendChild($modify_domain_assembly_node);
            }

        }
        elsif ( $old_domain_assembly_type{$old_domain_assembly} eq 'derived' ) {
            if ($DEBUG) {
                print "DERIVED_DOMAIN_ASSEMBLY $old_domain_assembly $old_domain_assembly_uid{$old_domain_assembly}\n";
            }

            my $modify_domain_assembly_node = $merge_xml_doc->createElement('modify');
            $modify_domain_assembly_node->setAttribute( 'type',      'domain_assembly' );
            $modify_domain_assembly_node->setAttribute( 'mod_type',  'derived_shift' );
            $modify_domain_assembly_node->setAttribute( 'domain_id', $primary_domain_id );
            print "primary domain id: $primary_domain_id\n";
            $modify_domain_assembly_node->setAttribute( 'old_uid', $old_domain_assembly_uid{$old_domain_assembly} );
            print "old_uid $old_domain_assembly_uid{$old_domain_assembly}\n";
            $modify_domain_assembly_node->setAttribute( 'old_parent_uid',
                $old_domain_assembly_derivedFrom{$old_domain_assembly} );
            print "old_parent_uid $old_domain_assembly_derivedFrom{$old_domain_assembly}\n";
            my $rep_domain_assembly_uid = $old_domain_assembly_derivedFrom{$old_domain_assembly};
            my $parent_f_id             = $old_domain_assembly_uid_fid{$rep_domain_assembly_uid};
            my $MAX                     = 10;
            my $i                       = 0;

            while ( $old_domain_assembly_derivedFrom_uid{$rep_domain_assembly_uid} && $i < $MAX ) {
                $rep_domain_assembly_uid = $old_domain_assembly_derivedFrom_uid{$rep_domain_assembly_uid};
                $parent_f_id             = $old_domain_assembly_uid_fid{$rep_domain_assembly_uid};
                $i++;
            }
            $modify_domain_assembly_node->setAttribute( 'old_parent_f_id', $parent_f_id );
            print "old_parent_f_id $parent_f_id\n";
            $modify_node->appendChild($modify_domain_assembly_node);

        }
        else {
            if ($DEBUG) {
                print "OBSOLETE_DOMAIN_ASSEMBLY $old_domain_assembly $old_domain_assembly_uid{$old_domain_assembly}\n";
            }

            my $obsolete_domain_node = $merge_xml_doc->createElement('obsolete');
            $obsolete_domain_node->setAttribute( 'type', 'domain_assembly' );
            $obsolete_domain_node->setAttribute( 'uid',  $old_domain_assembly_uid{$old_domain_assembly} );
            my $real_domain_id = $old_domain_assembly_unhash{$old_domain_assembly};
            $obsolete_domain_node->setAttribute( 'domain_id', $real_domain_id );
            $obsolete_domain_node->setAttribute( 'parent_f_id',
                $old_domain_assembly_uid_fid{ $old_domain_assembly_uid{$old_domain_assembly} } );

            $obsolete_node->appendChild($obsolete_domain_node);
        }
    }

    #Domains
    my (
        %new_domain_seen, %new_domain_parent, %new_domain_range,
        %new_domain_uid,  %new_domain_unhash, %new_domain_parent_uid,
        %new_domain_fid,  %new_domain_type,   %new_domain_assembly
    );
    my %new_domain_da;
    print "Read new domains\n";
    foreach my $domain ( $new_ecod_xml_doc->findnodes('//domain')->get_nodelist() ) {

        my $domain_id;
        my $uid = $domain->findvalue('@uid');
        if ( $domain->exists('scop_domain') ) {
            $domain_id = $domain->findvalue('scop_domain/@scop_domain_id');
        }
        elsif ( $domain->exists('ecod_domain') ) {
            $domain_id = $domain->findvalue('ecod_domain/@ecod_domain_id');
        }
        else {
            die "new: $uid?\n";
        }

        my $f_id;
        if ( $domain->parentNode->exists('@f_id') ) {
            $f_id = $domain->parentNode->findvalue('@f_id');
        }
        elsif ( $domain->parentNode->parentNode->exists('@f_id') ) {
            $f_id = $domain->parentNode->parentNode->findvalue('@f_id');
        }
        my $isAssemblyMember = 0;
        my $da_uid;
        my $da_type;
        if ( $domain->parentNode->nodeName eq 'domain_assembly' ) {
            $isAssemblyMember = 1;
            $da_uid           = $domain->parentNode->findvalue('@uid');
            $da_type          = $domain->parentNode->findvalue('@assembly_type');
        }
        my $range;
        my $type;
        if ( $domain->findvalue('@scop_implicit_range') eq 'true' ) {
            $range = 'imp';
            $type  = 'implicit';
        }
        else {
            $range = $domain->findvalue('manual_range');
            $type  = 'explicit';
        }

        my $hash = $domain_id . $range;
        $new_domain_seen{$hash}++;
        $new_domain_fid{$hash}      = $f_id;
        $new_domain_range{$hash}    = $range;
        $new_domain_uid{$hash}      = $uid;
        $new_domain_unhash{$hash}   = $domain_id;
        $new_domain_da{$hash}       = $da_uid;
        $new_domain_assembly{$hash} = $da_type;
        $new_domain_type{$hash}     = $type;

    }

    my (
        %old_domain_seen,    %old_domain_parent, %old_domain_range,      %old_domain_uid,
        %old_domain_type,    %old_domain_unhash, %old_domain_parent_uid, %old_domain_fid,
        %old_domain_uid_fid, %old_domain_domainid2uid
    );
    my %old_domain_derivedFrom;
    my %old_domain_derivedFrom_uid;
    my %old_domain_da;
    print "Read new domains\n";
    foreach my $domain ( $old_ecod_xml_doc->findnodes('//domain') ) {

        my $domain_id;
        my $uid = $domain->findvalue('@uid');
        if ( $domain->exists('scop_domain') ) {
            $domain_id = $domain->findvalue('scop_domain/@scop_domain_id');
        }
        elsif ( $domain->exists('ecod_domain') ) {
            $domain_id = $domain->findvalue('ecod_domain/@ecod_domain_id');
            $old_domain_domainid2uid{$domain_id} = $uid;
        }
        else {
            die "old: $uid?\n";
        }

        my $isAssemblyMember = 0;
        my $da_uid;
        if ( $domain->parentNode->nodeName eq 'domain_assembly' ) {
            $isAssemblyMember = 1;
            $da_uid           = $domain->parentNode->findvalue('@uid');
        }

        my $range;
        my $type;
        if ( $domain->findvalue('@scop_implicit_range') eq 'true' ) {
            $range                               = 'imp';
            $type                                = 'implicit';
            $old_domain_domainid2uid{$domain_id} = $uid;
        }
        elsif ( $domain->findvalue('@manual_range') eq 'true' && !$domain->exists('@yx_side_load') ) {
            $range = $domain->findvalue('manual_range');
            $type  = 'explicit';
        }
        elsif ( $domain->findvalue('@yx_side_load') eq 'true' ) {
            $range = $domain->findvalue('manual_range');
            $type  = 'side_load';
        }
        elsif ( $domain->findvalue('@derived_range') eq 'true' ) {
            $range = $domain->findvalue('derived_range');
            $type  = 'derived';

            #my $derivedFrom = $domain->findvalue('derived_range/@derivedFrom');
            my $derivedFrom = $domain->findvalue('ecod_representative_domain/@uid');
            my $hash        = $domain_id . $range;
            $old_domain_derivedFrom{$hash}    = $derivedFrom;
            $old_domain_derivedFrom_uid{$uid} = $derivedFrom;
        }
        elsif ( $domain->findvalue('@curated_range') eq 'true' ) {

            #I don't believe this is currently working
            $range = $domain->findvalue('curated_range');
            $type  = 'curated';
            my $derivedFrom = $domain->findvalue('ecod_representative_domain/@uid');
            my $hash        = $domain_id . $range;
            $old_domain_derivedFrom{$hash}    = $derivedFrom;
            $old_domain_derivedFrom_uid{$uid} = $derivedFrom;
        }
        else {
            die "ERROR! Unable to discern domain type of $domain_id, $uid\n";
        }
        my $hash = $domain_id . $range;

        my $f_id;
        if ( $domain->parentNode->exists('@f_id') ) {
            $f_id = $domain->parentNode->findvalue('@f_id');
        }
        elsif ( $domain->parentNode->parentNode->exists('@f_id') ) {
            $f_id = $domain->parentNode->parentNode->findvalue('@f_id');
        }
        elsif ( $domain->parentNode->parentNode->parentNode->exists('@f_id') ) {
            $f_id = $domain->parentNode->parentNode->parentNode->findvalue('@f_id');
        }
        else {
            die "ERROR! No f_id for domain $uid\n";
        }

        if ( $domain->findvalue('@manual_rep') ne 'true' ) {
            my $parent_domain_uid;
            if ( $domain->exists('ecod_representative_domain/@uid') ) {
                $parent_domain_uid = $domain->findvalue('ecod_representative_domain/@uid');
            }
            elsif ( $domain->exists('ecod_representative_domain/@ecod_domain_id') ) {
                my $parent_domain_id = $domain->findvalue('ecod_representative_domain/@ecod_domain_id');
                $parent_domain_uid = $old_domain_domainid2uid{$parent_domain_id};
            }
            elsif ( $domain->exists('derived_range/@derivedFrom') ) {
                $parent_domain_uid = $domain->findvalue('derived_range/@derivedFrom');
            }
            elsif ( $domain->exists('scop_representative_domain/@scop_domain_id') ) {
                my $parent_domain_id = $domain->findvalue('scop_representative_domain/@scop_domain_id');
                $parent_domain_uid = $old_domain_domainid2uid{$parent_domain_id};

            }
            else {
                print "WARNING! No representative notation for $uid $domain_id\n";
            }

            $new_domain_parent_uid{$hash} = $parent_domain_uid;
        }

        $old_domain_seen{$hash}++;
        $old_domain_fid{$hash}    = $f_id;
        $old_domain_uid_fid{$uid} = $f_id;
        $old_domain_range{$hash}  = $range;
        $old_domain_uid{$hash}    = $uid;
        $old_domain_type{$hash}   = $type;
        $old_domain_unhash{$hash} = $domain_id;
        $old_domain_da{$hash}     = $da_uid;
    }

    print "Generate old domain rules\n";
    foreach my $old_domain ( sort { $a cmp $b } keys %old_domain_seen ) {
        if ( $old_domain_da{$old_domain} ) { next }    #da moves only with da now
        if ( $new_domain_seen{$old_domain} && $old_domain_range{$old_domain} eq $new_domain_range{$old_domain} ) {

#print "DEBUG: nds $old_domain $new_domain_seen{$old_domain} $old_domain_fid{$old_domain} $new_domain_fid{$old_domain}\n";
            if ( $old_domain_fid{$old_domain} eq $new_domain_fid{$old_domain} ) {
                if ( !$new_domain_da{$old_domain} ) {
                    if ( $old_domain_type{$old_domain} eq $new_domain_type{$old_domain} ) {

                        #This is not correctly distinguishing domains that are changed from non-rep to rep
                        if ($DEBUG) {
                            print
"DOMAIN_CORRESPONDENCE $old_domain $new_domain_uid{$old_domain} $old_domain_uid{$old_domain} $new_domain_range{$old_domain} $old_domain_range{$old_domain} $old_domain_type{$old_domain} $new_domain_type{$old_domain}\n";
                        }

                        my $modify_domain_node = $merge_xml_doc->createElement('modify');
                        $modify_domain_node->setAttribute( 'type',     'domain' );
                        $modify_domain_node->setAttribute( 'mod_type', 'inherit_annotation' );
                        my $real_domain_id = $old_domain_unhash{$old_domain};
                        $modify_domain_node->setAttribute( 'domain_id', $real_domain_id );
                        $modify_domain_node->setAttribute( 'new_uid',   $new_domain_uid{$old_domain} );
                        $modify_domain_node->setAttribute( 'old_uid',   $old_domain_uid{$old_domain} );

                        $modify_node->appendChild($modify_domain_node);
                    }
                    elsif ( $old_domain_type{$old_domain} eq 'derived' && $new_domain_type{$old_domain} ) {

                        #CAse is broken or unnecessary?
                        if ($DEBUG) {
                            print
"DOMAIN_CONVERTED_TO_REP $old_domain $new_domain_uid{$old_domain} $old_domain_uid{$old_domain} $new_domain_range{$old_domain} $old_domain_range{$old_domain}\n";
                        }

                        my $modify_domain_node = $merge_xml_doc->createElement('modify');
                        $modify_domain_node->setAttribute( 'type',     'domain' );
                        $modify_domain_node->setAttribute( 'mod_type', 'nonrep_to_rep_conversion' );
                        my $real_domain_id = $old_domain_unhash{$old_domain};
                        $modify_domain_node->setAttribute( 'domain_id', $real_domain_id );
                        $modify_domain_node->setAttribute( 'new_uid',   $new_domain_uid{$old_domain} );
                        $modify_domain_node->setAttribute( 'old_uid',   $old_domain_uid{$old_domain} );

                        $modify_domain_node->setAttribute( 'new_domain_type', $new_domain_type{$old_domain} );
                        $modify_domain_node->setAttribute( 'old_domain_type', $old_domain_type{$old_domain} );

                        $modify_node->appendChild($modify_domain_node);

                    }
                }
                else {
                    if ( $old_domain_type{$old_domain} eq $new_domain_type{$old_domain} ) {
                        if ($DEBUG) {
                            print
"DOMAIN_CONVERTER_TO_ASSEMBLY_MEMBER $old_domain $new_domain_uid{$old_domain} $old_domain_uid{$old_domain} $new_domain_range{$old_domain} $old_domain_range{$old_domain}\n";
                        }
                        my $modify_domain_node = $merge_xml_doc->createElement('modify');
                        $modify_domain_node->setAttribute( 'type',     'domain' );
                        $modify_domain_node->setAttribute( 'mod_type', 'convert_to_domain_assembly' );
                        my $real_domain_id = $old_domain_unhash{$old_domain};
                        $modify_domain_node->setAttribute( 'domain_id',       $real_domain_id );
                        $modify_domain_node->setAttribute( 'new_uid',         $new_domain_uid{$old_domain} );
                        $modify_domain_node->setAttribute( 'old_uid',         $old_domain_uid{$old_domain} );
                        $modify_domain_node->setAttribute( 'old_parent_f_id', $old_domain_fid{$old_domain} );
                        $modify_domain_node->setAttribute( 'da',              $new_domain_da{$old_domain} );
                        $modify_domain_node->setAttribute( 'assembly_type',   $new_domain_assembly{$old_domain} );

                        $modify_node->appendChild($modify_domain_node);
                    }
                }

            }
            else {
                if ($DEBUG) {
                    print
"DOMAIN_FSHIFT $old_domain $new_domain_uid{$old_domain} $old_domain_uid{$old_domain} $new_domain_range{$old_domain} $old_domain_range{$old_domain}\n";
                }

                my $modify_domain_node = $merge_xml_doc->createElement('modify');
                $modify_domain_node->setAttribute( 'type',     'domain' );
                $modify_domain_node->setAttribute( 'mod_type', 'f_shift' );
                my $real_domain_id = $old_domain_unhash{$old_domain};
                $modify_domain_node->setAttribute( 'domain_id', $real_domain_id );
                $modify_domain_node->setAttribute( 'new_uid',   $new_domain_uid{$old_domain} );
                $modify_domain_node->setAttribute( 'old_uid',   $old_domain_uid{$old_domain} );
                $modify_domain_node->setAttribute( 'new_fid',   $new_domain_fid{$old_domain} );
                $modify_domain_node->setAttribute( 'old_fid',   $old_domain_fid{$old_domain} );

                $modify_domain_node->setAttribute( 'old_domain_type', $old_domain_type{$old_domain} );
                $modify_domain_node->setAttribute( 'new_domain_type', $new_domain_type{$old_domain} );

                $modify_node->appendChild($modify_domain_node);
            }

        }
        elsif ( $old_domain_type{$old_domain} eq 'side_load' ) {

#Side load domains are weird, they are like derived domains in that they do not have a new domain uid explicitly, but they also have no parent.
            if ($DEBUG) {
                print "SIDE_LOAD $old_domain $old_domain_uid{$old_domain}\n";
            }

            my $modify_domain_node = $merge_xml_doc->createElement('modify');
            $modify_domain_node->setAttribute( 'type',     'domain' );
            $modify_domain_node->setAttribute( 'mod_type', 'side_load_shift' );
            if ( $new_domain_da{$old_domain} || $old_domain_da{$old_domain} ) { next }
            my $real_domain_id = $old_domain_unhash{$old_domain};
            $modify_domain_node->setAttribute( 'domain_id', $real_domain_id );
            $modify_domain_node->setAttribute( 'old_uid',   $old_domain_uid{$old_domain} );

            #why do side load domains have a parent f_id? they are not derived...
            #$modify_domain_node->setAttribute('old_parent_f_id', $old_domain_uid_fid{$old_domain_uid{$old_domain}});

            $modify_node->appendChild($modify_domain_node);

        }
        elsif ( $old_domain_type{$old_domain} eq 'derived' ) {

            #print "DERIVED_DOMAIN $old_domain $old_domain_uid{$old_domain}\n";

            my $modify_domain_node = $merge_xml_doc->createElement('modify');
            $modify_domain_node->setAttribute( 'type',     'domain' );
            $modify_domain_node->setAttribute( 'mod_type', 'derived_shift' );
            my $real_domain_id = $old_domain_unhash{$old_domain};
            $modify_domain_node->setAttribute( 'domain_id',      $real_domain_id );
            $modify_domain_node->setAttribute( 'old_uid',        $old_domain_uid{$old_domain} );
            $modify_domain_node->setAttribute( 'old_parent_uid', $old_domain_derivedFrom{$old_domain} );
            my $rep_domain_uid = $old_domain_derivedFrom{$old_domain};
            my $parent_f_id    = $old_domain_uid_fid{$rep_domain_uid};
            my $MAX            = 10;
            my $i              = 0;
            while ( $old_domain_derivedFrom_uid{$rep_domain_uid} && $i < $MAX ) {
                $rep_domain_uid = $old_domain_derivedFrom_uid{$rep_domain_uid};
                $parent_f_id    = $old_domain_uid_fid{$rep_domain_uid};

                #print "$rep_domain_uid $parent_f_id\n";
                $i++;
            }
            $modify_domain_node->setAttribute( 'old_parent_f_id', $parent_f_id );
            $modify_node->appendChild($modify_domain_node);

        }
        else {
            #print "OBSOLETE_DOMAIN $old_domain $old_domain_uid{$old_domain} $old_domain_parent{$old_domain}\n";

            my $obsolete_domain_node = $merge_xml_doc->createElement('obsolete');
            $obsolete_domain_node->setAttribute( 'type', 'domain' );
            $obsolete_domain_node->setAttribute( 'uid',  $old_domain_uid{$old_domain} );
            my $real_domain_id = $old_domain_unhash{$old_domain};
            $obsolete_domain_node->setAttribute( 'domain_id', $real_domain_id );
            if ( $old_domain_uid_fid{ $old_domain_uid{$old_domain} } ) {
                $obsolete_domain_node->setAttribute( 'parent_f_id',
                    $old_domain_uid_fid{ $old_domain_uid{$old_domain} } );
            }
            else {
                print "WARNING! No parent f_id for $old_domain\n";
            }

            $obsolete_node->appendChild($obsolete_domain_node);
        }
    }

    print "Generate new domain rules\n";
    foreach my $new_domain ( sort { $a cmp $b } keys %new_domain_seen ) {

        if ( $old_domain_seen{$new_domain} && $old_domain_range{$new_domain} eq $new_domain_range{$new_domain} ) {

            #print "DOMAIN_CORRESPONDENCE $new_domain $new_domain_uid{$new_domain} $old_domain_uid{$new_domain}\n";
        }
        else {
            if ($DEBUG) {
                print "NEW_DOMAIN $new_domain $new_domain_uid{$new_domain} $new_domain_fid{$new_domain}\n";
            }

            my $new_domain_node = $merge_xml_doc->createElement('new');
            $new_domain_node->setAttribute( 'type', 'domain' );
            my $real_domain_id = $new_domain_unhash{$new_domain};
            $new_domain_node->setAttribute( 'domain_id', $real_domain_id );
            $new_domain_node->setAttribute( 'uid',       $new_domain_uid{$new_domain} );

            $new_node->appendChild($new_domain_node);
        }
    }

    return $merge_xml_doc;
}

sub ecod_group_summary {
    my $sub = 'ecod_group_summary';

    my ($ecod_xml_doc) = @_;

    my $x_groups          = $ecod_xml_doc->findnodes('//x_group')->size();
    my $h_groups          = $ecod_xml_doc->findnodes('//h_group')->size();
    my $f_groups          = $ecod_xml_doc->findnodes('//f_group')->size();
    my $domains           = $ecod_xml_doc->findnodes('//domain')->size();
    my $domain_assemblies = $ecod_xml_doc->findnodes('//domain_assembly')->size();

    if ($DEBUG) {
        printf "DEBUG $sub: x:%i h:%i f:%i d:%i da;%i\n", $x_groups, $h_groups, $f_groups, $domains, $domain_assemblies;
    }
}

sub apply_merge_v2 {
    my $sub = 'apply_merge_v2';

    my ( $old_ecod_xml_fn, $new_ecod_xml_fn, $merge_xml_fn ) = @_;

    print "$sub: open files\n";
    my $old_ecod_xml_doc = xml_open($old_ecod_xml_fn);
    my $new_ecod_xml_doc = xml_open($new_ecod_xml_fn);
    my $merge_xml_doc    = xml_open($merge_xml_fn);

    #Index
    print "$sub: index generate old\n";
    my ($old_index) = ecod_index($old_ecod_xml_doc);
    print "$sub: index generate new\n";
    my ($new_index) = ecod_index($new_ecod_xml_doc);

    #Record new version
    record_new_version( $old_ecod_xml_doc, $new_ecod_xml_doc );

    my $old_domain_dictionary_node = get_domain_dictionary_node($old_ecod_xml_doc);
    my $max_uid                    = get_max_uid($old_domain_dictionary_node);

    #Obsolete old nodes
    print "$sub: Generate obsolete merge_index\n";
    my $obsolete_merge_index = obsolete_merge_index($merge_xml_doc);
    print "$sub: Obsolete old nodes\n";
    annotate_obsolete_ecod_nodes( $old_index, $obsolete_merge_index );    #This may be vestigial

    #Index modifications
    print "$sub: Generate modify merge index\n";
    my $mod_merge_index = modification_merge_index($merge_xml_doc);
    print "$sub: Generate obsolete merge _index\n";
    my $new_merge_index = new_merge_index( $merge_xml_doc, $new_index );

    print "$sub: Apply rules\n";
    add_new_mod_old_ecod_nodes( $old_index, $new_index, $new_merge_index, $mod_merge_index, $obsolete_merge_index,
        $old_ecod_xml_doc );

    print "$sub: Remove obsolete nodes\n";
    remove_obsolete_ecod_nodes( $old_index, $obsolete_merge_index );

}

sub remove_obsolete_ecod_nodes {
    my ( $old_index, $obs_merge_index ) = @_;

    foreach my $obs_d_uid ( @{ $$obs_merge_index{d} } ) {
        $$old_index{d}{$obs_d_uid}->setAttribute( 'isObsolete', 'true' );

    }

    foreach my $x_id ( @{ $$obs_merge_index{x} } ) {
        my $x_node  = $$old_index{x}{$x_id};
        my $warning = 0;
        foreach my $child_h_node ( find_h_group_nodes($x_node) ) {
            my $h_id = get_h_id($child_h_node);
            if ( !isObsolete($child_h_node) ) {
                print "WARNING! Obsoletion of X-group $x_id orphans non-obsolete H-group $h_id\n";
                $warning++;
            }
        }
        if ( !$warning ) {
            $x_node->unbindNode;
        }
    }

    foreach my $h_id ( @{ $$obs_merge_index{h} } ) {
        my $h_node  = $$old_index{h}{$h_id};
        my $warning = 0;
        foreach my $child_f_node ( find_f_group_nodes($h_node) ) {
            my $f_id = get_f_id($child_f_node);
            if ( !isObsolete($child_f_node) ) {
                print "WARNING! Obsoletion of H-group $h_id orphans non-obsolete F-group $f_id\n";
                $warning++;
            }
        }
        if ( !$warning ) {
            $h_node->unbindNode;
        }
    }

    foreach my $f_id ( @{ $$obs_merge_index{f} } ) {
        my $f_node  = $$old_index{f}{$f_id};
        my $warning = 0;
        foreach my $child_domain_node ( find_domain_nodes($f_node) ) {
            my $child_uid = get_uid($child_domain_node);
            if ( !isObsolete($child_domain_node) ) {
                if ( $child_domain_node->findvalue('derived_range') eq 'true' ) {
                    print "WARNING! Obsoletion of F-group $f_id orphans derived domain $child_uid\n";
                }
                else {
                    print "WARNING! Obsoletion of F-group $f_id orphans non-obsolete domain $child_uid\n";
                }
                $warning++;
            }
        }
        if ( !$warning ) {
            $f_node->unbindNode;
        }
    }
}

sub add_new_mod_old_ecod_nodes {
    my ( $old_index, $new_index, $new_merge_index, $mod_merge_index, $obs_merge_index, $old_ecod_xml_doc ) = @_;

    my $old_domain_dictionary_node = $$old_index{dd};

    #Architectures added first
    foreach my $new_arch ( @{ $$new_merge_index{a} } ) {

        my $new_arch_name = $$new_arch{arch_name};

        #Name collision check for architectuere is necessary because of lack of id in manual table
        if ( has_arch_name( $old_ecod_xml_doc, $new_arch_name ) ) {
            die "ERROR! old->new arch_name conflict on $new_arch_name\n";
        }

        my $new_arch_node = create_new_arch( $$new_arch{arch_id}, $$new_arch{arch_name}, $old_ecod_xml_doc );
        $old_domain_dictionary_node->appendChild($new_arch_node);

        $$old_index{arch}{ $$new_arch{arch_id} } = $new_arch_node;

    }

    #Add new X-groups
    foreach my $new_x ( @{ $$new_merge_index{x} } ) {
        my $x_id = $$new_x{x_id};

        #		if (has_x_id($old_ecod_xml_doc, $x_id)) {
        #			die "ERROR! old->new x_group id conflict on $x_id\n";
        #		}
        if ( $$old_index{x}{$x_id} ) {
            die "ERROR! old->new x_group id conflict on $x_id\n";
        }

        my $new_x_group_node = create_new_x_group( $$new_x{x_id}, $$new_x{x_name}, $old_ecod_xml_doc );

        my $x_arch_id = $$new_x{xarch};
        $$old_index{arch}{$x_arch_id}->appendChild($new_x_group_node);
        $$old_index{x}{$x_id} = $new_x_group_node;

    }

    #Modify old X-groups
    foreach my $mod_x ( @{ $$mod_merge_index{x} } ) {
        my $x_id   = $$mod_x{x_id};
        my $x_node = $$old_index{x}{$x_id};

        if ( $$mod_x{remove_name} ) {
            $x_node->removeAttribute('name');
            $x_node->setAttribute( 'tmp_name_remove', 'true' );    #Why?
            if ($DEBUG) {
                print "DEBUG: mod_x remove name $x_id\n";
            }
        }
        if ( $$mod_x{new_name} ) {
            $x_node->setAttribute( 'name',            $$mod_x{new_name} );
            $x_node->setAttribute( 'tmp_name_change', 'true' );              #Again why?
            if ($DEBUG) {
                print "DEBUG: mod_x new_name $x_id $$mod_x{new_name}\n";
            }
        }
        if ( $$mod_x{x_shift} ) {
            my $new_arch_name = $$mod_x{new_arch_name};
            my $new_arch_id   = $$mod_x{new_arch_id};
            my $old_arch_name = $$mod_x{old_arch_name};
            my $old_arch_id   = $$mod_x{old_arch_id};

            my $trans_arch_id = $$new_index{arch_lookup}{$new_arch_name};
            my $new_arch_node = $$old_index{arch}{$trans_arch_id};
            if ($DEBUG) {
                print "DEBUG: mod_x x_shift $x_id\n";
            }

            $new_arch_node->appendChild($x_node);
        }
        if ( $$mod_x{ordinal_shift} ) {
            $x_node->setAttribute( 'x_ordinal', $$mod_x{new_ordinal} );
        }
    }

    #Obsolete old X-groups, why does this happen twice?
    foreach my $obs_x ( @{ $$obs_merge_index{x} } ) {
        if ( !$$old_index{x}{ $$obs_x{x_id} } ) {
            die "ERROR! Could not find obs x_group for $$obs_x{x_id}\n";
        }

        my $obs_x_node = $$old_index{x}{ $$obs_x{x_id} };
        $obs_x_node->setAttribute( 'isObsolete', 'true' );
    }

    #New H-groups
    foreach my $new_h ( @{ $$new_merge_index{h} } ) {
        my $new_h_id   = $$new_h{h_id};
        my $new_h_name = $$new_index{h}{name};
        if ( exists $$new_index{h}{$new_h_id} ) {
            die "ERROR! old->new h_group id conflict on $new_h_id...\n";
        }

        my $new_h_group_node = create_new_h_group( $new_h_id, $new_h_name, $old_ecod_xml_doc );

        my $parent_x_id = get_parent_id($new_h_id);

        my $parent_x_group_node =
          exists $$old_index{x}{$parent_x_id}
          ? $$old_index{x}{$parent_x_id}
          : die "ERROR! Could not find parent X-group node $parent_x_id for H-group $new_h_id\n";

        $parent_x_group_node->appendChild($new_h_group_node);
        $$old_index{h}{$new_h_id} = $new_h_group_node;
    }

    #Modify h_groups;
    foreach my $mod_h ( @{ $$mod_merge_index{h} } ) {

        my $mod_h_id   = $$mod_h{h_id};
        my $mod_h_node = $$old_index{h}{$mod_h_id};

        if ( $$mod_h{remove_name} ) {
            $mod_h_node->removeAttribute('name');
            $mod_h_node->setAttribute( 'tmp_name_remove', 'true' );    #Again why?
        }
        elsif ( $$mod_h{new_name} ) {
            $mod_h_node->setAttribute( 'name',            $$mod_h{new_name} );
            $mod_h_node->setAttribute( 'tmp_name_change', 'true' );              #Again why?
        }

        if ( $$mod_h{ordinal_shift} ) {
            $mod_h_node->setAttribute( 'h_ordinal', $$mod_h{new_ordinal} );
        }
    }

    #Obsolete h_groups
    foreach my $obs_h ( @{ $$obs_merge_index{h} } ) {
        my $obs_h_id   = $$obs_h{h_id};
        my $obs_h_node = $$old_index{h}{$obs_h_id};

        if ( !exists $$old_index{h}{$obs_h_id} ) {
            die "ERROR! COuld not find obs h_group for $obs_h_id\n";
        }

        $obs_h_node->setAttribute( 'isObsolete', 'true' );    #Are these being annotated twice?
    }

    #New F-groups
    foreach my $new_f ( @{ $$new_merge_index{f} } ) {
        my $new_f_id      = $$new_f{f_id};
        my $new_f_name    = $$new_index{f}{$new_f_id}{f_name};
        my $new_f_ordinal = $$new_f{new_ordinal};

        if ( exists $$old_index{f}{$new_f_id} && $$old_index{f}{$new_f_id}->findvalue('@isNew') eq 'true' ) {
            die "ERROR! old->new f_group id conflict on $new_f_id...?\n";
        }                                                     #Has this ever triggered?

        my $new_f_group_node = create_new_f_group( $new_f_id, $new_f_name, $new_f_ordinal, $old_ecod_xml_doc );

        my $parent_h_id = get_parent_id($new_f_id);

        my $parent_h_group_node =
          exists $$old_index{x}{$parent_h_id}
          ? $$old_index{h}{$parent_h_id}
          : die "ERROR! Could not find parent H-group node $parent_h_id for F-group $new_f_id\n";

        $parent_h_group_node->appendChild($new_f_group_node);
        $$old_index{f}{$new_f_id} = $new_f_group_node;
    }

    #Modify F-groups
    foreach my $mod_f ( @{ $$mod_merge_index{f} } ) {
        my $mod_f_id   = $$mod_f{f_id};
        my $mod_f_node = $$old_index{f}{$mod_f_id};

        if ( $$mod_f{remove_name} ) {
            $mod_f_node->removeAttribute('name');
            $mod_f_node->setAttribute( 'tmp_name_remove', 'true' );
        }
        elsif ( $$mod_f{new_name} ) {
            $mod_f_node->setAttribute( 'name', $$mod_f{new_name} );
        }

        if ( $$mod_f{f_shift} ) {
            if ( $$old_index{f}{$mod_f_id} ) {
                print "WARNING! potential F-group Id collision at $mod_f_id\n";
                $mod_f_node->setAttribute( 'old_f_id', $mod_f_id );    #Ancient debug
            }
            $mod_f_node->setAttribute( 'f_id', $mod_f_id );

            my $parent_h_id = get_parent_id($mod_f_id);

            my $parent_h_group_node =
              exists $$old_index{h}{$parent_h_id}
              ? $$old_index{h}{$parent_h_id}
              : die "ERROR! Could not find parent H-group node $parent_h_id for F-group $mod_f_id\n";

            $parent_h_group_node->appendChild($mod_f_node);

            #$$old_index{f}{$mod_f_id} = $mod_f_node;
        }

        if ( $$mod_f{ordinal_shift} ) {
            $mod_f_node->setAttribute( 'f_ordinal', $$mod_f{new_ordinal} );
        }
    }

    #New domains
    my %domain_assembly_uid;
    foreach my $new_d ( @{ $$new_merge_index{d} } ) {
        my $new_d_uid = $$new_d{uid};
        my $sql_d_uid = sprintf "%09i", sql_uid_gen();

        my $new_d_node =
          exists $$new_index{d}{$new_d_uid}
          ? $$new_index{d}{$new_d_uid}
          : print "WARNING! Could not find new domain node for $new_d_uid, skipping...\n";

        if ( isDomainAssembly($new_d_node) ) {
            my $parent_f_id = get_f_id( $new_d_node->parentNode->parentNode );

            my $parent_f_node =
              exists $$old_index{f}{$parent_f_id}
              ? $$old_index{f}{parent_f_id}
              : die "ERROR! Domain assembly parent F-group $parent_f_id not found\n";

            my $old_domain_assembly_node;
            if ( $domain_assembly_uid{$new_d_uid} ) {
                my $da_uid = $domain_assembly_uid{$new_d_uid};

                $old_domain_assembly_node = $$old_index{da}{$da_uid};
            }
            else {
                $old_domain_assembly_node = create_domain_assembly_node($old_ecod_xml_doc);
                $parent_f_node->appendChild($old_domain_assembly_node);
                my $da_uid = get_uid($old_domain_assembly_node);

                foreach my $new_sib_node ( find_domain_nodes( $new_d_node->parentNode ) ) {
                    if ( $new_sib_node->isSameNode($new_d_node) ) { next }
                    else {
                        my $new_sib_uid = get_uid($new_sib_node);
                        $domain_assembly_uid{$new_sib_uid} = $da_uid;
                    }
                }
                $$old_index{da}{$da_uid} = $old_domain_assembly_node;
            }
            $old_domain_assembly_node->appendChild($new_d_node);
            my $sql_uid = sprintf "%09i", sql_uid_gen();
            $new_d_node->setAttribute( 'uid', $sql_uid );
            $$old_index{d}{$sql_uid} = $new_d_node;
        }
        else {
            my $parent_f_id = get_f_id( $new_d_node->parentNode );
            my $parent_f_node =
              exists $$old_index{f}{$parent_f_id}
              ? $$old_index{f}{parent_f_id}
              : die "ERROR! domain parent F-group $parent_f_id not found\n";

            $parent_f_node->appendChild($new_d_node);
            my $sql_uid = sprintf "%09i", sql_uid_gen();
            $new_d_node->setAttribute( 'uid', $sql_uid );
            $$old_index{d}{$sql_uid} = $new_d_node;
        }
    }

    #Modify domains
    my %current_da_uids;
    foreach my $mod_d ( @{ $$mod_merge_index{d} } ) {
        my $mod_d_uid   = $$mod_d{uid};
        my $mod_d_node  = $$old_index{d}{$mod_d_uid};
        my $parent_f_id = $$mod_d{parent_f_id};
        my $domain_id   = $$mod_d{domain_id};
        my $parent_uid  = $$mod_d{parent_uid};

        if ( $$mod_d{derived_shift} ) {

            #my $parent_uid 	= $$mod_d{parent_uid};
            #my $parent_f_id = $$mod_d{parent_f_id};
            #my $domain_id	= $$mod_d{domain_id};

            my $parent_domain_node = exists $$old_index{d}{$parent_uid} ? $$old_index{d}{$parent_uid} : next;

            my $parent_f_node = find_parent_f_group_node($parent_domain_node);

            my $current_f_node = find_parent_f_group_node($mod_d_node);

            my $new_parent_f_id  = get_f_id($parent_f_node);
            my $new_current_f_id = get_f_id($current_f_node);

            if ( $new_parent_f_id ne $parent_f_id || $new_current_f_id ne $parent_f_id ) {
                $parent_f_node->appendChild($mod_d_node);
            }

        }
        elsif ( $$mod_d{side_load} ) {

            #my $parent_f_id = $$mod_d{parent_f_id};
            #my $domain_id 	= $$mod_d{domain_id};

            my $parent_f_node =
              exists $$old_index{f}{$parent_f_id}
              ? $$old_index{f}{$parent_f_id}
              : next;

            #What is this case for?
        }
        elsif ( $$mod_d{inherit_annotation} ) {
            my $new_d_uid = $$mod_d{new_uid};
            my $new_d_node =
              exists $$new_index{d}{$new_d_uid}
              ? $$new_index{d}{$new_d_uid}
              : die "ERROR! No new domain node found for UID $new_d_uid\n";

            transfer_annotation( $new_d_node, $mod_d_node );
        }
        elsif ( $$mod_d{domain_fshift} ) {
            my $new_d_uid = $$mod_d{new_uid};
            my $new_d_node =
              exists $$new_index{d}{$new_d_uid}
              ? $$new_index{d}{$new_d_uid}
              : die "ERROR! No new domain node found for UID $new_d_uid\n";
            my $new_f_id   = $$mod_d{new_fid};
            my $old_f_node = $$old_index{f}{$new_f_id};
            $old_f_node->appendChild($mod_d_node);

            my $old_domain_type = $$mod_d{old_domain_type};
            my $new_domain_type = $$mod_d{new_domain_type};

            if ( !is_manual_representative($new_d_node) ) {
                if ( has_scop_representative($new_d_node) && !has_scop_representative($mod_d_node) ) {
                    $mod_d_node->appendChild( get_scop_rep_node($new_d_node) );
                }
            }

            if ( $old_domain_type eq 'derived' && $new_domain_type eq 'explicit' ) {
                $mod_d_node->removeAttribute('derived_range');

                foreach my $old_derived_range_node ( find_derived_range_nodes($mod_d_node) ) {
                    $old_derived_range_node->unbindNode;
                }
                foreach my $old_derived_seqid_range_node ( find_derived_seqid_range_nodes($mod_d_node) ) {
                    $old_derived_seqid_range_node->unbindNode;
                }
                if ( $mod_d_node->exists('manual_range') ) {
                    foreach my $m_node ( $mod_d_node->findnodes('manual_range') ) {
                        $m_node->unbindNode;
                    }
                }
                $mod_d_node->appendChild( $new_d_node->findnodes('manual_range')->get_node(1)->cloneNode(1) );
                $mod_d_node->setAttribute( 'manual_range', 'true' );
                $mod_d_node->setAttribute( 'manual_rep',   'true' );

            }
            elsif ( $old_domain_type eq 'derived' && $new_domain_type eq 'implicit' ) {
                print "ERROR! Impossible upconvert on $mod_d_uid -> $new_d_uid\n";
            }

            transfer_annotation( $new_d_node, $mod_d_node );
        }
        elsif ( $$mod_d{nonrep_conversion} ) {
            my $new_d_uid = $$mod_d{new_uid};
            my $new_d_node =
              exists $$new_index{d}{$new_d_uid}
              ? $$new_index{d}{$new_d_uid}
              : die "ERROR! No new domain node found for uid $new_d_uid\n";

            my $old_domain_type = $$mod_d{old_domain_type};
            my $new_domain_type = $$mod_d{new_domain_type};

            if ( $old_domain_type eq 'derived' && $new_domain_type eq 'explicit' ) {
                $mod_d_node->removeAttribute('derived_range');

                foreach my $old_derived_range_node ( find_derived_range_nodes($mod_d_node) ) {
                    $old_derived_range_node->unbindNode;
                }
                foreach my $old_derived_seqid_range_node ( find_derived_seqid_range_nodes($mod_d_node) ) {
                    $old_derived_seqid_range_node->unbindNode;
                }
                if ( $mod_d_node->exists('manual_range') ) {
                    foreach my $m_node ( $mod_d_node->findnodes('manual_range') ) {
                        $m_node->unbindNode;
                    }
                }
                $mod_d_node->appendChild( $new_d_node->findnodes('manual_range')->get_node(1)->cloneNode(1) );
                $mod_d_node->setAttribute( 'manual_range', 'true' );
                $mod_d_node->setAttribute( 'manual_rep',   'true' );

            }
            elsif ( $old_domain_type eq 'derived' && $new_domain_type eq 'implicit' ) {
                print "ERROR! Impossible upconvert on $mod_d_uid -> $new_d_uid\n";
            }
            else {
                die "ERROR! Confusion case 1\n";
            }

        }
        elsif ( $$mod_d{domain_assembly_conversion} ) {
            my $new_d_uid = $$mod_d{new_uid};

            my $new_d_node =
              exists $$new_index{d}{$new_d_uid}
              ? $$new_index{d}{$new_d_uid}
              : die "ERROR! No new domain node for uid $new_d_uid\n";

            my $mod_d_node =
              exists $$old_index{c}{$mod_d_uid}
              ? $$old_index{d}{$mod_d_uid}
              : die "ERROR! No old domain node for ui $mod_d_uid\n";

            my $da_uid        = $$mod_d{da_uid};
            my $parent_f_id   = $$mod_d{parent_f_id};
            my $parent_f_node = $$old_index{f}{$parent_f_id};

            #my $assembly_type = $new_assembly_types{$da_uid}; #Where is this set?
            my $assembly_type = $$new_index{da}{$da_uid}{assembly_type};

            if ( $current_da_uids{$da_uid} ) {

                #$old_domain_assembly_nodes{$current_da_uids{$da_uid}}->appendChild($mod_d_node);
                $$old_index{da}{ $current_da_uids{$da_uid} }->appendChild($mod_d_node);
            }
            else {
                my $domain_assembly_node = create_domain_assembly_node( 'true', $assembly_type, $old_ecod_xml_doc );
                my $uid = get_uid($domain_assembly_node);

                #$old_domain_assembly_nodes{$uid} = $domain_assembly_node;
                $$old_index{da}{$uid} = $domain_assembly_node;
                $parent_f_node->appendChild($domain_assembly_node);
                $domain_assembly_node->appendChild($mod_d_node);
                $mod_d_node->setAttribute( 'primary', 'true' );
                $current_da_uids{$da_uid} = $uid;

            }
        }
    }
}

sub create_domain_assembly_node {
    my ( $isAssemblyRep, $assembly_type, $xml_doc ) = @_;

    my $assembly_node = $xml_doc->createElement('domain_assembly');
    my $sql_uid = sprintf "%09i", sql_uid_gen();
    $assembly_node->setAttribute( 'uid', $sql_uid );
    if ($isAssemblyRep) {
        $assembly_node->setAttribute( 'domain_assembly_representative', 'true' );
    }
    $assembly_node->setAttribute( 'assembly_type', $assembly_type );
    return $assembly_node;

}

sub transfer_annotation {
    my ( $from_domain_node, $to_domain_node ) = @_;

    #manual rep status;
    if ( $from_domain_node->exists('@manual_rep') ) {
        $to_domain_node->setAtttribute( 'manual_rep', $from_domain_node->findvalue('@manual_rep') );
    }

    #has SCOP representative
    if ( !is_manual_representative($from_domain_node) ) {
        if ( has_scop_representative($from_domain_node) && !has_scop_representative($to_domain_node) ) {
            $to_domain_node->appendChild( get_scop_rep_node($from_domain_node) );
        }
    }

    #comments
    foreach my $comment_type (qw(alert_comment scop_comment comment)) {
        if ( has_comment($from_domain_node), $comment_type ) {
            if ( !has_comment( $to_domain_node, $comment_type ) ) {
                $to_domain_node->appendChild( get_comment_node( $from_domain_node, $comment_type ) );
            }
            elsif ( get_comment( $to_domain_node, $comment_type ) eq get_comment( $from_domain_node, $comment_type ) ) {

                #do nothing
            }
            else {
                foreach my $comment_node ( get_comment_nodes( $to_domain_node, $comment_type ) ) {
                    $comment_node->unbindNode;
                }
                $to_domain_node->appendChild( get_comment_node( $from_domain_node, $comment_type ) );
            }
        }
    }

}

sub has_comment {
    if ( $_[1] eq 'alert_comment' ) {
        $_[0]->exists('alert_comment');
    }
    elsif ( $_[1] eq 'scop_comment' ) {
        $_[0]->exists('scop_comment');
    }
    elsif ( $_[1] eq 'comment' ) {
        $_[0]->exists('comment');
    }
}

sub get_comment_node {
    if ( $_[1] eq 'alert_comment' ) {
        $_[0]->findnodes('alert_comment')->get_node(1);
    }
    elsif ( $_[1] eq 'scop_comment' ) {
        $_[0]->findnodes('scop_comment')->get_node(1);
    }
    elsif ( $_[1] eq 'comment' ) {
        $_[0]->findnodes('comment')->get_node(1);
    }
}

sub is_manual_representative {
    $_[0]->findvalue('@manual_rep') eq 'true';
}

sub has_scop_representative {
    $_[0]->exists('scop_representative_domain');
}

sub get_scop_rep_node {
    $_[0]->findnodes('scop_representative_domain')->get_node(1);
}

sub find_parent_f_group_node {
    $_[0]->nodeName eq 'domain_dictionary' && return 0;
    $_[0]->nodeName eq 'f_group' ? return $_[0] : find_parent_f_group_node( $_[0]->parentNode );
}

sub get_parent_id {
    $_[0] =~ /(\d+)\.\d+/ or die "ERROR! Weird id $_[0]\n";
    return $1;
}

#sub create_domain_assembly_node {
#	my ($xml_doc) = @_;
#	my $domain_assembly_node = $xml_doc->createElement('domain_assembly');
#	$domain_assembly_node->setAttribute('uid', sprintf "%09i", sql_uid_gen());
#	return $domain_assembly_node;
#}

sub create_new_f_group {
    my ( $f_id, $f_name, $f_ordinal, $xml_doc ) = @_;
    my $new_f_group_node = $xml_doc->createElement('f_group');
    $new_f_group_node->setAttribute( 'f_id',      $f_id );
    $new_f_group_node->setAttribute( 'name',      $f_name );
    $new_f_group_node->setAttribute( 'f_ordinal', $f_ordinal );
    $new_f_group_node->setAttribute( 'isNew',     'true' );       #Is this still necessary?
    return $new_f_group_node;
}

sub create_new_h_group {
    my ( $h_id, $h_name, $xml_doc ) = @_;
    my $new_h_group_node = $xml_doc->createElement('h_group');
    $new_h_group_node->setAttribute( 'h_id', $h_id );
    $new_h_group_node->setAttribute( 'name', $h_name );
    return $new_h_group_node;
}

sub create_new_x_group {
    my ( $x_id, $x_name, $xml_doc ) = @_;
    my $new_x_group_node = $xml_doc->createElement('x_group');
    $new_x_group_node->setAttribute( 'x_id',   $x_id );
    $new_x_group_node->setAttribute( 'x_name', $x_name );
    return $new_x_group_node;

}

sub create_new_arch {
    my ( $arch_id, $arch_name, $xml_doc ) = @_;
    my $new_arch_node = $xml_doc->createElement('architectuure');
    $new_arch_node->setAttribute( 'arch_id',   $arch_id );
    $new_arch_node->setAttribute( 'arch_name', $arch_name );
    return $new_arch_node;
}

sub has_arch_name {
    $_[0]->exists(qq{.//architecture[\@arch_name="$_[1]"]});
}

sub new_merge_index {
    my ( $merge_xml_doc, $new_index ) = @_;

    my %index;
    foreach my $new_node ( find_new_merge_nodes($merge_xml_doc) ) {
        my $ecod_node_type = $new_node->findvalue('@type');

        if ( $ecod_node_type eq 'x_group' ) {
            my $x_id   = get_x_id($new_node);
            my $x_name = get_name($new_node);

            my ( $x_arch_id, $x_arch_name ) = get_arch_attrs( $$new_index{x}{$x_id}->parentNode );
            my $x_ordinal = get_x_ordinal( $$new_index{x}{$x_id} );

            my %x_index;
            $x_index{xarch}     = $x_arch_id;
            $x_index{x_name}    = $x_name;
            $x_index{arch_name} = $x_arch_name;
            $x_index{x_ordinal} = $x_ordinal;
            $x_index{x_id}      = $x_id;
            push( @{ $index{x} }, \%x_index );

        }
        elsif ( $ecod_node_type eq 'arch' ) {

            my ( $arch_id, $arch_name ) = get_arch_attrs($new_node);
            my $arch_ordinal = get_arch_ordinal($new_node);

            my %arch_index;
            $arch_index{arch_id}      = $arch_id;
            $arch_index{arch_name}    = $arch_name;
            $arch_index{arch_ordinal} = $arch_ordinal;
            push( @{ $index{a} }, \%arch_index );

        }
        elsif ( $ecod_node_type eq 'h_group' ) {
            my $h_id      = get_h_id($new_node);
            my $h_name    = get_name($new_node);
            my $h_ordinal = get_h_ordinal($new_node);

            my %h_index;
            $h_index{h_id}      = $h_id;
            $h_index{h_name}    = $h_name;
            $h_index{h_ordinal} = $h_ordinal;

            push( @{ $index{h} }, \%h_index );

        }
        elsif ( $ecod_node_type eq 'f_group' ) {
            my $f_id      = get_f_id($new_node);
            my $f_name    = get_name($new_node);
            my $f_ordinal = get_f_ordinal($new_node);

            my %f_index;
            $f_index{f_id}      = $f_id;
            $f_index{f_name}    = $f_name;
            $f_index{f_ordinal} = $f_ordinal;

            push( @{ $index{f} }, \%f_index );

        }
        elsif ( $ecod_node_type eq 'domain' ) {
            my $domain_id  = get_domain_id($new_node);
            my $domain_uid = get_uid($new_node);

            my %domain_index;
            $domain_index{domain_id} = $domain_id;
            $domain_index{uid}       = $domain_uid;

            push( @{ $index{d} }, \%domain_index );

        }
        else {
            print "WARNING! No rules for ecod_node_type $ecod_node_type, skipping...\n";
        }
    }
    return \%index;
}

sub get_arch_ordinal {
    $_[0]->findvalue('@arch_ordinal');
}

sub get_arch_attrs {
    ( $_[0]->findvalue('@arch_id'), $_[0]->findvalue('@arch_name') );
}

sub annotate_obsolete_ecod_nodes {
    my ( $index, $obs_merge_index ) = @_;

    foreach my $type (qw{a x h f d da}) {
        if ( exists $$obs_merge_index{$type} ) {
            foreach my $id ( @{ $$obs_merge_index{$type} } ) {
                setObsolete( $$index{$type}{$id} );
            }
        }
    }
}

sub setObsolete {
    $_[0]->setAttribute( 'obsolete', 'true' );
}

sub modification_merge_index {
    my ( $merge_xml_doc, $old_version, $new_version ) = @_;

    my %index;
    foreach my $mod_node ( find_modification_merge_nodes($merge_xml_doc) ) {
        my $ecod_node_type    = $mod_node->findvalue('@type');
        my $modification_type = $mod_node->findvalue('@mod_type');

        if ( $ecod_node_type eq 'x_group' ) {
            my $x_id = get_x_id($mod_node);
            if ( $modification_type eq 'name_change' ) {

                #my $old_x_name = get_name_by_version($mod_node, $old_version);
                my $new_x_name = get_name_by_version( $mod_node, $new_version );
                $index{x}{$x_id}{name_change} = 1;

                set_name( $index{x}{$x_id}, $new_x_name );
            }
            elsif ( $modification_type eq 'arch_shift' ) {

                my ( $old_arch_id, $old_arch_name ) = get_arch_attrs_by_version( $mod_node, $old_version );
                my ( $new_arch_id, $new_arch_name ) = get_arch_attrs_by_version( $mod_node, $old_version );

                $index{x}{$x_id}{x_shift} = 1;

                $index{x}{$x_id}{old_arch_id}   = $old_arch_id;
                $index{x}{$x_id}{old_arch_name} = $old_arch_name;
                $index{x}{$x_id}{new_arch_id}   = $new_arch_id;
                $index{x}{$x_id}{new_arch_name} = $new_arch_name;

            }
            elsif ( $modification_type eq 'ordinal_shift' ) {

                $index{x}{$x_id}{old_ordinal} = get_old_ordinal($mod_node);
                $index{x}{$x_id}{new_ordinal} = get_new_ordinal($mod_node);

                $index{x}{$x_id}{ordinal_shift} = 1;
            }
            else {
                print "WARNING! mod_type $modification_type not supported for node $ecod_node_type, skipping...\n";
            }
        }
        elsif ( $ecod_node_type eq 'h_group' ) {
            my $h_id = get_h_id($mod_node);
            if ( $modification_type eq 'name_change' ) {
                my $h_id = get_h_id($mod_node);    #schema issues here

                #my $old_h_name = get_name_by_version($mod_node, $old_version);
                my $new_h_name = get_name_by_version( $mod_node, $new_version );

                set_name( $index{h}{$h_id}, $new_h_name );
            }
            elsif ( $modification_type eq 'h_shift' ) {
                my ( $old_h_id, $new_h_id ) = get_mod_h_ids($mod_node);

                my $new_h_name = get_name_by_version( $mod_node, $new_version );

                $index{h}{$old_h_id}{h_shift} = 1;

                set_name( $index{h}{$h_id}, $new_h_name );
                $index{x}{$old_h_id}{new_h_id} = $new_h_id;
            }
            elsif ( $modification_type eq 'ordinal_shift' ) {
                my $h_id = get_h_id($mod_node);

                $index{h}{$h_id}{old_ordinal} = get_old_ordinal($mod_node);
                $index{h}{$h_id}{new_ordinal} = get_new_ordinal($mod_node);

                $index{h}{$h_id}{ordinal_shift} = 1;

            }
            else {
                print "WARNING! mod type $modification_type not supported for $ecod_node_type, skipping...\n";
            }
        }
        elsif ( $ecod_node_type eq 'f_group' ) {
            my $f_id = get_f_id($mod_node);

            my $old_f_name = get_name_by_version( $mod_node, $old_version );
            my $new_f_name = get_name_by_version( $mod_node, $new_version );

            if ( $modification_type eq 'name_change' ) {

                $index{f}{$f_id}{name_change} = 1;

                set_name( $index{f}{$f_id}, $new_f_name );
            }
            elsif ( $modification_type eq 'f_shift' ) {
                my ( $old_f_id, $new_f_id ) = get_mod_f_ids($mod_node);

                $index{f}{$old_f_id}{old_ordinal} = get_old_ordinal { $mod_node };
                $index{f}{$new_f_id}{new_ordinal} = get_new_ordinal { $mod_node };

                $index{f}{$f_id}{ordinal_shift} = 1;
            }
            else {
                print "WARNING! mod type $modification_type not supported for $ecod_node_type, skipping...\n";
            }
        }
        elsif ( $ecod_node_type eq 'domain' ) {
            my $old_domain_id = get_domain_id($mod_node);

            my $old_domain_uid = get_old_uid($mod_node);
            my $new_domain_uid = get_new_uid($mod_node);

            my $old_domain_parent_f_id = get_old_domain_parent_f_id($mod_node);
            my $old_domain_parent_uid  = get_old_domain_parent_uid($mod_node);

            my ( $old_domain_type, $new_domain_type ) = get_domain_type_from_mod_node($mod_node);

            if ( $modification_type eq 'derived_shift' ) {

                $index{d}{$old_domain_uid}{derived_shift} = 1;
                $index{d}{$old_domain_uid}{domain_id}     = $old_domain_id;
                $index{d}{$old_domain_uid}{parent_uid}    = $old_domain_parent_uid;
                $index{d}{$old_domain_uid}{parent_f_id}   = $old_domain_parent_f_id;

            }
            elsif ( $modification_type eq 'inherit_annotation' ) {

                $index{d}{$old_domain_uid}{inherit_annotation} = 1;
                $index{d}{$old_domain_uid}{domain_id}          = $old_domain_id;
                $index{d}{$old_domain_uid}{new_uid}            = $new_domain_uid;
            }
            elsif ( $modification_type eq 'f_shift' ) {

                my $new_domain_fid = get_new_domain_fid($mod_node);

                $index{d}{$old_domain_uid}{domain_fshift} = 1;
                $index{d}{$old_domain_uid}{domain_id}     = $old_domain_id;
                $index{d}{$old_domain_uid}{new_uid}       = $new_domain_uid;
                $index{d}{$old_domain_uid}{new_fid}       = $new_domain_fid;

                $index{d}{$old_domain_uid}{old_domain_type} = $old_domain_type;
                $index{d}{$old_domain_uid}{new_domain_type} = $new_domain_type;

            }
            elsif ( $modification_type eq 'nonrep_to_rep_conversion' ) {

                $index{d}{$old_domain_uid}{nonrep_conversion} = 1;
                $index{d}{$old_domain_uid}{domain_id}         = $old_domain_id;
                $index{d}{$old_domain_uid}{new_uid}           = $new_domain_uid;
                $index{d}{$old_domain_uid}{old_domain_type}   = $old_domain_type;
                $index{d}{$old_domain_uid}{new_domain_type}   = $new_domain_type;
            }
            elsif ( $modification_type eq 'convert_to_domain_assembly' ) {

                my $new_da_uid    = get_mod_da_uid($mod_node);
                my $assembly_type = get_assembly_type($mod_node);

                $index{da}{$new_da_uid}{assembly_type} = $assembly_type;

                $index{d}{$old_domain_uid}{domain_assembly_conversion} = 1;
                $index{d}{$old_domain_uid}{domain_id}                  = $old_domain_id;
                $index{d}{$old_domain_uid}{new_uid}                    = $new_domain_uid;

                $index{d}{$old_domain_id}{da_uid}       = $new_da_uid;
                $index{d}{$old_domain_uid}{parent_f_id} = $old_domain_parent_f_id;
            }
            else {
                print "WARNING! modification type $modification_type not supported for $ecod_node_type, skipping...\n";
            }
        }
        elsif ( $ecod_node_type eq 'domain_assembly' ) {
            my $old_domain_id = get_domain_id($mod_node);
            my $old_da_uid    = get_old_uid($mod_node);
            my $new_da_uid    = get_new_uid($mod_node);

            my $old_da_parent_uid  = get_old_domain_parent_uid($mod_node);
            my $old_da_parent_f_id = get_old_domain_parent_f_id($mod_node);

            my $new_da_fid = get_new_fid($mod_node);

            if ( $modification_type eq 'derived_shift' ) {
                $index{da}{$old_da_uid}{derived_shift} = 1;

                $index{da}{$old_da_uid}{domain_id}   = $old_domain_id;
                $index{da}{$old_da_uid}{parent_uid}  = $old_da_parent_uid;
                $index{da}{$old_da_uid}{parent_f_id} = $old_da_parent_f_id;

            }
            elsif ( $modification_type eq 'inherit_annotation' ) {

                $index{da}{$old_da_uid}{inherit_annotation} = 1;
                $index{da}{$old_da_uid}{domain_id}          = $old_domain_id;
                $index{da}{$old_da_uid}{new_uid}            = $new_da_uid;

            }
            elsif ( $modification_type eq 'fshift' ) {

                $index{da}{$old_da_uid}{domain_fshift} = 1;
                $index{da}{$old_da_uid}{domain_id}     = $old_domain_id;
                $index{da}{$old_da_uid}{new_uid}       = $new_da_uid;
                $index{da}{$old_da_uid}{new_fid}       = $new_da_fid;

            }
            else {
                print
"WARNING! modification_type $modification_type not supported for node type $ecod_node_type, skipping...\n";

            }
        }
    }
    return \%index;
}

sub get_new_fid {
    $_[0]->findvalue('@new_fid');
}

sub get_assembly_type {
    $_[0]->findvalue('@assembly_type');
}

sub get_mod_da_uid {
    $_[0]->findvalue('@da');
}

sub get_domain_type_from_mod_node {
    ( $_[0]->findvalue('@old_domain_type'), $_[0]->findvalue('@new_domain_type') ) = @_;
}

sub get_new_domain_fid {
    $_[0]->findvalue('@new_fid');
}

sub get_old_uid {
    $_[0]->findvalue('@old_uid');
}

sub get_new_uid {
    $_[0]->findvalue('@new_uid');
}

sub get_old_domain_parent_f_id {
    $_[0]->findvalue('@old_parent_f_id');
}

sub get_old_domain_parent_uid {
    $_[0]->findvalue('@old_parent_uid');
}

sub get_mod_h_ids {
    ( $_[0]->findvalue('@old_h_id'), $_[0]->findvalue('@new_h_id') );
}

sub set_name {
    my ( $href, $name );
    if ( $name eq 'NA' ) {
        $$href{remove_name} = 1;
    }
    else {
        $$href{new_name} = $name;
    }
}

sub get_arch_attrs_by_version {
    ( $_[0]->findvalue(qq{arch_name[\@version="$_[1]"]/\@arch_id}),
        $_[0]->findvalue(qq{arch_name[\@version="$_[1]"]}) );
}

sub get_name_by_version {
    $_[0]->findvalue(qq{name[\@version="$_[1]"]});
}

sub obsolete_merge_index {
    my ($merge_xml_doc) = @_;

    my %index;
    foreach my $obs_node ( find_obsolete_merge_nodes($merge_xml_doc) ) {

        my $type = get_type($obs_node);

        if ( $type eq 'x_group' ) {
            push( @{ $index{x} }, $obs_node->findvalue('@x_id') );
        }
        elsif ( $type eq 'h_group' ) {
            push( @{ $index{h} }, $obs_node->findvalue('@h_id') );
        }
        elsif ( $type eq 'f_group' ) {
            push( @{ $index{f} }, $obs_node->findvalue('@f_id') );
        }
        elsif ( $type eq 'domain' ) {
            push( @{ $index{d} }, $obs_node->findvalue('@uid') );
        }
        elsif ( $type eq 'domain_assembly' ) {
            push( @{ $index{da} }, $obs_node->findvalue('@uid') );
        }
        elsif ( $type eq 'arch' ) {
            push( @{ $index{a} }, $obs_node->findvalue('@arch_id') );
        }
    }

    return \%index;
}

sub get_type {
    $_[0]->findvalue('@type');
}

sub find_obsolete_merge_nodes {
    $_[0]->findnodes('//obsolete_elements/obsolete');
}

sub get_max_uid {
    $_[0]->findvalue('//@maxUID');
}

sub get_domain_dictionary_node {
    $_[0]->findnodes('//domain_dictionary')->get_node(1);
}

sub record_new_version {
    my ( $old_xml_doc, $new_xml_doc ) = @_;

    foreach my $old_version_node ( $old_xml_doc->findnodes('//createdOn') ) {
        $old_version_node->unbindNode;
    }
    my $new_version_node = $new_xml_doc->findnodes('//createdOn')->get_node(1)->cloneNode(1);
    my $old_doc_node     = $old_xml_doc->findnodes('//ecod_document')->get_node(1);
    $old_doc_node->appendChild($new_version_node);
}

sub ecod_index {
    my $sub = 'ecod_index';
    my ($ecod_xml_doc) = @_;

    my ( %arch_lookup, %arch_index, %x_index, %h_index, %f_index, %domain_index );
    my $domain_dictionary_node = get_domain_dictionary_node($ecod_xml_doc);
    foreach my $arch_node ( find_architecture_nodes($ecod_xml_doc) ) {
        my $arch_id   = get_arch_id($arch_node);
        my $arch_name = get_arch_name($arch_node);
        print "arch: $arch_id $arch_name\n";
        $arch_lookup{$arch_name} = $arch_id;
        $arch_index{$arch_id}    = $arch_node;

        foreach my $x_node ( find_x_group_nodes($arch_node) ) {
            my $x_id = get_x_id($x_node);
            $x_index{$x_id} = $x_node;
            print "x: $x_id\n";

            foreach my $h_node ( find_h_group_nodes($x_node) ) {
                my $h_id = get_h_id($h_node);
                print "h: $h_id\n";
                $h_index{$h_id} = $h_node;

                foreach my $f_node ( find_f_group_nodes($h_node) ) {

                    my $f_id = get_f_id($f_node);
                    print "f: $f_id\n";
                    $f_index{$f_id} = $f_node;

                    foreach my $domain_node ( find_domain_nodes($f_node) ) {

                        my $uid = get_uid($domain_node);

                        #print "d: $uid\n";
                        $domain_index{$uid} = $domain_node;

                    }
                }
            }
        }
    }
    print "$sub: finish up\n";
    my $version = $ecod_xml_doc->findvalue('//createdOn/@version');

    my %ecod_index;
    $ecod_index{arch_lookup} = \%arch_lookup;
    $ecod_index{arch}        = \%arch_index;
    $ecod_index{x}           = \%x_index;
    $ecod_index{h}           = \%h_index;
    $ecod_index{f}           = \%f_index;
    $ecod_index{d}           = \%domain_index;
    $ecod_index{version}     = $version;
    $ecod_index{dd}          = $domain_dictionary_node;
    return \%ecod_index;

}

sub apply_merge {
    my $sub = 'apply_merge';

    my ( $old_ecod_xml_fn, $new_ecod_xml_fn, $merge_xml_fn ) = @_;

    my $old_ecod_xml_doc = xml_open($old_ecod_xml_fn);
    my $new_ecod_xml_doc = xml_open($new_ecod_xml_fn);
    my $merge_xml_doc    = xml_open($merge_xml_fn);

    my %arch_index;
    foreach my $arch_node ( $old_ecod_xml_doc->findnodes('//architecture') ) {
        my $arch_id   = $arch_node->findvalue('@arch_id');
        my $arch_name = $arch_node->findvalue('@arch_name');
        $arch_index{$arch_name} = $arch_id;
    }

    #Index domain nodes
    my %old_domain_nodes;
    foreach my $domain_node ( $old_ecod_xml_doc->findnodes('//domain') ) {
        my $uid = $domain_node->findvalue('@uid');
        $old_domain_nodes{$uid} = $domain_node;
    }

    my %new_domain_nodes;
    foreach my $domain_node ( $new_ecod_xml_doc->findnodes('//domain') ) {
        my $uid = $domain_node->findvalue('@uid');
        $new_domain_nodes{$uid} = $domain_node;
    }

    my $domain_assembly_XPath = '//domain_assembly';
    my %old_domain_assembly_nodes;
    foreach my $domain_assembly_node ( $old_ecod_xml_doc->findnodes($domain_assembly_XPath) ) {
        my $uid = $domain_assembly_node->findvalue('@uid');
        $old_domain_assembly_nodes{$uid} = $domain_assembly_node;

    }
    my %new_domain_assembly_nodes;
    foreach my $domain_assembly_node ( $new_ecod_xml_doc->findnodes($domain_assembly_XPath) ) {
        my $uid = $domain_assembly_node->findvalue('@uid');
        $new_domain_assembly_nodes{$uid} = $domain_assembly_node;
    }

    #Index f_nodes;
    my $f_group_XPath = '//f_group';
    my %old_f_nodes;
    foreach my $f_node ( $old_ecod_xml_doc->findnodes($f_group_XPath) ) {
        my $f_id = $f_node->findvalue('@f_id');
        $old_f_nodes{$f_id} = $f_node;

    }
    my %new_f_nodes;
    foreach my $f_node ( $new_ecod_xml_doc->findnodes($f_group_XPath) ) {
        my $f_id = $f_node->findvalue('@f_id');
        $new_f_nodes{$f_id} = $f_node;
    }

    #Index h_nodes
    my $h_group_XPath = '//h_group';
    my %old_h_nodes;
    foreach my $h_node ( $old_ecod_xml_doc->findnodes($h_group_XPath) ) {
        my $h_id = $h_node->findvalue('@h_id');
        $old_h_nodes{$h_id} = $h_node;
    }
    my %new_h_nodes;
    foreach my $h_node ( $new_ecod_xml_doc->findnodes($h_group_XPath) ) {
        my $h_id = $h_node->findvalue('@h_id');
        $new_h_nodes{$h_id} = $h_node;
    }

    #Index x_nodes
    my $x_group_XPath = '//x_group';
    my %old_x_nodes;
    foreach my $x_node ( $old_ecod_xml_doc->findnodes($x_group_XPath) ) {
        my $x_id = $x_node->findvalue('@x_id');
        $old_x_nodes{$x_id} = $x_node;
    }
    my %new_x_nodes;
    foreach my $x_node ( $new_ecod_xml_doc->findnodes($x_group_XPath) ) {
        my $x_id = $x_node->findvalue('@x_id');
        $new_x_nodes{$x_id} = $x_node;
    }

    #Index arch nodes;
    my %old_arch_nodes;
    my $arch_XPath = '//architecture';
    foreach my $arch_node ( $old_ecod_xml_doc->findnodes($arch_XPath) ) {
        my $arch_id = $arch_node->findvalue('@arch_id');
        $old_arch_nodes{$arch_id} = $arch_node;
    }
    my %new_arch_nodes;
    foreach my $arch_node ( $new_ecod_xml_doc->findnodes($arch_XPath) ) {
        my $arch_id = $arch_node->findvalue('@arch_id');
        $new_arch_nodes{$arch_id} = $arch_node;
    }

    my $old_version = $old_ecod_xml_doc->findvalue('//createdOn/@version');
    my $new_version = $new_ecod_xml_doc->findvalue('//createdOn/@version');

    print "DEBUG old_v $old_version new_v $new_version\n";

    foreach my $old_version_node ( $old_ecod_xml_doc->findnodes('//createdOn')->get_nodelist() ) {
        $old_version_node->unbindNode;
    }

    my $new_version_node       = $new_ecod_xml_doc->findnodes('//createdOn')->get_node(1);
    my $new_clone_version_node = $new_version_node->cloneNode(1);

    my $old_doc_node = $old_ecod_xml_doc->findnodes('//ecod_document')->get_node(1);
    $old_doc_node->appendChild($new_clone_version_node);

    my $domain_dictionary_node = $old_ecod_xml_doc->findnodes('//domain_dictionary')->get_node(1);
    my $maxUID                 = $domain_dictionary_node->findvalue('@maxUID');

#Obseletion happens first, obsolete domains, f_groups, h_groups, then x_groups, then look for orphan nodes, where possible move those domains to new f_groups using modify nodes

    my $obsolete_nodes = $merge_xml_doc->findnodes('//obsolete_elements/obsolete');

    #What is this doing?
    my ( @obs_arch, @obs_x, @obs_h, @obs_f, @obs_d, @obs_da );
    foreach my $obs_node ( $obsolete_nodes->get_nodelist() ) {

        my $type = $obs_node->findvalue('@type');

        if ( $type eq 'x_group' ) {
            my $x_id = $obs_node->findvalue('@x_id');
            push( @obs_x, $x_id );
        }
        elsif ( $type eq 'h_group' ) {
            my $h_id = $obs_node->findvalue('@h_id');
            push( @obs_h, $h_id );
        }
        elsif ( $type eq 'f_group' ) {
            my $f_id = $obs_node->findvalue('@f_id');
            push( @obs_f, $f_id );
        }
        elsif ( $type eq 'domain' ) {
            my $uid = $obs_node->findvalue('@uid');
            push( @obs_d, $uid );
        }
        elsif ( $type eq 'domain_assembly' ) {
            my $uid = $obs_node->findvalue('@uid');
            push( @obs_da, $uid );
        }
        elsif ( $type eq 'arch' ) {
            my $arch_id = $obs_node->findvalue('@arch_id');
            push( @obs_arch, $arch_id );
        }
        else {
            print "WARNING! Obsolete type $type not recognized\n";
        }
    }
    if ($DEBUG) {
        printf "DEBUG: Found %i arch, %i x, %i h, %i f, %i d, obsolete nodes\n", scalar(@obs_arch), scalar(@obs_x),
          scalar(@obs_h), scalar(@obs_f), scalar(@obs_d);
    }

    foreach my $obs_arch_id (@obs_arch) {
        $old_arch_nodes{$obs_arch_id}->setAttribute( 'obsolete', 'true' );
    }

    foreach my $obs_x_id (@obs_x) {
        $old_x_nodes{$obs_x_id}->setAttribute( 'obsolete', 'true' );
    }

    foreach my $obs_h_id (@obs_h) {
        $old_h_nodes{$obs_h_id}->setAttribute( 'obsolete', 'true' );
    }

    foreach my $obs_f_id (@obs_f) {
        $old_f_nodes{$obs_f_id}->setAttribute( 'obsolete', 'true' );
    }

    foreach my $obs_d_uid (@obs_d) {
        $old_domain_nodes{$obs_d_uid}->setAttribute( 'obsolete', 'true' );
    }

    foreach my $obs_da_uid (@obs_da) {
        $old_domain_assembly_nodes{$obs_da_uid}->setAttribute( 'obsolete', 'true' );
    }

    #Apply modifications.

    #F-groups
    my %mod_x;
    my %mod_h;
    my %mod_f;
    my %mod_d;
    my %mod_da;

    my $modification_nodes = $merge_xml_doc->findnodes('//modify_elements/modify');

    my %new_assembly_types;
    foreach my $modify_node ( $modification_nodes->get_nodelist() ) {

        my $modify_type     = $modify_node->findvalue('@type');
        my $modify_mod_type = $modify_node->findvalue('@mod_type');

        if ( $modify_type eq 'x_group' ) {
            if ( $modify_mod_type eq 'name_change' ) {
                my $x_id = $modify_node->findvalue('@x_id');

                my $old_x_name = $modify_node->findvalue(qq{name[\@version="$old_version"]});
                my $new_x_name = $modify_node->findvalue(qq{name[\@version="$new_version"]});

                $mod_x{$x_id}{name_change} = 1;

                if ($DEBUG) {
                    print "DEBUG: $modify_type $x_id $modify_mod_type $new_x_name $old_version $new_version\n";
                }

                if ( $new_x_name eq 'NA' ) {

                    #Remove name case;
                    $mod_x{$x_id}{remove_name}++;
                }
                else {
                    $mod_x{$x_id}{new_name} = $new_x_name;
                }
            }
            elsif ( $modify_mod_type eq 'arch_shift' ) {
                my $x_id = $modify_node->findvalue('@x_id');

                my $old_arch_id = $modify_node->findvalue(qq{arch_name[\@version="$old_version"]/\@arch_id});
                my $new_arch_id = $modify_node->findvalue(qq{arch_name[\@version="$new_version"]/\@arch_id});

                $mod_x{$x_id}{x_shift} = 1;

                my $old_arch_name = $modify_node->findvalue(qq{arch_name[\@version="$old_version"]});
                my $new_arch_name = $modify_node->findvalue(qq{arch_name[\@version="$new_version"]});

                if ($DEBUG) {
                    print
"DEBUG: $modify_type $x_id $modify_mod_type $old_arch_id:$old_arch_name $new_arch_id:$new_arch_name $old_version $new_version\n";
                }

                $mod_x{$x_id}{old_arch_id}   = $old_arch_id;
                $mod_x{$x_id}{old_arch_name} = $old_arch_name;
                $mod_x{$x_id}{new_arch_id}   = $new_arch_id;
                $mod_x{$x_id}{new_arch_name} = $new_arch_name;

            }
            elsif ( $modify_mod_type eq 'ordinal_shift' ) {
                my $x_id = $modify_node->findvalue('@x_id');

                my $old_x_ordinal = $modify_node->findvalue('@old_ordinal');
                my $new_x_ordinal = $modify_node->findvalue('@new_ordinal');

                $mod_x{$x_id}{old_ordinal} = $old_x_ordinal;
                $mod_x{$x_id}{new_ordinal} = $new_x_ordinal;

                $mod_x{$x_id}{ordinal_shift}++;

            }
            else {
                print "WARNING! mod_type $modify_mod_type not supported for $modify_type, skipping...\n";
            }
        }
        elsif ( $modify_type eq 'h_group' ) {
            if ( $modify_mod_type eq 'name_change' ) {
                my $h_id = $modify_node->findvalue('@h_id');

                $mod_h{$h_id}{name_change} = 1;

                my $old_h_name = $modify_node->findvalue(qq{name[\@version="$old_version"]});
                my $new_h_name = $modify_node->findvalue(qq{name[\@version="$new_version"]});

                if ( $new_h_name eq 'NA' ) {
                    $mod_h{$h_id}{remove_name}++;
                }
                else {
                    $mod_h{$h_id}{new_name} = $new_h_name;
                }

            }
            elsif ( $modify_mod_type eq 'h_shift' ) {

                my $old_h_id = $modify_node->findvalue('@old_h_id');
                my $new_h_id = $modify_node->findvalue('@new_h_id');

                my $old_h_name = $modify_node->findvalue(qq{name[\@version="$old_version"]});
                my $new_h_name = $modify_node->findvalue(qq{name[\@version="$new_version"]});

                $mod_h{$old_h_id}{h_shift} = 1;

                if ( $new_h_name eq 'NA' ) {
                    $mod_h{$old_h_id}{remove_name}++;
                }
                elsif ( $old_h_name ne $new_h_name ) {
                    $mod_h{$old_h_id}{new_name} = $new_h_name;
                }

                $mod_h{$old_h_id}{new_h_id} = $new_h_id;
            }
            elsif ( $modify_mod_type eq 'ordinal_shift' ) {
                my $h_id = $modify_node->findvalue('@h_id');

                my $old_h_ordinal = $modify_node->findvalue('@old_ordinal');
                my $new_h_ordinal = $modify_node->findvalue('@new_ordinal');

                $mod_h{$h_id}{old_ordinal} = $old_h_ordinal;
                $mod_h{$h_id}{new_ordinal} = $new_h_ordinal;

                $mod_h{$h_id}{ordinal_shift}++;

            }
            else {
                print "WARNING! mod type $modify_mod_type not supported for $modify_type, skipping..\n";
            }
        }
        elsif ( $modify_type eq 'f_group' ) {
            if ( $modify_mod_type eq 'name_change' ) {
                my $f_id = $modify_node->findvalue('@f_id');

                my $old_f_name = $modify_node->findvalue(qq{name[\@version="$old_version"]});
                my $new_f_name = $modify_node->findvalue(qq{name[\@version="$new_version"]});

                $mod_f{$f_id}{name_change} = 1;

                if ( $new_f_name eq 'NA' ) {
                    $mod_f{$f_id}{remove_name}++;
                }
                else {
                    $mod_f{$f_id}{new_name} = $new_f_name;
                }
            }
            elsif ( $modify_mod_type eq 'f_shift' ) {

                my $old_f_id = $modify_node->findvalue('@old_f_id');
                my $new_f_id = $modify_node->findvalue('@new_f_id');

                my $old_f_name = $modify_node->findvalue(qq{name[\@version="$old_version"]});
                my $new_f_name = $modify_node->findvalue(qq{name[\@version="$new_version"]});

                $mod_f{$old_f_id}{f_shift} = 1;

                if ( $new_f_name eq 'NA' ) {
                    $mod_f{$old_f_id}{remove_name}++;
                }
                elsif ( $new_f_name ne $old_f_name ) {
                    $mod_f{$old_f_id}{new_name} = $new_f_name;
                }

                $mod_f{$old_f_id}{new_f_id} = $new_f_id;
            }
            elsif ( $modify_mod_type eq 'ordinal_shift' ) {
                my $f_id = $modify_node->findvalue('@f_id');

                my $old_f_ordinal = $modify_node->findvalue('@old_ordinal');
                my $new_f_ordinal = $modify_node->findvalue('@new_ordinal');

                $mod_f{$f_id}{old_ordinal} = $old_f_ordinal;
                $mod_f{$f_id}{new_ordinal} = $new_f_ordinal;

                $mod_f{$f_id}{ordinal_shift}++;

            }
            else {
                print "WARNING! mod type $modify_mod_type not supported for $modify_type, skipping...\n";
            }
        }
        elsif ( $modify_type eq 'domain' ) {
            if ( $modify_mod_type eq 'derived_shift' ) {

                my $old_domain_id          = $modify_node->findvalue('@domain_id');
                my $old_domain_uid         = $modify_node->findvalue('@old_uid');
                my $old_domain_parent_uid  = $modify_node->findvalue('@old_parent_uid');
                my $old_domain_parent_f_id = $modify_node->findvalue('@old_parent_f_id');

                $mod_d{$old_domain_uid}{derived_shift}++;
                $mod_d{$old_domain_uid}{domain_id}   = $old_domain_id;
                $mod_d{$old_domain_uid}{parent_uid}  = $old_domain_parent_uid;
                $mod_d{$old_domain_uid}{parent_f_id} = $old_domain_parent_f_id;
            }
            elsif ( $modify_mod_type eq 'inherit_annotation' ) {

                my $old_domain_id  = $modify_node->findvalue('@domain_id');
                my $old_domain_uid = $modify_node->findvalue('@old_uid');
                my $new_domain_uid = $modify_node->findvalue('@new_uid');

                $mod_d{$old_domain_uid}{inherit_annotation}++;
                $mod_d{$old_domain_uid}{domain_id} = $old_domain_id;
                $mod_d{$old_domain_uid}{new_uid}   = $new_domain_uid;
            }
            elsif ( $modify_mod_type eq 'f_shift' ) {

                my $old_domain_id  = $modify_node->findvalue('@domain_id');
                my $old_domain_uid = $modify_node->findvalue('@old_uid');
                my $new_domain_uid = $modify_node->findvalue('@new_uid');
                my $new_domain_fid = $modify_node->findvalue('@new_fid');

                my $old_domain_type = $modify_node->findvalue('@old_domain_type');
                my $new_domain_type = $modify_node->findvalue('@new_domain_type');

                $mod_d{$old_domain_uid}{domain_fshift}++;
                $mod_d{$old_domain_uid}{domain_id} = $old_domain_id;
                $mod_d{$old_domain_uid}{new_uid}   = $new_domain_uid;
                $mod_d{$old_domain_uid}{new_fid}   = $new_domain_fid;

                $mod_d{$old_domain_uid}{old_domain_type} = $old_domain_type;
                $mod_d{$old_domain_uid}{new_domain_type} = $new_domain_type;

            }
            elsif ( $modify_mod_type eq 'side_load_shift' ) {

                my $old_domain_id          = $modify_node->findvalue('@domain_id');
                my $old_domain_uid         = $modify_node->findvalue('@old_uid');
                my $old_domain_parent_f_id = $modify_node->findvalue('@old_parent_f_id');

                $mod_d{$old_domain_uid}{side_load}++;
                $mod_d{$old_domain_uid}{domain_id}   = $old_domain_id;
                $mod_d{$old_domain_uid}{parent_f_id} = $old_domain_parent_f_id;

            }
            elsif ( $modify_mod_type eq 'nonrep_to_rep_conversion' ) {
                my $old_domain_id  = $modify_node->findvalue('@domain_id');
                my $old_domain_uid = $modify_node->findvalue('@old_uid');
                my $new_domain_uid = $modify_node->findvalue('@new_uid');

                my $old_domain_type = $modify_node->findvalue('@old_domain_type');
                my $new_domain_type = $modify_node->findvalue('@new_domain_type');

                $mod_d{$old_domain_uid}{nonrep_conversion}++;
                $mod_d{$old_domain_uid}{domain_id} = $old_domain_id;
                $mod_d{$old_domain_uid}{new_uid}   = $new_domain_uid;

                $mod_d{$old_domain_uid}{old_domain_type} = $old_domain_type;
                $mod_d{$old_domain_uid}{new_domain_type} = $new_domain_type;

            }
            elsif ( $modify_mod_type eq 'convert_to_domain_assembly' ) {
                my $old_domain_id  = $modify_node->findvalue('@domain_id');
                my $old_domain_uid = $modify_node->findvalue('@old_uid');
                my $new_domain_uid = $modify_node->findvalue('@new_uid');

                my $old_domain_parent_f_id = $modify_node->findvalue('@old_parent_f_id');
                my $new_da_uid             = $modify_node->findvalue('@da');
                my $assembly_type          = $modify_node->findvalue('@assembly_type');
                $new_assembly_types{$new_da_uid} = $assembly_type;

                $mod_d{$old_domain_uid}{domain_assembly_conversion}++;
                $mod_d{$old_domain_uid}{domain_id} = $old_domain_id;
                $mod_d{$old_domain_uid}{new_uid}   = $new_domain_uid;

                $mod_d{$old_domain_uid}{da_uid}      = $new_da_uid;
                $mod_d{$old_domain_uid}{parent_f_id} = $old_domain_parent_f_id;

            }
            else {
                print "WARNING! mod type $modify_mod_type not supported for $modify_type, skipping...\n";
            }
        }
        elsif ( $modify_type eq 'domain_assembly' ) {
            if ( $modify_mod_type eq 'derived_shift' ) {
                my $old_domain_id          = $modify_node->findvalue('@domain_id');
                my $old_domain_uid         = $modify_node->findvalue('@old_uid');
                my $old_domain_parent_uid  = $modify_node->findvalue('@old_parent_uid');
                my $old_domain_parent_f_id = $modify_node->findvalue('@old_parent_f_id');

                $mod_da{$old_domain_uid}{derived_shift}++;
                $mod_da{$old_domain_uid}{domain_id}   = $old_domain_id;
                $mod_da{$old_domain_uid}{parent_uid}  = $old_domain_parent_uid;
                $mod_da{$old_domain_uid}{parent_f_id} = $old_domain_parent_f_id;
            }
            elsif ( $modify_mod_type eq 'inherit_annotation' ) {

                my $old_domain_id  = $modify_node->findvalue('@domain_id');
                my $old_domain_uid = $modify_node->findvalue('@old_uid');
                my $new_domain_uid = $modify_node->findvalue('@new_uid');

                $mod_da{$old_domain_uid}{inherit_annotation}++;
                $mod_da{$old_domain_uid}{domain_id} = $old_domain_id;
                $mod_da{$old_domain_uid}{new_uid}   = $new_domain_uid;

            }
            elsif ( $modify_mod_type eq 'f_shift' ) {
                my $old_domain_id  = $modify_node->findvalue('@domain_id');
                my $old_domain_uid = $modify_node->findvalue('@old_uid');
                my $new_domain_uid = $modify_node->findvalue('@new_uid');
                my $new_domain_fid = $modify_node->findvalue('@new_fid');

                $mod_da{$old_domain_uid}{domain_fshift}++;
                $mod_da{$old_domain_uid}{domain_id} = $old_domain_id;
                $mod_da{$old_domain_uid}{new_uid}   = $new_domain_uid;
                $mod_da{$old_domain_uid}{new_fid}   = $new_domain_fid;

            }
            else {
                print "WARNING! mod type $modify_mod_type not supported for $modify_type, skipping...\n";
            }
        }
    }
    if ($DEBUG) {
        printf "DEBUG: Found %i x, %i h, %i f, %i d modify nodes\n", scalar( keys %mod_x ), scalar( keys %mod_h ),
          scalar keys %mod_f, scalar keys %mod_d;
    }

    #Read new nodes;
    my @new_x;
    my %new_xname;
    my @new_arch;
    my %new_xarch;
    my %new_xarchname;
    my %new_archname;
    my @new_h;
    my %new_hname;
    my @new_f;
    my %new_fname;
    my @new_d;
    my %new_domain_id;

    my %new_arch_ordinal;
    my %new_x_ordinal;
    my %new_h_ordinal;
    my %new_f_ordinal;

    my $new_nodes = $merge_xml_doc->findnodes('//new_elements/new');

    foreach my $new_node ( $new_nodes->get_nodelist() ) {

        my $new_type = $new_node->findvalue('@type');

        if ( $new_type eq 'x_group' ) {
            my $x_id   = $new_node->findvalue('@x_id');
            my $x_name = $new_node->findvalue('@name');

            #my $x_arch	= $new_node->findvalue('@architecture');
            #my $new_ecod_node	= $new_ecod_xml_doc->findnodes(qq{//x_group[\@x_id="$x_id"]})->get_node(1);
            my $new_ecod_node = $new_x_nodes{$x_id};
            my $x_arch_id     = $new_ecod_node->parentNode->findvalue('@arch_id');
            my $x_arch_name   = $new_ecod_node->parentNode->findvalue('@arch_name');

            my $x_ordinal = $new_node->findvalue('@x_ordinal');

            push( @new_x, $x_id );
            $new_xname{$x_id}     = $x_name;
            $new_xarch{$x_id}     = $x_arch_id;
            $new_xarchname{$x_id} = $x_arch_name;
            $new_x_ordinal{$x_id} = $x_ordinal;
        }
        elsif ( $new_type eq 'arch' ) {

            my $arch_id      = $new_node->findvalue('@arch_id');
            my $arch_name    = $new_node->findvalue('@name');
            my $arch_ordinal = $new_node->findvalue('@arch_ordinal');

            push( @new_arch, $arch_id );
            $new_archname{$arch_id}     = $arch_name;
            $new_arch_ordinal{$arch_id} = $arch_ordinal;

        }
        elsif ( $new_type eq 'h_group' ) {
            my $h_id      = $new_node->findvalue('@h_id');
            my $h_name    = $new_node->findvalue('@name');
            my $h_ordinal = $new_node->findvalue('@h_ordinal');

            push( @new_h, $h_id );
            $new_hname{$h_id}     = $h_name;
            $new_h_ordinal{$h_id} = $h_ordinal;

        }
        elsif ( $new_type eq 'f_group' ) {
            my $f_id      = $new_node->findvalue('@f_id');
            my $f_name    = $new_node->findvalue('@name');
            my $f_ordinal = $new_node->findvalue('@f_ordinal');

            push( @new_f, $f_id );
            $new_fname{$f_id}     = $f_name;
            $new_f_ordinal{$f_id} = $f_ordinal;

        }
        elsif ( $new_type eq 'domain' ) {
            my $domain_id = $new_node->findvalue('@domain_id');
            my $domain_uid =
              $new_node->findvalue('@uid');    #this UID needs to be dropped/converted to be non-colliding with old UIDs

            push( @new_d, $domain_uid );
            $new_domain_id{$domain_uid} =
              $domain_id;    #this is new_uid, only good for finding checking the new and old f_group id;
        }
        else {
            print "WARNING! No rules for new type $new_type, skipping...\n";
        }
    }

    if ($DEBUG) {
        printf "DEBUG: Found %i x, %i h, %i f, %i d new nodes\n", scalar(@new_x), scalar(@new_h), scalar(@new_f),
          scalar(@new_d);
    }

    my $old_top_node = $old_ecod_xml_doc->findnodes('//domain_dictionary')->get_node(1);

    #Architecture
    foreach my $new_arch_id (@new_arch) {
        my $new_arch_name = $new_archname{$new_arch_id};

        if ( $old_ecod_xml_doc->exists(qq{//architecture[\@arch_name="$new_arch_name"]}) ) {
            die "ERROR! old->new arch name conflict on $new_arch_name... improper modify rules?\n";
        }

        my $new_arch_node = $old_ecod_xml_doc->createElement('architecture');

        $new_arch_node->setAttribute( 'arch_id',   $new_arch_id );
        $new_arch_node->setAttribute( 'arch_name', $new_arch_name );

        $old_top_node->appendChild($new_arch_node);

    }

    #X-groups;
    foreach my $new_x_id (@new_x) {

        if ( $old_ecod_xml_doc->exists(qq{//x_group[\@x_id="$new_x_id"]}) ) {
            die "ERROR! old->new x_group id conflict on $new_x_id.... improper modify rules?\n";
        }

        my $new_x_group_node = $old_ecod_xml_doc->createElement('x_group');

        $new_x_group_node->setAttribute( 'x_id',      $new_x_id );
        $new_x_group_node->setAttribute( 'x_ordinal', $new_x_ordinal{$new_x_id} );
        if ( $new_xname{$new_x_id} ne 'NA' ) {
            if ($DEBUG) {
                print "DEBUG: x new name $new_x_id $new_xname{$new_x_id}\n";
            }
            $new_x_group_node->setAttribute( 'name', $new_xname{$new_x_id} );
        }

        #my $x_arch_id = $new_xarch{$new_x_id}; # NO
        my $x_arch_name = $new_xarchname{$new_x_id};
        my $x_arch_id   = $arch_index{$x_arch_name};

        if ($DEBUG) {
            print "DEBUG: x_arch $x_arch_name $x_arch_id\n";
        }

        #$new_x_group_node->setAttribute('architecture', $new_arch{$new_x_id});
        #$new_x_group_node->setAttribute('tmp_new', 'true');
        #my $arch_node = $old_ecod_xml_doc->findnodes(qq{//architecture[\@arch_id="$x_arch_id"]})->get_node(1);
        my $arch_node = $old_arch_nodes{$x_arch_id};

        #$old_top_node->appendChild($new_x_group_node);
        $arch_node->appendChild($new_x_group_node);
        $old_x_nodes{$new_x_id} = $new_x_group_node;
    }

    foreach my $mod_x_id ( sort { $a cmp $b } keys %mod_x ) {

        #my $old_x_group_node	= $old_ecod_xml_doc->findnodes(qq{//x_group[\@x_id="$mod_x_id"]})->get_node(1);
        my $old_x_group_node = $old_x_nodes{$mod_x_id};

        if ( $mod_x{$mod_x_id}{remove_name} ) {
            $old_x_group_node->removeAttribute('name');
            $old_x_group_node->setAttribute( 'tmp_name_remove', 'true' );

            if ($DEBUG) {
                print "DEBUG: mod_x remove name $mod_x_id\n";
            }
        }
        if ( $mod_x{$mod_x_id}{new_name} ) {

            $old_x_group_node->setAttribute( 'name',            $mod_x{$mod_x_id}{new_name} );
            $old_x_group_node->setAttribute( 'tmp_name_change', 'true' );

            if ($DEBUG) {
                print "DEBUG: mod_x new_name $mod_x_id $mod_x{$mod_x_id}{new_name}\n";
            }
        }
        if ( $mod_x{$mod_x_id}{x_shift} ) {

            my $new_arch_id   = $mod_x{$mod_x_id}{new_arch_id};
            my $old_arch_name = $mod_x{$mod_x_id}{old_arch_name};
            my $old_arch_id   = $mod_x{$mod_x_id}{old_arch_id};
            my $new_arch_name = $mod_x{$mod_x_id}{new_arch_name};

           #my $new_arch_node = $old_ecod_xml_doc->findnodes(qq{//architecture[\@arch_id="$new_arch_id"]})->get_node(1);
            my $trans_arch_id = $arch_index{$new_arch_name};
            my $new_arch_node =
              $old_ecod_xml_doc->findnodes(qq{//architecture[\@arch_id="$trans_arch_id"]})->get_node(1);

            if ($DEBUG) {
                print
"DEBUG: mod_x $mod_x_id old arch $old_arch_id new_arch $new_arch_id trans_arch $trans_arch_id $old_arch_name=>$new_arch_name\n";
            }

            $new_arch_node->appendChild($old_x_group_node);
        }
        if ( $mod_x{$mod_x_id}{ordinal_shift} ) {
            my $old_x_ordinal = $mod_x{$mod_x_id}{old_ordinal};
            my $new_x_ordinal = $mod_x{$mod_x_id}{new_ordinal};

            $old_x_group_node->setAttribute( 'x_ordinal', $new_x_ordinal );
        }

    }

    foreach my $obs_x_id (@obs_x) {

        if ( !$old_ecod_xml_doc->exists(qq{//x_group[\@x_id="$obs_x_id"]}) ) {
            die "ERROR! Could not find obs x_group for $obs_x_id\n";
        }

        #my $obs_x_group_node	= $old_ecod_xml_doc->findnodes(qq{//x_group[\@x_id="$obs_x_id"]})->get_node(1);
        my $obs_x_group_node = $old_x_nodes{$obs_x_id};

        if ($DEBUG) {
            print "DEBUG obs x: $obs_x_id\n";
        }

        $obs_x_group_node->setAttribute( 'isObsolete', 'true' );
    }

    #H-groups
    foreach my $new_h_id (@new_h) {

        if ( $old_ecod_xml_doc->exists(qq{//h_group[\@h_id="$new_h_id"]}) ) {
            die "ERROR! old->new x_group id conflict on $new_h_id.... improper modify rules?\n";
        }

        my $new_h_group_node = $old_ecod_xml_doc->createElement('h_group');

        $new_h_group_node->setAttribute( 'h_id',      $new_h_id );
        $new_h_group_node->setAttribute( 'h_ordinal', $new_h_ordinal{$new_h_id} );
        if ( $new_hname{$new_h_id} ne 'NA' ) {
            $new_h_group_node->setAttribute( 'name', $new_hname{$new_h_id} );
        }

        #$new_h_group_node->setAttribute('tmp_new', 'true');

        my $parent_x_id;
        if ( $new_h_id =~ /(\d+)\.(\d+)/ ) {
            $parent_x_id = $1;
        }
        else {
            die "ERROR! Weird H-group id $new_h_id\n";
        }

        my $parent_x_group_node;
        if ( $old_ecod_xml_doc->exists(qq{//x_group[\@x_id="$parent_x_id"]}) ) {

            #$parent_x_group_node = $old_ecod_xml_doc->findnodes(qq{//x_group[\@x_id="$parent_x_id"]})->get_node(1);
            $parent_x_group_node = $old_x_nodes{$parent_x_id};
        }
        else {
            die "ERROR! Could not find parent X-group node $parent_x_id for H-group $new_h_id\n";
        }

        if ($DEBUG) {
            print "DEBUG new_h: $new_h_id $new_hname{$new_h_id} px:$parent_x_id\n";
        }

        $parent_x_group_node->appendChild($new_h_group_node);
        $old_h_nodes{$new_h_id} = $new_h_group_node;
    }

    foreach my $mod_h_id ( sort { $a cmp $b } keys %mod_h ) {

        #my $old_h_group_node	= $old_ecod_xml_doc->findnodes(qq{//h_group[\@h_id="$mod_h_id"]})->get_node(1);
        my $old_h_group_node = $old_h_nodes{$mod_h_id};

        if ( $mod_h{$mod_h_id}{remove_name} ) {
            $old_h_group_node->removeAttribute('name');
            $old_h_group_node->setAttribute( 'tmp_name_remove', 'true' );

            if ($DEBUG) {
                print "DEBUG: mod_h remove name $mod_h_id\n";
            }
        }
        elsif ( $mod_h{$mod_h_id}{new_name} ) {

            $old_h_group_node->setAttribute( 'name',            $mod_h{$mod_h_id}{new_name} );
            $old_h_group_node->setAttribute( 'tmp_name_change', 'true' );

            if ($DEBUG) {
                print "DEBUG: mod_h new_name $mod_h_id $mod_h{$mod_h_id}{new_name}\n";
            }
        }

        if ( $mod_h{$mod_h_id}{ordinal_shift} ) {
            my $old_h_ordinal = $mod_h{$mod_h_id}{old_ordinal};
            my $new_h_ordinal = $mod_h{$mod_h_id}{new_ordinal};

            $old_h_group_node->setAttribute( 'h_ordinal', $new_h_ordinal );
        }

    }

    foreach my $obs_h_id (@obs_h) {

        if ( !$old_ecod_xml_doc->exists(qq{//h_group[\@h_id="$obs_h_id"]}) ) {
            die "ERROR! Could not find obs h_group for $obs_h_id\n";
        }

        #my $obs_h_group_node	= $old_ecod_xml_doc->findnodes(qq{//h_group[\@h_id="$obs_h_id"]})->get_node(1);
        my $obs_h_group_node = $old_h_nodes{$obs_h_id};

        $obs_h_group_node->setAttribute( 'isObsolete', 'true' );

        if ($DEBUG) {
            print "DEBUG obs_h: $obs_h_id\n";
        }
    }

    #F-groups

    my %old_f_group_nodes;
    foreach my $new_f_id (@new_f) {

        if ( $old_ecod_xml_doc->exists(qq{//f_group[\@f_id="$new_f_id"][\@isNew="true"]}) ) {
            die "ERROR! old->new f_group id conflict on $new_f_id... improper modify rules?\n";
        }

        my $new_f_group_node = $old_ecod_xml_doc->createElement('f_group');

        $new_f_group_node->setAttribute( 'f_id',      $new_f_id );
        $new_f_group_node->setAttribute( 'name',      $new_fname{$new_f_id} );
        $new_f_group_node->setAttribute( 'f_ordinal', $new_f_ordinal{$new_f_id} );
        $new_f_group_node->setAttribute( 'isNew',     'true' );

        my $parent_h_id;
        if ( $new_f_id =~ /(\d+\.\d+)\.\d+/ ) {
            $parent_h_id = $1;
        }
        else {
            die "ERROR! Weird F-group id $new_f_id\n";
        }

        my $parent_h_group_node;
        if ( $old_ecod_xml_doc->exists(qq{//h_group[\@h_id="$parent_h_id"]}) ) {

            #$parent_h_group_node	= $old_ecod_xml_doc->findnodes(qq{//h_group[\@h_id="$parent_h_id"]})->get_node(1);
            $parent_h_group_node = $old_h_nodes{$parent_h_id};
        }
        else {
            die "ERROR! Could not find parent H-group node $parent_h_id for F-group $new_f_id\n";
        }

        $parent_h_group_node->appendChild($new_f_group_node);
        $old_f_group_nodes{$new_f_id} = $new_f_group_node;
        $old_f_nodes{$new_f_id}       = $new_f_group_node;

        if ($DEBUG) {
            print "DEBUG new_f: $new_f_id $new_fname{$new_f_id} pf: $parent_h_id\n";
        }
    }

    foreach my $mod_f_id ( sort { $a cmp $b } keys %mod_f ) {

        my $old_f_group_node = $old_f_nodes{$mod_f_id};

        #	my $old_f_group_nodes = $old_ecod_xml_doc->findnodes(qq{//f_group[\@f_id="$mod_f_id"]});
        #	if ($old_f_group_nodes->size() > 1) {
        #		print "WARNING! Found more than one f_group node with $mod_f_id\n";
        #		foreach my $node ($old_f_group_nodes->get_nodelist() ) {
        #			if ($node->exists('@old_f_id')) {
        #				next;
        #			}elsif($node->exists('@isNew')) {
        #				next;
        #			}else{
        #				$old_f_group_node = $node;
        #			}
        #		}
        #	}else{
        #		$old_f_group_node = $old_f_group_nodes->get_node(1);
        #	}

        if ( $mod_f{$mod_f_id}{remove_name} ) {
            $old_f_group_node->removeAttribute('name');
            $old_f_group_node->setAttribute( 'tmp_name_remove', 'true' );

            if ($DEBUG) {
                print "DEBUG: mod_f remove name $mod_f_id\n";
            }
        }
        elsif ( $mod_f{$mod_f_id}{new_name} ) {

            $old_f_group_node->setAttribute( 'name',            $mod_f{$mod_f_id}{new_name} );
            $old_f_group_node->setAttribute( 'tmp_name_change', 'true' );

            if ($DEBUG) {
                print "DEBUG: mod_f new name $mod_f_id $mod_f{$mod_f_id}{new_name}\n";
            }
        }

        if ( $mod_f{$mod_f_id}{f_shift} ) {

            my $new_f_id = $mod_f{$mod_f_id}{new_f_id};

            if ( $old_ecod_xml_doc->exists(qq{//f_group[\@f_id="$new_f_id"]}) ) {
                print "WARNING! potential F-group ID collision at $new_f_id\n";
                $old_f_group_node->setAttribute( 'old_f_id', $mod_f_id );
            }
            $old_f_group_node->setAttribute( 'f_id', $new_f_id );

            my $parent_h_id;
            if ( $new_f_id =~ /(\d+\.\d+)\.\d+/ ) {
                $parent_h_id = $1;
            }
            else {
                die "ERROR! Weird F-group id $new_f_id\n";
            }

            my $parent_h_group_node;
            if ( $old_ecod_xml_doc->exists(qq{//h_group[\@h_id="$parent_h_id"]}) ) {

                #$parent_h_group_node	= $old_ecod_xml_doc->findnodes(qq{//h_group[\@h_id="$parent_h_id"]})->get_node(1);
                $parent_h_group_node = $old_h_nodes{$parent_h_id};
            }
            else {
                die "ERROR! Could not find parent H-group node $parent_h_id for F-group $mod_f_id\n";
            }

            $parent_h_group_node->appendChild($old_f_group_node);

            if ($DEBUG) {
                print "DEBUG: mod_f f_shift: $mod_f_id->$new_f_id name? px: $parent_h_id\n";
            }
        }

        if ( $mod_f{$mod_f_id}{ordinal_shift} ) {
            my $new_f_ordinal = $mod_f{$mod_f_id}{new_ordinal};
            $old_f_group_node->setAttribute( 'f_ordinal', $new_f_ordinal );
        }
    }

    #Domains.
    my %domain_assembly_uid;
    foreach my $new_domain_uid (@new_d) {

        #new domain uid needs to be replaced by whatever is current in old domain file
        #my $old_domain_uid = sprintf "%09i", $maxUID++;
        my $sql_uid = sprintf "%09i", sql_uid_gen();

        if ($DEBUG) {
            print "DEBUG new_d: $new_domain_uid->$sql_uid";
        }

        my $new_ecod_domain_node;

        if ( !$new_domain_nodes{$new_domain_uid} ) {
            print "WARNING! Could not find new domain node for $new_domain_uid, skipping...\n";
            next;
        }
        else {
            $new_ecod_domain_node = $new_domain_nodes{$new_domain_uid};
        }

        #Domain assembly case
        my $parent_f_id;
        if ( $new_ecod_domain_node->parentNode->nodeName eq 'domain_assembly' ) {

            #$isDomainAssembly = 1;

            $parent_f_id = $new_ecod_domain_node->parentNode->parentNode->findvalue('@f_id');

            my $parent_f_node;
            if ( !$old_ecod_xml_doc->exists(qq{//f_group[\@f_id="$parent_f_id"]}) ) {
                die "ERROR! domain_assembly parent F-group $parent_f_id not found\n";
            }
            else {
                #$parent_f_node	= $old_ecod_xml_doc->findnodes(qq{//f_group[\@f_id="$parent_f_id"]})->get_node(1);
                $parent_f_node = $old_f_nodes{$parent_f_id};
            }

            my $old_domain_assembly_node;
            if ( $domain_assembly_uid{$new_domain_uid} ) {
                my $da_uid = $domain_assembly_uid{$new_domain_uid};
                $old_domain_assembly_node =
                  $old_ecod_xml_doc->findnodes(qq{//domain_assembly[\@uid="$da_uid"]})->get_node(1);

                #$old_domain_assembly_node = $old_domain_assembly_nodes{$da_uid};

                if ($DEBUG) {
                    print "DEBUG new_domain existing da: $new_domain_uid in $da_uid\n";
                }
            }
            else {
                $old_domain_assembly_node = $old_ecod_xml_doc->createElement('domain_assembly');

                #my $old_domain_assembly_uid	= sprintf "%09i", $maxUID++;
                my $sql_uid = sprintf "%09i", sql_uid_gen();
                $old_domain_assembly_node->setAttribute( 'uid', $sql_uid );    #WHAT
                $parent_f_node->appendChild($old_domain_assembly_node);

                foreach my $sib_node ( $new_ecod_domain_node->parentNode->findnodes('domain')->get_nodelist() ) {
                    if ( $sib_node->isSameNode($new_ecod_domain_node) ) { next }
                    else {
                        my $sib_new_uid = $sib_node->findvalue('@uid');
                        $domain_assembly_uid{$sib_new_uid} = $sql_uid;
                    }
                }

                if ($DEBUG) {
                    print "DEBUG new_domain new da: $new_domain_uid $sql_uid\n";
                }
            }

            $old_domain_assembly_node->appendChild($new_ecod_domain_node);

            #$new_ecod_domain_node->setAttribute('uid', $old_domain_uid);
            $new_ecod_domain_node->setAttribute( 'uid', $sql_uid );

        }
        else {
            $parent_f_id = $new_ecod_domain_node->parentNode->findvalue('@f_id');

            my $parent_f_node;
            if ( !$old_ecod_xml_doc->exists(qq{//f_group[\@f_id="$parent_f_id"]}) ) {
                die "ERROR! domain parent F-group $parent_f_id not found\n";
            }
            else {
                #$parent_f_node	= $old_ecod_xml_doc->findnodes(qq{//f_group[\@f_id="$parent_f_id"]})->get_node(1);
                $parent_f_node = $old_f_nodes{$parent_f_id};
            }

            $parent_f_node->appendChild($new_ecod_domain_node);
            $new_ecod_domain_node->setAttribute( 'uid', $sql_uid );

            if ($DEBUG) {
                print "DEBUG new_d: $sql_uid pf:$parent_f_id\n";
            }

        }
    }

    #mod domain
    my %current_da_uids;
    foreach my $mod_domain_uid ( sort { $a cmp $b } keys %mod_d ) {

        if ( $mod_d{$mod_domain_uid}{derived_shift} ) {

            my $parent_uid  = $mod_d{$mod_domain_uid}{parent_uid};
            my $parent_f_id = $mod_d{$mod_domain_uid}{parent_f_id};    #This is just F-id, not parent F-id
            my $domain_id   = $mod_d{$mod_domain_uid}{domain_id};

            my $old_domain_node = $old_domain_nodes{$mod_domain_uid};
            my $parent_domain_node;

            if ( !$old_domain_nodes{$parent_uid} ) {

                #print "WARNING! Parent domain not found for derived domain $domain_id, obsolete by assocation?\n";
                next;
            }
            else {
                $parent_domain_node = $old_domain_nodes{$parent_uid};
            }

            my $parent_f_node;
            if ( $parent_domain_node->parentNode->nodeName eq 'domain_assembly' ) {
                $parent_f_node = $parent_domain_node->parentNode->parentNode;
            }
            elsif ( $parent_domain_node->parentNode->parentNode->nodeName eq 'f_group' ) {
                $parent_f_node =
                  $parent_domain_node->parentNode->parentNode;    #hack for mixed pfgroup fgroup versions, remove ASAP
            }
            else {
                $parent_f_node = $parent_domain_node->parentNode;
            }

            my $current_f_node;
            if ( $old_domain_node->parentNode->nodeName eq 'domain_assembly' ) {
                $current_f_node = $old_domain_node->parentNode->parentNode;
            }
            elsif ( $old_domain_node->parentNode->parentNode->nodeName eq 'f_group' ) {
                $current_f_node =
                  $old_domain_node->parentNode->parentNode;       #hack for mixed pfgroup fgroup versions, remove ASAP
            }
            else {
                $current_f_node = $old_domain_node->parentNode;
            }

            my $new_parent_f_id  = $parent_f_node->findvalue('@f_id');
            my $new_current_f_id = $current_f_node->findvalue('@f_id');

          #print "DEBUG $domain_id\t$new_parent_f_id\t$parent_f_id\t$new_current_f_id\t$mod_domain_uid\t$parent_uid\n";
          #Shouldn't this happen much earlier, why even make an derived_shift mod in the first place if it is not needed
            if (   $parent_f_node->findvalue('@f_id') ne $parent_f_id
                || $current_f_node->findvalue('@f_id') ne $parent_f_id )
            {
                #print "WARNING! F-group shift for derived domain $domain_id\n";
                $parent_f_node->appendChild($old_domain_node);
            }

            if ($DEBUG) {
                print "DEBUG mod_d derived_shift: $mod_domain_uid $domain_id child of $parent_uid $parent_f_id\n";
            }

        }
        elsif ( $mod_d{$mod_domain_uid}{side_load} ) {

            my $parent_f_id = $mod_d{$mod_domain_uid}{parent_f_id};
            my $domain_id   = $mod_d{$mod_domain_uid}{domain_id};

            #my $old_domain_node	= $old_ecod_xml_doc->findnodes(qq{//domain[\@uid="$mod_domain_uid"]})->get_node(1);
            my $old_domain_node = $old_domain_nodes{$mod_domain_uid};

            my $parent_f_node;
            if ( $old_ecod_xml_doc->exists(qq{//f_group[\@f_id="$parent_f_id"]}) ) {

                #$parent_f_node = $old_ecod_xml_doc->findnodes(qq{//f_group[\@f_id="$parent_f_id"]})->get_node(1);
                $parent_f_node = $old_f_nodes{$parent_f_id};
            }

            #Then nothing happens?

        }
        elsif ( $mod_d{$mod_domain_uid}{inherit_annotation} ) {

            my $new_domain_uid = $mod_d{$mod_domain_uid}{new_uid};

            if ( !$new_domain_nodes{$new_domain_uid} ) {
                die "ERROR! No new domain node found for UID $new_domain_uid\n";
            }
            my $new_domain_node = $new_domain_nodes{$new_domain_uid};

            if ( !$old_domain_nodes{$mod_domain_uid} ) {
                die "ERROR! No old domain node found for UID $mod_domain_uid\n";
            }
            my $old_domain_node = $old_domain_nodes{$mod_domain_uid};

            if ($DEBUG) {
                print "DEBUG mod_d inherit_annotation: $mod_domain_uid $new_domain_uid\n";
            }

            #Annotation types that should be transferred

            #Is the domain a manual representative.. i.e. appears in Hua's file and not as a fluke domain

            if ( $new_domain_node->exists('@manual_rep') ) {
                $old_domain_node->setAttribute( 'manual_rep', $new_domain_node->findvalue('@manual_rep') );
            }

            if ( $new_domain_node->findvalue('@manual_rep') eq 'false' ) {
                if ( $new_domain_node->exists('scop_representative_domain')
                    && !$old_domain_node->exists('scop_representative_domain') )
                {
                    my $scop_rep_node = $new_domain_node->findnodes('scop_representative_domain')->get_node(1);
                    $old_domain_node->appendChild($scop_rep_node);
                }
                elsif ( $new_domain_node->exists('ecod_representative_domain')
                    && !$old_domain_node->exists('ecod_representative_domain') )
                {
                    my $ecod_rep_node = $new_domain_node->findnodes('ecod_representative_domain')->get_node(1);
                    $old_domain_node->appendChild($ecod_rep_node);
                }
            }

            #alert comment
            if ( $new_domain_node->exists('alert_comment') ) {
                my $new_alert_comment      = $new_domain_node->findvalue('alert_comment');
                my $new_alert_comment_node = $new_domain_node->findnodes('alert_comment')->get_node(1);
                if ( $old_domain_node->exists('alert_comment') ) {
                    my $old_alert_comment      = $old_domain_node->findvalue('alert_comment');
                    my $old_alert_comment_node = $old_domain_node->findnodes('alert_comment')->get_node(1);
                    if ( $old_alert_comment ne $new_alert_comment ) {
                        print "WARNING! overwriting '$old_alert_comment' with '$new_alert_comment'\n";
                        $old_domain_node->replaceChild( $new_alert_comment_node, $old_alert_comment_node );
                    }
                }
                else {
                    $old_domain_node->appendChild($new_alert_comment_node);
                }
            }

            #SCOP comment
            if ( $new_domain_node->exists('scop_comment') ) {
                my $new_scop_comment      = $new_domain_node->findvalue('scop_comment');
                my $new_scop_comment_node = $new_domain_node->findnodes('scop_comment')->get_node(1);
                if ( $old_domain_node->exists('scop_comment') ) {
                    my $old_scop_comment      = $old_domain_node->findvalue('scop_comment');
                    my $old_scop_comment_node = $old_domain_node->findnodes('scop_comment')->get_node(1);
                    if ( $old_scop_comment ne $new_scop_comment ) {
                        print "WARNING! overwriting '$old_scop_comment' with '$new_scop_comment'\n";
                        $old_domain_node->replaceChild( $new_scop_comment_node, $old_scop_comment_node );
                    }
                }
                else {
                    $old_domain_node->appendChild($new_scop_comment_node);
                }
            }

            #comment
            if ( $new_domain_node->exists('comment') ) {
                my $new_comment      = $new_domain_node->findvalue('comment');
                my $new_comment_node = $new_domain_node->findnodes('comment')->get_node(1);
                if ( $old_domain_node->exists('comment') ) {
                    my $old_comment      = $old_domain_node->findvalue('comment');
                    my $old_comment_node = $old_domain_node->findnodes('comment')->get_node(1);
                    if ( $old_comment ne $new_comment ) {
                        print "WARNING! overwriting '$old_comment' with '$new_comment'\n";
                        $old_domain_node->replaceChild( $new_comment_node, $old_comment_node );
                    }
                }
                else {
                    $old_domain_node->appendChild($new_comment_node);
                }
            }
        }
        elsif ( $mod_d{$mod_domain_uid}{domain_fshift} ) {

            my $new_domain_uid = $mod_d{$mod_domain_uid}{new_uid};

            if ( !$new_domain_nodes{$new_domain_uid} ) {
                die "ERROR! No new domain node found for UID $new_domain_uid\n";
            }
            my $new_domain_node = $new_domain_nodes{$new_domain_uid};

            if ( !$old_domain_nodes{$mod_domain_uid} ) {
                die "ERROR! No old domain node found for UID $mod_domain_uid\n";
            }
            my $old_domain_node = $old_domain_nodes{$mod_domain_uid};

            my $new_f_id   = $mod_d{$mod_domain_uid}{new_fid};
            my $old_f_node = $old_ecod_xml_doc->findnodes(qq{//f_group[\@f_id="$new_f_id"]})->get_node(1);
            $old_f_node->appendChild($old_domain_node);

            my $old_domain_type = $mod_d{$mod_domain_uid}{old_domain_type};
            my $new_domain_type = $mod_d{$mod_domain_uid}{new_domain_type};

            if ($DEBUG) {
                print "DEBUG mod_d inherit_annotation: $mod_domain_uid $new_domain_uid\n";
            }

            #Annotation types that should be transferred

            #Is the domain a manual representative.. i.e. appears in Hua's file and not as a fluke domain

            if ( $new_domain_node->exists('@manual_rep') ) {
                $old_domain_node->setAttribute( 'manual_rep', $new_domain_node->findvalue('@manual_rep') );
            }

            if ( $new_domain_node->findvalue('@manual_rep') eq 'false' ) {
                if ( $new_domain_node->exists('scop_representative_domain')
                    && !$old_domain_node->exists('scop_representative_domain') )
                {
                    my $scop_rep_node = $new_domain_node->findnodes('scop_representative_domain')->get_node(1);
                    $old_domain_node->appendChild($scop_rep_node);
                }
                elsif ( $new_domain_node->exists('ecod_representative_domain')
                    && !$old_domain_node->exists('ecod_representative_domain') )
                {
                    my $ecod_rep_node = $new_domain_node->findnodes('ecod_representative_domain')->get_node(1);
                    $old_domain_node->appendChild($ecod_rep_node);
                }
            }

            #Upconvert
            if ( $old_domain_type eq 'derived' && ( $new_domain_type eq 'explicit' ) ) {

                $old_domain_node->removeAttribute('derived_range');

                foreach my $old_derived_range_node ( $old_domain_node->findnodes('derived_range')->get_nodelist() ) {
                    $old_derived_range_node->unbindNode;
                }
                foreach my $old_derived_seqid_range_node (
                    $old_domain_node->findnodes('derived_seqid_range')->get_nodelist() )
                {
                    $old_derived_seqid_range_node->unbindNode;
                }

                if ( $old_domain_node->exists('manual_range') ) {
                    foreach my $m_node ( $old_domain_node->findnodes('manual_range')->get_nodelist() ) {
                        $m_node->unbindNode;
                    }
                }

                my $manual_range_node = $new_domain_node->findnodes('manual_range')->get_node(1)->cloneNode(1);

                $old_domain_node->appendChild($manual_range_node);

                $old_domain_node->setAttribute( 'manual_range', 'true' );
                $old_domain_node->setAttribute( 'manual_rep',   'true' );

            }
            elsif ( $old_domain_type eq 'derived' && $new_domain_type eq 'implicit' ) {
                print "ERROR! Impossible upconvert on $mod_domain_uid -> $new_domain_uid\n";
            }

            #alert comment
            if ( $new_domain_node->exists('alert_comment') ) {
                my $new_alert_comment      = $new_domain_node->findvalue('alert_comment');
                my $new_alert_comment_node = $new_domain_node->findnodes('alert_comment')->get_node(1);
                if ( $old_domain_node->exists('alert_comment') ) {
                    my $old_alert_comment      = $old_domain_node->findvalue('alert_comment');
                    my $old_alert_comment_node = $old_domain_node->findnodes('alert_comment')->get_node(1);
                    if ( $old_alert_comment ne $new_alert_comment ) {
                        print "WARNING! overwriting '$old_alert_comment' with '$new_alert_comment'\n";
                        $old_domain_node->replaceChild( $new_alert_comment_node, $old_alert_comment_node );
                    }
                }
                else {
                    $old_domain_node->appendChild($new_alert_comment_node);
                }
            }

            #SCOP comment
            if ( $new_domain_node->exists('scop_comment') ) {
                my $new_scop_comment      = $new_domain_node->findvalue('scop_comment');
                my $new_scop_comment_node = $new_domain_node->findnodes('scop_comment')->get_node(1);
                if ( $old_domain_node->exists('scop_comment') ) {
                    my $old_scop_comment      = $old_domain_node->findvalue('scop_comment');
                    my $old_scop_comment_node = $old_domain_node->findnodes('scop_comment')->get_node(1);
                    if ( $old_scop_comment ne $new_scop_comment ) {
                        print "WARNING! overwriting '$old_scop_comment' with '$new_scop_comment'\n";
                        $old_domain_node->replaceChild( $new_scop_comment_node, $old_scop_comment_node );
                    }
                }
                else {
                    $old_domain_node->appendChild($new_scop_comment_node);
                }
            }

            #comment
            if ( $new_domain_node->exists('comment') ) {
                my $new_comment      = $new_domain_node->findvalue('comment');
                my $new_comment_node = $new_domain_node->findnodes('comment')->get_node(1);
                if ( $old_domain_node->exists('comment') ) {
                    my $old_comment      = $old_domain_node->findvalue('comment');
                    my $old_comment_node = $old_domain_node->findnodes('comment')->get_node(1);
                    if ( $old_comment ne $new_comment ) {
                        print "WARNING! overwriting '$old_comment' with '$new_comment'\n";
                        $old_domain_node->replaceChild( $new_comment_node, $old_comment_node );
                    }
                }
                else {
                    $old_domain_node->appendChild($new_comment_node);
                }
            }
        }
        elsif ( $mod_d{$mod_domain_uid}{nonrep_conversion} ) {

            my $new_domain_uid = $mod_d{$mod_domain_uid}{new_uid};

            if ( !$new_domain_nodes{$new_domain_uid} ) {
                die "ERROR! No new domain node found for UID $new_domain_uid\n";
            }
            my $new_domain_node = $new_domain_nodes{$new_domain_uid};

            if ( !$old_domain_nodes{$mod_domain_uid} ) {
                die "ERROR! No old domain node found for UID $mod_domain_uid\n";
            }
            my $old_domain_node = $old_domain_nodes{$mod_domain_uid};

            my $old_domain_type = $mod_d{$mod_domain_uid}{old_domain_type};
            my $new_domain_type = $mod_d{$mod_domain_uid}{new_domain_type};
            if ( $old_domain_type eq 'derived' && ( $new_domain_type eq 'explicit' ) ) {

                $old_domain_node->removeAttribute('derived_range');

                foreach my $old_derived_seqid_range_node (
                    $old_domain_node->findnodes('derived_seqid_range')->get_nodelist() )
                {
                    $old_derived_seqid_range_node->unbindNode;
                }
                foreach my $old_derived_range_node ( $old_domain_node->findnodes('derived_range')->get_nodelist() ) {
                    $old_derived_range_node->unbindNode;
                }
                if ( $old_domain_node->exists('manual_range') ) {
                    foreach my $m_node ( $old_domain_node->findnodes('manual_range')->get_nodelist() ) {
                        $m_node->unbindNode;
                    }
                }
                my $manual_range_node = $new_domain_node->findnodes('manual_range')->get_node(1)->cloneNode(1);

                $old_domain_node->appendChild($manual_range_node);

                $old_domain_node->setAttribute( 'manual_range', 'true' );
                $old_domain_node->setAttribute( 'manual_rep',   'true' );

            }
            elsif ( $old_domain_type eq 'derived' && $new_domain_type eq 'implicit' ) {
                die "ERROR! Impossible upconvert on $mod_domain_uid -> $new_domain_uid\n";
            }
            else {
                die "ERROR! Confusion case 1\n";
            }

        }
        elsif ( $mod_d{$mod_domain_uid}{domain_assembly_conversion} ) {
            my $new_domain_uid = $mod_d{$mod_domain_uid}{new_uid};
            if ( !$new_domain_nodes{$new_domain_uid} ) {
                die "ERROR! No new domain node found for UID $new_domain_uid\n";
            }
            my $new_domain_node = $new_domain_nodes{$new_domain_uid};
            if ( !$old_domain_nodes{$mod_domain_uid} ) {
                die "ERROR! No old domain node found for $mod_domain_uid\n";
            }
            my $old_domain_node = $old_domain_nodes{$mod_domain_uid};

            my $da_uid        = $mod_d{$mod_domain_uid}{da_uid};
            my $parent_f_id   = $mod_d{$mod_domain_uid}{parent_f_id};
            my $parent_f_node = $old_f_nodes{$parent_f_id};

            my $assembly_type = $new_assembly_types{$da_uid};

   #check if any partner already has been assigned a current domain assembly uid, if so use that and the associated node
            if ( $current_da_uids{$da_uid} ) {
                $old_domain_assembly_nodes{ $current_da_uids{$da_uid} }->appendChild($old_domain_node);
            }
            else {
                my $old_domain_assembly_node = $old_ecod_xml_doc->createElement('domain_assembly');

                #my $old_domain_assembly_uid = sprintf "%09i", $maxUID++;
                my $sql_uid = sprintf "%09i", sql_uid_gen();
                $old_domain_assembly_node->setAttribute( 'uid',                            $sql_uid );
                $old_domain_assembly_node->setAttribute( 'domain_assembly_representative', 'true' );
                $old_domain_assembly_node->setAttribute( 'assembly_type',                  $assembly_type );
                $old_domain_assembly_nodes{$sql_uid} = $old_domain_assembly_node;
                $parent_f_node->appendChild($old_domain_assembly_node);
                $old_domain_assembly_node->appendChild($old_domain_node);
                $old_domain_node->setAttribute( 'primary', 'true' );
                $current_da_uids{$da_uid} = $sql_uid;
            }

            #else assign a new domain assembly node to the correct current f_id

        }

    }

    #Obsolete domains
    foreach my $obs_domain_uid (@obs_d) {

        my $old_domain_node = $old_ecod_xml_doc->findnodes(qq{//domain[\@uid="$obs_domain_uid"]})->get_node(1);

        $old_domain_node->setAttribute( 'isObsolete', 'true' );

    }

#Orphaned nodes, for each obsolete group are all child groups already obsoleted? Note that derived domains are necessarily orphans at this point

    foreach my $obs_x_id (@obs_x) {
        my $XPath = qq{//x_group[\@x_id="$obs_x_id"]};
        if ( $old_ecod_xml_doc->exists($XPath) ) {
            my $x_node  = $old_ecod_xml_doc->findnodes($XPath)->get_node(1);
            my $warning = 0;
            foreach my $child_h_node ( $x_node->findnodes('h_group')->get_nodelist() ) {
                my $child_h_id = $child_h_node->findvalue('@h_id');
                if ( $child_h_node->findvalue('@isObsolete') ne 'true' ) {
                    print "WARNING! Obsoletion of X-group $obs_x_id orphans non-obsolete H-group $child_h_id\n";
                    $warning++;
                }
            }
            if ( !$warning ) {
                $x_node->parentNode->removeChild($x_node);
            }
        }
    }

    foreach my $obs_h_id (@obs_h) {
        my $XPath = qq{//h_group[\@h_id="$obs_h_id"]};
        if ( $old_ecod_xml_doc->exists($XPath) ) {
            my $h_node  = $old_ecod_xml_doc->findnodes($XPath)->get_node(1);
            my $warning = 0;
            foreach my $child_f_node ( $h_node->findnodes('f_group')->get_nodelist() ) {
                my $child_f_id = $child_f_node->findvalue('@f_id');
                if ( $child_f_node->findvalue('@isObsolete') ne 'true' ) {
                    print "WARNING! Obsoletion of H-group $obs_h_id orphans non-obsolete F-group $child_f_id\n";
                    $warning++;
                }
            }
            if ( !$warning ) {
                $h_node->parentNode->removeChild($h_node);
            }
        }
    }

    foreach my $obs_f_id (@obs_f) {
        my $XPath = qq{//f_group[\@f_id="$obs_f_id"]};
        if ( $old_ecod_xml_doc->exists($XPath) ) {
            my $f_node  = $old_ecod_xml_doc->findnodes($XPath)->get_node(1);
            my $warning = 0;
            foreach my $child_domain_node ( $f_node->findnodes('domain')->get_nodelist() ) {
                my $child_domain_id = $child_domain_node->findvalue('@uid');
                if ( $child_domain_node->findvalue('@isObsolete') ne 'true' ) {
                    if ( $child_domain_node->findvalue('derived_range') eq 'true' ) {
                        print "WARNING! Obsoletion of F-group $obs_f_id orphans derived domain $child_domain_id\n";
                    }
                    else {
                        print "WARNING! Obsoletion of F-group $obs_f_id orphans non-obsolete domain $child_domain_id\n";
                    }
                    $warning++;
                }
            }
            if ( !$warning ) {
                $f_node->parentNode->removeChild($f_node);
            }
        }
    }

    $domain_dictionary_node->setAttribute( 'maxUID', $maxUID );

    return $old_ecod_xml_doc;
}

sub fix_manual_range_nodes {
    my $sub = 'fix_manual_range_nodes';

    my ($ecod_xml_doc) = @_;

    #Check current obsolete domains
    print "DEBUG: check obsolete\n";
    my $obsolete_status = '/usr2/pdb/data/status/obsolete.dat';
    open( IN, $obsolete_status ) or die "ERROR! Could not open $obsolete_status for reading:$!\n";
    my %obsolete;
    while ( my $ln = <IN> ) {
        if ( $ln =~ /OBSLTE/ ) {
            my @F = split( /\s+/, $ln );
            $obsolete{ lc( $F[2] ) }++;
        }
    }
    close IN;

    print "DEBUG: check all manual range nodes\n";
    my $domain_XPath = '//domain[@manual_range="true"]';
    foreach my $domain_node ( $ecod_xml_doc->findnodes($domain_XPath) ) {

        my $uid          = $domain_node->findvalue('@uid');
        my $manual_range = $domain_node->findvalue('manual_range');
        print "$uid $manual_range\n";

        if ( $domain_node->exists('scop_domain') ) {
            my $scop_domain_id = $domain_node->findvalue('scop_domain/@scop_domain_id');
            $scop_domain_id =~ /d(\w{4})/;
            my $pdb_id = $1;

            if ( !$domain_node->exists('structure') ) {
                my ( $pdb_id, $chain_id ) = get_scop_domain_pdb_chain($scop_domain_id);
                my $structure_node = create_structure_node( $pdb_id, $chain_id, 'scop', $ecod_xml_doc );
                $domain_node->appendChild($structure_node);
            }
        }
        elsif ( $domain_node->exists('ecod_domain') ) {

            my $ecod_domain_id = $domain_node->findvalue('ecod_domain/@ecod_domain_id');
            $ecod_domain_id =~ /e(\w{4})/;
            my $pdb_id = $1;

            my %chains;
            my @segs = split( /,/, $manual_range );
            foreach my $seg (@segs) {
                while ( $seg =~ /(.+):/g ) {
                    $chains{$1}++;
                }
            }
            my $chain_id;
            if ( scalar keys %chains == 1 ) {
                my @c = keys %chains;
                $chain_id = $c[0];
            }
            else {
                $chain_id = '.';
            }

            #This is deprecated, cannot infer pdb/chain from ECOD domain id
            if ( !$domain_node->exists('structure') ) {
                my $structure_node = create_structure_node( $pdb_id, $chain_id, 'ecod', $ecod_xml_doc );
                $domain_node->appendChild($structure_node);
            }
        }
        my ( $pdb_id, $chain_id ) = get_pdb_chain($domain_node);
        print "$pdb_id $chain_id\n";

        $pdb_id =~ /\w(\w{2})\w/;
        my $two = $1;

        if ( $obsolete{$pdb_id} ) {
            print "WARNING! $pdb_id obsolete, skipping...\n";
            $domain_node->findnodes('structure')->get_node(1)->setAttribute( 'structure_obsolete', 'true' );
            if ( !$USE_OBSOLETE ) { next }
        }
        if ( $domain_node->exists('seqid_range') && $domain_node->findvalue('seqid_range') !~ /\d+/ ) {
            print "unbind $uid!\n";
            $domain_node->findnodes('seqid_range')->get_node(1)->unbindNode;
        }
        if ( $domain_node->exists('range') && $domain_node->findvalue('range') !~ /\d+/ ) {
            print "unbind $uid!\n";
            $domain_node->findnodes('range')->get_node(1)->unbindNode;
        }

        if ( !$domain_node->exists('seqid_range') ) {
            if ( !isMultiChainDomain($domain_node) ) {
                my ( $seqid_range, $pdb_range ) =
                  convert_manual_range_to_seqid_range( $manual_range, $pdb_id, $chain_id );
                print "DEBUG: rangeconvert $pdb_id ($manual_range) => $seqid_range $pdb_range\n";
                $domain_node->appendTextChild( 'seqid_range', $seqid_range );
                if ( !$domain_node->exists('range') ) {
                    $domain_node->appendTextChild( 'range', $pdb_range );
                }
            }
            else {
                #I think this is where the mc_domain bug is.
                my ( $mc_seqid_range, $mc_pdb_range ) =
                  convert_mc_manual_range_to_seqid_range( $manual_range, $pdb_id );
                $domain_node->appendTextChild( 'seqid_range', $mc_seqid_range );
                if ( !$domain_node->exists('range') ) {
                    $domain_node->appendTextChild( 'range', $mc_pdb_range );
                }
            }
        }

    }
}

sub get_scop_domain_pdb_chain { 
	if ($_[0] =~ /d(\w{4})(\w)/) { 
		return ($1, uc($2));
	}else{
		die "ERROR! $_[0] is not a SCOP domain id\n";
	}
}

sub convert_manual_range_to_seqid_range {
    my ( $manual_range, $pdb_id, $chain_id ) = @_;

    my ( $seqid_aref, $struct_seqid_aref, $pdbnum_aref, $asym_id ) = pdbml_seq_parse( $pdb_id, $chain_id );

    my $m_seqid_aref = pdb_range_expand( $manual_range, $pdbnum_aref, $pdb_id, $chain_id );

    my $seqid_range = rangify(@$m_seqid_aref);
    $seqid_range = scopify_range( $seqid_range, $chain_id );

    my $pdb_range = pdb_rangify( $m_seqid_aref, $pdbnum_aref );
    $pdb_range = scopify_range( $pdb_range, $chain_id );

    return ( $seqid_range, $pdb_range );

}

sub convert_mc_manual_range_to_seqid_range {
    my ( $manual_range, $pdb_id ) = @_;

    my @segs = split( /\,/, $manual_range );
    my %chains;
    my @c;
    foreach my $seg (@segs) {
        if ( $seg =~ /(.+):/ ) {

            #$chains{$1}++;
            push( @c, $1 );
        }
        else {
            die "regexp error $seg\n";
        }
    }

    my ( $seqid_aref, $struct_seqid_aref, $pdbnum_href, $asym_id, $chain_aref, $struct_chain_aref ) =
      pdbml_mc_seq_parse( $pdb_id, \@c );

    my ( $m_range, $m_chain ) = multi_chain_pdb_range_expand( $manual_range, $seqid_aref, $chain_aref, $pdbnum_href );
    my $seqid_range = multi_chain_rangify( $m_range, $m_chain );
    my $pdb_range = multi_chain_pdb_rangify( $m_range, $pdbnum_href, $m_chain );

    return ( $seqid_range, $pdb_range );

}

sub isMultiChainDomain {
    if ( $_[0]->findvalue('structure/@chain_id') eq '.' ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub create_structure_node {
    my ( $pdb_id, $chain_id, $type, $xml_doc ) = @_;
    my $structure_node = $xml_doc->createElement('structure');
    if ( $type ne 'ecod' ) {
        $structure_node->setAttribute( 'chain_case_ambiguous', 'true' );
    }
    $structure_node->setAttribute( 'db',       'PDB' );
    $structure_node->setAttribute( 'pdb_id',   $pdb_id );
    $structure_node->setAttribute( 'chain_id', $chain_id );

    return $structure_node;

}

sub fix_implicit_range_nodes {
    my $sub = 'fix_implicit_range_nodes';

    my ($ecod_xml_doc) = @_;

    my $scop_cla_file = '/home/rschaeff/data/scop/v1.75/dir.cla.scop.txt_1.75.txt';
    if ( !-f $scop_cla_file ) { die "ERROR! Could not find scop cla file...\n"; }

    open( IN, $scop_cla_file ) or die "ERROR! Could not open $scop_cla_file for reading:$!\n";

    my %scop_cla_range;
    my %scop_cla_pdb_range;
    while ( my $ln = <IN> ) {
        $ln =~ /^#/ and next;
        my @F              = split( /\s+/, $ln );
        my $scop_domain_id = $F[0];
        my $scop_pdb_range = $F[2];

        $scop_cla_range{$scop_domain_id} = $scop_pdb_range;

    }
    close IN;

    my $domain_XPath = '//domain[@scop_implicit_range="true"]';

    foreach my $domain_node ( $ecod_xml_doc->findnodes($domain_XPath)->get_nodelist() ) {

        if ( $domain_node->exists('scop_domain') ) {
            my $scop_domain_id = $domain_node->findvalue('scop_domain/@scop_domain_id');
            if ( !$domain_node->exists('structure') ) {

                my ( $pdb, $chain );
                if ( $scop_domain_id =~ /d(\w{4})(.)/ ) {
                    $pdb   = $1;
                    $chain = uc($2);
                    print "pdb: $pdb chain: $chain\n";
                }
                else {
                    die "$scop_domain_id doesn't regexp as a SCOP domain...\n";
                }
                my $structure_node = $ecod_xml_doc->createElement('structure');
                $structure_node->setAttribute( 'chain_case_ambiguous', 'true' );
                $structure_node->setAttribute( 'chain_id',             $chain );
                $structure_node->setAttribute( 'pdb_id',               $pdb );
                $structure_node->setAttribute( 'DB',                   'PDB' );

                $domain_node->appendChild($structure_node);
            }

            if ( $scop_cla_range{$scop_domain_id} ) {

                my $pdb        = $domain_node->findvalue('structure/@pdb_id');
                my $chain      = $domain_node->findvalue('structure/@chain_id');
                my $scop_range = $scop_cla_range{$scop_domain_id};

                #my $scop_seqid_range	= $scop_cla_range{$scop_domain_id};

                if ( !$domain_node->exists('range') || !$domain_node->exists('seqid_range') ) {

                    print "pc: $pdb $chain\n";

                    #my ($scop_pdb_range, $scop_seqid_range) = scop_range_convert($pdb, $scop_range);
                    my ( $scop_seqid_range, $scop_pdb_range ) = scop_range_convert( $pdb, $scop_range );
                    if ( !$domain_node->exists('range') ) {
                        my $range_node = $ecod_xml_doc->createElement('range');
                        print "r:$scop_domain_id $scop_pdb_range\n";
                        $range_node->appendTextNode($scop_pdb_range);
                        $domain_node->appendChild($range_node);
                    }

                    if ( !$domain_node->exists('seqid_range') ) {
                        my $seqid_range_node = $ecod_xml_doc->createElement('seqid_range');

                        #print "sr:$scop_domain_id $scop_seqid_range\n";
                        $seqid_range_node->appendTextNode($scop_seqid_range);
                        $domain_node->appendChild($seqid_range_node);
                    }
                }
            }
            else {
                print "WARNING! $scop_domain_id not in cla file\n";
            }
        }
    }
}

sub scop_range_convert {

    my ( $scop_pdb, $range_str ) = @_;

    my $sub = "scop_range_convert";

    if ($DEBUG) {
        print "DEBUG $sub: $scop_pdb $range_str\n";
    }

    #implicit range case
    my %chains;
    if ( $range_str =~ /,/ ) {
        my @rsegs = split( /\,/, $range_str );
        my ( $range_str, $pdb_range_str );

        my ( @nsegs, @npsegs );
        foreach my $rseg (@rsegs) {
            if ( $rseg =~ /^([A-Za-z0-9]):$/ ) {
                my $chain = $1;
                my ( $seqid_aref, $struct_seqid_aref, $pdbnum_aref ) = pdbml_seq_parse( $scop_pdb, $chain );

                if ( !$seqid_aref ) { return 0 }

                my $domain_range = rangify(@$struct_seqid_aref);
                my $domain_pdb_range = pdb_rangify( $struct_seqid_aref, $pdbnum_aref );

                my @segs  = split( ",", $domain_range );
                my @psegs = split( ",", $domain_pdb_range );

                foreach my $seg (@segs) {
                    push( @nsegs, "$chain:$seg" );
                }
                foreach my $pseg (@psegs) {
                    push( @npsegs, "$chain:$pseg" );
                }
            }
            elsif ( $rseg =~ /^([A-Za-z0-9]):(\-?\d+[A-Z]?)\-(\-?\d+[A-Z]?)/ ) {

                my $chain = $1;
                my $start = $2;
                my $end   = $3;

                my ( $seqid_aref, $stuct_eqid_aref, $pdbnum_aref ) = pdbml_seq_parse( $scop_pdb, $chain );

                if ( !$seqid_aref ) {
                    return 0;
                }

                my %pdbnum_lookup;
                for ( my $i = 0 ; $i < scalar(@$pdbnum_aref) ; $i++ ) {
                    if ( defined( $$pdbnum_aref[$i] ) ) {
                        $pdbnum_lookup{ $$pdbnum_aref[$i] } = $i;
                    }
                }

                my $start_seqid = $pdbnum_lookup{$start};
                my $end_seqid   = $pdbnum_lookup{$end};

                push( @nsegs,  "$chain:$start_seqid-$end_seqid" );
                push( @npsegs, "$chain:$$pdbnum_aref[$start_seqid]-$$pdbnum_aref[$end_seqid]" );
            }
        }
        $range_str     = join( ",", @nsegs );
        $pdb_range_str = join( ",", @npsegs );

        return ( $range_str, $pdb_range_str );

    }
    else {
        if ( $range_str =~ /^([A-Za-z0-9]):$/ ) {

            my $chain = $1;
            my ( $seqid_aref, $struct_seqid_aref, $pdbnum_aref ) = pdbml_seq_parse( $scop_pdb, $chain );
            if ( !$seqid_aref ) { return 0 }

            my $domain_range = rangify(@$struct_seqid_aref);
            my $domain_pdb_range = pdb_rangify( $struct_seqid_aref, $pdbnum_aref );

            my @segs  = split( ",", $domain_range );
            my @psegs = split( ",", $domain_pdb_range );
            my @nsegs;
            my @npsegs;
            foreach my $seg (@segs) {
                push( @nsegs, "$chain:$seg" );
            }
            foreach my $pseg (@psegs) {
                push( @npsegs, "$chain:$pseg" );
            }

            my $range_str     = join( ",", @nsegs );
            my $pdb_range_str = join( ",", @npsegs );
            return ( $range_str, $pdb_range_str );
        }
        elsif ( $range_str =~ /^[A-Za-z0-9]:\-?\d+[A-Z]?\-\-?\d+[A-Z]?/ ) {
            my @segs;
            my $i = 0;

            my %chains;

            while ( $range_str =~ /([A-Za-z0-9]):(\-?\d+)\-(\-?\d+)/g ) {
                my $chain = $1;
                my $start = $2;
                my $end   = $3;

                $segs[$i]{chain} = $chain;
                $segs[$i]{start} = $start;
                $segs[$i]{end}   = $end;

                $chains{$chain}++;

                $i++;
            }
            if ( scalar( keys %chains ) > 1 ) {
                my @nsegs;
                my @npsegs;
                for ( my $i = 0 ; $i < scalar(@segs) ; $i++ ) {
                    my $chain = $segs[$i]{chain};
                    my ( $seqid_aref, $struct_seqid_aref, $pdbnum_aref ) = pdbml_seq_parse( $scop_pdb, $chain );
                    if ( !$seqid_aref ) {
                        return 0;
                    }
                    my %pdbnum_lookup;
                    for ( my $i = 0 ; $i < scalar(@$pdbnum_aref) ; $i++ ) {
                        if ( $$pdbnum_aref[$i] ) {
                            $pdbnum_lookup{ $$pdbnum_aref[$i] } = $i;
                        }
                    }
                    my $start = $segs[$i]{start};
                    my $end   = $segs[$i]{end};

                    my $start_seqid = $pdbnum_lookup{$start};
                    my $end_seqid   = $pdbnum_lookup{$end};

                    push( @nsegs,  "$chain:$start_seqid-$end_seqid" );
                    push( @npsegs, "$chain:$$pdbnum_aref[$start_seqid]-$$pdbnum_aref[$end_seqid]" );

                }
                my $range_str     = join( ",", @nsegs );
                my $pdb_range_str = join( ",", @npsegs );
                return ( $range_str, $pdb_range_str );
            }
            else {
                my ( $domain_pdb_range, $domain_chain ) = scop_range_split($range_str);
                my ( $seqid_aref, $struct_seqid_aref, $pdbnum_aref ) = pdbml_seq_parse( $scop_pdb, $domain_chain );
                if ( !$seqid_aref ) {
                    return 0;
                }

                my $domain_seqid_aref = pdb_range_expand( $domain_pdb_range, $pdbnum_aref, $scop_pdb, $domain_chain );

                my $debug2 = rangify(@$domain_seqid_aref);
                my %u;
                my $domain_struct_seqid_aref;
                foreach my $s (@$domain_seqid_aref) {
                    $u{$s}++;
                }
                foreach my $s (@$struct_seqid_aref) {
                    if ( $u{$s} ) {
                        push( @$domain_struct_seqid_aref, $s );
                    }
                }
                my $domain_range         = rangify(@$domain_struct_seqid_aref);
                my $new_domain_pdb_range = pdb_rangify( $domain_struct_seqid_aref, $pdbnum_aref );
                my @segs                 = split( /,/, $domain_range );
                my @psegs                = split( /,/, $new_domain_pdb_range );
                my @nsegs;
                my @npsegs;
                foreach my $seg (@segs) {
                    push( @nsegs, "$domain_chain:$seg" );
                }
                foreach my $pseg (@psegs) {
                    push( @npsegs, "$domain_chain:$pseg" );
                }
                my $range_str     = join( ',', @nsegs );
                my $pdb_range_str = join( ',', @npsegs );
                return ( $range_str, $pdb_range_str );
            }
        }
    }
}

sub add_ecod_domain_ids {
    my $sub = 'add_ecod_domain_ids';

    my ($ecod_xml_doc) = @_;

    my %known_ids;
    foreach my $d_node ( find_domain_nodes($ecod_xml_doc) ) {
        if ( $d_node->exists('@ecod_domain_id') ) {
            my $ecod_domain_id = $d_node->findvalue('@ecod_domain_id');
            $known_ids{$ecod_domain_id}++;
        }
    }

    my %pdb_chains;
    my %ranges;
    foreach my $d_node ( find_domain_nodes($ecod_xml_doc) ) {

        my $uid = $d_node->findvalue('@uid');

        my $ecod_domain_id = 'NA';
        if ( $d_node->exists('@ecod_domain_id') ) {
            $ecod_domain_id = $d_node->findvalue('@ecod_domain_id');
        }
        elsif ( $d_node->exists('ecod_domain/@ecod_domain_id') && !$FORCE_OVERWRITE ) {
            $ecod_domain_id = $d_node->findvalue('ecod_domain/@ecod_domain_id');
        }
        elsif ($d_node->exists('ecod_domain/@ecod_domain_id')
            && !$d_node->exists('scop_domain/@scop_domain_id')
            && $FORCE_OVERWRITE )
        {
            $ecod_domain_id = $d_node->findvalue('ecod_domain/@ecod_domain_id');
        }

        my $scop_domain_id = 'NA';
        if ( $d_node->exists('scop_domain/@scop_domain_id') ) {
            $scop_domain_id = $d_node->findvalue('scop_domain/@scop_domain_id');
        }

        #print "$uid $ecod_domain_id $scop_domain_id\n";

        my ( $pdb, $chain );

        if ( $ecod_domain_id ne 'NA' ) {
            $ecod_domain_id =~ /e(\w{4})(.+)/;
            $pdb   = $1;
            $chain = $2;
        }
        elsif ( $scop_domain_id ne 'NA' ) {
            $scop_domain_id =~ /d(\w{4})+(.)/;
            $pdb   = $1;
            $chain = $2;
        }
        else {
            print "WARNING! No scop or ecod domain id for $uid? orphaned?\n";
            next;
        }

        $chain = uc($chain);
        if ( $d_node->exists('structure/@chain_id') ) {
            $chain = $d_node->findvalue('structure/@chain_id');
        }
        if ( $d_node->findvalue('structure/@structure_obsolete') eq 'true' ) {

            #print "WARNNING! $pdb obsolete\n";
        }

        my $range = 'NA';
        if ( $d_node->exists('range') ) {
            $range = $d_node->findvalue('range');
        }
        elsif ( $d_node->exists('derived_range') ) {
            $range = $d_node->findvalue('derived_range');
        }
        elsif ( $d_node->exists('manual_range') ) {
            $range = $d_node->findvalue('manual_range');
        }
        else {
            print "WARNING! No range for $uid\n";
            next;
        }

        if ( $ecod_domain_id eq 'NA' && $scop_domain_id ne 'NA' && !$ranges{$pdb}{$chain}{$range} ) {
            $pdb_chains{$pdb}{$chain}++;
            $ecod_domain_id = "e$pdb$chain$pdb_chains{$pdb}{$chain}";
            while ( $known_ids{$ecod_domain_id} ) {
                $pdb_chains{$pdb}{$chain}++;
                $ecod_domain_id = "e$pdb$chain$pdb_chains{$pdb}{$chain}";
            }
            $ranges{$pdb}{$chain}{$range} = $ecod_domain_id;
        }
        elsif ( $range ne 'NA' && $ecod_domain_id eq 'NA' ) {

            #What is this case for?
            $ecod_domain_id = $ranges{$pdb}{$chain}{$range};
        }

        if ( !$d_node->exists('ecod_domain') ) {
            my $ed_node = $ecod_xml_doc->createElement('ecod_domain');
            $ed_node->setAttribute( 'ecod_domain_id',       $ecod_domain_id );
            $ed_node->setAttribute( 'chain_case_ambiguous', 'true' );
            $d_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
            $d_node->appendChild($ed_node);
        }
        elsif ( $d_node->exists('ecod_domain') && $FORCE_OVERWRITE ) {
            my $ed_node = $d_node->findnodes('ecod_domain')->get_node(1);
            $d_node->removeChild($ed_node);
            $ed_node = $ecod_xml_doc->createElement('ecod_domain');
            $ed_node->setAttribute( 'ecod_domain_id',       $ecod_domain_id );
            $ed_node->setAttribute( 'chain_case_ambiguous', 'true' );
            $d_node->appendChild($ed_node);
            $d_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
        }
        if ( !$d_node->exists('structure') ) {
            my $s_node = $ecod_xml_doc->createElement('structure');
            $s_node->setAttribute( 'pdb_id',               $pdb );
            $s_node->setAttribute( 'chain_id',             $chain );
            $s_node->setAttribute( 'chain_case_ambiguous', 'true' );
            $s_node->setAttribute( 'db',                   'PDB' );
            $d_node->appendChild($s_node);
        }

        if ( !$d_node->exists('@ecod_domain_id') || !$d_node->findvalue('@ecod_domain_id') =~ /\w+/ ) {
            $d_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
        }
        elsif ( $d_node->exists('@ecod_domain_id') && $FORCE_OVERWRITE ) {
            $d_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
        }

        #print "$uid $ecod_domain_id\n";
    }
}

#This should no longer be used
sub pf_reorder {
    my $sub = 'pf_reorder';
    my ($ecod_xml_doc) = @_;
    foreach my $f_group ( $ecod_xml_doc->findnodes('//f_group')->get_nodelist() ) {

        my @order;
        my %pf_ids;
        foreach my $pf_group ( $f_group->findnodes('pf_group[@pfam_cluster="true"]')->get_nodelist() ) {
            push( @order, $pf_group );
            my $pf_id = $pf_group->findvalue('@pf_id');
            $pf_ids{$pf_id}++;

            if ( $pf_ids{$pf_id} > 1 ) {
                print "WARNING! Non unique $pf_id\n";
                $pf_group->setAttribute( 'tmp_reid', 'true' );
            }
        }

        foreach my $pf_group ( $f_group->findnodes('pf_group[@hh_cluster="true"]')->get_nodelist() ) {
            push( @order, $pf_group );
            my $pf_id = $pf_group->findvalue('@pf_id');
            $pf_ids{$pf_id}++;
            if ( $pf_ids{$pf_id} > 1 ) {
                print "WARNING! Non unique $pf_id\n";
                $pf_group->setAttribute( 'tmp_reid', 'true' );
            }
        }

        foreach my $pf_group ( $f_group->findnodes('pf_group[@jb_cluster="true"]')->get_nodelist() ) {
            push( @order, $pf_group );
            my $pf_id = $pf_group->findvalue('@pf_id');
            $pf_ids{$pf_id}++;
            if ( $pf_ids{$pf_id} > 1 ) {
                print "WARNING! Non unique $pf_id\n";
                $pf_group->setAttribute( 'tmp_reid', 'true' );
            }
        }
        my $id_max = 0;
        foreach my $pf_id ( keys %pf_ids ) {
            if ( $pf_id > $id_max ) {
                $id_max = $pf_id + 1;
            }
        }

        for ( my $i = 0 ; $i < scalar(@order) ; $i++ ) {
            $order[$i]->setAttribute( 'pf_ordinal', $i );
            if ( $order[$i]->exists('@tmp_reid') ) {
                $order[$i]->setAttribute( 'pf_id', $id_max++ );
                $order[$i]->removeAttribute('tmp_reid');
            }
        }

    }
}

sub rep_check_ecodf {
    my $sub = 'rep_check_ecodf';

    my ( $ecod_xml_doc, $UNBIND_EMPTY_PFGROUP ) = @_;

    my $warning      = 0;
    my $total        = 0;
    my $external_rep = 0;

    my %ecod_manual_reps;
    my %scop_compkey_to_ecod_uid;
    my %rep_lookup;
    foreach my $manual_domain_node ( $ecod_xml_doc->findnodes('//domain[@manual_rep="true"]') ) {
        my $uid            = $manual_domain_node->findvalue('@uid');
        my $ecod_domain_id = $manual_domain_node->findvalue('@ecod_domain_id');
        $rep_lookup{$uid}++;
        if ( !$ecod_manual_reps{$ecod_domain_id} ) {
            $ecod_manual_reps{$ecod_domain_id} = $uid;
        }
        if ( $manual_domain_node->findvalue('@manual_range') eq 'true' ) {
            my $manual_range = $manual_domain_node->findvalue('manual_range');
            if ( $manual_domain_node->exists('scop_domain/@scop_domain_id') ) {
                my $scop_domain_id = $manual_domain_node->findvalue('scop_domain/@scop_domain_id');
                my $comp_key       = $scop_domain_id . $manual_range;
                $scop_compkey_to_ecod_uid{$comp_key} = $uid;
            }
        }
        elsif ( $manual_domain_node->findvalue('@scop_implicit_range') eq 'true' ) {
            if ( $manual_domain_node->exists('scop_domain/@scop_domain_id') ) {
                my $scop_domain_id = $manual_domain_node->findvalue('scop_domain/@scop_domain_id');
                $scop_compkey_to_ecod_uid{$scop_domain_id} = $uid;
            }
        }
    }
    foreach my $prov_domain_node ( $ecod_xml_doc->findnodes('//domain[@provisional_manual_rep="true"]') ) {
        my ( $uid, $ecod_domain_id ) = get_ids($prov_domain_node);
        $ecod_manual_reps{$ecod_domain_id} = $uid;
        $rep_lookup{$uid}++;

    }
    foreach my $f_group ( $ecod_xml_doc->findnodes('//f_group') ) {
        my $pf_iter = 1;
        if ( $f_group->exists('pf_group') ) {
            foreach my $pf_group ( $f_group->findnodes('pf_group') ) {
                my $pf_id = $pf_group->findvalue('@pf_id');
                if ( $pf_id =~ /\d+\.\d+\.\d+\.(\d+)/ ) {
                    my $local_pf_iter = $1;
                    if ( $local_pf_iter >= $pf_iter ) {
                        $pf_iter = $local_pf_iter + 1;
                    }
                }
            }
        }

        foreach my $pf_group ( $f_group->findnodes('pf_group') ) {

            #Can no longer remove pf_groups by emptying domains, must remove ecodf family.
            if (   $pf_group->findnodes('domain')->size() == 0
                && $pf_group->findnodes('domain_assembly')->size() == 0
                && $UNBIND_EMPTY_PFGROUP )
            {    # DOMAIN ASSEMBLIES
                print "WARNING empty pfgroup removing\n";
                $pf_group->unbindNode;
                next;
            }

            my $pf_id = $pf_group->findvalue('@pf_id');

            #my $f_id = $pf_group->parentNode->findvalue('@f_id');
            #print "$f_id $pf_id\n";
            #my $pf_id = "$f_id.$pf_iter";
            #$pf_group->setAttribute('pf_id', $pf_id);

            #$pf_id++;

            my $rep_domain_node;
            my $rep_uid;
            my $rep_domain_id;
            my %manual_rep_uids;
            my %manual_domain_ids;
            my %side_load_uids;
            my %side_load_domain_ids;

            if ( $pf_group->exists('domain[@manual_rep="true"]') ) {
                $rep_domain_node = $pf_group->findnodes('domain[@manual_rep="true"]')->get_node(1);
                $rep_uid         = $rep_domain_node->findvalue('@uid');
                $rep_domain_id   = $rep_domain_node->findvalue('@ecod_domain_id');
                foreach my $manual_rep_node ( $pf_group->findnodes('.//domain[@manual_rep="true"]') ) {

                    my $uid       = $manual_rep_node->findvalue('@uid');
                    my $domain_id = $manual_rep_node->findvalue('@ecod_domain_id');

                    $manual_rep_uids{$uid}         = $domain_id;
                    $manual_domain_ids{$domain_id} = $uid;

                    if ( $manual_rep_node->findvalue('@yx_side_load') eq 'true' ) {
                        $side_load_uids{$uid}             = $domain_id;
                        $side_load_domain_ids{$domain_id} = $uid;
                    }
                }
            }
            elsif ( $pf_group->exists('domain_assembly/domain[@manual_rep="true"]') ) {
                $rep_domain_node = $pf_group->findnodes('domain_assembly/domain[@manual_rep="true"]')->get_node(1);
                $rep_uid         = $rep_domain_node->findvalue('@uid');
                $rep_domain_id   = $rep_domain_node->findvalue('@ecod_domain_id');

                foreach my $manual_rep_node ( $pf_group->findnodes('.//domain[@manual_rep="true"]')->get_nodelist() ) {
                    my $uid       = $manual_rep_node->findvalue('@uid');
                    my $domain_id = $manual_rep_node->findvalue('@ecod_domain_id');

                    $manual_rep_uids{$uid}         = $domain_id;
                    $manual_domain_ids{$domain_id} = $uid;
                    if ( $manual_rep_node->findvalue('@yx_side_load') eq 'true' ) {
                        $side_load_uids{$uid}             = $domain_id;
                        $side_load_domain_ids{$domain_id} = $uid;
                    }
                }
            }
            elsif ( $pf_group->exists('domain[@provisional_manual_rep="true"]') ) {

                #print "WARNING! $pf_id lacks manual rep (has prov)\n";

                $rep_domain_node = $pf_group->findnodes('domain[@provisional_manual_rep="true"]')->get_node(1);
                $rep_uid         = $rep_domain_node->findvalue('@uid');
                $rep_domain_id   = $rep_domain_node->findvalue('@ecod_domain_id');

                $manual_rep_uids{$rep_uid}         = $rep_domain_id;
                $manual_domain_ids{$rep_domain_id} = $rep_uid;
                $warning++;
            }
            else {
                print "WARNING! $pf_id lacks manual rep\n";

                #cheat

                my $rep_domain_node;
                if ( $pf_group->exists('domain') ) {
                    $rep_domain_node = $pf_group->findnodes('domain')->get_node(1);
                }
                elsif ( $pf_group->exists('domain_assembly/domain') ) {
                    $rep_domain_node = $pf_group->findnodes('domain_assembly/domain')->get_node(1);
                }
                else {
                    print "WARNING! empty pf_group\n";
                    next;
                }
                $rep_domain_node->setAttribute( 'provisional_manual_rep', 'true' );
                $rep_uid       = $rep_domain_node->findvalue('@uid');
                $rep_domain_id = $rep_domain_node->findvalue('@ecod_domain_id');

                $manual_rep_uids{$rep_uid}         = $rep_domain_id;
                $manual_domain_ids{$rep_domain_id} = $rep_uid;
                $warning++;
            }
            $total++;

            if ( $rep_domain_id =~ /\w+/ && $rep_uid =~ /\d+/ ) {
                foreach my $domain ( $pf_group->findnodes('domain')->get_nodelist() ) {
                    my $uid       = $domain->findvalue('@uid');
                    my $domain_id = $domain->findvalue('@ecod_domain_id');

                    #If a domain is a known representative, don't alter it.

                    #if ($manual_rep_uids{$uid} && !$side_load_uids{$uid}) { next }
                    if ( $rep_lookup{$uid} ) { next }
                    if ( $manual_domain_ids{$domain_id} && !$side_load_domain_ids{$domain_id} ) { next }
                    my $rep_node;
                    if ( !$domain->exists('ecod_representative_domain') ) {

                        #print "WARNING $uid no rep domain node and not manual? wtf?\n";
                        $rep_node = $ecod_xml_doc->createElement('ecod_representative_domain');
                        $domain->appendChild($rep_node);
                    }
                    else {
                        $rep_node = $domain->findnodes('ecod_representative_domain')->get_node(1);
                    }
                    my $current_rep_uid;
                    if (   $domain->exists('ecod_representative_domain/@uid')
                        && $domain->findvalue('ecod_representative_domain/@uid') =~ /\d+/ )
                    {
                        $current_rep_uid = sprintf "%09i", $domain->findvalue('ecod_representative_domain/@uid');
                    }
                    elsif ($domain->exists('derived_range/@derivedFrom')
                        && $domain->findvalue('derived_range/@derivedFrom') =~ /\d+/ )
                    {
                        $current_rep_uid = sprintf "%09i", $domain->findvalue('derived_range/@derivedFrom');
                    }
                    elsif ( $domain->exists('ecod_representative_domain/@ecod_domain_id') ) {
                        my $rep_domain_id = $domain->findvalue('ecod_representative_domain/@ecod_domain_id');
                        if ( $ecod_manual_reps{$rep_domain_id} ) {
                            $current_rep_uid = $ecod_manual_reps{$rep_domain_id};
                            $rep_node->setAttribute( 'uid', $current_rep_uid );
                        }
                        else {
                            die "ERROR! Rep recorded for $domain_id is not a manual_rep\n";
                        }
                    }
                    elsif ( $domain->exists('scop_representative_domain/@scop_domain_id') ) {

                        #my $manual_range = $domain->findvalue('manual_range');
                        my $scop_domain_id = $domain->findvalue('scop_representative_domain/@scop_domain_id');
                        if ( $scop_compkey_to_ecod_uid{$scop_domain_id} ) {
                            $current_rep_uid = $scop_compkey_to_ecod_uid{$scop_domain_id};
                        }
                        else {
                            if ( $domain->exists('scop_representative_domain/rep_manual_range') ) {
                                my $manual_range = $domain->findvalue('scop_representative_domain/rep_manual_range');
                                my $compkey      = $scop_domain_id . $manual_range;
                                $current_rep_uid = $scop_compkey_to_ecod_uid{$compkey};
                            }
                            else {
                                die
"ERROR! Could not deference manual range and scop domain id to an ecod uid for $domain_id\n";
                            }

                        }
                    }
                    elsif ( $side_load_uids{$uid} && $side_load_domain_ids{$domain_id} ) {
                        print "WARNING! $uid $domain_id is a side_load domain with no rep, skipping...\n";
                        next;

                    }
                    else {
                        die "ERROR! No rep recorded for $uid $domain_id\n";
                    }
                    if (   $rep_node->findvalue('@ecod_domain_id') !~ /\w+/
                        && $manual_rep_uids{$rep_uid}
                        && $manual_rep_uids{$rep_uid} =~ /\w+/ )
                    {
                        $rep_node->setAttribute( 'ecod_domain_id', $manual_rep_uids{$rep_uid} );
                    }

                    #if (!$current_rep_uid) {
                    #	die "$uid $domain_id\n";
                    #}
                    my $current_rep_domain_id = $domain->findvalue('ecod_representative_domain/@ecod_domain_id');

                    #If the current manual rep is not known either by recorded uid or domain id

                    if ( $manual_rep_uids{$current_rep_uid} && $manual_domain_ids{$current_rep_domain_id} ) {
                        next;
                    }
                    elsif ( !$manual_rep_uids{$current_rep_uid} && !$manual_domain_ids{$current_rep_domain_id} ) {

                        #print "WARNING! $uid has external rep $rep_uid\n";
                        $external_rep++;
                        $rep_node->setAttribute( 'uid',            $rep_uid );
                        $rep_node->setAttribute( 'ecod_domain_id', $rep_domain_id );
                    }
                    elsif ( !$manual_rep_uids{$current_rep_uid} && $manual_domain_ids{$current_rep_domain_id} ) {
                        $rep_node->setAttribute( 'uid', $manual_domain_ids{$current_rep_domain_id} );
                    }
                    elsif ( $manual_rep_uids{$rep_uid} && !$manual_domain_ids{$current_rep_domain_id} ) {

                        #print "DEBUG: set rep $rep_uid $manual_rep_uids{$rep_uid}\n";
                        $rep_node->setAttribute( 'ecod_domain_id', $manual_rep_uids{$rep_uid} );
                        $rep_node->setAttribute( 'uid',            $rep_uid );
                    }
                    else {
                        die "$domain_id $uid?";
                    }

                }
            }
            else {
                die "$rep_uid $rep_domain_id unk?\n";
            }
        }
    }
}
###
sub rep_check {

    my $sub            = 'rep_check';
    my ($ecod_xml_doc) = @_;
    my $warning        = 0;
    my $total          = 0;
    my $external_rep   = 0;

    foreach my $f_group ( $ecod_xml_doc->findnodes('//f_group')->get_nodelist() ) {
        my $pf_id = 1;
        foreach my $pf_group ( $f_group->findnodes('pf_group')->get_nodelist() ) {

            if ( $pf_group->findnodes('domain')->size() == 0 && $pf_group->findnodes('domain_assembly')->size() == 0 )
            {    # DOMAIN ASSEMBLIES
                print "WARNING empty pfgroup removing\n";
                $pf_group->unbindNode;
                next;
            }

            #my $pf_id = $pf_group->findvalue('@pf_id');
            my $f_id = $pf_group->parentNode->findvalue('@f_id');

            #print "$f_id $pf_id\n";
            $pf_group->setAttribute( 'pf_id', $pf_id );
            $pf_id++;

            my $rep_domain_node;
            my $rep_uid;
            my $rep_domain_id;
            my %manual_rep_uids;
            my %manual_domain_ids;
            if ( $pf_group->exists('domain[@manual_rep="true"]') ) {
                $rep_domain_node = $pf_group->findnodes('domain[@manual_rep="true"]')->get_node(1);
                $rep_uid         = $rep_domain_node->findvalue('@uid');
                $rep_domain_id   = $rep_domain_node->findvalue('@ecod_domain_id');
                foreach my $manual_rep_node ( $pf_group->findnodes('domain[@manual_rep="true"]')->get_nodelist() ) {
                    my $uid       = $manual_rep_node->findvalue('@uid');
                    my $domain_id = $manual_rep_node->findvalue('@ecod_domain_id');
                    $manual_rep_uids{$uid}         = $domain_id;
                    $manual_domain_ids{$domain_id} = $uid;
                }
            }
            elsif ( $pf_group->exists('domain_assembly/domain[@manual_rep="true"]') ) {
                $rep_domain_node = $pf_group->findnodes('domain_assembly/domain[@manual_rep="true"]')->get_node(1);
                $rep_uid         = $rep_domain_node->findvalue('@uid');
                $rep_domain_id   = $rep_domain_node->findvalue('@ecod_domain_id');
                foreach my $manual_rep_node ( $pf_group->findnodes('domain[@manual_rep="true"]')->get_nodelist() ) {
                    my $uid       = $manual_rep_node->findvalue('@uid');
                    my $domain_id = $manual_rep_node->findvalue('@ecod_domain_id');
                    $manual_rep_uids{$uid}         = $domain_id;
                    $manual_domain_ids{$domain_id} = $uid;
                }
            }
            elsif ( $pf_group->exists('domain[@provisional_manual_rep="true"]') ) {

                #print "WARNING! $f_id $pf_id lacks manual rep (has prov)\n";
                $rep_domain_node = $pf_group->findnodes('domain[@provisional_manual_rep="true"]')->get_node(1);
                $rep_uid         = $rep_domain_node->findvalue('@uid');
                $rep_domain_id   = $rep_domain_node->findvalue('@ecod_domain_id');
                $manual_rep_uids{$rep_uid}         = $rep_domain_id;
                $manual_domain_ids{$rep_domain_id} = $rep_uid;
                $warning++;
            }
            else {
                print "WARNING! $f_id $pf_id lacks manual rep\n";

                #cheat
                my $rep_domain_node;
                if ( $pf_group->exists('domain') ) {
                    $rep_domain_node = $pf_group->findnodes('domain')->get_node(1);
                }
                elsif ( $pf_group->exists('domain_assembly/domain') ) {
                    $rep_domain_node = $pf_group->findnodes('domain_assembly/domain')->get_node(1);
                }
                else {
                    print "WARNING! ? next\n";
                    next;
                }
                $rep_domain_node->setAttribute( 'provisional_manual_rep', 'true' );
                $rep_uid                           = $rep_domain_node->findvalue('@uid');
                $rep_domain_id                     = $rep_domain_node->findvalue('@ecod_domain_id');
                $manual_rep_uids{$rep_uid}         = $rep_domain_id;
                $manual_domain_ids{$rep_domain_id} = $rep_uid;
                $warning++;
            }
            $total++;

            if ( $rep_domain_id =~ /\w+/ && $rep_uid =~ /\d+/ ) {
                foreach my $domain ( $pf_group->findnodes('domain')->get_nodelist() ) {
                    my $uid       = $domain->findvalue('@uid');
                    my $domain_id = $domain->findvalue('@ecod_domain_id');

                    #If a domain is a known representative, don't try to mess with it
                    if ( $manual_rep_uids{$uid} )         { next }
                    if ( $manual_domain_ids{$domain_id} ) { next }
                    my $rep_node;
                    if ( !$domain->exists('ecod_representative_domain') ) {

                        #print "WARNING $uid no rep domain node and not manual? wtf?\n";
                        $rep_node = $ecod_xml_doc->createElement('ecod_representative_domain');
                        $domain->appendChild($rep_node);
                    }
                    else {
                        $rep_node = $domain->findnodes('ecod_representative_domain')->get_node(1);
                    }
                    my $current_rep_uid;

                    my $apparent_rep_domain_id;
                    if ( $domain->findvalue('ecod_representative_domain/@ecod_domain_id') =~ /(e\w{4}.\d+)\[/ ) {
                        $apparent_rep_domain_id = $1;
                        $rep_node->setAttribute( 'ecod_domain_id', $apparent_rep_domain_id );
                    }

                    if (   $domain->exists('ecod_representative_domain/@uid')
                        && $domain->findvalue('ecod_representative_domain/@uid') =~ /\d+/ )
                    {
                        $current_rep_uid = sprintf "%09i", $domain->findvalue('ecod_representative_domain/@uid');
                    }
                    elsif ($domain->exists('derived_range/@derivedFrom')
                        && $domain->findvalue('derived_range/@derivedFrom') =~ /\d+/ )
                    {
                        $current_rep_uid = sprintf "%09i", $domain->findvalue('derived_range/@derivedFrom');
                    }
                    elsif ( $manual_domain_ids{$apparent_rep_domain_id} ) {    #What does this case solve?
                        $current_rep_uid = $manual_domain_ids{$apparent_rep_domain_id};
                    }
                    else {
                        print "WARNING! No rep recorded for $uid $domain_id\n";
                        next;
                    }
                    if (   $rep_node->findvalue('@ecod_domain_id') !~ /\w+/
                        && $manual_rep_uids{$rep_uid}
                        && $manual_rep_uids{$rep_uid} =~ /\w+/ )
                    {
                        $rep_node->setAttribute( 'ecod_domain_id', $manual_rep_uids{$rep_uid} );
                    }

                    #if (!$current_rep_uid) {
                    #	die "$uid $domain_id\n";
                    #}
                    my $current_rep_domain_id = $domain->findvalue('ecod_representative_domain/@ecod_domain_id');
                    if ( !$manual_rep_uids{$current_rep_uid} && !$manual_domain_ids{$current_rep_domain_id} ) {

                        #print "WARNING! $uid has external rep $rep_uid\n";
                        $external_rep++;
                        $rep_node->setAttribute( 'uid',            $rep_uid );
                        $rep_node->setAttribute( 'ecod_domain_id', $rep_domain_id );
                    }
                    elsif ( !$manual_rep_uids{$rep_uid} && $manual_domain_ids{$domain_id} ) {
                        $rep_node->setAttribute( 'uid', $manual_domain_ids{$domain_id} );
                    }
                    elsif ( $manual_rep_uids{$rep_uid} && !$manual_domain_ids{$domain_id} ) {

                        #print "DEBUG: set rep $rep_uid $manual_rep_uids{$rep_uid}\n";
                        $rep_node->setAttribute( 'ecod_domain_id', $manual_rep_uids{$rep_uid} );
                    }
                    else {
                        die "?";
                    }

                }
            }
            else {
                die "$rep_uid $rep_domain_id unk?\n";
            }
        }
    }
}

sub build_chainwise {
    my $sub = 'build_chainwise';

    #Strip down into pdb, chains, domains, and ranges of those domains.
    #Resolve case-ambiguity in later step.
    my ( $ecod_xml_doc, $USE_ECODF ) = @_;

    my %pdb_chain;
    my %domain_ranges;
    my %domain_pdb_ranges;

    my $domain_nodes = find_domain_nodes($ecod_xml_doc);

    my %ecod_domain_ids;
    my %scop_domain_ids;

    my %recut;
    my %range_type;
    my %domain_fids;
    my %domain_pfids;

    my %side_load;

    my %isAssembly;

    my %prov_man;
    my %obsolete;
    my %ecodf_annot;
    foreach my $d_node ( $domain_nodes->get_nodelist() ) {

        my ( $uid, $ecod_domain_id ) = get_ids($d_node);

        if ( $d_node->exists('scop_domain') ) {
            my $scop_domain_id = $d_node->findvalue('scop_domain/@scop_domain_id');
            $scop_domain_ids{$uid} = $scop_domain_id;
        }

        my $fid;
        my $pfid = 'Unk';
        if ( $d_node->parentNode->nodeName eq 'domain_assembly' ) {
            $isAssembly{$uid}++;
            if ( $d_node->parentNode->parentNode->parentNode->nodeName eq 'f_group' ) {
                $fid  = $d_node->parentNode->parentNode->parentNode->findvalue('@f_id');
                $pfid = $d_node->parentNode->parentNode->findvalue('@pf_id');
            }
            else {
                $fid  = $d_node->parentNode->parentNode->findvalue('@f_id');
                $pfid = $d_node->parentNode->findvalue('@pf_id');
            }
        }
        elsif ( $d_node->parentNode->parentNode->nodeName eq 'f_group' ) {
            $fid  = $d_node->parentNode->parentNode->findvalue('@f_id');
            $pfid = $d_node->parentNode->findvalue('@pf_id');
        }
        else {
            $fid = $d_node->parentNode->findvalue('@f_id');
        }

        my ( $pdb, $chain ) = get_pdb_chain($d_node);

        push( @{ $pdb_chain{$pdb}{$chain} }, $uid );
        $ecod_domain_ids{$uid} = $ecod_domain_id;

        my $seqid_range = get_seqid_range($d_node);
        my $pdb_range   = get_pdb_range($d_node);
        if ( $d_node->findvalue('@derived_range') eq 'true' ) {
            $range_type{$uid} = 'derived';
        }
        elsif ( $d_node->findvalue('@scop_implicit_range') eq 'true' ) {
            $range_type{$uid} = 'implicit';
        }
        elsif ( $d_node->findvalue('@manual_range') eq 'true' ) {
            $range_type{$uid} = 'manual';
        }
        elsif ( $d_node->findvalue('@curated_range') eq 'true' ) {
            $range_type{$uid} = 'curated';
        }
        else {
            print "WARNING! Did not recognize domain type for $ecod_domain_id, skipping...\n";
            next;
        }

        if ( $d_node->findvalue('@provisional_manual_rep') eq 'true' ) {
            $prov_man{$uid}++;
        }

        if ( $d_node->findvalue('@isObsolete') eq 'true' ) {
            $obsolete{$uid}++;
        }

        #Domains derived by "recutting" ASTRAL representative domains
        if ( $d_node->findvalue('@recut_derived') eq 'true' ) {
            $recut{$uid} = 1;
        }

        #Domains loaded "sideways" from yuxing's manual curation table rather than from Hua's
        if ( $d_node->findvalue('@yx_side_load') eq 'true' ) {
            $side_load{$uid} = 1;
        }

        if ($USE_ECODF) {
            if ( $d_node->exists('ecodf_annot') ) {
                foreach my $ecodf_annot ( $d_node->findnodes('ecodf_annot') ) {
                    my $ecodf_acc = $ecodf_annot->findvalue('@ecodf_acc');
                    push @{ $ecodf_annot{$uid} }, $ecodf_acc;
                }
            }
        }

        #print "DEBUG: $ecod_domain_id $pdb $chain $seqid_range\n";
        $domain_ranges{$uid}     = $seqid_range;
        $domain_pdb_ranges{$uid} = $pdb_range;
        $domain_fids{$uid}       = $fid;
        $domain_pfids{$uid}      = $pfid;

    }

    my ( $ecod_chain_xml_doc, $ecod_chain_root ) = xml_create('ecod_chain_doc');

    my $pdb_chain_list_node = $ecod_chain_xml_doc->createElement('pdb_chain_list');
    $ecod_chain_root->appendChild($pdb_chain_list_node);

    foreach my $pdb_id ( keys %pdb_chain ) {
        my $pdb_node = $ecod_chain_xml_doc->createElement('pdb');

        $pdb_node->setAttribute( 'pdb_id', $pdb_id );
        $pdb_chain_list_node->appendChild($pdb_node);

        foreach my $chain_id ( sort { $a cmp $b } keys %{ $pdb_chain{$pdb_id} } ) {
            my $chain_node;
            if ( $pdb_node->exists(qq{chain[\@chain_id="$chain_id"]}) ) {
                $chain_node = $pdb_node->findnodes(qq{chain[\@chain_id="$chain_id"]})->get_node(1);
            }
            else {
                $chain_node = $ecod_chain_xml_doc->createElement('chain');
                $chain_node->setAttribute( 'chain_id', $chain_id );
                $pdb_node->appendChild($chain_node);
            }

            #print "$key1\t$key2:\n";
            foreach my $dom_uid ( @{ $pdb_chain{$pdb_id}{$chain_id} } ) {
                my $domain_node = $ecod_chain_xml_doc->createElement('domain');
                my $dom         = $ecod_domain_ids{$dom_uid};
                $domain_node->setAttribute( 'ecod_domain_id', $dom );

                if ( $scop_domain_ids{$dom_uid} ) {
                    $domain_node->setAttribute( 'scop_domain_id', $scop_domain_ids{$dom_uid} );
                }

                $domain_node->setAttribute( 'uid',        $dom_uid );
                $domain_node->setAttribute( 'range_type', $range_type{$dom_uid} );

                if ( $recut{$dom_uid} ) {
                    $domain_node->setAttribute( 'recut', 'true' );
                }
                if ( $side_load{$dom_uid} ) {
                    $domain_node->setAttribute( 'side_load', 'true' );
                }
                if ( $isAssembly{$dom_uid} ) {
                    $domain_node->setAttribute( 'domain_assembly', 'true' );
                }
                if ( $prov_man{$dom_uid} ) {
                    $domain_node->setAttribute( 'provisional', 'true' );
                }
                if ( $obsolete{$dom_uid} ) {
                    $domain_node->setAttribute( 'obsolete', 'true' );
                }
                if ($USE_ECODF) {
                    if ( $ecodf_annot{$dom_uid} ) {
                        foreach my $ecodf_acc ( @{ $ecodf_annot{$dom_uid} } ) {
                            my $ecodf_annot_node = $ecod_chain_xml_doc->createElement('ecodf_annot');
                            $ecodf_annot_node->setAttribute( 'ecodf_acc', $ecodf_acc );
                            $domain_node->appendChild($ecodf_annot_node);
                        }
                    }
                }
                $domain_node->setAttribute( 'f_id',  $domain_fids{$dom_uid} );
                $domain_node->setAttribute( 'pf_id', $domain_pfids{$dom_uid} );

                #print "\t$dom\t$domain_ranges{$dom_uid}\n";

                $domain_node->appendTextChild( 'seqid_range', $domain_ranges{$dom_uid} );
                $domain_node->appendTextChild( 'pdb_range',   $domain_pdb_ranges{$dom_uid} );

                $chain_node->appendChild($domain_node);

                if ( $chain_id eq '.' ) {

                    my ( $uid, $ecod_domain_id ) = get_ids($domain_node);

                    my $seqid_range = $domain_node->findvalue('seqid_range');
                    my $pdb_range   = $domain_node->findvalue('pdb_range');
                    my @segs        = split( /,/, $seqid_range );
                    my @pdb_segs    = split( /,/, $pdb_range );

                    my $seg_count = 1;

                    #print "e: $ecod_domain_id s: $seqid_range p: $pdb_range\n";

                    for ( my $i = 0 ; $i < scalar(@segs) ; $i++ ) {
                        my $seg = $segs[$i];
                        $seg =~ /((.):.*)/ or die "ERROR! $seg\n";
                        my $seg_chain       = $2;
                        my $seg_seqid_range = $1;
                        my $pdb_seg         = $pdb_segs[$i];
                        $pdb_seg =~ /((.):.*)/ or die "ERROR! $pdb_seg\n";
                        my $seg_pdb_range = $1;
                        if ( $pdb_node->exists(qq{chain[\@chain_id='$seg_chain']}) ) {
                            my $pseudo_chain_node =
                              $pdb_node->findnodes(qq{chain[\@chain_id='$seg_chain']})->get_node(1);
                            if ( $pseudo_chain_node->exists(qq{domain[\@ecod_domain_id='$ecod_domain_id']}) ) {
                                $pseudo_chain_node->removeChild(
                                    $pseudo_chain_node->findnodes(qq{domain[\@ecod_domain_id='$ecod_domain_id']})
                                      ->get_node(1) );
                            }
                            my $seg_domain_node = $ecod_chain_xml_doc->createElement('domain');

                            $seg_domain_node->setAttribute( 'domain_fragment', 'true' );
                            $seg_domain_node->setAttribute( 'range_type',      $range_type{$uid} );
                            $seg_domain_node->setAttribute( 'ecod_domain_id',  $ecod_domain_id );
                            $seg_domain_node->setAttribute( 'uid',             "$uid.$seg_count" );
                            $seg_count++;

                            my $seg_range_node = $ecod_chain_xml_doc->createElement('seqid_range');
                            $seg_range_node->appendTextNode($seg_seqid_range);
                            my $seg_pdb_range_node = $ecod_chain_xml_doc->createElement('pdb_range');
                            $seg_pdb_range_node->appendTextNode($seg_pdb_range);

                            $seg_domain_node->appendChild($seg_range_node);
                            $seg_domain_node->appendChild($seg_pdb_range_node);

                            if ( $recut{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'recut', 'true' );
                            }
                            if ( $side_load{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'side_load', 'true' );
                            }
                            if ( $isAssembly{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'domain_assembly', 'true' );
                            }
                            if ( $prov_man{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'provisional', 'true' );
                            }
                            if ( $obsolete{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'obsolete', 'true' );
                            }

                            if ($USE_ECODF) {
                                if ( $ecodf_annot{$dom_uid} ) {
                                    foreach my $ecodf_acc ( @{ $ecodf_annot{$dom_uid} } ) {
                                        my $ecodf_annot_node = $ecod_chain_xml_doc->createElement('ecodf_annot');
                                        $ecodf_annot_node->setAttribute( 'ecodf_acc', $ecodf_acc );
                                        $domain_node->appendChild($ecodf_annot_node);
                                    }
                                }
                            }
                            $pseudo_chain_node->appendChild($seg_domain_node);
                        }
                        else {
                            my $pseudo_chain_node = $ecod_chain_xml_doc->createElement('chain');
                            $pseudo_chain_node->setAttribute( 'chain_id', $seg_chain );
                            $pdb_node->appendChild($pseudo_chain_node);

                            #my ($seqid_aref, $struct_seqid_aref, $pdbnum_aref) = pdbml_seq_parse($pdb_id, $seg_chain);
                            #$pseudo_chain_node->setAttribute('chain_case_ambiguous', 'true');
                            #my $seqid_range = rangify(@$seqid_aref);
                            #$pseudo_chain_node->appendTextChild('seqid_range', $seqid_range);

                            my $seg_domain_node = $ecod_chain_xml_doc->createElement('domain');
                            $seg_domain_node->setAttribute( 'domain_fragment', 'true' );
                            $seg_domain_node->setAttribute( 'ecod_domain_id',  $ecod_domain_id );
                            $seg_domain_node->setAttribute( 'uid',             "$uid.$seg_count" );
                            $seg_count++;

                            $uid =~ /(\d{9}).\d+/;
                            $seg_domain_node->setAttribute( 'range_type', $range_type{$uid} );

                            $pseudo_chain_node->appendChild($seg_domain_node);
                            my $seg_range_node = $ecod_chain_xml_doc->createElement('seqid_range');
                            $seg_range_node->appendTextNode($seg_seqid_range);
                            my $seg_pdb_range_node = $ecod_chain_xml_doc->createElement('pdb_range');
                            $seg_pdb_range_node->appendTextNode($seg_pdb_range);
                            $seg_domain_node->appendChild($seg_range_node);
                            $seg_domain_node->appendChild($seg_pdb_range_node);

                            if ( $recut{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'recut', 'true' );
                            }
                            if ( $side_load{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'side_load', 'true' );
                            }
                            if ( $isAssembly{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'domain_assembly', 'true' );
                            }
                            if ( $prov_man{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'provisional', 'true' );
                            }
                            if ( $obsolete{$dom_uid} ) {
                                $seg_domain_node->setAttribute( 'obsolete', 'true' );
                            }
                            $pdb_node->appendChild($pseudo_chain_node);
                        }
                    }

                    #nothing for now

                }
            }
        }
    }
    return $ecod_chain_xml_doc;
}

sub find_overlaps {
    my $sub                 = 'find_overlaps';
    my ($chainwise_xml_doc) = @_;
    my $chain_xpath         = '//chain';

    foreach my $chain_node ( $chainwise_xml_doc->findnodes($chain_xpath)->get_nodelist() ) {

        my $pdb_id   = $chain_node->parentNode->findvalue('@pdb_id');
        my $chain_id = $chain_node->findvalue('@chain_id');
        my $f_id     = $chain_node->findvalue('@fid');

        #if ($chain_id eq '.') { next } #No mc

        my %domains;
        my %domain_lookup;
        my %ranges;
        my %range_type;
        my %domain_nodes;
        my %domain_fragment;
        my %side_load;
        my %provisional;

        foreach my $domain_node ( find_domain_nodes($chain_node) ) {

            my ( $domain_uid, $domain_id ) = get_ids($domain_node);

            my $range      = $domain_node->findvalue('seqid_range');
            my $range_type = $domain_node->findvalue('@range_type');
            my $side_load  = 0;
            if ( $domain_node->findvalue('@side_load') eq 'true' ) {
                $side_load = 1;
            }
            $side_load{$domain_uid} = $side_load;

            my $provisional = 0;
            if ( $domain_node->findvalue('@provisional') eq 'true' ) {
                $provisional = 1;

            }
            $provisional{$domain_uid} = $provisional;

            $domains{$domain_id}++;
            $domain_lookup{$domain_uid} = $domain_id;
            if ($range) {
                $ranges{$domain_uid} = $range;
            }
            $range_type{$domain_uid}   = $range_type;
            $domain_nodes{$domain_uid} = $domain_node;
            if ( $domain_node->findvalue('@domain_fragment') eq 'true' ) {
                $domain_fragment{$domain_uid}++;
            }

        }

        my %conflicts;
        my %duplicate_domain_ids;
        foreach my $domain_id ( keys %domains ) {
            if ( $domains{$domain_id} > 1 ) {
                print "WARNING! Duplicate domain_id $domain_id in $pdb_id, $chain_id\n";
                $duplicate_domain_ids{$domain_id}++;

                foreach
                  my $domain_node ( $chain_node->findnodes(qq{domain[\@ecod_domain_id='$domain_id']})->get_nodelist() )
                {
                    $domain_node->setAttribute( 'duplicate', 'true' );
                    my $uid = $domain_node->findvalue('@uid');
                    print "$domain_id $uid\n";

                    #$conflicts{$uid}++;
                }
            }
        }

        my %done;
        foreach my $uid1 ( sort { $a cmp $b } keys %ranges ) {
            my $range1       = $ranges{$uid1};
            my $range_type1  = $range_type{$uid1};
            my $provisional1 = $provisional{$uid1};

            #my ($nc_range1, $chain_str1) = scop_range_split($range1);
            my ( $range1_aref, $chain1_aref ) = multi_chain_range_expand($range1);

            foreach my $uid2 ( sort { $a cmp $b } keys %ranges ) {
                my $range2       = $ranges{$uid2};
                my $range_type2  = $range_type{$uid2};
                my $provisional2 = $provisional{$uid2};

                if ( $uid1 == $uid2 ) { next }
                if ( $done{$uid1}{$uid2} || $done{$uid2}{$uid1} ) { next }

                if ( $range1 eq $range2 ) {
                    print "WARNING! duplicate range $uid1-$uid2\n";
                    $done{$uid1}{$uid2}++;
                    $done{$uid2}{$uid1}++;
                    if ( $range_type1 eq 'manual' && $range_type2 ne 'manual' ) {
                        $conflicts{$uid2}++;
                    }
                    elsif ( $range_type2 eq 'manual' && $range_type1 ne 'manual' ) {
                        $conflicts{$uid1}++;
                    }
                    else {
                        if ( $provisional1 && !$provisional2 ) {
                            $conflicts{$uid2}++;
                        }
                        elsif ( $provisional2 && !$provisional1 ) {
                            $conflicts{$uid1}++;
                        }
                    }

                    if ( $domain_fragment{$uid2} && !$domain_fragment{$uid1} ) {
                        print "FRAG1 $uid1\n";
                        $conflicts{$uid1}++;
                    }
                    elsif ( $domain_fragment{$uid1} && !$domain_fragment{$uid2} ) {
                        print "FRAG2 $uid2\n";
                        $conflicts{$uid2}++;
                    }
                }

                #my ($nc_range2, $chain_str2) = scop_range_split($range2);
                my ( $range2_aref, $chain2_aref ) = multi_chain_range_expand($range2);

                #if ($chain_str1 ne $chain_str2) {
                #	die "Something terrible: $uid1:$chain_str1 $uid2:$chain_str2\n";
                #}

                my $r = multi_chain_residue_coverage( $range1_aref, $chain1_aref, $range2_aref, $chain2_aref );

                if ( $r > 20 ) {
                    print
"WARNING! non-zero coverage between $uid1:$domain_lookup{$uid1} $uid2:$domain_lookup{$uid2} => $r\n";
                    if ( $range_type1 eq 'derived' && $range_type2 eq 'derived' ) {
                        if ( $domain_fragment{$uid2} && !$domain_fragment{$uid1} ) {
                            print "FRAG1 $uid1\n";
                            $conflicts{$uid1}++;
                        }
                        elsif ( $domain_fragment{$uid1} && !$domain_fragment{$uid2} ) {
                            print "FRAG2 $uid2\n";
                            $conflicts{$uid2}++;
                        }
                        else {
                            if ( $uid1 > $uid2 ) {
                                $conflicts{$uid2}++;
                            }
                            else {
                                $conflicts{$uid1}++;
                            }
                        }
                    }
                    else {
                        $conflicts{$uid1}++;
                        $conflicts{$uid2}++;
                    }
                }

                $done{$uid1}{$uid2}++;
                $done{$uid2}{$uid1}++;
            }
        }

        foreach my $uid ( keys %conflicts ) {
            print "CONFLICT_SUMM $f_id: $uid $conflicts{$uid} $range_type{$uid}\n";
            if (   $conflicts{$uid} > 0
                && $range_type{$uid} ne "manual"
                && $range_type{$uid} ne "implicit"
                && $range_type{$uid} ne 'curated' )
            {
                $domain_nodes{$uid}->setAttribute( 'delete', 'true' );
            }
            elsif ( $conflicts{$uid} > 0 && $side_load{$uid} ) {
                $domain_nodes{$uid}->setAttribute( 'delete', 'true' );
            }

        }

    }
}

#Search a chainwise XML structure for duplicate domain nodes, remove them from ECOD xml structure {
sub resolve_duplicate_domains {
    my ( $chainwise_xml_doc, $ecod_xml_doc ) = @_;
    my $sub = 'resolve_duplicate_domains';

    my %domain_nodes;
    my %domain_rep_uid;
    my %domain_children;
    foreach my $domain_node ( $ecod_xml_doc->findnodes('//domain') ) {
        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
        my $uid            = $domain_node->findvalue('@uid');
        $domain_nodes{$uid} = $domain_node;

        if ( $domain_node->exists('ecod_representative_domain') ) {
            my $rep_uid = $domain_node->findvalue('ecod_representative_domain/@uid');
            $domain_rep_uid{$uid} = $rep_uid;
            push( @{ $domain_children{$rep_uid} }, $uid );
        }
    }

    foreach my $chain_node ( $chainwise_xml_doc->findnodes('//chain') ) {

        if ( $chain_node->exists('./domain[@duplicate="true"]') ) {
            my %domain_hash;

            my $pdb_id   = $chain_node->findvalue('@pdb_id');
            my $chain_id = $chain_node->findvalue('@chain_id');

            #hash names;
            my %domain_ids;
            my %uid2domain_id;
            foreach my $domain_node ( $chain_node->findnodes('domain') ) {
                my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
                my $uid            = $domain_node->findvalue('@uid');
                $uid2domain_id{$uid} = $ecod_domain_id;
                $domain_ids{$ecod_domain_id}++;
            }

            my $duplicate_seen = 0;
            foreach my $uid ( keys %uid2domain_id ) {
                if ( $domain_ids{ $uid2domain_id{$uid} } > 1 ) {
                    my $old_domain_id = $uid2domain_id{$uid};
                    if ($duplicate_seen) {
                        my $num           = 1;
                        my $new_domain_id = "e$pdb_id$chain_id$num";
                        while ( $domain_ids{$$new_domain_id} ) {
                            $num++;
                            $new_domain_id = "e$pdb_id$chain_id$num";
                        }
                        print "$sub: Changing $old_domain_id to $new_domain_id\n";
                        $domain_nodes{$uid}->setAttribute( 'ecod_domain_id', $new_domain_id );
                        $domain_nodes{$uid}->findnodes('ecod_domain')->get_node(1)
                          ->setAttribute( 'ecod_domain_id', $new_domain_id );
                        $domain_nodes{$uid}->findnodes('ecod_domain')->get_node(1)->setAttribute( 'relabel', 'true' )
                          ;    #This is temporary

                    }
                    else {
                        $duplicate_seen++;
                    }
                }
            }
        }

    }
}

#Search a chainwise XMl structure for deleted domain nodes, remove them from ecod xml structure
sub delete_chainwise_domains {
    my ( $chainwise_xml_doc, $ecod_xml_doc ) = @_;
    my $sub = 'delete_chainwise_domains';

    my %domain_nodes;
    my %domain_derivedFrom;
    my %domain_children;

    foreach my $domain_node ( $ecod_xml_doc->findnodes('//domain')->get_nodelist() ) {

        #if ($domain_node->findvalue('@manual_rep') eq 'true') { next }
        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
        my $uid            = $domain_node->findvalue('@uid');
        $domain_nodes{$uid} = $domain_node;

        if ( $domain_node->exists('ecod_representative_domain') ) {
            my $derivedFrom = $domain_node->findvalue('ecod_representative_domain/@uid');
            $domain_derivedFrom{$uid} = $derivedFrom;
            push( @{ $domain_children{$derivedFrom} }, $uid );
        }
    }

    foreach my $chain_domain_node ( $chainwise_xml_doc->findnodes('//domain')->get_nodelist() ) {

        my $ecod_domain_id = $chain_domain_node->findvalue('@ecod_domain_id');
        my $uid            = $chain_domain_node->findvalue('@uid');

        if (
               $chain_domain_node->parentNode->findvalue('@delete') eq 'true'
            || $chain_domain_node->findvalue('@delete') eq 'true'
            || $chain_domain_node->findvalue('@obsolete') eq 'true'
            || (   $chain_domain_node->parentNode->findvalue('@purge_derived') eq 'true'
                && $chain_domain_node->findvalue('@range_type') eq 'derived' )
          )
        {
            #if ($chain_domain_node->findvalue('@domain_assembly') eq 'true') {
            #	print "WARNING! No deleting out of domain assemblies automatically, skipping $ecod_domain_id $uid\n";
            #	next;
            #}
            if ($DEBUG) {
                print "DELETE $ecod_domain_id $uid\n";
            }
            if ( $domain_nodes{$uid} ) {
                my $domain_node = $domain_nodes{$uid};
                $domain_node->unbindNode;

                if ( $domain_children{$uid} ) {
                    if ($DEBUG) {
                        printf "%i children of $ecod_domain_id $uid\n", scalar( @{ $domain_children{$uid} } );
                    }
                    foreach my $child_uid ( @{ $domain_children{$uid} } ) {
                        my $child_domain_node = $domain_nodes{$child_uid};
                        $child_domain_node->unbindNode;
                    }
                }
            }
        }
    }
    return $ecod_xml_doc;
}

sub generate_domain_seqinput_file {
    my $sub = 'generate_domain_seqinput_file';

    my ( $ecod_xml_doc, $out_fn ) = @_;

    my $SCOUNT   = 0;
    my $COUNT    = 0;
    my $PC_COUNT = 0;

    my $FILE_TEST = 1;

    my $DOMAIN_DATA_DIR = '/data/ecod/domain_data/';

    open( OUT, ">$out_fn" ) or die "Could not open $out_fn for writing:$!\n";
    my %pc_seen;
  DOMAIN:
    foreach my $d_node ( find_domain_nodes($ecod_xml_doc) ) {
        my $domain_uid = $d_node->findvalue('@uid');
        my $short_uid = substr( $domain_uid, 2, 5 );
        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid");
        }
        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/$domain_uid" ) {
            mkdir("$DOMAIN_DATA_DIR/$short_uid/$domain_uid");
        }

        if (   $FILE_TEST
            && -f "$DOMAIN_DATA_DIR/$short_uid/$domain_uid/$domain_uid.fa"
            && -s "$DOMAIN_DATA_DIR/$short_uid/$domain_uid/$domain_uid.fa" > 0 )
        {
            next;
        }

        if ($DEBUG) {
            print "\t$domain_uid\n";
        }

        my $seqid_range;

        if ( $d_node->findvalue('@derived_range') eq 'true' ) {
            $seqid_range = $d_node->findvalue('derived_seqid_range');
        }
        elsif ( $d_node->findvalue('@scop_implicit_range') eq 'true' ) {
            $seqid_range = $d_node->findvalue('seqid_range');
        }
        elsif ( $d_node->findvalue('@manual_range') eq 'true' ) {
            $seqid_range = $d_node->findvalue('seqid_range');
        }
        elsif ( $d_node->findvalue('@curated_range') eq 'true' ) {
            $seqid_range = $d_node->findvalue('seqid_range');
        }
        else {
            print "WARNING! unknown domain type for $domain_uid\n";
            $COUNT++;
            next;
        }

        my ( $pdb, $chain ) = get_pdb_chain($d_node);
        if ( $d_node->exists('structure') ) {
            my $pdb_chain = $pdb . "_" . $chain;
            if ( !$pc_seen{$pdb_chain} ) {
                $PC_COUNT++;
                $pc_seen{$pdb_chain}++;
            }
        }
        else {
            print "WARNING! No pdb for $domain_uid\n";
            $COUNT++;
            next;
        }
        my $ecod_domain_id = $d_node->findvalue('@ecod_domain_id');

        print OUT "$ecod_domain_id\t$seqid_range\t$domain_uid\t$pdb\n";
        $COUNT++;

    }
    if ($DEBUG) {
        print "DEBUG $sub: $COUNT $PC_COUNT\n";
    }
    close OUT;
    return ( $COUNT, $PC_COUNT );
}

sub generate_domain_range_cache_file {
    my $sub = 'generate_domain_range_cache_file';

    my ( $ecod_xml_doc, $out_fn ) = @_;

    open( OUT, ">$out_fn" ) or die "ERROR! Could not open $out_fn for writing:$!\n";

    my $xpath = '//domain';

    #my $manual_xpath = '//domain[@manual_rep="true"]';
    foreach my $manual_domain ( $ecod_xml_doc->findnodes($xpath)->get_nodelist() ) {

        my $seqid_range;
        if ( $manual_domain->exists('seqid_range') ) {
            $seqid_range = $manual_domain->findvalue('seqid_range');
        }
        else {
            $seqid_range = $manual_domain->findvalue('derived_seqid_range');
        }
        my $ecod_domain_id = $manual_domain->findvalue('@ecod_domain_id');
        my $uid            = $manual_domain->findvalue('@uid');

        my $pdb   = $manual_domain->findvalue('structure/@pdb_id');
        my $chain = $manual_domain->findvalue('structure/@chain_id');

        print OUT "$uid\t$ecod_domain_id\t$seqid_range\t$pdb\t$chain\n";
    }
    close OUT;

}

sub generate_domain_fasta_jobs {
    my $sub = 'generate_domain_fasta_jobs';

    my ( $seq_input_fn, $dir, $N ) = @_;

    print "DEBUG: $sub: $seq_input_fn $dir $N\n";
    my @jobs;
    $INC = 10;
    my $max = int( $N / $INC );
    for my $i ( 0 .. ( $max + 1 ) ) {

        open( OUT, ">$dir/para.$i.job" ) or die "ERROR! $sub: Could not open para.$i.job for writing:$!\n";

        print OUT "#!/bin/bash\n";
        print OUT "#\$ -cwd\n";
        print OUT "#\$ -j y \n";
        print OUT "#\$ -S /bin/bash\n";
        print OUT "#\$ -M dustin.schaeffer\@gmail.com\n";

        #print OUT "generate_seqdb_from_ecod_para.pl ecod.twig.xml $i\n";
        print OUT "$GENERATE_SEQDB_EXE $dir/$seq_input_fn $i\n";

        push( @jobs, "$dir/para.$i.job" );

    }
    return \@jobs;

}

sub generate_chainwise_fasta_jobs {
    my $sub = 'generate_chainwise_fasta_jobs';

    my ( $chainwise_xml_fn, $dir, $N ) = @_;

    my @jobs;
    for my $i ( 0 .. $N ) {

        open( OUT, ">$dir/chpara.$i.job" ) or die "ERROR! $sub: Could not open chpara.$i.job for writing:$!\n";

        print OUT "#!/bin/bash\n";
        print OUT "#\$ -cwd\n";
        print OUT "#\$ -j y \n";
        print OUT "#\$ -S /bin/bash\n";
        print OUT "#\$ -M dustin.schaeffer\@gmail.com\n";
        print OUT "$GENERATE_CHAINWISE_SEQDB_EXE $dir/$chainwise_xml_fn $i\n";

        push( @jobs, "$dir/chpara.$i.job" );

    }
    return \@jobs;
}

sub generate_hhblits_profile_jobs {
    my $sub = 'generate_hhblits_profile_jobs';

    my ($new_ecod_xml_doc) = @_;

    my %new_manual_reps;
    my %new_hhm;
    foreach my $domain_node ( find_manual_domain_nodes($new_ecod_xml_doc) ) {
        my ( $domain_uid, $domain_id ) = get_ids($domain_node);
        my $short_uid = substr( $domain_uid, 2, 5 );
        $new_manual_reps{$domain_uid} = $domain_node;
        if ( -f "$DOMAIN_DATA_DIR/$short_uid/$domain_uid/$domain_uid.hhm" ) {
            $new_hhm{$domain_uid}++;
        }
    }

  DOMAIN:
    foreach my $new_uid ( keys %new_manual_reps ) {
        my @jobs;
        my $ecod_domain_id = get_ecod_domain_id( $new_manual_reps{$new_uid} );
        my $short_uid = substr( $new_uid, 2, 5 );

        my $uid_path = "$DOMAIN_DATA_DIR/$short_uid/$new_uid";
        my $hhm_path = "$uid_path/$new_uid.hhm";
        my $fa_path  = "$uid_path/$new_uid.fa";
        next if -f $hhm_path;

        my $d_node = $new_manual_reps{$new_uid};

        my $seqid_range = get_seqid_range($d_node);
        my ( $pdb, undef ) = get_pdb_chain($d_node);

        my @chains;
        my @segs;

        if ( $seqid_range =~ /\w:\-?\d+\-\d+/ ) {
            while ( $seqid_range =~ /(\w):(\-?\d+\-\d+)/g ) {
                my $chain = $1;
                my $seg   = $2;
                push( @chains, $chain );
                push( @segs,   $seg );
            }
        }
        else {
            print "WARNING! range regexp failure on $seqid_range\n";
            next DOMAIN;
        }

        if ( !-f $fa_path ) {
            warn "WARNING! This fasta creation for $new_uid is deprecated, please generate it elsewhere!\n";
            my $fasta_string;
            if ( !-d "$DOMAIN_DATA_DIR/$short_uid" ) {
                mkdir("$DOMAIN_DATA_DIR/$short_uid");
            }
            chown( $UID, $GID, "$DOMAIN_DATA_DIR/$short_uid" );
            if ( !-d $uid_path ) {
                mkdir($uid_path);
            }
            chown( $UID, $GID, $uid_path );

            for ( my $i = 0 ; $i < scalar(@chains) ; $i++ ) {
                my $chain = $chains[$i];
                my ( $seqid_aref, $struct_seqid_aref, $pdbnum_aref, $asym_id ) = pdbml_seq_parse( $pdb, $chain );
                if ( !$seqid_aref ) {
                    print "WARNING! $pdb, $chain did not seq parse\n";
                    next DOMAIN;
                }
                my $range = range_expand( $segs[$i] );
                if ( scalar(@$range) == 0 ) {
                    print "WARNING $pdb, $chain, $segs[$i] is empty\n";
                    next DOMAIN;
                }
                my $fasta_seg = pdbml_fasta_fetch( $pdb, $asym_id, $chain, $range );
                $fasta_string .= $fasta_seg;
            }

            open( OUT, ">", $fa_path ) or die "ERROR! Could not open $fa_path for writing:$!\n";
            print OUT ">$ecod_domain_id $seqid_range $new_uid\n$fasta_string\n";
            close OUT;
        }

        #HHBLITS it
        #my $hh_lib = "export HHLIB=$HHLIB\n"; #Does this even work?
        #push (@jobs, $hhlib);

        my $a3m_path = "$uid_path/$new_uid.a3m";
        if ( -f $fa_path & !-f $a3m_path ) {
            my $ln =
"$HHBLITS_EXE -i $fa_path -d $NR_LIB -oa3m $a3m_path -ohhm $hhm_path -n 3 -cpu $HH_BLITS_CPU -addss -psipred $PSIPRED_DIR -psipred_data $PSIPRED_DATA_DIR\n";
            push( @jobs, $ln );

        }
        if ( scalar @jobs > 0 ) {
            job_create( "q.$new_uid.hhblits.job", \@jobs );
            qsub("-pe orte 8 q.$new_uid.hhblits.job");
        }
    }
}

sub generate_hmmer_jobs {
    my $sub = 'generate_hmmer_jobs';

    my ($ecod_xml_doc) = @_;

    my @jobs;
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {

        my $domain_uid = get_uid($domain_node);
        my $short_uid = substr( $domain_uid, 2, 5 );

        my @job_lns;
        if ( !-f "$DOMAIN_DATA_DIR/$short_uid/$domain_uid/$domain_uid.ecodf.hmm_results" ) {
            my $job_ln = sprintf
"$HMMSCAN_EXE --acc --noali -o $DOMAIN_DATA_DIR/$short_uid/$domain_uid/$domain_uid.ecodf.hmm_results $HMM_LIB $DOMAIN_DATA_DIR/$short_uid/$domain_uid/$domain_uid.fa\n";
            push( @job_lns, $job_ln );

            job_create( "hmmer.$domain_uid.job", \@job_lns );
            push( @jobs, "hmmer.$domain_uid.job" );
        }
    }
    return \@jobs;
}

sub is_rep_only {
    if ( $_[0]->findvalue('@manual_rep') eq 'true' ) {
        return 1;
    }
    elsif ( $_[0]->findvalue('@provsional_manual_rep') eq 'true' ) {
        if ( $_[0]->findvalue('@curated_range') eq 'true' ) {
            return 1;
        }
        elsif ( $_[0]->findvalue('@canonical_by_curation') eq 'true' ) {
            return 1;
        }
    }
}

#assemble_domain_fasta_library - given an ECOD XML, domain data direcotry, and a fasta output, dump either all-domains or reps_only domains to a file
sub assemble_domain_fasta_library {
    my $sub = 'assemble_domain_fasta_library';

    my ( $ecod_xml_doc, $fa_db_fn, $reps_only ) = @_;

    if ( !-d $DOMAIN_DATA_DIR ) { die "ERROR! Could not find domain data dir $DOMAIN_DATA_DIR\n"; }

    my @global_fa;
    open( OUT, ">$fa_db_fn" ) or die "ERROR $sub: Could not open $fa_db_fn for writing:$!\n";

    my $i = 0;
  DOMAIN:
    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {

        if ( !is_rep_only($domain_node) && $reps_only ) { next }

        my $uid = get_uid($domain_node);

        my $short_uid = substr( $uid, 2, 5 );
        my $path = "$DOMAIN_DATA_DIR/$short_uid/$uid";

        if ( !-f "$path/$uid.fa" ) { print "WARNING! $path/$uid.fa not found\n"; next }
        open( IN, "$path/$uid.fa" ) or die "ERROR! $path/$uid.fa file not found\n";

        my @fa_lns;
        while ( my $ln = <IN> ) {
            if ( $ln =~ /^0/ ) { next DOMAIN }    #Skip empty sequences from obsolete PDB
            push( @fa_lns, $ln );
        }
        close IN;

        foreach my $ln (@fa_lns) {
            print OUT $ln;
        }
        $i++;

    }
    close OUT;
}

sub assemble_chainwise_fasta_library {
    my $sub = 'assemble_chainwise_fasta_library';
    my ( $chainwise_xml_doc, $fa_db_fn, $reps_only ) = @_;

    open( OUT, ">$fa_db_fn" ) or die "ERROR $sub! Could not open $fa_db_fn for writing:$!\n";
  CHAIN:
    foreach my $chain_node ( $chainwise_xml_doc->findnodes('//chain')->get_nodelist() ) {

        my $domain_count = $chain_node->findnodes('domain')->size();
        my $rep_count    = $chain_node->findnodes('domain[@range_type="manual"]')->size();
        $rep_count += $chain_node->findnodes('domain[@range_type="implicit"]')->size();

        #For chains, if a chain is not completely built of reps (range_type{manual,implicit}) skip
        if ( $rep_count != $domain_count && $reps_only ) { next }

        my $pdb   = $chain_node->parentNode->findvalue('@pdb_id');
        my $two   = substr( $pdb, 1, 2 );
        my $chain = $chain_node->findvalue('@chain_id');
        my $pc    = $pdb . "_" . $chain;

        #Check for presence of chain_data directory
        if ( !-d "$CHAIN_DATA_DIR/$two/$pc" ) {
            print "WARNING! $sub: dir $CHAIN_DATA_DIR/$two/$pc not found\n";
            next;
        }
        if ( !-f "$CHAIN_DATA_DIR/$two/$pc/$pc.fa" ) {
            print "WARNING! $sub: dir $CHAIN_DATA_DIR/$two/$pc/$pc.fa not found\n";
            next;
        }
        open( IN, "$CHAIN_DATA_DIR/$two/$pc/$pc.fa" )
          or die "ERROR! $sub:  $CHAIN_DATA_DIR/$two/$pc/$pc.fa file not found\n";

        my @fa_lns;
        while ( my $ln = <IN> ) {
            if ( $ln =~ /^0/ ) { next CHAIN }    #This is hacky fix for 0-residue chains
            push( @fa_lns, $ln );
        }

        foreach my $ln (@fa_lns) {
            print OUT $ln;
        }
        close IN;

    }
    close OUT

}

sub assemble_hh_profile_library {
    my $sub = 'assemble_hh_profile_library';

    my ( $ecod_xml_doc, $hhm_db_fn, $rep_flag ) = @_;

    !$rep_flag and $rep_flag = 'manual';

    my @global_fa;

    my %nodes = (
        manual => \&find_manual_domain_nodes,
        F40    => \&find_domain_reps,
        F70    => \&find_domain_reps,
        F99    => \&find_domain_reps,
    );

    open my $out_fh, ">", $hhm_db_fn or die "ERROR! $sub: Could not open $hhm_db_fn for writing $!\n";
	
	my $method = "blastclust";
    foreach my $domain_node ( $nodes{$rep_flag}->( $ecod_xml_doc, $rep_flag, $method ) ) {

        my $uid = $domain_node->findvalue('@uid');
        my $short_uid = substr( $uid, 2, 5 );

        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/$uid" ) {
            print "WARNING! dir $DOMAIN_DATA_DIR/$short_uid/$uid not found\n";
            next;
        }

        my $suffix = 'hhm';
        my $file   = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.$suffix";
        if ( !-f $file ) {
            print "WARNING! hmm results $file not found\n";
            next;
        }
        if ( !-r $file ) {
            print "WARNING! hmm results $file not readable (perms)\n";
            next;
        }

        open my $fh, "<", $file or die "ERROR! $file file not found:$!\n";

        while ( my $ln = <$fh> ) {
            print $out_fh $ln;
        }
        close $fh;

    }
    close $out_fh;
}

sub assemble_hmmer_profile_library {
    my $sub = "assemble_hmmer_profile_library";

    my ( $ecod_xml_doc, $hmmer_db_fn, $pfam_type ) = @_;
    my @global_fa;
    my $suffix;
    if ( !$pfam_type ) {    #Default is ECODf
        $pfam_type = 'ECODf';
        $suffix    = 'ecodf.hmm_results';
    }
    elsif ( $pfam_type eq 'PfamA' ) {
        $pfam_type = 'PfamA';
        $suffix    = 'fa.hmm_results';
    }
    elsif ( $pfam_type eq 'PfamB' ) {
        $pfam_type = 'PfamB';
        $suffix    = 'fa.pfam_b_hmm_results';
    }
    elsif ( $pfam_type eq 'ECODf' ) {
        $pfam_type = 'ECODf';
        $suffix    = 'ecodf.hmm_results';
    }
    else {
        die "ERROR! $sub: unknown pfam_type $pfam_type\n";
    }

    open( OUT, ">$hmmer_db_fn" ) or die "ERROR! $sub: Could not open $hmmer_db_fn for writing:$!\n";
    my $XPath = '//f_group/domain';

    #my $XPath = '//pf_group[@hh_cluster="true"]/domain';
    #my $XPath = '//domain';
    foreach my $domain_node ( $ecod_xml_doc->findnodes($XPath)->get_nodelist() ) {

        my $uid = $domain_node->findvalue('@uid');

        #print "$uid\n";

        my $short_uid = substr( $uid, 2, 5 );

        if ( !-d "$DOMAIN_DATA_DIR/$short_uid/$uid" ) {
            print "WARNING! dir $DOMAIN_DATA_DIR/$short_uid/$uid not found\n";
            next;
        }
        if ( !-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.fa" ) {
            print "WARNING! No fasta for $uid\n";
            next;
        }
        if ( !-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.$suffix" ) {
            print "WARNING! hmm results $DOMAIN_DATA_DIR/$short_uid/$uid/$uid.$suffix not found\n";
            next;
        }

        open( IN, "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.$suffix" )
          or die "ERROR! $DOMAIN_DATA_DIR/$short_uid/$uid/$uid.$suffix file not found\n";

        while ( my $ln = <IN> ) {
            print OUT $ln;
        }
        close IN;

    }
    close OUT;
}

sub convert_hmmer_txt_to_xml {
    my $sub = 'convert_hmmer_txt_to_xml';

    my ( $txt_input_fn, $xml_output_fn, $ecodf_fn ) = @_;

    my $ecodf_xml_doc;
    if ( !$ecodf_fn ) {
        $ecodf_fn = $ECODF_XML_FN;
    }
    $ecodf_xml_doc = xml_open($ecodf_fn);

    my %ecodf_acc2ecodf_name;
    my %ecodf_acc2pfam_acc;

    foreach my $family_node ( $ecodf_xml_doc->findnodes('//family') ) {
        my $ecodf_acc  = get_short_ecodf_acc($family_node);
        my $pfam_acc   = $family_node->findvalue('@pfam_acc');
        my $ecodf_name = $family_node->findvalue('@ecodf_id');

        $ecodf_acc2ecodf_name{$ecodf_acc} = $ecodf_name;
        $ecodf_acc2pfam_acc{$ecodf_acc}   = $pfam_acc;

    }

    open( IN, $txt_input_fn ) or die "ERROR! Could not open $txt_input_fn for reading:$!\n";

    my ( $hmm_xml_doc, $hmm_root_node ) = xml_create('hmm_parse_document');

    my $hmm_domain_list_node = $hmm_xml_doc->createElement('hmm_domain_list');
    $hmm_root_node->appendChild($hmm_domain_list_node);

    my ( $query, $length, $range, $uid );
    my ($query_node);
    my $hmm_type = 'ecodf';
    while ( my $ln = <IN> ) {

        $ln =~ /^#/ and next;

        if ( $ln =~ /Query:\s+([\w\.]+)\s+\[L=(\d+)\]/ ) {
            $query  = $1;
            $length = $2;

            $query_node = $hmm_xml_doc->createElement('hmm_query');
            $query_node->setAttribute( 'ecod_domain_id', $query );
            $query_node->setAttribute( 'length',         $length );

        }
        if ( $ln =~ /Description:\s+(.+:.*?)\s+(\d+)/ ) {
            $range = $1;
            $uid   = $2;
            $query_node->setAttribute( 'range', $range );
            $query_node->setAttribute( 'uid',   $uid );
        }

        if ( $ln =~ /^Scores for complete sequence/ ) {
            while ( my $ln = <IN> ) {
                if ( $ln =~ /No hits detected/ )    { last }
                if ( $ln =~ /inclusion threshold/ ) { last }
                my $sci_regexp = qr/[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?/;
                if ( $ln =~
/^\s+($sci_regexp)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+($sci_regexp)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d)\s+(\d+)\s+([\w\-\_]+)/
                  )
                {
                    my $full_eval  = $1;
                    my $full_score = $2;
                    my $full_bias  = $3;

                    my $best_eval  = $4;
                    my $best_score = $5;
                    my $best_bias  = $6;

                    my $exp   = $7;
                    my $N     = $8;
                    my $model = $9;

                  #print "$full_eval\t$full_score\t$full_bias\t$best_eval\t$best_score\t$best_bias\t$exp\t$N\t$model\n";
                  #exit();

                    my $hmm_summ_node = $hmm_xml_doc->createElement('hmm_summ');

                    $hmm_summ_node->setAttribute( 'full_eval',  $full_eval );
                    $hmm_summ_node->setAttribute( 'full_score', $full_score );
                    $hmm_summ_node->setAttribute( 'full_bias',  $full_bias );

                    $hmm_summ_node->setAttribute( 'best_eval',  $best_eval );
                    $hmm_summ_node->setAttribute( 'best_score', $best_score );
                    $hmm_summ_node->setAttribute( 'best_bias',  $best_bias );

                    if ( $model =~ /(EF\d+)/ ) {
                        my $ecodf_acc = $1;
                        if ( $ecodf_acc2ecodf_name{$ecodf_acc} ) {
                            my $ecodf_name = $ecodf_acc2ecodf_name{$ecodf_acc};
                            my $pfam_acc   = $ecodf_acc2pfam_acc{$ecodf_acc};

                            $hmm_summ_node->setAttribute( 'ecodf_acc',  $ecodf_acc );
                            $hmm_summ_node->setAttribute( 'ecodf_name', $ecodf_name );
                            if ($pfam_acc) {
                                $hmm_summ_node->setAttribute( 'pfam_acc', $pfam_acc );
                            }
                        }
                        else {
                            warn "WARNING! ecodf_acc $ecodf_acc not in $ecodf_fn, reclassify $query\n";
                            next;
                        }
                    }

                    $query_node->appendChild($hmm_summ_node);

                }
                if ( $ln =~ /^\s+$/ ) { last }
            }
        }

        if ( $ln =~ /^Domain annotation for each model/ ) {
            my ( $model, $desc );
            while ( my $ln = <IN> ) {
                if ( $ln =~ /No targets/ ) { last }
                if ( $ln =~ /^\>\> ([\w\.]+) (.*)/ ) {
                    $model = $1;
                    $desc  = $2;
                    my $junk_ln = <IN>;
                    if ( $junk_ln =~ /No individual/ ) { last }
                    $junk_ln = <IN>;
                    while ( my $good_ln = <IN> ) {
                        if ( $good_ln =~ /^\s+$/ ) { last }
                        $good_ln =~ s/^\s+//;
                        my @F       = split( /\s+/, $good_ln );
                        my $n       = $F[0];
                        my $ex      = $F[1];
                        my $score   = $F[2];
                        my $bias    = $F[3];
                        my $ceval   = $F[4];
                        my $ieval   = $F[5];
                        my $hmmfrom = $F[6];
                        my $hmmto   = $F[7];
                        my $hmmenc  = $F[8];
                        my $alifrom = $F[9];
                        my $alito   = $F[10];
                        my $alienc  = $F[11];
                        my $envfrom = $F[12];
                        my $envto   = $F[13];
                        my $envenc  = $F[14];
                        my $acc     = $F[15];

                        my $hmm_domannot_node = $hmm_xml_doc->createElement('hmm_domannot');

                        $hmm_domannot_node->setAttribute( 'acc',      $model );
                        $hmm_domannot_node->setAttribute( 'desc',     $desc );
                        $hmm_domannot_node->setAttribute( 'N',        $n );
                        $hmm_domannot_node->setAttribute( 'c_eval',   $ceval );
                        $hmm_domannot_node->setAttribute( 'i_eval',   $ieval );
                        $hmm_domannot_node->setAttribute( 'hmm_from', $hmmfrom );
                        $hmm_domannot_node->setAttribute( 'hmm_to',   $hmmto );
                        $hmm_domannot_node->setAttribute( 'ali_from', $alifrom );
                        $hmm_domannot_node->setAttribute( 'ali_to',   $alito );
                        $hmm_domannot_node->setAttribute( 'env_from', $envfrom );
                        $hmm_domannot_node->setAttribute( 'env_to',   $envto );

                        $query_node->appendChild($hmm_domannot_node);

                    }

                }

                if ( $ln =~ /Internal pipeline statistics summary/ ) { last }
            }
        }

        if ( $ln =~ /\/\// ) {    #all done
                                  #print "$query\t$length\t$range\t$uid\n";
            $hmm_domain_list_node->appendChild($query_node);
        }
    }
    xml_write( $hmm_xml_doc, $xml_output_fn );

}

sub cluster_pfgroups_by_ecodf {
    my $sub = 'cluster_pf_groups_by_ecodf';
    my ($ecod_xml_doc) = @_;

    foreach my $f_group ( $ecod_xml_doc->findnodes('//f_group')->get_nodelist() ) {

        my $f_id = $f_group->findvalue('@f_id');

        my %hmmer_ecodf;
        my %hmmer_ecodf_name;
        my $domain_nodes          = $f_group->findnodes('domain');
        my $domain_assembly_nodes = $f_group->findnodes('domain_assembly');

        foreach my $domain_node ( $f_group->findnodes('domain')->get_nodelist() ) {
            my $ecodf_string;
            if ( $domain_node->findnodes('ecodf_annot')->size() == 1 ) {
                my $ecodf_acc  = $domain_node->findvalue('ecodf_annot/@ecodf_acc');
                my $ecodf_name = $domain_node->findvalue('ecodf_annot/@ecodf_name');
                $ecodf_string = $ecodf_acc;

                #print "ecf_str $ecodf_string\n";
                $hmmer_ecodf{$ecodf_string}++;
                $hmmer_ecodf_name{$ecodf_string} = $ecodf_name;
            }
            elsif ( $domain_node->findnodes('ecodf_annot')->size() > 1 ) {
                my @ecodfs;
                my %ecodf_names;
                foreach my $ecodf_node ( $domain_node->findnodes('ecodf_annot')->get_nodelist() ) {
                    my $ecodf_acc = $ecodf_node->findvalue('@ecodf_acc');
                    push( @ecodfs, $ecodf_acc );
                    my $ecodf_name = $ecodf_node->findvalue('@ecodf_name');
                    $ecodf_names{$ecodf_acc} = $ecodf_name;
                }
                @ecodfs = sort { $a cmp $b } @ecodfs;
                $ecodf_string = join( ",", @ecodfs );

                #print "ecf_str_cat $ecodf_string\n";
                my @ecodf_names;
                foreach my $ecodf_acc (@ecodfs) {
                    push( @ecodf_names, $ecodf_names{$ecodf_acc} );
                }
                my $ecodf_name_string = join( ",", @ecodf_names );
                $hmmer_ecodf{$ecodf_string}++;
                $hmmer_ecodf_name{$ecodf_string} = $ecodf_name_string;
            }

        }

        foreach my $domain_assembly_node ( $f_group->findnodes('domain_assembly') ) {
            my $ecodf_string;
            my %assembly_ecodfs;
            my %assembly_ecodf_names;

            foreach my $domain_node ( $domain_assembly_node->findnodes('domain') ) {
                my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');

                #print "DA? $ecod_domain_id\n";
                if ( $domain_node->findnodes('ecodf_annot')->size() == 1 ) {
                    my $ecodf_acc  = $domain_node->findvalue('ecodf_annot/@ecodf_acc');
                    my $ecodf_name = $domain_node->findvalue('ecodf_annot/@ecodf_name');
                    $assembly_ecodfs{$ecodf_acc}++;
                    $assembly_ecodf_names{$ecodf_acc} = $ecodf_name;
                }
            }
            if ( keys %assembly_ecodfs == 1 ) {
                my @strings = sort { $a cmp $b } keys %assembly_ecodfs;

                #my @names = sort {$a cmp $b} keys %assembly_ecodf_names;
                my @names;
                foreach my $ecodf_acc (@strings) {
                    push( @names, $assembly_ecodf_names{$ecodf_acc} );
                }
                my $ecodf_string = $strings[0];
                my $ecodf_name   = $names[0];

                #print "ecf_str_da $ecodf_string\n";
                #print "DA1: $ecodf_string $ecodf_name\n";
                $hmmer_ecodf{$ecodf_string}++;
                $hmmer_ecodf_name{$ecodf_string} = $ecodf_name;
            }
        }

        #my $pf_iter = 1;
        my $pf_iter = 1;
        if ( $f_group->exists('pf_group') ) {
            foreach my $pf_group ( $f_group->findnodes('pf_group') ) {
                my $pf_id = $pf_group->findvalue('@pf_id');
                if ( $pf_id =~ /\d+\.\d+\.\d+\.(\d+)/ ) {
                    my $local_pf_iter = $1;
                    if ( $local_pf_iter >= $pf_iter ) {
                        $pf_iter = $local_pf_iter + 1;
                    }
                }
            }
        }

        my %pf_nodes;
        my @acc;
        my %acc2;
        push( @acc, keys %hmmer_ecodf );
        foreach my $acc (@acc) { $acc2{$acc}++ }
        @acc = sort { $a cmp $b } keys %acc2;
        foreach my $good_ecodf_acc (@acc) {
            my $ecodf_group_node;
            if ( $f_group->exists(qq{pf_group/ecodf[\@ecodf_acc="$good_ecodf_acc"]}) ) {
                $ecodf_group_node =
                  $f_group->findnodes(qq{pf_group/ecodf[\@ecodf_acc="$good_ecodf_acc"]})->get_node(1)->parentNode;
            }
            else {
                $ecodf_group_node = $ecod_xml_doc->createElement('pf_group');
                $ecodf_group_node->setAttribute( 'pf_id',         "$f_id.$pf_iter" );
                $ecodf_group_node->setAttribute( 'ecodf_cluster', 'true' );
                my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf');
                $ecodf_annot_node->setAttribute( 'ecodf_acc', $good_ecodf_acc );
                if ( $hmmer_ecodf_name{$good_ecodf_acc} ) {
                    $ecodf_annot_node->setAttribute( 'ecodf_name', $hmmer_ecodf_name{$good_ecodf_acc} );
                }
                $ecodf_group_node->appendChild($ecodf_annot_node);
                $f_group->appendChild($ecodf_group_node);
            }
            $pf_nodes{$good_ecodf_acc} = $ecodf_group_node;
            $pf_iter++;

        }

        #print "pf: $pf_iter?\n";

        #This needs to be fixed

        my $append = 0;
        foreach my $domain ( $f_group->findnodes('domain')->get_nodelist() ) {
            if ( $domain->findnodes('ecodf_annot')->size() == 1 ) {
                my $ecodf_acc = $domain->findvalue('ecodf_annot/@ecodf_acc');
                $pf_nodes{$ecodf_acc}->appendChild($domain);
                $append++;
            }
            else {
                my %check;
                foreach my $ecodf_annot_node ( $domain->findnodes('ecodf_annot')->get_nodelist() ) {
                    my $ecodf_acc = $ecodf_annot_node->findvalue('@ecodf_acc');

                    #print "2:$ecodf_acc\n";
                    $check{$ecodf_acc}++;
                }
                if ( scalar keys %check == 1 ) {
                    my @ecodf_accs = keys %check;
                    my $ecodf_acc  = $ecodf_accs[0];
                    $pf_nodes{$ecodf_acc}->appendChild($domain);
                    $append++;
                }
                elsif ( scalar keys %check > 1 ) {
                    my @ecodf_accs = sort { $a cmp $b } keys %check;
                    my $ecodf_acc = join( ",", @ecodf_accs );

                    #print "4:$ecodf_acc\n";
                    if ( !$pf_nodes{$ecodf_acc} ) {
                        my $ecodf_group_node = $ecod_xml_doc->createElement('pf_group');
                        $ecodf_group_node->setAttribute( 'pf_id',         "$f_id.$pf_iter" );
                        $ecodf_group_node->setAttribute( 'ecodf_cluster', 'true' );
                        my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf');
                        $ecodf_annot_node->setAttribute( 'ecodf_acc', $ecodf_acc );
                        if ( $hmmer_ecodf_name{$ecodf_acc} ) {
                            $ecodf_annot_node->setAttribute( 'ecodf_name', $hmmer_ecodf_name{$ecodf_acc} );
                        }
                        $ecodf_group_node->appendChild($ecodf_annot_node);
                        $f_group->appendChild($ecodf_group_node);
                        $pf_nodes{$ecodf_acc} = $ecodf_group_node;
                        $pf_iter++;
                    }
                    $pf_nodes{$ecodf_acc}->appendChild($domain);
                    $append++;
                    $pf_iter++;
                }
            }
        }

        #print "pf: $pf_iter ap: $append\n";

      DA:
        foreach my $domain_assembly_node ( $f_group->findnodes('domain_assembly') ) {
            my $uid = $domain_assembly_node->findvalue('@uid');
            my %da_ecodfs;
            foreach my $domain_node ( $domain_assembly_node->findnodes('domain[@primary="true"]') ) {
                if ( $domain_node->exists('ecodf_annot') ) {
                    if ( $domain_node->findnodes('ecodf_annot')->size() == 1 ) {
                        my $ecodf_acc = $domain_node->findvalue('ecodf_annot/@ecodf_acc');
                        print "1:$ecodf_acc\n";
                        $pf_nodes{$ecodf_acc}->appendChild($domain_assembly_node);
                        next DA;

                        #$da_ecodfs{$ecodf_acc}++
                        #$append++;
                    }
                    else {
                        my %check;
                        foreach my $ecodf_annot_node ( $domain_node->findnodes('ecodf_annot')->get_nodelist() ) {
                            my $ecodf_acc = $ecodf_annot_node->findvalue('@ecodf_acc');

                            #print "2:$ecodf_acc\n";
                            $check{$ecodf_acc}++;
                        }
                        if ( scalar keys %check == 1 ) {
                            my @ecodf_accs = keys %check;
                            my $ecodf_acc  = $ecodf_accs[0];

                            #print "3:$ecodf_acc\n";
                            $pf_nodes{$ecodf_acc}->appendChild($domain_assembly_node);
                            next DA;

                            #$da_ecodfs{$ecodf_acc}++;
                            #$append++;
                        }
                        elsif ( scalar keys %check > 1 ) {
                            my @ecodf_accs = sort { $a cmp $b } keys %check;
                            my $ecodf_acc = join( ",", @ecodf_accs );

                            #print "4:$ecodf_acc\n";
                            if ( !$pf_nodes{$ecodf_acc} ) {
                                my $ecodf_group_node = $ecod_xml_doc->createElement('pf_group');
                                $ecodf_group_node->setAttribute( 'pf_id',         "$f_id.$pf_iter" );
                                $ecodf_group_node->setAttribute( 'ecodf_cluster', 'true' );
                                my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf');
                                $ecodf_annot_node->setAttribute( 'ecodf_acc', $ecodf_acc );
                                if ( $hmmer_ecodf_name{$ecodf_acc} ) {
                                    $ecodf_annot_node->setAttribute( 'ecodf_name', $hmmer_ecodf_name{$ecodf_acc} );
                                }
                                $ecodf_group_node->appendChild($ecodf_annot_node);
                                $f_group->appendChild($ecodf_group_node);
                                $pf_nodes{$ecodf_acc} = $ecodf_group_node;
                                $pf_iter++;
                            }
                            $pf_nodes{$ecodf_acc}->appendChild($domain_assembly_node);
                            $append++;
                            next DA;
                            $pf_iter++;
                        }
                    }

                }
            }
        }

#Catch those representative domains
#		foreach my $domain ($f_group->findnodes('domain')->get_nodelist() ) {
#			my $ecod_domain_id = $domain->findvalue('@ecod_domain_id');
#			my $ecod_domain_uid = $domain->findvalue('@uid');
#			if ($domain->exists('ecod_representative_domain/@uid') && $domain->findvalue('@provisional_manual_rep') ne 'true') {
#				my $ecod_rep_uid = $domain->findvalue('ecod_representative_domain/@uid');
#				if ($f_group->exists(qq{pf_group/domain[\@uid="$ecod_rep_uid"]})) {
#					my $ecod_rep_node = $f_group->findnodes(qq{pf_group/domain[\@uid="$ecod_rep_uid"]})->get_node(1);
#					$ecod_rep_node->parentNode->appendChild($domain); #This is terrible and causing a lot of problems.
#				}elsif($f_group->exists(qq{domain[\@uid="$ecod_rep_uid"]})) {
#					#do nothing the rep node has not moved
#				}else{
#					print "WARNING! No rep node found for $ecod_domain_id $ecod_domain_uid?\n";
#				}
#			}
#		}
#

        if ( $domain_nodes->size() > 0 || scalar keys %hmmer_ecodf > 0 ) {
            printf "F: $f_id - found %i hmmer ecodf hits, %i domains\n", scalar keys %hmmer_ecodf,
              $domain_nodes->size();
        }
    }
}

#This is now obsolete, we no longer build pfgroups using pfam, only ecodf
sub cluster_pfgroups_by_pfam {
    my $sub = 'cluster_pf_groups_by_pfam';
    my ($ecod_xml_doc) = @_;
    foreach my $f_group ( $ecod_xml_doc->findnodes('//f_group')->get_nodelist() ) {

        my $f_id = $f_group->findvalue('@f_id');

        my %hmmer_ecodf;
        my %hmmer_ecodf_name;
        my $domain_nodes          = $f_group->findnodes('domain');
        my $domain_assembly_nodes = $f_group->findnodes('domain_assembly');

        foreach my $domain_node ( $f_group->findnodes('domain')->get_nodelist() ) {
            my $ecodf_string;
            if ( $domain_node->findnodes('ecodf_annot')->size() == 1 ) {
                my $ecodf_acc  = $domain_node->findvalue('ecodf_annot/@ecodf_acc');
                my $ecodf_name = $domain_node->findvalue('ecodf_annot/@ecodf_name');
                $ecodf_string = $ecodf_acc;

                #print "ecf_str $ecodf_string\n";
                $hmmer_ecodf{$ecodf_string}++;
                $hmmer_ecodf_name{$ecodf_string} = $ecodf_name;
            }
            elsif ( $domain_node->findnodes('ecodf_annot')->size() > 1 ) {
                my @ecodfs;
                my %ecodf_names;
                foreach my $ecodf_node ( $domain_node->findnodes('ecodf_annot')->get_nodelist() ) {
                    my $ecodf_acc = $ecodf_node->findvalue('@ecodf_acc');
                    push( @ecodfs, $ecodf_acc );
                    my $ecodf_name = $ecodf_node->findvalue('@ecodf_name');
                    $ecodf_names{$ecodf_acc} = $ecodf_name;
                }
                @ecodfs = sort { $a cmp $b } @ecodfs;
                $ecodf_string = join( ",", @ecodfs );

                #print "ecf_str_cat $ecodf_string\n";
                my @ecodf_names;
                foreach my $ecodf_acc (@ecodfs) {
                    push( @ecodf_names, $ecodf_names{$ecodf_acc} );
                }
                my $ecodf_name_string = join( ",", @ecodf_names );
                $hmmer_ecodf{$ecodf_string}++;
                $hmmer_ecodf_name{$ecodf_string} = $ecodf_name_string;
            }

        }

        foreach my $domain_assembly_node ( $f_group->findnodes('domain_assembly') ) {
            my $ecodf_string;
            my %assembly_ecodfs;
            my %assembly_ecodf_names;

            foreach my $domain_node ( $domain_assembly_node->findnodes('domain') ) {
                my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
                print "DA? $ecod_domain_id\n";
                if ( $domain_node->findnodes('ecodf_annot')->size() == 1 ) {
                    my $ecodf_acc  = $domain_node->findvalue('ecodf_annot/@ecodf_acc');
                    my $ecodf_name = $domain_node->findvalue('ecodf_annot/@ecodf_name');
                    $hmmer_ecodf{$ecodf_acc}++;
                    $hmmer_ecodf_name{$ecodf_acc} = $ecodf_name;

                }
                else {
                    my @ecodfs;
                    my %ecodf_names;
                    foreach my $ecodf_annot_node ( $domain_node->findnodes('ecodf_annot') ) {
                        my $ecodf_acc  = $ecodf_annot_node->findvalue('@ecodf_acc');
                        my $ecodf_name = $ecodf_annot_node->findvalue('@ecodf_name');
                        push( @ecodfs, $ecodf_acc );

                        $ecodf_names{$ecodf_acc} = $ecodf_name;
                    }
                    @ecodfs = sort { $a cmp $b } @ecodfs;
                    $ecodf_string = join( ",", @ecodfs );

                    #print "ecf_str_cat $ecodf_string\n";
                    my @ecodf_names;
                    foreach my $ecodf_acc (@ecodfs) {
                        push( @ecodf_names, $ecodf_names{$ecodf_acc} );
                    }
                    my $ecodf_name_string = join( ",", @ecodf_names );
                    $hmmer_ecodf{$ecodf_string}++;
                    $hmmer_ecodf_name{$ecodf_string} = $ecodf_name_string;

                }
            }
        }

        my $pf_iter = 1;
        my %pf_nodes;
        my @acc;
        my %acc2;
        push( @acc, keys %hmmer_ecodf );
        foreach my $acc (@acc) { $acc2{$acc}++ }
        @acc = sort { $a cmp $b } keys %acc2;
        foreach my $good_ecodf_acc (@acc) {
            my $ecodf_group_node;
            if ( $f_group->exists(qq{pf_group/ecodf[\@ecodf_acc="$good_ecodf_acc"]}) ) {
                $ecodf_group_node =
                  $f_group->findnodes(qq{pf_group/ecodf[\@ecodf_acc="$good_ecodf_acc"]})->get_node(1)->parentNode;
            }
            else {
                $ecodf_group_node = $ecod_xml_doc->createElement('pf_group');
                $ecodf_group_node->setAttribute( 'pf_id',         "$f_id.$pf_iter" );
                $ecodf_group_node->setAttribute( 'ecodf_cluster', 'true' );
                my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf');
                $ecodf_annot_node->setAttribute( 'ecodf_acc', $good_ecodf_acc );
                if ( $hmmer_ecodf_name{$good_ecodf_acc} ) {
                    $ecodf_annot_node->setAttribute( 'ecodf_name', $hmmer_ecodf_name{$good_ecodf_acc} );
                }
                $ecodf_group_node->appendChild($ecodf_annot_node);
                $f_group->appendChild($ecodf_group_node);
            }
            $pf_nodes{$good_ecodf_acc} = $ecodf_group_node;
            $pf_iter++;

        }

        print "pf: $pf_iter?\n";

        #This needs to be fixed

        my $append = 0;
        foreach my $domain ( $f_group->findnodes('domain')->get_nodelist() ) {
            if ( $domain->findnodes('ecodf_annot')->size() == 1 ) {
                my $ecodf_acc = $domain->findvalue('ecodf_annot/@ecodf_acc');

                #print "1:$ecodf_acc\n";
                $pf_nodes{$ecodf_acc}->appendChild($domain);
                $append++;
            }
            else {
                my %check;
                foreach my $ecodf_annot_node ( $domain->findnodes('ecodf_annot')->get_nodelist() ) {
                    my $ecodf_acc = $ecodf_annot_node->findvalue('@ecodf_acc');

                    #print "2:$ecodf_acc\n";
                    $check{$ecodf_acc}++;
                }
                if ( scalar keys %check == 1 ) {
                    my @ecodf_accs = keys %check;
                    my $ecodf_acc  = $ecodf_accs[0];

                    #print "3:$ecodf_acc\n";
                    $pf_nodes{$ecodf_acc}->appendChild($domain);
                    $append++;
                }
                elsif ( scalar keys %check > 1 ) {
                    my @ecodf_accs = sort { $a cmp $b } keys %check;
                    my $ecodf_acc = join( ",", @ecodf_accs );

                    #print "4:$ecodf_acc\n";
                    if ( !$pf_nodes{$ecodf_acc} ) {
                        my $ecodf_group_node = $ecod_xml_doc->createElement('pf_group');
                        $ecodf_group_node->setAttribute( 'pf_id',         "$f_id.$pf_iter" );
                        $ecodf_group_node->setAttribute( 'ecodf_cluster', 'true' );
                        my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf');
                        $ecodf_annot_node->setAttribute( 'ecodf_acc', $ecodf_acc );
                        if ( $hmmer_ecodf_name{$ecodf_acc} ) {
                            $ecodf_annot_node->setAttribute( 'ecodf_name', $hmmer_ecodf_name{$ecodf_acc} );
                        }
                        $ecodf_group_node->appendChild($ecodf_annot_node);
                        $f_group->appendChild($ecodf_group_node);
                        $pf_nodes{$ecodf_acc} = $ecodf_group_node;
                        $pf_iter++;
                    }
                    $pf_nodes{$ecodf_acc}->appendChild($domain);
                    $append++;
                    $pf_iter++;
                }
            }
        }
        print "pf: $pf_iter ap: $append\n";

      DA:
        foreach my $domain_assembly_node ( $f_group->findnodes('domain_assembly') ) {
            my $uid = $domain_assembly_node->findvalue('@uid');
            my %da_ecodfs;
            foreach my $domain_node ( $domain_assembly_node->findnodes('domain[@primary="true"]') ) {
                if ( $domain_node->exists('ecodf_annot') ) {
                    if ( $domain_node->findnodes('ecodf_annot')->size() == 1 ) {
                        my $ecodf_acc = $domain_node->findvalue('ecodf_annot/@ecodf_acc');
                        print "1:$ecodf_acc\n";
                        $pf_nodes{$ecodf_acc}->appendChild($domain_assembly_node);
                        next DA;

                        #$da_ecodfs{$ecodf_acc}++
                        #$append++;
                    }
                    else {
                        my %check;
                        foreach my $ecodf_annot_node ( $domain_node->findnodes('ecodf_annot')->get_nodelist() ) {
                            my $ecodf_acc = $ecodf_annot_node->findvalue('@ecodf_acc');

                            #print "2:$ecodf_acc\n";
                            $check{$ecodf_acc}++;
                        }
                        if ( scalar keys %check == 1 ) {
                            my @ecodf_accs = keys %check;
                            my $ecodf_acc  = $ecodf_accs[0];

                            #print "3:$ecodf_acc\n";
                            $pf_nodes{$ecodf_acc}->appendChild($domain_assembly_node);
                            next DA;

                            #$da_ecodfs{$ecodf_acc}++;
                            #$append++;
                        }
                        elsif ( scalar keys %check > 1 ) {
                            my @ecodf_accs = sort { $a cmp $b } keys %check;
                            my $ecodf_acc = join( ",", @ecodf_accs );

                            #print "4:$ecodf_acc\n";
                            if ( !$pf_nodes{$ecodf_acc} ) {
                                my $ecodf_group_node = $ecod_xml_doc->createElement('pf_group');
                                $ecodf_group_node->setAttribute( 'pf_id',         "$f_id.$pf_iter" );
                                $ecodf_group_node->setAttribute( 'ecodf_cluster', 'true' );
                                my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf');
                                $ecodf_annot_node->setAttribute( 'ecodf_acc', $ecodf_acc );
                                if ( $hmmer_ecodf_name{$ecodf_acc} ) {
                                    $ecodf_annot_node->setAttribute( 'ecodf_name', $hmmer_ecodf_name{$ecodf_acc} );
                                }
                                $ecodf_group_node->appendChild($ecodf_annot_node);
                                $f_group->appendChild($ecodf_group_node);
                                $pf_nodes{$ecodf_acc} = $ecodf_group_node;
                                $pf_iter++;
                            }
                            $pf_nodes{$ecodf_acc}->appendChild($domain_assembly_node);
                            $append++;
                            next DA;
                            $pf_iter++;
                        }
                    }

                }
            }
        }

        #Catch those representative domains
        foreach my $domain ( $f_group->findnodes('domain')->get_nodelist() ) {
            my $ecod_domain_id  = $domain->findvalue('@ecod_domain_id');
            my $ecod_domain_uid = $domain->findvalue('@uid');
            if (   $domain->exists('ecod_representative_domain/@uid')
                && $domain->findvalue('@provisional_manual_rep') ne 'true' )
            {
                my $ecod_rep_uid = $domain->findvalue('ecod_representative_domain/@uid');
                if ( $f_group->exists(qq{pf_group/domain[\@uid="$ecod_rep_uid"]}) ) {
                    my $ecod_rep_node = $f_group->findnodes(qq{pf_group/domain[\@uid="$ecod_rep_uid"]})->get_node(1);
                    $ecod_rep_node->parentNode->appendChild($domain);
                }
                elsif ( $f_group->exists(qq{domain[\@uid="$ecod_rep_uid"]}) ) {

                    #do nothing the rep node has not moved
                }
                else {
                    print "WARNING! No rep node found for $ecod_domain_id $ecod_domain_uid?\n";
                }
            }
        }

        printf "F: $f_id - found %i hmmer ecodf hits, %i domains\n", scalar keys %hmmer_ecodf, $domain_nodes->size();
    }

}

sub hmmer_annotation {
    my $sub = 'hmmer_annotation';
    my ( $ecod_xml_doc, $hmm_dump_xml_doc ) = @_;

    my %hmmer_map;
    foreach my $hmmer_node ( $hmm_dump_xml_doc->findnodes('//hmm_query')->get_nodelist() ) {
        my $ecod_domain_id = $hmmer_node->findvalue('@ecod_domain_id');
        my $uid            = $hmmer_node->findvalue('@uid');
        $hmmer_map{$uid} = $hmmer_node;
    }

    foreach my $domain_node ( $ecod_xml_doc->findnodes('//f_group/domain')->get_nodelist() ) {

        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
        my $uid            = $domain_node->findvalue('@uid');

        if ( $hmmer_map{$uid} ) {
            my $found_domain_id = $hmmer_map{$uid}->findvalue('@ecod_domain_id');
            if ( $hmmer_map{$uid}->findvalue('@ecod_domain_id') eq $ecod_domain_id ) {
                if ( !$domain_node->exists('hmm_query') ) {
                    $domain_node->appendChild( $hmmer_map{$uid} );
                }
                else {
                    print "WARNING! $ecod_domain_id $uid already has hmmer node\n";
                }
            }
            else {
                print "$uid $ecod_domain_id/$found_domain_id mismatch\n";
            }
        }
    }

}

sub convert_hmmer_to_ecodf_annotation {
    my $sub = 'convert_hmmer_to_ecodf_annotation';

    my ($ecod_xml_doc) = @_;

    my $EVAL_THRESHOLD = '1e-2';    #This shouldn't be here, move to update_ecod.pl or top of file

    foreach my $domain_node ( $ecod_xml_doc->findnodes('//domain')->get_nodelist() ) {

        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
        my $uid            = $domain_node->findvalue('@uid');

        if ( $domain_node->exists('ecodf_annot') ) { next }

        my $pfam_type = $domain_node->findvalue('hmm_query/@ecod_type');

        my %ranges;
        if ( $domain_node->exists('hmm_query/hmm_domannot') ) {
            foreach my $hmmer_rangenode ( $domain_node->findnodes('hmm_query/hmm_domannot')->get_nodelist ) {

                my $acc_full = $hmmer_rangenode->findvalue('@acc');
                unless ( $acc_full =~ /(EF\d+)/ ) {
                    warn "WARNING! $acc_full doesn't regexp, skipping...\n";
                    next;
                }
                my $acc = $1;
                my $n   = $hmmer_rangenode->findvalue('@N');
                if ( $n !~ /\d+/ ) { die "$ecod_domain_id $uid $acc_full $acc $n\n" }
                my $hmm_from = $hmmer_rangenode->findvalue('@hmm_from');
                my $hmm_to   = $hmmer_rangenode->findvalue('@hmm_to');
                my $ali_from = $hmmer_rangenode->findvalue('@ali_from');
                my $ali_to   = $hmmer_rangenode->findvalue('@ali_to');
                my $env_from = $hmmer_rangenode->findvalue('@env_from');
                my $env_to   = $hmmer_rangenode->findvalue('@env_to');

                my $ieval = $hmmer_rangenode->findvalue('@i_eval');
                my $ceval = $hmmer_rangenode->findvalue('@c_eval');

                if ( $ceval > $EVAL_THRESHOLD ) {

                    #print "FAIL domannot ceval: $ecod_domain_id $uid $acc_full $ceval\n";
                    next;
                }

                $ranges{$acc}{$n}{hmm_from} = $hmm_from;
                $ranges{$acc}{$n}{hmm_to}   = $hmm_to;
                $ranges{$acc}{$n}{ali_from} = $ali_from;
                $ranges{$acc}{$n}{ali_to}   = $ali_to;
                $ranges{$acc}{$n}{env_from} = $env_from;
                $ranges{$acc}{$n}{env_to}   = $env_to;

            }
        }

        if ( $domain_node->exists('hmm_query/hmm_summ') ) {
            my @global_hmm_used;
            my @global_ali_used;
            my @global_env_used;
            foreach my $hmmer_node ( $domain_node->findnodes('hmm_query/hmm_summ')->get_nodelist() ) {
                my $best_score = $hmmer_node->findvalue('@best_score');
                my $best_eval  = $hmmer_node->findvalue('@best_eval');

                if ( $best_eval > $EVAL_THRESHOLD ) {

                    #print "Fail summeval: $ecod_domain_id $uid $best_eval\n";
                    next;
                }

                #my $ecodf_name;
                #my $ecodf_acc;

                my $ecodf_acc = $hmmer_node->findvalue('@ecodf_acc');
                $ecodf_acc =~ /(EF\d+)/;
                $ecodf_acc = $1;

                my $ecodf_name = $hmmer_node->findvalue('@ecodf_name');

                if ( $domain_node->exists(qq{ecodf_annot[\@ecodf_acc="$ecodf_acc"]}) ) { next }

                my $ecodf_annot_node = $ecod_xml_doc->createElement('ecodf_annot');
                $ecodf_annot_node->setAttribute( 'ecodf_name', $ecodf_name );
                $ecodf_annot_node->setAttribute( 'ecodf_acc',  $ecodf_acc );

                if ( $hmmer_node->exists('@pfam_acc') ) {
                    my $pfam_acc = $hmmer_node->findvalue('@pfam_acc');
                    $ecodf_annot_node->setAttribute( 'pfam_acc', $pfam_acc );
                }
                my @hmm_used;
                my @ali_used;
                my @env_used;
                if ( $ranges{$ecodf_acc} ) {
                    foreach my $n ( sort { $a <=> $b } keys %{ $ranges{$ecodf_acc} } ) {
                        my @sub_hmm_used = ( $ranges{$ecodf_acc}{$n}{hmm_from} .. $ranges{$ecodf_acc}{$n}{hmm_to} );
                        my @sub_ali_used = ( $ranges{$ecodf_acc}{$n}{ali_from} .. $ranges{$ecodf_acc}{$n}{ali_to} );
                        my @sub_env_used = ( $ranges{$ecodf_acc}{$n}{env_from} .. $ranges{$ecodf_acc}{$n}{env_to} );

                        range_include( \@sub_hmm_used, \@hmm_used );
                        range_include( \@sub_ali_used, \@ali_used );
                        range_include( \@sub_env_used, \@env_used );
                    }

             #if hmm used range global previously used hmm range, or ali used range overlaps global  used ali range SKIP
             #Do you care about hmm used range? probably not, not relative, only want for comparison at later step
                    if ( residue_coverage( \@ali_used, \@global_ali_used ) > 10 ) {
                        next;
                    }
                    else {
                        range_include( \@ali_used, \@global_ali_used );
                    }

                    my $hmm_used_range = rangify(@hmm_used);
                    my $ali_used_range = rangify(@ali_used);
                    my $env_used_range = rangify(@env_used);

                    my $pfam_hmm_range_node = $ecod_xml_doc->createElement('hmm_range');
                    $pfam_hmm_range_node->appendTextNode($hmm_used_range);
                    my $pfam_ali_range_node = $ecod_xml_doc->createElement('ali_range');
                    $pfam_ali_range_node->appendTextNode($ali_used_range);
                    my $pfam_env_range_node = $ecod_xml_doc->createElement('env_range');
                    $pfam_env_range_node->appendTextNode($env_used_range);

                    $ecodf_annot_node->appendChild($pfam_hmm_range_node);
                    $ecodf_annot_node->appendChild($pfam_ali_range_node);
                    $ecodf_annot_node->appendChild($pfam_env_range_node);

                    $domain_node->appendChild($ecodf_annot_node);
                }
                else {
                    my $str = join( ",", keys %ranges );
                    print "WARNING! No ranges for ecodf_acc $ecodf_acc?! $ecod_domain_id $uid  ($str)\n";

                }
            }
        }
    }
}

#This is obsolete now
sub convert_hmmer_to_pfam_annotation {
    my $sub = 'convert_hmmer_to_pfam_annotation';

    my ($ecod_xml_doc) = @_;
    my $Pfam_clan_file = '/home/rschaeff/data/pfam/Pfam-A.clans.tsv';
    if ( !-f $Pfam_clan_file ) { die }
    my ( $acc_lookup, $name_lookup ) = clans_read($Pfam_clan_file);

    my $EVAL_THRESHOLD = '1e-3';

    foreach my $domain_node ( $ecod_xml_doc->findnodes('//domain')->get_nodelist() ) {

        my $ecod_domain_id = $domain_node->findvalue('@ecod_domain_id');
        my $uid            = $domain_node->findvalue('@uid');

        if ( $domain_node->exists('pfam_annot') ) { next }

        my %ranges;
        if ( $domain_node->exists('hmm_query/hmm_domannot') ) {
            foreach my $hmmer_rangenode ( $domain_node->findnodes('hmm_query/hmm_domannot')->get_nodelist ) {

                my $acc_full = $hmmer_rangenode->findvalue('@acc');
                $acc_full =~ /(PF\d+)/;
                my $acc = $1;
                my $n   = $hmmer_rangenode->findvalue('@N');
                if ( $n !~ /\d+/ ) { die "$ecod_domain_id $uid $acc_full $acc $n\n" }
                my $hmm_from = $hmmer_rangenode->findvalue('@hmm_from');
                my $hmm_to   = $hmmer_rangenode->findvalue('@hmm_to');
                my $ali_from = $hmmer_rangenode->findvalue('@ali_from');
                my $ali_to   = $hmmer_rangenode->findvalue('@ali_to');
                my $env_from = $hmmer_rangenode->findvalue('@env_from');
                my $env_to   = $hmmer_rangenode->findvalue('@env_to');

                my $ieval = $hmmer_rangenode->findvalue('@i_eval');
                my $ceval = $hmmer_rangenode->findvalue('@c_eval');

                if ( $ceval > $EVAL_THRESHOLD ) {

                    #print "FAIL domannot ceval: $ecod_domain_id $uid $acc_full $ceval\n";
                    next;
                }

                $ranges{$acc}{$n}{hmm_from} = $hmm_from;
                $ranges{$acc}{$n}{hmm_to}   = $hmm_to;
                $ranges{$acc}{$n}{ali_from} = $ali_from;
                $ranges{$acc}{$n}{ali_to}   = $ali_to;
                $ranges{$acc}{$n}{env_from} = $env_from;
                $ranges{$acc}{$n}{env_to}   = $env_to;

            }
        }

        if ( $domain_node->exists('hmm_query/hmm_summ') ) {
            my @global_hmm_used;
            my @global_ali_used;
            my @global_env_used;
            foreach my $hmmer_node ( $domain_node->findnodes('hmm_query/hmm_summ')->get_nodelist() ) {
                my $best_score = $hmmer_node->findvalue('@best_score');
                my $best_eval  = $hmmer_node->findvalue('@best_eval');

                #			if ($best_score < $BITS_THRESHOLD) {
                #				next;
                #			}

                if ( $best_eval > $EVAL_THRESHOLD ) {

                    #print "Fail summeval: $ecod_domain_id $uid $best_eval\n";
                    next;
                }

                my $pfam_name;
                my $pfam_acc;
                if ( $hmmer_node->exists('@model') ) {
                    $pfam_name = $hmmer_node->findvalue('@model');
                    $pfam_acc  = $$acc_lookup{$pfam_name};
                }
                else {
                    $pfam_acc  = $hmmer_node->findvalue('@acc');
                    $pfam_name = $$name_lookup{$pfam_acc};
                }
                if ( $domain_node->exists(qq{pfam_annot[\@pfam_name="$pfam_name"]}) ) { next }

                my $pfam_annot_node = $ecod_xml_doc->createElement('pfam_annot');
                $pfam_annot_node->setAttribute( 'match_type', 'hmmer' );
                $pfam_annot_node->setAttribute( 'pfam_name',  $pfam_name );
                $pfam_annot_node->setAttribute( 'pfam_acc',   $pfam_acc );
                $pfam_annot_node->setAttribute( 'pfam_type',  $HMM_LABEL );    #Usually PfamA?

                my @hmm_used;
                my @ali_used;
                my @env_used;
                if ( $ranges{$pfam_acc} ) {
                    foreach my $n ( sort { $a <=> $b } keys %{ $ranges{$pfam_acc} } ) {
                        my @sub_hmm_used = ( $ranges{$pfam_acc}{$n}{hmm_from} .. $ranges{$pfam_acc}{$n}{hmm_to} );
                        my @sub_ali_used = ( $ranges{$pfam_acc}{$n}{ali_from} .. $ranges{$pfam_acc}{$n}{ali_to} );
                        my @sub_env_used = ( $ranges{$pfam_acc}{$n}{env_from} .. $ranges{$pfam_acc}{$n}{env_to} );

                        range_include( \@sub_hmm_used, \@hmm_used );
                        range_include( \@sub_ali_used, \@ali_used );
                        range_include( \@sub_env_used, \@env_used );
                    }

             #if hmm used range global previously used hmm range, or ali used range overlaps global  used ali range SKIP
             #Do you care about hmm used range? probably not, not relative, only want for comparison at later step
                    if ( residue_coverage( \@ali_used, \@global_ali_used ) > 10 ) {
                        next;
                    }
                    else {
                        range_include( \@ali_used, \@global_ali_used );
                    }

                    my $hmm_used_range = rangify(@hmm_used);
                    my $ali_used_range = rangify(@ali_used);
                    my $env_used_range = rangify(@env_used);

                    my $pfam_hmm_range_node = $ecod_xml_doc->createElement('hmm_range');
                    $pfam_hmm_range_node->appendTextNode($hmm_used_range);
                    my $pfam_ali_range_node = $ecod_xml_doc->createElement('ali_range');
                    $pfam_ali_range_node->appendTextNode($ali_used_range);
                    my $pfam_env_range_node = $ecod_xml_doc->createElement('env_range');
                    $pfam_env_range_node->appendTextNode($env_used_range);

                    $pfam_annot_node->appendChild($pfam_hmm_range_node);
                    $pfam_annot_node->appendChild($pfam_ali_range_node);
                    $pfam_annot_node->appendChild($pfam_env_range_node);

                    $domain_node->appendChild($pfam_annot_node);
                }
                else {
                    my $str = join( ",", keys %ranges );
                    print "WARNING! No ranges for pfam_acc $pfam_acc?! $ecod_domain_id $uid  ($str)\n";

                }
            }
        }
    }
}

sub read_hhsearch_dir {
    my $sub = 'read_hhsearch_dir';
    my ($dir) = @_;

    opendir( DIR, $dir ) or die "ERROR! $sub: Could not opendir $dir for reading:$!\n";
    my @files = grep { $_ =~ /.*\.hhm/ } readdir(DIR);
    close DIR;

    printf "Found %i hhsearch profiles\n", scalar(@files);

    my %files;
    foreach my $file (@files) {
        $files{$file} = "$dir/$file";
    }

    return \%files;
}

###
sub read_fasta_file {
    my $sub = 'read_fasta_file';
    my ($file) = @_;

    open( FA, $file ) or die "ERROR! $sub: Could not open $file for reading:$!\n";
    my %fasta;
    my $ecod_domain_id;
    while ( my $ln = <FA> ) {
        chomp $ln;
        if ( $ln =~ /^>(e\w{4}.\d+)/ ) {
            $ecod_domain_id = $1;
        }
        else {
            $fasta{$ecod_domain_id} .= $ln;
        }
    }
    return \%fasta;
}

sub generate_blast_library {
    my $sub = "generate_blast_library";
    my ( $fa_file, $blast_prefix ) = @_;

    if ( !-f $fa_file ) { die "ERROR! $sub: $fa_file not found\n"; }

    my $command = "$MAKEBLAST_EXE -in $fa_file -out $blast_prefix -title $blast_prefix\n";
    system($command);

}

sub clans_read {
    my $sub = 'clans_read';
    my ($file) = @_;

    open( IN, $file ) or die;
    my %acc_lookup;
    my %name_lookup;
    while ( my $ln = <IN> ) {
        my @F = split( /\t/, $ln );

        my $acc      = $F[0];
        my $clan_acc = $F[1];
        my $cla_name = $F[2];
        my $pf_name  = $F[3];
        my $pf_desc  = $F[4];
        $acc_lookup{$pf_name} = $acc;
        $name_lookup{$acc}    = $pf_name;
    }
    return ( \%acc_lookup, \%name_lookup );
}
1;
