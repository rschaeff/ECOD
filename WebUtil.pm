package ECOD::WebUtil;
require Exporter;

use warnings;
use strict;

use LWP::UserAgent;
use XML::Grishin;
use Domains::Partition;
use Domains::PDBml;
use Domains::Range;
use SQL::Util;

my $WEB_SPECIAL_URL = "http://morpho.swmed.edu:5001/add";

our @ISA = qw(Exporter);
our @EXPORT = (
	"&web_register_peptide",
	"&run_list_web_register_special",
	"&transfer_special_annotation",
	);
					

sub transfer_special_annotation { 

	my $ecod_curation_dbh = db_connect('ecod_curation');
	my $ecod_production_dbh = db_connect('ecod');

#Gather non-nucleotide special ranges

	my $special_prod_sth = $ecod_production_dbh->prepare(qq{SELECT uid, pdb, chain, seqid_range, type FROM special WHERE type != 'nonpeptide_poly'});
	$special_prod_sth->execute();
	my $insert_sth = $ecod_curation_dbh->prepare(sql_seq_insert() );
	my $filter_sth = $ecod_curation_dbh->prepare("SELECT uid FROM special_seq");
	$filter_sth->execute();
	my %filter;
	while (my $row_aref = $filter_sth->fetchrow_arrayref() ) { $filter{$$row_aref[0]}++ }
	while (my $row_aref = $special_prod_sth->fetchrow_arrayref() ) { 
		next if $filter{$$row_aref[0]};
		my $uid	= sprintf "%09i", $$row_aref[0];
		my $pdb_id 		= $$row_aref[1];
		my $chain_id 	= $$row_aref[2];
		my $seqid_range = $$row_aref[3];
		my $type		= $$row_aref[4];

		my ($seqid_aref, $struct_seqid_aref, $pdbnum_href, $asym_id) = pdbml_seq_parse($pdb_id, $chain_id);
		
		if (!$seqid_aref || !$asym_id) { 
			warn "WARNING! No PDB for $pdb_id, $chain_id\n";
			next;
		}

		my $fa;
		if ($seqid_range) { 
			my ($range_str, $chain_str) = scop_range_split($seqid_range);
			$fa = pdbml_fasta_fetch($pdb_id, $asym_id, $chain_id, range_expand($range_str));
		}else{
			$fa = pdbml_fasta_fetch($pdb_id, $asym_id, $chain_id, $seqid_aref);
		}
		$insert_sth->execute($uid, $pdb_id, $chain_id, $seqid_range,$fa, $type);
	}


}

sub web_register_peptide  { 
	my $job_xml_doc = $_[0];
	my $type = $_[1];


	my $seqid_range;
	my $ua = LWP::UserAgent->new();
	foreach my $job_node (find_job_nodes($job_xml_doc)) { 
		my ($pdb_id, $chain_id) = get_pdb_chain($job_node);
		print "$pdb_id $chain_id $type\n";

		if ($job_node->findvalue('peptide_filter/@apply') eq 'true') { 
			my $resp = $ua->post($WEB_SPECIAL_URL, 
						{ pdb_id 	=> $pdb_id,
						  chain_id 	=> $chain_id,
						  seqid_range => $seqid_range,
						  type => $type		});
			if ($resp->is_error()) { 
				printf "%s\n", $resp->status_line;
			}
		}
	}
}

sub run_list_web_register_special { 

	my $sub = 'run_list_web_register_special';

	my $run_list_file = $_[0];
	my $type = $_[1];

	print "$sub $run_list_file $type\n";

	my $run_list_xml_doc = xml_open($run_list_file);

	my $run_list_dir = get_run_list_dir($run_list_xml_doc);

	my %type2action = ( emuq => '@to_emuq',
						coil => '@to_coil'  );

	my $action = $type2action{$type};

	foreach my $job_list_node (find_job_list_nodes($run_list_xml_doc)) { 
		my $job_list_fn = get_job_list_fn($job_list_node, "seq_iter");
		my $job_xml_doc = xml_open("$run_list_dir/$job_list_fn");
		job_list_register_special($job_xml_doc, $action, $type);
	}
}

sub job_list_register_special { 
	my $job_xml_doc = $_[0];
	my $action = $_[1];
	my $type = $_[2];

	
	my $ua = LWP::UserAgent->new();
	foreach my $job_node (find_job_nodes($job_xml_doc)) { 
		my ($pdb_id, $chain_id) = get_pdb_chain($job_node);
		print "$pdb_id $chain_id $type\n";
		if ($job_node->findvalue($action) eq 'true') { 
			my $resp = $ua->post($WEB_SPECIAL_URL,
				{ 	pdb_id 		=> $pdb_id,
					chain_id	=> $chain_id,
					type 		=> $type});
			if ($resp->is_error()) { 
				printf "ERROR! %s\n", $resp->status_line
			}
		}
	}
}

sub sql_seq_insert { 
my $q = <<EOF;	
	INSERT INTO special_seq ( 
	uid,
	pdb_id,
	chain_id,
	seqid_range,
	fasta_seq,
	type
	) VALUES ( ?, ?, ?, ?, ?, ? );
EOF
}
