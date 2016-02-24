package ECOD::Reference;
require Exporter;

use warnings;
use strict;

use XML::LibXML;

our @ISA = qw(Exporter);

our @EXPORT = (	'%DOMAIN_REF',
		'%PREVIOUS_VERSIONS', 
		'%LOCAL_DIR',
		'%PDB_XML_DOC',
		'%CHAIN_REF',
		'%CHAIN_REPSONLY_REF',
		'%DOMAIN_REPSONLY_REF', 
		'%HH_REF',
		'%REF_DALI_CACHE',
		'%REF_XML',
		'%REF_RANGE_CACHE',
		'$LATEST_REFERENCE',
		'$DOMAIN_DATA_DIR',
		'&load_references',
		'&register_version',
		'&update_current_version',
		'&reference_cache_load',
		'&reference_chainwise_transform',
		'&register_pdb_doc');

our %CHAIN_REF;
our %DOMAIN_REF; 

our %CHAIN_REPSONLY_REF;
our %DOMAIN_REPSONLY_REF;

our %HH_REF;
our %REF_RANGE_CACHE;
our %REF_DALI_CACHE;
our %REF_XML;
our %LOCAL_DIR;
our %PDB_XML_DOC;
our %PREVIOUS_VERSIONS;

our $DOMAIN_DATA_DIR;


#my $REFERENCE_INDEX = '/home/rschaeff/lib/ECOD/ecod.tmp.ref.xml';
my $REFERENCE_INDEX = '/data/ecod/database_versions/ecod.master_ref.xml';

#Loads the current set of known ECOD references and populates global data structures
sub load_references { 

	if (!-f $REFERENCE_INDEX) { 
		die "ERROR! Could not find reference index: $REFERENCE_INDEX\n";
	}
	$DOMAIN_DATA_DIR = '/data/ecod/domain_data';
	if (! -d $DOMAIN_DATA_DIR) { 
		die "ERROR! Could not find ECOD domain data directory: $DOMAIN_DATA_DIR\n";
	}

	open (my $xml_fh, $REFERENCE_INDEX) or die "ERROR! Could not open $REFERENCE_INDEX for reading:$!\n";
	my $ref_xml_doc = XML::LibXML->load_xml(IO => $xml_fh);

	foreach my $ref_node ($ref_xml_doc->findnodes('//reference_list/reference') ) { 
		
		my $version		= $ref_node->findvalue('@version');

		if ($ref_node->exists('chainwise_fa_ref')) { 
			my $chainwise_fa_ref 	= $ref_node->findvalue('chainwise_fa_ref');
			$CHAIN_REF{$version}	= $chainwise_fa_ref;
		}
		if ($ref_node->exists('chainwise_reps_only_fa_ref')) { 
			my $chainwise_reps_only_fa_ref 	= $ref_node->findvalue('chainwise_reps_only_fa_ref');
			$CHAIN_REPSONLY_REF{$version}	= $chainwise_reps_only_fa_ref;
		}
		if ($ref_node->exists('domain_fa_ref')) { 
			my $domain_fa_ref	= $ref_node->findvalue('domain_fa_ref');
			$DOMAIN_REF{$version}	= $domain_fa_ref;
		}
		if ($ref_node->exists('domain_reps_only_fa_ref')) { 
			my $domain_reps_only_fa_ref	= $ref_node->findvalue('domain_reps_only_fa_ref');
			$DOMAIN_REPSONLY_REF{$version}	= $domain_reps_only_fa_ref;
		}
		if ($ref_node->exists('hhm_ref')) { 
			my $hhm_ref		= $ref_node->findvalue('hhm_ref');
			$HH_REF{$version} = $hhm_ref;
		}
		if ($ref_node->exists('ref_xml')) { 
			my $ref_xml		= $ref_node->findvalue('ref_xml');
			$REF_XML{$version}	= $ref_xml
		}
		if ($ref_node->exists('ref_range_cache')) { 
			my $range_cache		= $ref_node->findvalue('ref_range_cache');
			$REF_RANGE_CACHE{$version}	= $range_cache;
		}
		if ($ref_node->exists('pdb_doc')) { 
			my $pdb_doc 		= $ref_node->findvalue('pdb_doc');
			$PDB_XML_DOC{$version}	= $pdb_doc;
		}
		if ($ref_node->exists('local_dir')) { 
			my $local_dir		= $ref_node->findvalue('local_dir');
			$LOCAL_DIR{$version}	= $local_dir;
		}
		if ($ref_node->exists('@previous_version')) { 
			my $previous_version 	= $ref_node->findvalue('@previous_version');
			$PREVIOUS_VERSIONS{$version}	= $previous_version;
		}
		if ($ref_node->exists('dali_cache')) { 
			my $dali_cache	= $ref_node->findvalue('dali_cache');
			$REF_DALI_CACHE{$version} = $dali_cache;
		}
	}
	our $LATEST_REFERENCE = $ref_xml_doc->findvalue('//@currentVersion');

}
sub reference_cache_load { 
	my $sub = 'reference_cache_load';

	my ($reference) = @_;

	my $ref_loc;
	if ($REF_RANGE_CACHE{$reference} && -f $REF_RANGE_CACHE{$reference}) { 
		$ref_loc = $REF_RANGE_CACHE{$reference};
	}else{
		die "ERROR ref range cahce for $reference $ref_loc not found\n";
	}

	open (IN, $ref_loc) or die "ERROR! $sub	: Could not open $ref_loc for reading:$!\n";

	my %ref_range;
	my %uid_lookup;
	while (my $ln = <IN>) { 
		if ($ln =~ /^#/) { next } 
		my @F = split(/\s+/, $ln);

		my $uid			= $F[0];
		my $ecod_domain_id	= $F[1];
		my $seqid_range		= $F[2];
		my $pdb			= $F[3];
		my $chain		= $F[4];

		$ref_range{$uid}{ecod_domain_id}	= $ecod_domain_id;
		$ref_range{$uid}{seqid_range}		= $seqid_range;
		$ref_range{$uid}{pdb}			= $pdb;
		$ref_range{$uid}{chain}			= $chain;

		#print "CACHE LOAD u:$uid ed:$ecod_domain_id sr:$seqid_range p:$pdb c:$chain\n";

		$uid_lookup{$ecod_domain_id} 		= $uid;
	}
	close IN;

	return (\%ref_range, \%uid_lookup);
}

sub reference_chainwise_transform { 
	my $sub = 'reference_chainwise_transform';

	my ($ref_range_href) = @_;

	my %chain_domains;
	foreach my $domain_uid (keys %$ref_range_href) { 

		my $pdb 		= $$ref_range_href{$domain_uid}{pdb};
		my $chain 		= $$ref_range_href{$domain_uid}{chain};
		my $seqid_range		= $$ref_range_href{$domain_uid}{seqid_range};
		my $ecod_domain_id 	= $$ref_range_href{$domain_uid}{ecod_domain_id};

		push (@{$chain_domains{$pdb}{$chain}}, $domain_uid); 

	}

	return \%chain_domains;
}
#Registers a new version of ECOD (libraries included) 
sub register_version { 
	my $sub = 'register_version';
	my ($update_node, $force_overwrite) = @_;

	if (!-f $REFERENCE_INDEX) { 
		die "ERROR! Could not find reference index: $REFERENCE_INDEX\n";
	}

	open (my $xml_fh, $REFERENCE_INDEX) or die "ERROR! Could not open $REFERENCE_INDEX for reading:$!\n";
	my $ref_xml_doc = XML::LibXML->load_xml(IO => $xml_fh);
	close $xml_fh;

	my $version = $update_node->findvalue('@version');
	if ($ref_xml_doc->exists(qq{//reference_list/reference[\@version="$version"]})) { 
		print "WARNING! $sub: version $version already registered\n";
		return 0;
	}
	my $previous_version = $update_node->findvalue('@previous_version');

	my $top_dir = $update_node->findvalue('directory');

	my $chainwise_fa_ref 	= $update_node->findvalue('derived_libraries/chainwise_blast_lib');
	my $domain_fa_ref 	= $update_node->findvalue('derived_libraries/domain_blast_lib');
	
	my $chainwise_reps_only_fa_ref = $update_node->findvalue('derived_libraries/chainwise_blast_reps_only_lib');
	my $domain_reps_only_fa_ref = $update_node->findvalue('derived_libraries/domain_blast_reps_only_lib');

	my $domain_fa	= $update_node->findvalue('derived_libraries/domain_fasta');

	my $hhm_ref		= $update_node->findvalue('derived_libraries/hhm_db');
	my $ref_xml		= $update_node->findvalue('ecod_xml');

	my $range_cache		= $update_node->findvalue('derived_libraries/domain_range_cache');
	my $dali_cache		= $update_node->findvalue('derived_libraries/dali_cache');

	my $reference_node 	= $ref_xml_doc->createElement('reference');
	$reference_node->setAttribute('version', $version);
	$reference_node->setAttribute('previous_version', $previous_version);
	$reference_node->appendTextChild('chainwise_fa_ref', "$top_dir/$chainwise_fa_ref");
	$reference_node->appendTextChild('domain_fa_ref', "$top_dir/$domain_fa_ref");

	$reference_node->appendTextChild('chainwise_reps_only_fa_ref', "$top_dir/$chainwise_reps_only_fa_ref");
	$reference_node->appendTextChild('domain_reps_only_fa_ref', "$top_dir/$domain_reps_only_fa_ref");

	$reference_node->appendTextChild('domain_fa', "$top_dir/$domain_fa");
	$reference_node->appendTextChild('hhm_ref', "$top_dir/$hhm_ref");
	$reference_node->appendTextChild('ref_xml', "$top_dir/$ref_xml");
	$reference_node->appendTextChild('ref_range_cache', "$top_dir/$range_cache");
	$reference_node->appendTextChild('dali_cache', "$top_dir/$dali_cache");
	$reference_node->appendTextChild('local_dir', $top_dir);

	my $reference_list_node	= $ref_xml_doc->findnodes('//reference_list')->get_node(1);
	$reference_list_node->appendChild($reference_node);

	my $doc_string = $ref_xml_doc->toString(1);	
	open (OUT, ">$REFERENCE_INDEX") or die "ERROR! Could not open $REFERENCE_INDEX for writing:$!\n";
	print OUT $doc_string;
	close OUT;

	return 1;
}

sub register_pdb_doc { 
	my $sub = 'register_pdb_doc';

	my ($pdb_fn, $version) = @_;

	if (!-f $REFERENCE_INDEX) { 
		die "ERROR! Could not find reference index; $REFERENCE_INDEX\n";;
	}

	open (my $xml_fh, $REFERENCE_INDEX) or die "ERROR! Could not open $REFERENCE_INDEX for reading:$!\n";
	my $ref_xml_doc = XML::LibXML->load_xml(IO => $xml_fh);
	close $xml_fh;

	if (!$ref_xml_doc->exists(qq{//reference_list/reference[\@version="$version"]})) { 
		print "WARNING! $sub: version $version not created\n";
		return 0;
	}

	my $reference_node = $ref_xml_doc->findnodes(qq{//reference_list/reference[\@version="$version"]})->get_node(1);
	if ($reference_node->exists('pdb_doc')) { 
		print "WARNING! PDB doc already registered for $version\n";
		return 0;
	}else{
		$reference_node->appendTextChild('pdb_doc', $pdb_fn);
	}
	my $doc_string = $ref_xml_doc->toString(1);	

	open (OUT, ">$REFERENCE_INDEX") or die "ERROR! Could not open $REFERENCE_INDEX for writing:$!\n";
	print OUT $doc_string;
	close OUT;

	return 1;
}
#updates the current_version attribute of the reference index
sub update_current_version { 
	my $sub = 'update_current_version';
	my ($version) = @_;
	if (!-f $REFERENCE_INDEX) { 
		die "ERROR! Could not find reference index: $REFERENCE_INDEX\n";
	}

	open (my $xml_fh, $REFERENCE_INDEX) or die "ERROR! Could not open $REFERENCE_INDEX for reading:$!\n";
	my $ref_xml_doc = XML::LibXML->load_xml(IO => $xml_fh);
	close $xml_fh;

	#if ($ref_xml_doc->exists(qq{//reference_list/reference[\@version="$version"]})) { 
#		print "WARNING! $sub: version $version already registered\n";
#		return 0;
#	}

	my $doc_node = $ref_xml_doc->findnodes('/ecod_reference_document')->get_node(1);
	my $current_version = $doc_node->findvalue("\@currentVersion");
	if ($current_version ne $version) { 
		$doc_node->setAttribute('currentVersion', $version);
		our $LATEST_REFERENCE = $version;
	}else{
		print "WARNING! $version already currentVersion\n";
		return 0;
	}
	my $doc_string = $ref_xml_doc->toString(1);	
	open (OUT, ">$REFERENCE_INDEX") or die "ERROR! Could not open $REFERENCE_INDEX for writing:$!\n";
	print OUT $doc_string;
	close OUT;



}

#$CHAIN_REF{'scop175'} = '/home/rschaeff/data/scop/v1.75/chainwise.scop.v1.75';  
#$CHAIN_REF{'scop169'} = '/home/rschaeff/data/scop/v1.69/chainwise.scop.v1.69';  
#$CHAIN_REF{'develop2'} = '/home/rschaeff/data/ecod/v2/chainwise_ecod100.develop2';  
#$CHAIN_REF{'develop8'} = '/home/rschaeff/data/ecod/v8/chainwise_ecod100.develop8'; 

#$CHAIN_REF{'develop10'} = '/home/rschaeff/data/ecod/v10/chainwise100.develop10';  
#$CHAIN_REF{'develop12'} = '/home/rschaeff/data/ecod/v12/chainwise100.develop12'; 
#$CHAIN_REF{'develop13'} = '/home/rschaeff/data/ecod/v13/chainwise100.develop13';
#$CHAIN_REF{'develop14'} = '/home/rschaeff/data/ecod/v14/chainwise100.develop14';
#$CHAIN_REF{'develop17'} = '/home/rschaeff/data/ecod/v17/chainwise100.develop17';
#$CHAIN_REF{'develop18'} = '/home/rschaeff/data/ecod/v18/chainwise100.develop18c';
#$CHAIN_REF{'develop19'} = '/home/rschaeff/data/ecod/v19/chainwise100.develop19';
#$CHAIN_REF{'develop20'} = '/home/rschaeff/data/ecod/v20/chainwise100.develop20';
#$CHAIN_REF{'develop21'} = '/home/rschaeff/data/ecod/v21/chainwise100.develop21';
#$CHAIN_REF{'develop24'} = '/home/rschaeff/data/ecod/v24/chainwise100.develop24';
#$CHAIN_REF{'develop25'} = '/home/rschaeff/data/ecod/v25/chainwise100.develop25e';
#$CHAIN_REF{'develop26'} = '/home/rschaeff/data/ecod/v26/chainwise100.develop26b';
#$CHAIN_REF{'develop27'} = '/home/rschaeff/data/ecod/v27/chainwise100.develop27e';
#$CHAIN_REF{'develop28'} = '/home/rschaeff/data/ecod/v28/chainwise100.develop28b';
#$CHAIN_REF{'develop29'} = '/home/rschaeff/data/ecod/v29/chainwise100.develop29';
#$CHAIN_REF{'develop29.reps_only'} = '/home/rschaeff/data/ecod/v29/chainwise100.reps_only.develop29';
#$CHAIN_REF{'develop30'} = '/home/rschaeff/data/ecod/v30/chainwise100.develop30';
#$CHAIN_REF{'develop30.reps_only'} = '/home/rschaeff/data/ecod/v30/chainwise100.develop30.reps_only';
#$CHAIN_REF{'develop31'} = '/home/rschaeff/data/ecod/v31/chainwise100.develop31a';
#$CHAIN_REF{'develop32'} = '/home/rschaeff/data/ecod/v32/chainwise100.develop32';
#$CHAIN_REF{'develop33'} = '/home/rschaeff/data/ecod/v33/chainwise100.develop33';
#
##
#$DOMAIN_REF{'develop1'} = '/home/rschaeff/data/ecod/v1/ecod100.struct';
#$DOMAIN_REF{'develop2'} = '/home/rschaeff/data/ecod/v2/ecod100.develop2';
#$DOMAIN_REF{'develop8'}  = '/home/rschaeff/data/ecod/v8/ecod100.develop8';
#$DOMAIN_REF{'develop10'} = '/home/rschaeff/data/ecod/v10/ecod100.develop10';
#$DOMAIN_REF{'develop12'} = '/home/rschaeff/data/ecod/v12/ecod100.develop12';
#$DOMAIN_REF{'develop13'} = '/home/rschaeff/data/ecod/v13a/ecod100.develop13';
#$DOMAIN_REF{'develop14'} = '/home/rschaeff/data/ecod/v14/ecod100.develop14';
#$DOMAIN_REF{'develop17'} = '/home/rschaeff/data/ecod/v17/ecod100.develop17';
#$DOMAIN_REF{'develop18'} ='/home/rschaeff/data/ecod/v18/ecod100.develop18c';
#$DOMAIN_REF{'develop19'} ='/home/rschaeff/data/ecod/v19/ecod100.develop19';
#$DOMAIN_REF{'develop20'} ='/home/rschaeff/data/ecod/v20/ecod100.develop20';
#$DOMAIN_REF{'develop21'} ='/home/rschaeff/data/ecod/v21/ecod100.develop21';
#$DOMAIN_REF{'develop24'} ='/home/rschaeff/data/ecod/v24/ecod100.develop24';
#$DOMAIN_REF{'develop26'} ='/home/rschaeff/data/ecod/v26/ecod100.develop26b';
#$DOMAIN_REF{'develop27'} ='/home/rschaeff/data/ecod/v27/ecod100.develop27e';
#$DOMAIN_REF{'develop28'} ='/home/rschaeff/data/ecod/v28/ecod100.develop28b';
#$DOMAIN_REF{'develop29'} ='/home/rschaeff/data/ecod/v29/ecod100.develop29';
#$DOMAIN_REF{'develop29.reps_only'} ='/home/rschaeff/data/ecod/v29/ecod100.reps_only.develop29';
#$DOMAIN_REF{'develop30'} ='/home/rschaeff/data/ecod/v30/ecod100.develop30';
#$DOMAIN_REF{'develop31'} ='/home/rschaeff/data/ecod/v31/ecod100.develop31a';
#$DOMAIN_REF{'develop32'} ='/home/rschaeff/data/ecod/v32/ecod100.develop32';
#$DOMAIN_REF{'develop33'} ='/home/rschaeff/data/ecod/v33/ecod100.develop33';
#$DOMAIN_REF{'scop175'} ='/home/rschaeff/data/scop/v1.75/scop100.v1.75';
#$DOMAIN_REF{'scop169'} ='/home/rschaeff/data/scop/v1.69/scop100.v1.69';
#
#$HH_REF{'develop8'}                = '/home/rschaeff/bin/evdb/local_dbs/ecod40.develop8.new/ecod.develop8.rebuild_hhm_db';
#$HH_REF{'develop10'}               = '/home/rschaeff/data/ecod/v10/ecod.develop10.rebuild_hhm_db';
#$HH_REF{'develop12'}               = '/home/rschaeff/data/ecod/v12/ecod.develop12_all.hhm_db';
#$HH_REF{'develop13'}               = '/home/rschaeff/data/ecod/v13a/ecod.develop13a.hhm_db';
#$HH_REF{'develop14'}               = '/home/rschaeff/data/ecod/v14/ecod.develop14.hhm_db';
#$HH_REF{'develop18'}               = '/home/rschaeff/data/ecod/v18/ecod.develop18.hhm_db';
#$HH_REF{'develop19'}               = '/home/rschaeff/data/ecod/v19/ecod.develop19.hhm_db';
#$HH_REF{'develop21'}               = '/home/rschaeff/data/ecod/v21/ecod.develop21.hhm_db';
#$HH_REF{'develop24'}               = '/home/rschaeff/data/ecod/v24/ecod.develop24.hhm_db';
#$HH_REF{'develop26'}               = '/home/rschaeff/data/ecod/v26/ecod.develop26b.hhm_db';
#$HH_REF{'develop27'}               = '/home/rschaeff/data/ecod/v27/ecod.develop27g.hhm_db';
#$HH_REF{'develop28'}               = '/home/rschaeff/data/ecod/v28/ecod.develop28b.hhm_db';
#$HH_REF{'develop29'}               = '/home/rschaeff/data/ecod/v29/ecod.develop29.hhm_db';
#$HH_REF{'develop29.reps_only'}               = '/home/rschaeff/data/ecod/v29/ecod.develop29.hhm_db';
#$HH_REF{'develop30'}               = '/home/rschaeff/data/ecod/v30/ecod.develop30.hhm_db';
#$HH_REF{'develop30.reps_only'}               = '/home/rschaeff/data/ecod/v30/ecod.develop30.hhm_db';
#$HH_REF{'develop31'}               = '/home/rschaeff/data/ecod/v31/ecod.develop31a.hhm_db';
#$HH_REF{'develop32'}               = '/home/rschaeff/data/ecod/v32/ecod.develop32.hhm_db';
#$HH_REF{'develop33'}               = '/home/rschaeff/data/ecod/v33/ecod.develop33.hhm_db';
#$HH_REF{'scop175'}		   = '/home/rschaeff/data/scop/v1.75/scop175.hhm_db';
#$HH_REF{'scop169'}		   = '/home/rschaeff/data/scop/v1.69/scop169.hhm_db';
#
#$REF_XML{'scop175'}	= '/home/rschaeff/data/scop/v1.75/scop.v1.75.xml';
#$REF_XML{'scop169'}	= '/home/rschaeff/data/scop/v1.69/scop.v1.69.xml';
#$REF_XML{'develop2'}    = '/home/rschaeff/data/ecod/v2/ecod.develop2.xml';
#$REF_XML{'develop8'}    = '/home/rschaeff/data/ecod/v8/ecod.develop8.xml';
#$REF_XML{'develop10'}   = '/home/rschaeff/data/ecod/v10/ecod.develop10.xml';
#$REF_XML{'develop12'}   = '/home/rschaeff/data/ecod/v12/ecod.develop12.xml';
#$REF_XML{'develop13'}   = '/home/rschaeff/data/ecod/v13a/ecod.develop13a.xml';
#$REF_XML{'develop14'}   = '/home/rschaeff/data/ecod/v14/ecod.develop14.xml';
#$REF_XML{'develop18'}   = '/home/rschaeff/data/ecod/v18/ecod.develop18c.xml';
#$REF_XML{'develop19'}   = '/home/rschaeff/data/ecod/v19/ecod.develop19.xml';
#$REF_XML{'develop21'}   = '/home/rschaeff/data/ecod/v21/ecod.develop21.xml';  
#$REF_XML{'develop24'}   = '/home/rschaeff/data/ecod/v24/ecod.develop24.xml';  
#$REF_XML{'develop25'} 	= '/home/rschaeff/data/ecod/v25/ecod.develop25e.xml';
#$REF_XML{'develop26'} 	= '/home/rschaeff/data/ecod/v26/ecod.develop26b.xml';
#$REF_XML{'develop27'} 	= '/home/rschaeff/data/ecod/v27/ecod.develop27g.xml';
#$REF_XML{'develop28'} 	= '/home/rschaeff/data/ecod/v28/ecod.develop28b.xml';
#$REF_XML{'develop29'} 	= '/home/rschaeff/data/ecod/v29/ecod.develop29.xml';
#$REF_XML{'develop29.reps_only'} 	= '/home/rschaeff/data/ecod/v29/ecod.develop29.xml';
#$REF_XML{'develop30'} 	= '/home/rschaeff/data/ecod/v30/ecod.develop30.xml';
#$REF_XML{'develop30.reps_only'} 	= '/home/rschaeff/data/ecod/v30/ecod.develop30.xml';
#$REF_XML{'develop31'} 	= '/home/rschaeff/data/ecod/v31/ecod.develop31a.xml';
#$REF_XML{'develop32'} 	= '/home/rschaeff/data/ecod/v32/ecod.develop32.xml';
#$REF_XML{'develop33'} 	= '/home/rschaeff/data/ecod/v33/ecod.develop33.xml';
#
#$REF_RANGE_CACHE{'develop29.reps_only'} = '/home/rschaeff/data/ecod/v29/ecod.develop29.manual_ranges.txt';
#$REF_RANGE_CACHE{'develop29'} = '/home/rschaeff/data/ecod/v29/ecod.develop29.manual_ranges.txt';
#$REF_RANGE_CACHE{'develop30.reps_only'} = '/home/rschaeff/data/ecod/v30/ecod.develop30.manual_ranges.txt';
#$REF_RANGE_CACHE{'develop30'} = '/home/rschaeff/data/ecod/v30/ecod.develop30.manual_ranges.txt';
#$REF_RANGE_CACHE{'develop31'} = '/home/rschaeff/data/ecod/v31/ecod.develop31a.range_cache.txt';
#$REF_RANGE_CACHE{'develop32'} = '/home/rschaeff/data/ecod/v32/ecod.develop32.range_cache.txt';
#$REF_RANGE_CACHE{'develop33'} = '/home/rschaeff/data/ecod/v33/ecod.develop33.range_cache.txt';
#$REF_RANGE_CACHE{'scop175'} = '/home/rschaeff/data/scop/v1.75/scop.v1.75.range_cache.txt';
#$REF_RANGE_CACHE{'scop169'} = '/home/rschaeff/data/scop/v1.69/scop.v1.69.range_cache.txt';
#
1;
