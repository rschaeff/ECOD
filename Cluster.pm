package ECOD::Cluster;
require Exporter;

use warnings;
use strict;

our @ISA = qw(Exporter);

our @EXPORT = ( 
				"&blastclust_tag_domains",
				"&cdhit_tag_domains",
				"&assemble_pf_hhm_db",
				"&build_rep_hhms",
				"&generate_pf_group_blastclust_jobs",
				"&generate_pf_group_cdhit_jobs",
				"&generate_hhsearch_jobs_for_domain_reps",
				"&generate_hh_matrix_jobs",
				"&generate_domain_parse_jobs",
				"&generate_pf_uids_inx",
				"&select_domain_representatives",
				"&read_reference_cache",
				"&cluster_strip",
				);

use ECOD::Update;
use XML::Grishin;

use Carp;
use XML::LibXML;
use File::Path qw(make_path);
use File::Copy;
use File::Remove qw(remove);
use Getopt::Long;

my $HHBLITS_EXE  = '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhblits';
my $HHSEARCH_EXE = '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhsearch';
my $HHMAKE_EXE   = '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhmake';
my $HH_MAT_EXE   = '/data/ecod/database_versions/bin/hh_matrix_build.pl';
my $UNIPROT_LIB  = '/home/rschaeff_1/side/wtf_hh/uniprot20_2012_10_klust20_dc_2012_12_10/uniprot20_2012_10_klust20_dc_2012_12_10';
my $NR_LIB 	     = '/home/rschaeff_1/side/wtf_hh/nr20_12Aug11';
my $HH_BLITS_CPU 		= 8;
my $PSIPRED_DIR 		= '/usr1/psipred/bin';
my $PSIPRED_DATA_DIR	= '/usr1/psipred/data';

my $mode         = "F_ANALYZE";

my $DOMAIN_DATA_DIR = '/data/ecod/domain_data';
if ( !-d $DOMAIN_DATA_DIR ) { die }

my $DOMAIN_PARSE_SCRIPT = '/data/ecod/database_versions/bin/single_parse_domain_rep_hhsearch.pl';
if ( !-f $DOMAIN_PARSE_SCRIPT ) { die }

sub has_cluster_tag {
    $_[0]->exists(qq{cluster[\@level="$_[1]"][\@method="$_[2]"]});
}

sub create_cluster_tag {
    my ( $xml_doc, $d_node, $tag, $tag_id, $method ) = @_;
    my $c_node = $xml_doc->createElement('cluster');
    $c_node->setAttribute( 'method', $method );
    $c_node->setAttribute( 'level',  $tag );
    $c_node->setAttribute( 'id',     $tag_id );
    $d_node->appendChild($c_node);
}

sub build_rep_hhms { 
	my ($ecod_xml_doc, $tag, $method) = @_;
	my @cmds;
	foreach my $domain_node (find_domain_reps($ecod_xml_doc, $tag, $method)) { 
		my ($uid, $ecod_domain_id) = get_ids($domain_node);
		my $short_uid = substr($uid, 2, 5);
		if (!-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.hhm" && 
			 -f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.fa") { 
			push(@cmds, hh_build($uid));
		}
	}
	return \@cmds;
}

sub hh_build { 
	my ($uid) = @_;
	my $short_uid = substr($uid, 2, 5);

	my $uid_path = "$DOMAIN_DATA_DIR/$short_uid/$uid";
	my $fa_path  = "$uid_path/$uid.fa";
	my $a3m_path = "$uid_path/$uid.a3m";
	my $hhm_path = "$uid_path/$uid.hhm";

	my $cmd = "$HHBLITS_EXE -i $fa_path -d $NR_LIB -oa3m $a3m_path -ohhm $hhm_path ".
			  "-n 3 -cpu 8 -addss -psipred $PSIPRED_DIR -psipred_data $PSIPRED_DATA_DIR";

	return $cmd;

}

sub select_domain_representatives { 
	my ($ecod_xml_doc, $tag, $pf_fasta_dir, $prev, $method) = @_;
	my %clusters;
	my $missing = 0;
	my %cluster_has_rep;
	foreach my $domain_node (find_domain_nodes($ecod_xml_doc)) { 
		my ($uid, $ecod_domain_id) 	= get_ids($domain_node);
		if (has_cluster_tag($domain_node, $tag, $method)) { 
			my $cluster_id		= get_cluster_id($domain_node, $tag, $method);
			$cluster_has_rep{$cluster_id}++ if is_domain_rep($domain_node, $tag, $method);
			push (@{$clusters{$cluster_id}}, $domain_node);	
		}else{
			$missing++;
		}
	}
	my $prev_reps;
	if ($prev) { 
		my $prev_ecod_xml_doc = xml_open($prev);
		$prev_reps = get_prev_reps($prev_ecod_xml_doc, $tag, $method);
	}

	printf "#select_domain_representatives: Found %i clusters, %i domains not clustered on this tag\n", scalar keys %clusters, $missing;
	#Select domain reps
	foreach my $cluster_id (keys %clusters) { 
		next if $cluster_has_rep{$cluster_id};	

		my @domains;
		foreach my $domain_node (@{$clusters{$cluster_id}}) { 

			my %domain;
			my ($uid, $ecod_domain_id) = get_ids($domain_node);
			my $isXray			= 0;
			my $xray_res		= '99.99';

			my $isPrevRep 		= exists $$prev_reps{$uid} ? 1 : 0;
			my $isManualRep 	= is_manual_domain_rep($domain_node);
			my $isProvisionalRep = is_provisional_domain_rep($domain_node);
			my $isCurrentRep 	= is_domain_rep($domain_node, $tag, $method);
			#print "man: $isManualRep prov: $isProvisionalRep cure: $isCurrentRep\n";
			if (!defined $isCurrentRep) { die "curr undefined: $uid\n" }
			if (!defined $isProvisionalRep) { die "prov undefined; $uid\n" }
			if (!defined $isManualRep) { die "man undefined: $uid\n"}
			
			if ($domain_node->findvalue('structure/@method') eq 'X-RAY DIFFRACTION') { 
				$isXray = 1;
				if ($domain_node->exists('structure/@resolution')) { 
					$xray_res = $domain_node->findvalue('structure/@resolution');
				}
			}
			my $deposition_date = '9999-99-99';
			if ($domain_node->exists('structure/@deposition_date') && $domain_node->findvalue('structure/@deposition_data') =~ /\d{4}\-\d+\-\d+/) { 
				$deposition_date = $domain_node->findvalue('structure/@deposition_date');
			}
			my $cluster_node = $domain_node->findnodes(qq{cluster[\@level="$tag"][\@method="$method"]})->get_node(1);

			$deposition_date =~ s/\-//g;

			$domain{ecod_domain_id}	= $ecod_domain_id;
			$domain{uid}			= $uid;
			$domain{isManualRep} 	= $isManualRep;
			$domain{isCurrentRep} 	= $isCurrentRep;
			$domain{isProvRep}		= $isProvisionalRep;
			$domain{isPrevRep}		= $isPrevRep;
			$domain{isXray}			= $isXray;
			$domain{xray_res}		= $xray_res;
			$domain{deposition_date}	= $deposition_date;
			$domain{node}			= $cluster_node;

			push (@domains, \%domain);
			
		}
		@domains = sort { 
				$b->{isCurrentRep} <=> $a->{isCurrentRep} || 
				$b->{isManualRep} <=> $a->{isManualRep} ||
				$b->{isProvRep} <=> $a->{isProvrep} ||
				$b->{isPrevRep} <=> $a->{isPrevRep} ||
				$b->{isXray} <=> $a->{isXray} ||
				$a->{xray_res} <=> $b->{xray_res} || 
				$b->{deposition_date} <=> $a->{deposition_date}
		} @domains;
		$domains[0]{node}->setAttribute('domain_rep', 'true');
		
	}

}

sub get_prev_reps { 
	my $ecod_xml_doc = $_[0];
	my $tag 		= $_[1];
	my $method		= $_[2];
	my %reps;
	foreach my $domain (find_domain_nodes($ecod_xml_doc))  { 
		my ($uid, $ecod_domain_id) = get_ids($domain);
		if (is_domain_rep($domain, $tag, $method)) { 
			$reps{$uid}++;
		}
	}
	return \%reps;
}


sub blastclust_tag_domains {

    my ( $ecod_xml_doc, $tag, $pf_fasta_dir ) = @_;

    my $clust_i = 1;
	my $method ='blastclust';

    foreach my $pf_node ( find_nodes( $ecod_xml_doc, "//pf_group" ) ) {

        my $pf_id     = get_pf_id($pf_node);
        my $ecodf_acc = get_ecodf_acc($pf_node);

        my $path        = "$pf_fasta_dir/$pf_id";
        my $pf_clust_fn = "$path/$pf_id.$tag.clust";

        open my $fh, "<", $pf_clust_fn
          or die "ERROR! Could not open $pf_clust_fn for reading:$!\n";
        my %tags;
        while ( my $ln = <$fh> ) {
            next if $ln =~ /^#/;
            my @F = split /\s+/, $ln;
            foreach my $f (@F) { $tags{$f} = $clust_i }
            $clust_i++;
        }
        close $fh;

        foreach my $domain ( find_domain_nodes($pf_node) ) {
            my ( $uid, $ecod_domain_id ) = get_ids($domain);

            if ( has_cluster_tag( $domain, $tag, $method ) ) {
                foreach my $cluster_node ( find_nodes( $domain, "cluster[\@level='$tag'][\@method='$method']" ) ) {
                    $cluster_node->unbindNode;
                }
            }

            if ( $tags{$uid} ) {
                create_cluster_tag( $ecod_xml_doc, $domain, $tag, $tags{$uid}, $method );
            }
        }
    }
}

sub cdhit_tag_domains { 
	my ($ecod_xml_doc, $tag, $pf_fasta_dir ) = @_;
	my $clust_i = 1;
	my $method = 'cdhit';
	foreach my $pf_node (find_nodes( $ecod_xml_doc, '//pf_group') ) { 
		next unless $pf_node->findnodes('.//domain')->size() > 0;
        my $pf_id     = get_pf_id($pf_node);
        my $ecodf_acc = get_ecodf_acc($pf_node);

        my $path        = "$pf_fasta_dir/$pf_id";
        my $pf_clust_fn = "$path/$pf_id.$tag.cdhit.clust.clstr";

        open my $fh, "<", $pf_clust_fn
          or die "ERROR! Could not open $pf_clust_fn for reading:$!\n";
        my %tags;
        while ( my $ln = <$fh> ) {
			if ($ln =~ /^\>Cluster\s+(\d+)/) { 
				$clust_i++;
			}elsif ($ln =~ /^\d+/) {
				my @F = split( /\s+/, $ln);
				my $n 	= $F[0];
				my $aa 	= $F[1];
				my $uid = $F[2];
				$uid =~ s/>//;
				$uid =~ s/\.\.\.//;
				$tags{$uid} = $clust_i;
			}else{
				warn "Eh?: $ln\n";
			}	
        }
        close $fh;

        foreach my $domain ( find_domain_nodes($pf_node) ) {
            my ( $uid, $ecod_domain_id ) = get_ids($domain);

            if ( has_cluster_tag( $domain, $tag, $method ) ) {
                foreach my $cluster_node ( find_nodes( $domain, "cluster[\@level='$tag'][\@method='$method']" ) ) {
                    $cluster_node->unbindNode;
                }
            }

            if ( $tags{$uid} ) {
                create_cluster_tag( $ecod_xml_doc, $domain, $tag, $tags{$uid}, $method );
            }
        }
	}
}

sub build_hh_matrix {
    my ( $ecod_pf_uid_summ_fn, $pf_id, $tag, $pf_fasta_fn ) = @_;

    my $pf_id_uids = read_pf_uids($ecod_pf_uid_summ_fn);

    my %domain_reps;
    my %mat;
    my %reps;

    if ( exists $$pf_id_uids{$pf_id} ) {
        foreach my $uid ( @{ $$pf_id_uids{$pf_id} } ) {
            $reps{$uid}++;
            my $path = "$pf_fasta_fn/$pf_id/$uid.$tag.hh.xml";
            if ( -f $path ) {
                parse_hh_xml( $uid, $path, \%mat );
            }
            else {
                print "WARNING! Missing $path\n";
            }
        }
    }
    else {
        die "ERROR! pf_id not found in cache: $pf_id\n";
    }

    open( my $out_fh, ">", "$pf_fasta_fn/$pf_id/$pf_id.$tag.hh.mat" )
      or die "ERROR! On write for $pf_id $tag\n";

    foreach my $uid1 ( keys %reps ) {
        foreach my $uid2 ( keys %reps ) {
            if ( $uid1 eq $uid2 ) { next }

            if ( $uid1 < $uid2 ) {
                if ( $mat{$uid1}{$uid2} ) {
                    printf $out_fh "%09i\t%09i\t%.2f\t%.2f\n", $uid1, $uid2,
                      $mat{$uid1}{$uid2}{hh_prob}, $mat{$uid1}{$uid2}{coverage};
                }
                else {
                    printf $out_fh "%09i\t%09i\t%.2f\t%.2f\n", $uid1, $uid2, 0, 0;
                }
            }
        }
    }
    close $out_fh;

}

sub generate_pf_uids_inx { 
	my ($ecod_xml_doc, $tag, $out_fn, $method) = @_;
	open my $out_fh, ">", $out_fn or die "ERROR! Couldn't open $out_fn for writing:$!\n";
	foreach my $pf_group (find_pf_group_nodes($ecod_xml_doc)) { 
		my $pf_id = get_pf_id($pf_group);
		#Is this working correctly?
		foreach my $domain (find_domain_reps($pf_group, $tag, $method)) { 
			my ($uid, $ecod_domain_id) = get_ids($domain);
			print $out_fh "$pf_id\t$uid\n";
		}
	}
}

sub generate_hh_matrix_jobs {

    my ( $ecod_xml_doc, $tag, $pf_uids_inx_fn, $pf_fasta_dir, $method ) = @_;
    my ( %reps_analyzed, %total_reps );
    my @cmds;
	my $local_force_overwrite = 1;
	foreach my $pf_group (find_pf_group_nodes($ecod_xml_doc)) { 
		my $pf_id = get_pf_id($pf_group);

		foreach my $domain_rep ( find_domain_reps( $pf_group, $tag, $method ) ) {
			my ($uid, $ecod_domain_id) = get_ids($domain_rep);
			#print "$pf_id $uid\n";
			my $result_fn = "$uid.$tag.hhr";

			if ( -f "$pf_fasta_dir/$pf_id/$result_fn" ) {
				$reps_analyzed{$pf_id}++;
			}
			
			$total_reps{$pf_id}++;
		}
	}

    foreach my $pf_id ( keys %reps_analyzed ) {
        if ( $total_reps{$pf_id} == $reps_analyzed{$pf_id} ) {
            my $mat_fn = "$pf_id.$tag.hh.mat";
            if ( !-f "$pf_fasta_dir/$pf_id/$mat_fn" || $local_force_overwrite ) {
                my $cmd = "$HH_MAT_EXE $pf_uids_inx_fn $pf_id $tag";
                push( @cmds, $cmd );
            }
        }
    }
    return \@cmds;
}

sub parse_hh_xml {
    my $sub = 'parse_hh_xml';
    my ( $query_uid, $fn, $mat_href ) = @_;

    my $hh_xml_doc = xml_open($fn);

    foreach my $hh_hit_node ( $hh_xml_doc->findnodes('//hh_hit') ) {
        my $hit_uid   = $hh_hit_node->findvalue('@uid');
        my $hh_prob   = $hh_hit_node->findvalue('@hh_prob');
        my $hit_cover = sprintf "%.2f", $hh_hit_node->findvalue('template_seqid_range/@coverage');
        if ( $query_uid < $hit_uid ) {
            if ( !$$mat_href{$query_uid}{$hit_uid} ) {
                $$mat_href{$query_uid}{$hit_uid}{hh_prob}  = $hh_prob;
                $$mat_href{$query_uid}{$hit_uid}{coverage} = $hit_cover;
            }
        }
        else {
            if ( !$$mat_href{$hit_uid}{$query_uid} ) {
                $$mat_href{$hit_uid}{$query_uid}{hh_prob}  = $hh_prob;
                $$mat_href{$hit_uid}{$query_uid}{coverage} = $hit_cover;
            }
        }
    }
}
sub generate_pf_group_cdhit_jobs { 
	my ($ecod_xml_doc, $tag, $pf_fasta_dir, $global_opt) = @_;

    my ( $L, $S, $n );
    if ( $tag eq "F99" ) {
        $L = 0.9;
        $S = 0.99;
		$n = 5;
    }
    elsif ( $tag eq 'F40' ) {
        $L = 0.7;
        $S = 0.40;
		$n = 2;
    }
	elsif ($tag eq 'F70') { 
		$L = 0.7;
		$S = 0.70;
		$n = 4;
	}
    else {
        die "ERROR! tag $tag unknown\n";
    }

	my @cmds;
    foreach my $pf_node ( find_nodes( $ecod_xml_doc, '//pf_group' ) ) {
		next if $pf_node->findnodes('.//domain')->size() == 0;

       # my ( $pf_id, $ecodf_acc ) = get_pf_id($pf_node);
	   my $pf_id = get_pf_id($pf_node);
	   print "#$pf_id\n";

        my $path = "$pf_fasta_dir/$pf_id";
        if ( !-d $path ) {
            if ( make_path($path) == 0 ) {
                croak "ERROR! make_path failed on $path:$!\n";
            }
        }

        #my $pf_fasta_fn = "$path/$pf_id.$tag.fa"; #Why tagged?
        my $pf_fasta_fn = "$path/$pf_id.fa";
        if ( -f $pf_fasta_fn ) {
            if ( $$global_opt{verify} ) {
                unless ( verify( $pf_node, $pf_fasta_fn ) ) {
                    remove($pf_fasta_fn);
                    gen_pf_fasta( $pf_node, $pf_fasta_fn );
                }
            }
            elsif ( $$global_opt{force_overwrite} ) {
                remove($pf_fasta_fn);
                gen_pf_fasta( $pf_node, $pf_fasta_fn );
            }
            else {
                #warn "WARNING! $pf_fasta_fn exists, skipping...\n";
            }
        }
        else {
            gen_pf_fasta( $pf_node, $pf_fasta_fn );
        }
        ( my $pf_clust_fn = $pf_fasta_fn ) =~ s/fa$/$tag.cdhit.clust/;
        if ( !-f $pf_clust_fn || $$global_opt{force_overwrite} ) {
            my $cmd = gen_cdhit_job( $pf_fasta_fn, $pf_clust_fn, $S, $L, $n );
			push (@cmds, $cmd);
        }
    }
	return \@cmds;
}




sub generate_pf_group_blastclust_jobs {

    my ( $ecod_xml_doc, $tag, $pf_fasta_dir, $global_opt ) = @_;

    my ( $L, $S );
    if ( $tag eq "F99" ) {
        $L = 0.9;
        $S = 99;
    }
    elsif ( $tag eq 'F40' ) {
        $L = 0.4;
        $S = 40;
    }
	elsif ($tag eq 'F70') { 
		$L = 0.7;
		$S = 70;
	}
    else {
        die "ERROR! tag $tag unknown\n";
    }

	my @cmds;
    foreach my $pf_node ( find_nodes( $ecod_xml_doc, '//pf_group' ) ) {

       # my ( $pf_id, $ecodf_acc ) = get_pf_id($pf_node);
	   my $pf_id = get_pf_id($pf_node);
	   print "#$pf_id\n";

        my $path = "$pf_fasta_dir/$pf_id";
        if ( !-d $path ) {
            if ( make_path($path) == 0 ) {
                croak "ERROR! make_path failed on $path:$!\n";
            }
        }

        #my $pf_fasta_fn = "$path/$pf_id.$tag.fa"; #Why tagged?
        my $pf_fasta_fn = "$path/$pf_id.fa";
        if ( -f $pf_fasta_fn ) {
            if ( $$global_opt{verify} ) {
                unless ( verify( $pf_node, $pf_fasta_fn ) ) {
                    remove($pf_fasta_fn);
                    gen_pf_fasta( $pf_node, $pf_fasta_fn );
                }
            }
            elsif ( $$global_opt{force_overwrite} ) {
                remove($pf_fasta_fn);
                gen_pf_fasta( $pf_node, $pf_fasta_fn );
            }
            else {
                #warn "WARNING! $pf_fasta_fn exists, skipping...\n";
            }
        }
        else {
            gen_pf_fasta( $pf_node, $pf_fasta_fn );
        }
        ( my $pf_clust_fn = $pf_fasta_fn ) =~ s/fa$/$tag.clust/;
        if ( !-f $pf_clust_fn || $$global_opt{force_overwrite} ) {
            my $cmd = gen_blastclust_job( $pf_fasta_fn, $pf_clust_fn, $S, $L );
			push (@cmds, $cmd);
        }
    }
	return \@cmds;
}

sub gen_blastclust_job { 
	my ($in, $out, $S, $L) = @_;
	my $BLASTCLUST = '/data/usr1/seals/bin/blastclust';
	return "$BLASTCLUST -i $in -o $out -S $S -L $L -e F -b T -p T"; 
}

sub gen_cdhit_job { 
	my ($in, $out, $S, $L, $n) = @_;
	my $CDHIT = '/usr1/cd-hit-v4.5.6/cd-hit';
	!-f $CDHIT and die "ERROR! CD-hit executable not found: $CDHIT\n";
	return "$CDHIT -i $in -o $out -aL $L -c $S -n $n -G 0";
}

sub generate_domain_parse_jobs {
    my ( $ecod_xml_doc, $tag, $pf_fasta_dir, $reference, $global_opt, $method ) = @_;

	my @cmds;
    foreach my $domain_rep ( find_domain_reps( $ecod_xml_doc, $tag, $method ) ) {

        my $pf_id = get_pf_id($domain_rep);

        my ( $uid, $ecod_domain_id ) = get_ids($domain_rep);
        my $short_uid = substr( $uid, 2, 5 );

        my $hhm_file = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.hhm";

        if ( !-f $hhm_file ) {
            print "WARNING! $hhm_file not found, skipping...\n";
            next;
        }

        my $result_fn     = "$uid.$tag.hhr";
        my $result_xml_fn = "$uid.$tag.hh.xml";

        my $pf_ref_db = "$pf_fasta_dir/$pf_id/$pf_id.$tag.hhm_db";
        if ( !-f $pf_ref_db ) {
            next;
        }
        if ( -f "$pf_fasta_dir/$pf_id/$result_fn" ) {
            if ( !-f "$pf_fasta_dir/$pf_id/$result_xml_fn"
                || $$global_opt{force_overwrite} )
            {
				my $cmd = $$global_opt{force_overwrite} ? 
				"$DOMAIN_PARSE_SCRIPT $pf_id $uid $reference $tag --force_overwrite" : 
				"$DOMAIN_PARSE_SCRIPT $pf_id $uid $reference $tag";
				push (@cmds, $cmd);
            }
        } else {
            print "WARNING! Missing hhr for $pf_id $uid\n";
        }
    }
	return \@cmds;
}

sub hh_parse {
    my $sub = 'hh_parse';
    my ( $hhr_in_fn, $xml_out_fn, $reference, $chain_href, $pdb_href, $global_opt ) = @_;

    open my $fh, "<", $hhr_in_fn
      or die "ERROR! $sub: Could not open $hhr_in_fn for reading:$!\n";

    my $input = 'struct_seqid';

    my ( $hh_xml_doc, $hh_root_node ) = xml_create('hh_summ_doc');

    my $hh_hit_list_node = $hh_xml_doc->createElement('hh_hit_list');
    $hh_root_node->appendChild($hh_hit_list_node);

    my ( $uid_href, $uid_lookup ) = read_reference_cache($reference);

    my ( $query, $short_query, $match_columns, $matched_seqs, $total_seqs, $neff, $searched_hmms, $date, $command );

    my ( $query_seqid_aref, $query_struct_seqid_aref, $query_asym_id,
        $query_chain_seqid_aref, $query_chain_struct_seqid_aref );

    my %query;

    my %template_ends;
    my @aligned_cols;

    while ( my $ln = <$fh> ) {
        if ( $ln =~ /^Query\s+(.*)/ ) {
            $query{long_name} = $1;

            my ( $ecod_domain_id, $query_uid );
            if ( $query =~ /(e\w{4}[\w\.]+\d+)\s+.*(\d{9})$/ ) {    #Does it look like an ecod_domain_id and a uid?
                $ecod_domain_id = $1;
                $query{uid} = $2;
            }
            else {
                die "ERROR! $query doesn't regexp as a uid and ecod_domain_id\n";
            }

            my $pdb       = $$pdb_href{ $query{uid} };
            my $chain_str = $$chain_href{ $query{uid} };

            my @chain = split( '', $chain_str );

            $query{long_name} =~ /([\w\,\.]+_)/;    #What is this being used for?
            $query{short_name} = $1;                #CANDIDATE TO DATA STRUCTURE

            (
                $query{seqid_aref}, $query{struct_seqid_aref},
                $query{pdbnum_href}, $query{asym_id},
                $query{chain_seqid_aref},
                $query{chain_struct_seqid_aref}
            ) = pdbml_mc_seq_parse( $pdb, \@chain );

            if ( !$query{seqid_aref} ) {
                print "ERROR! $sub: $pdb did not parse, skipping...\n";
                return 0;
            }

            $hh_root_node->appendTextChild( 'hh_query', $query{long_name} );
        }
        elsif ( $ln =~ /^Match_columns\s+(\d+)/ ) {

            $query{match_columns} = $1;    #CANDIDATE_TO_DATA_STRUCTURE

            $hh_root_node->appendTextChild( 'hh_match_cols', $query{match_columns} );
        }
        elsif ( $ln =~ /^No_of_seqs\s+(\d+)\sout\sof\s(\d+)/ ) {
            $query{matched_seqs} = $1;
            $query{total_seqs}   = $2;

            my $hh_matched_seqs_node = $hh_xml_doc->createElement('hh_match_seqs');
            $hh_matched_seqs_node->setAttribute( 'matched', $query{matched_seqs} );
            $hh_matched_seqs_node->setAttribute( 'total',   $query{total_seqs} );
            $hh_root_node->appendChild($hh_matched_seqs_node);
        }
        elsif ( $ln =~ /^Neff\s+([\d\.]+)/ ) {
            $query{neff} = $1;
            $hh_root_node->appendTextChild( 'hh_neff', $query{neff} );
        }
        elsif ( $ln =~ /^Searched_HMMs\s+(\d+)/ ) {
            $query{searched_hmms} = $1;
            $hh_root_node->appendTextChild( 'searched_hhms', $query{searched_hmms} );
        }
        elsif ( $ln =~ /^Date\s+(.*)/ ) {
            $query{hh_date} = $1;
            $hh_root_node->appendTextChild( 'hh_date', $query{hh_date} );
        }
        elsif ( $ln =~ /Command\s+(\/.*)/ ) {
            $query{command} = $1;
            $hh_root_node->appendTextChild( 'hh_command', $query{command} );
        }

        if ( $ln =~ /^\s+No\s+Hit\s+Prob/ ) {
            my $hh_hit_summ_list_node = $hh_xml_doc->createElement('hh_hit_summ_list');

            while ( my $ln = <IN> ) {
                unless ( $ln =~ /\d+\s+\w{4}/ ) { last }
                $ln =~ s/^\s+//;
                my @F = split( /\s+/, $ln );
                my @G = split( /\s+/, substr( $ln, 36, 94 ) );

                my $hh_hit_summ_node = $hh_xml_doc->createElement('hh_hit_summ');

                my $hit_num = $F[0];
                $hh_hit_summ_node->setAttribute( 'hit_num', $hit_num );

                my $hit_ecod_domain_id = $F[1];
                $hh_hit_summ_node->setAttribute( 'ecod_domain_id', $hit_ecod_domain_id );
                my $hh_prob = $G[0];
                $hh_hit_summ_node->setAttribute( 'hh_prob', $hh_prob );
                my $hh_eval = $G[1];
                $hh_hit_summ_node->setAttribute( 'hh_eval', $hh_eval );
                my $hh_pval = $G[2];
                $hh_hit_summ_node->setAttribute( 'hh_pval', $hh_pval );
                my $hh_score = $G[3];
                $hh_hit_summ_node->setAttribute( 'hh_score', $hh_score );
                my $hh_ss = $G[4];
                $hh_hit_summ_node->setAttribute( 'hh_ss', $hh_ss );
                my $hh_cols = $G[5];
                $hh_hit_summ_node->setAttribute( 'hh_cols', $hh_cols );
                $aligned_cols[$hit_num] = $hh_cols;
                my $query_hmm_range = $G[6];
                $hh_hit_summ_node->appendTextChild( 'hh_query_hmm_range', $query_hmm_range );

                my $temp_hmm_range = $G[7];
                $temp_hmm_range =~ /-(\d+)/;
                my $end = $1;
                $template_ends{$hit_num}{$hit_ecod_domain_id} = $end;
                $hh_hit_summ_node->appendTextChild( 'hh_template_hmm_range', $temp_hmm_range );

                my $temp_hmm_num = $G[8];
                $hh_hit_summ_node->appendTextChild( 'hh_temp_hmm_num', $temp_hmm_num );

                $hh_hit_summ_list_node->appendChild($hh_hit_summ_node);
            }
        }

        my $last_template_block_seen = 0;
        if ( $ln =~ /^No (\d+)\s+$/ ) {
            my $hit_num = $1;
            my $uid;
            my $ecod_domain_id;
            my $ecod_pdb;
            my ( $hh_prob, $hh_eval, $hh_score, $aligned_cols, $idents, $sum_probs, $similarities );
            my $observed_cols = 0;
            my $query_regexp  = qr/$query{short_name}/;
            my $ecod_domain_id_regexp;
            my ( $template_ln, $query_ln );
            my ( $query_start, $template_start, $query_end, $template_end );

            if ( !$aligned_cols[$hit_num] ) {
                die "No aligned_cols found for $hit_num $aligned_cols[$hit_num]\n";
            }
            if ( $aligned_cols[$hit_num] < 10 ) {
                while ( my $ln = <IN> ) {
                    if ( $ln =~ /T ss_conf/ ) {
                        last;
                    }
                }
            }
            else {

                my $hh_hit_node = $hh_xml_doc->createElement('hh_hit');
                $hh_hit_list_node->appendChild($hh_hit_node);

                while ( my $ln = <IN> ) {
                    if ( $ln =~ /^\s+$/ && $last_template_block_seen ) {

                        my @query_align_chars    = split( '', $query_ln );
                        my @template_align_chars = split( '', $template_ln );

                        my @query_aligned_positions;
                        my @template_aligned_positions;

                        my $query_pos    = $query_start;
                        my $template_pos = $template_start;

                        for ( my $i = 0 ; $i < scalar(@query_align_chars) ; $i++ ) {
                            if ( $query_align_chars[$i] eq '-' ) {
                                $template_pos++;
                            }
                            elsif ( $template_align_chars[$i] eq '-' ) {
                                $query_pos++;
                            }
                            elsif ($template_align_chars[$i] =~ /[A-Z]/
                                && $query_align_chars[$i] =~ /[A-Z]/ )
                            {
                                push( @query_aligned_positions,    $query_pos );
                                push( @template_aligned_positions, $template_pos );
                                $template_pos++;
                                $query_pos++;
                            }
                        }
                        my $template_range = 'NA';
                        if ( scalar(@template_aligned_positions) > 2 ) {
                            $template_range = rangify(@template_aligned_positions);
                        }
                        my $query_range = 'NA';
                        if ( scalar(@query_aligned_positions) > 2 ) {
                            $query_range = rangify(@query_aligned_positions);
                        }

             #3/21/2012 -- modification from struct_seqid to seqid.
             #Whether you use struct_seqid or seqid as FA input, buildali.pl will expand to full sequence.
             #2/26/2013
             #But does hhblits?
             #my $query_seqid_aref = struct_region(\@query_aligned_positions, $query_seqid_aref);
             #multi_chain_struct_region -- returns 3 values, correct use of pos array. What's going on, profile change??
                        my ( $query_structregion_seqid_aref, $query_structregion_chain_seqid_aref, $warning );
                        if ( $input eq 'seqid' ) {
                            ( $query_structregion_seqid_aref, $query_structregion_chain_seqid_aref, $warning ) =
                              multi_chain_struct_region( \@query_aligned_positions, $query{seqid_aref},
                                $query{chain_seqid_aref} );
                        }
                        else {
                            ( $query_structregion_seqid_aref, $query_structregion_chain_seqid_aref, $warning ) =
                              multi_chain_struct_region(
                                \@query_aligned_positions,
                                $query{struct_seqid_aref},
                                $query{chain_struct_seqid_aref}
                              );
                        }

                        if ( $query_structregion_seqid_aref == 0 ) {
                            printf "Fail on mc struct region: %i %i %i\n",
                              scalar(@query_aligned_positions),
                              scalar(@$query_seqid_aref),
                              scalar(@$query_chain_seqid_aref);
                            die;
                        }
                        my $query_seqid_range = 'NA';
                        if ( scalar(@$query_structregion_seqid_aref) > 2 ) {
                            $query_seqid_range = multi_chain_rangify( $query_structregion_seqid_aref,
                                $query_structregion_chain_seqid_aref );
                        }
                        my $query_pdb_range = multi_chain_pdb_rangify( $query_structregion_seqid_aref,
                            $query{pdbnum_href}, $query_structregion_chain_seqid_aref );

                        my $ecod_range;

                        my @ecod_chains;
                        my @ecod_segs;
                        my $imp_segs = 0;

                        #TEMPLATE RANGE MAY BE MULTI CHAIN
                        if ( $$uid{$uid}{seqid_range} ) {
                            my @ecod_chunks =
                              split( /,/, $$uid{$uid}{seqid_range} );
                            foreach my $chunk (@ecod_chunks) {
                                if ( $chunk =~ /(\w+):\s*$/ ) {
                                    my $ecod_chain = $1;
                                    push( @ecod_chains, $ecod_chain );
                                    push( @ecod_segs,   'FULL' );
                                    $imp_segs++;
                                }
                                elsif ( $chunk =~ /(\w+):(\-?\d+\-\d+)/ ) {
                                    my $ecod_chain = $1;
                                    my $ecod_seg   = $2;
                                    push( @ecod_chains, $ecod_chain );
                                    push( @ecod_segs,   $ecod_seg );
                                }
                                else {
                                    print "WARNING! $sub: Bad range chunk $chunk\n";
                                    next;
                                }
                            }
                        }
                        else {
                            print "$ecod_domain_id deprecated, skipping...\n";
                            $hh_hit_node->setAttribute( 'structure_obsolete', 'true' );
                            last;
                        }
                        my $template_struct_seqid_range;
                        my $template_struct_pdb_range;
                        my $template_coverage          = 'Unk';
                        my $ungapped_template_coverage = 'Unk';
                        if ( $imp_segs == scalar(@ecod_segs) ) {

                            #print "DEBUG $sub: case 1\n";
                            my ( $template_seqid_aref, $template_struct_seqid_aref, $ecod_pdbnum_href, $asym_id,
                                $template_chain_aref )
                              = pdbml_mc_seq_parse( $ecod_pdb, \@ecod_chains );
                            if ( !$template_seqid_aref ) {
                                print "WARNING! $ecod_pdb obsolete!\n";
                                last;
                            }

                            #1.Change pos-inx into seqid inx
                            my $template_range_aref = range_expand($template_range);    #alignment pos-inx
                            my ( $template_aligned_seqid_range_aref, $template_aligned_chain_aref, $warning ) =
                              multi_chain_struct_region( $template_range_aref, $template_struct_seqid_aref,
                                $template_chain_aref );

                            #Output template struct_seqid and struct_pdb ranges
                            $template_struct_seqid_range =
                              multi_chain_rangify( $template_aligned_seqid_range_aref, $template_aligned_chain_aref );
                            $template_struct_pdb_range = multi_chain_pdb_rangify( $template_aligned_seqid_range_aref,
                                $ecod_pdbnum_href, $template_aligned_chain_aref );

                           #Template coverage should be number of aligned positions / total number of tmeplate positions

                            #Generate coverage
                            $template_coverage = multi_chain_region_coverage( $template_aligned_seqid_range_aref,
                                $template_aligned_chain_aref, $template_seqid_aref, $template_chain_aref );

                            #Generate ungapped coverage
                            my $ungapped_template_struct_seqid_range =
                              multi_chain_ungap_range( $template_struct_seqid_range, $$global_opt{gap_tol} );
                            my ( $ungapped_template_struct_seqid_aref, $ungapped_template_struct_chain_aref ) =
                              multi_chain_range_expand($ungapped_template_struct_seqid_range);
                            $ungapped_template_coverage =
                              multi_chain_region_coverage( $template_aligned_seqid_range_aref,
                                $template_aligned_chain_aref, $template_seqid_aref, $template_chain_aref );

                        }
                        else {
                            my @template_seqid;
                            my @template_chain;
                            for ( my $i = 0 ; $i < scalar(@ecod_chains) ; $i++ ) {
                                my $ecod_chain = $ecod_chains[$i];
                                my $ecod_seg   = $ecod_segs[$i];
                                if ( $ecod_seg eq 'FULL' ) {
                                    my ( $ecod_seqid_aref, $ecod_struct_seqid_aref, $ecod_pdbnum_aref, $asym_id ) =
                                      pdbml_seq_parse( $ecod_pdb, $ecod_chain );    #?
                                    if ( !$ecod_seqid_aref ) {
                                        print "WARNING! $ecod_pdb obsolete\n";
                                        last;
                                    }
                                    $ecod_seg = rangify($ecod_struct_seqid_aref);                  #This is problematic
                                    $ecod_seg = ungap_range( $ecod_seg, $$global_opt{gap_tol} );
                                }

                                #print "i: $i seg: $ecod_seg\n";
                                my $ecod_seg_aref = range_expand($ecod_seg);
                                push( @template_seqid, @$ecod_seg_aref );
                                for ( my $i = 0 ; $i < scalar(@$ecod_seg_aref) ; $i++ ) {
                                    push( @template_chain, $ecod_chain );
                                }
                            }
                            if ( scalar(@template_seqid) == 0 ) {
                                print "WARNING! No structured residues in template $ecod_domain_id\n";
                                last;
                            }

                            my ( $ecod_seqid_aref, $ecod_struct_seqid_aref, $ecod_pdbnum_href, $asym_id, $chain_aref )
                              = pdbml_mc_seq_parse( $ecod_pdb, \@ecod_chains );
                            if ( !$ecod_seqid_aref ) {
                                print "WARNING! $ecod_pdb obsolete!\n";
                                last;
                            }

                            #Change pos-inx into seqid_inx
                            my $template_range_aref = range_expand($template_range);    #Alignment pos-inx
                            my ( $template_aligned_seqid_aref, $template_aligned_chain_aref, $warning ) =
                              multi_chain_struct_region( $template_range_aref, \@template_seqid, \@template_chain );

                            if ( $warning == scalar(@$template_range_aref) ) {
                                print "WARNING! $ecod_domain_id is empty, check for corrupt HHM!\n";
                                last;
                            }

                            #Output tmeplate struct_seqid and struct_pdb ranges
                            $template_struct_seqid_range =
                              multi_chain_rangify( $template_aligned_seqid_aref, $template_aligned_chain_aref );
                            $template_struct_pdb_range =
                              multi_chain_pdb_rangify( $template_aligned_seqid_aref, $ecod_pdbnum_href,
                                $template_aligned_chain_aref );

                            $template_coverage = multi_chain_region_coverage( \@template_seqid, \@template_chain,
                                $template_aligned_seqid_aref, $template_aligned_chain_aref );
                            my $ungapped_template_struct_seqid_range =
                              multi_chain_ungap_range( $template_struct_seqid_range, $$global_opt{gap_tol} );
                            my ( $ungapped_template_struct_seqid_aref, $ungapped_template_struct_chain_aref ) =
                              multi_chain_range_expand($ungapped_template_struct_seqid_range);

                            $ungapped_template_coverage = multi_chain_region_coverage(
                                \@template_seqid, \@template_chain,
                                $ungapped_template_struct_seqid_aref,
                                $ungapped_template_struct_chain_aref
                            );

                        }

                        if ( $template_coverage < 0.7 ) {

                            $hh_hit_node->setAttribute( 'low_coverage_hit', 'true' );

                        }

                        $hh_hit_node->appendTextChild( 'template_range', $template_range );
                        $hh_hit_node->appendTextChild( 'query_range',    $query_range );
                        $hh_hit_node->appendTextChlid( 'query_seqid_range', $query_seqid_range );
                        $hh_hit_node->appendTextChild( 'query_pdb_range', $query_pdb_range );

                        my $template_seqid_range_node = $hh_xml_doc->createElement('template_seqid_range');
                        $template_seqid_range_node->appendTextNode($template_struct_seqid_range);
                        $hh_hit_node->appendChild($template_seqid_range_node);
                        $template_seqid_range_node->setAttribute( 'coverage',          $template_coverage );
                        $template_seqid_range_node->setAttribute( 'ungapped_coverage', $ungapped_template_coverage );

                        my $template_pdb_range_node = $hh_xml_doc->createElement('template_pdb_range');
                        $template_pdb_range_node->appendTextNode($template_struct_pdb_range);
                        $hh_hit_node->appendChild($template_pdb_range_node);
                        $template_pdb_range_node->setAttribute( 'coverage', $template_coverage );
                        $template_seqid_range_node->setAttribute( 'coverage', $template_coverage );

                        last;

                    }

                    if ( $ln =~ /^>((d|e)\w{4}[\w\.]+)/ ) {

                        $ecod_domain_id = $1;
                        $uid            = $$uid_lookup{$ecod_domain_id};

                        if ( $ecod_domain_id =~ /(d|e)(\w{4})/ ) {
                            $ecod_pdb = $2;
                        }
                        else {
                            $ecod_pdb = 'Unk';
                        }

                        $ecod_domain_id_regexp = qr/$ecod_domain_id/;
                        $hh_hit_node->setAttribute( 'uid',            $uid );
                        $hh_hit_node->setAttribute( 'ecod_domain_id', $ecod_domain_id );
                        next;
                    }
                    if ( $ln =~
/Probab=(\d+\.\d{2})\s+E-value=([\w\.\-\+]+)\s+Score=([\-\d\.]+)\s+Aligned_cols=(\d+)\s+Identities=([\d\.]+%)\s+Similarity=(\-?[\d\.]+)\s+Sum_probs=([\d\.]+)/
                      )
                    {
                        $hh_prob      = $1;
                        $hh_eval      = $2;
                        $hh_score     = $3;
                        $aligned_cols = $4;
                        $idents       = $5;
                        $similarities = $6;
                        $sum_probs    = $7;

                        $hh_hit_node->setAttribute( 'hh_prob',      $hh_prob );
                        $hh_hit_node->setAttribute( 'hh_eval',      $hh_eval );
                        $hh_hit_node->setAttribute( 'hh_score',     $hh_score );
                        $hh_hit_node->setAttribute( 'aligned_cols', $aligned_cols );
                        $hh_hit_node->setAttribute( 'idents',       $idents );
                        $hh_hit_node->setAttribute( 'similarities', $similarities );
                        $hh_hit_node->setAttribute( 'sum_probs',    $sum_probs );

                        next;
                    }

                    if ( $ln =~ /^(Q|T) ss_pred/ ) { next }
                    if ( $ln =~ /^(Q|T) ss_conf/ ) { next }
                    if ( $ln =~ /^Confidence/ )    { next }

                    if ( $ln =~ /^Q\s+$query_regexp\s+(\d+)\s+([\w\-]+)\s+(\d+)\s+\((\d+)\)/ ) {
                        if ( !defined($query_start) ) {
                            $query_start = $1;
                        }
                        $query_ln .= $2;
                        $query_end = $3;

                        next;
                    }
                    if ( $ln =~ /^T\s+$ecod_domain_id_regexp\s+(\d+)\s+([\w\-]+)\s+(\d+)\s+\((\d+)\)/ ) {
                        if ( !defined($template_start) ) {
                            $template_start = $1;
                        }
                        $template_ln .= $2;
                        $template_end = $3;

                        if ( $template_end == $template_ends{$hit_num}{$ecod_domain_id} ) {
                            $last_template_block_seen = 1;
                        }
                        next;
                    }
                    if ( $ln =~ /^\s+$/ ) { next }
                }
            }
        }
    }
}

sub read_reference_cache {
    use ECOD::Reference;
    load_references();

    my ($reference) = @_;

    my %uid;
    my %uid_lookup;

    open my $fh, "<", $REF_RANGE_CACHE{$reference}
      or die "ERROR! Could not open $REF_RANGE_CACHE{$reference} for reading:$!\n";

    while ( my $ln = <$fh> ) {
        next if $ln =~ /^#/;

        my @F = split /\s+/, $ln;

        my $uid            = $F[0];
        my $ecod_domain_id = $F[1];
        my $seqid_range    = $F[2];
        my $pdb            = $F[3];
        my $chain          = $F[4];

        $uid{$uid}{seqid_range}      = $seqid_range;
        $uid{$uid}{ecod_domain_id}   = $ecod_domain_id;
        $uid{$uid}{pdb}              = $pdb;
        $uid_lookup{$ecod_domain_id} = $uid;

        if ( $chain eq '.' ) {
            my %lchains;
            my @split = split /,/, $seqid_range;
            foreach my $s (@split) {
                $s =~ /(.+):\d+/;
                $lchains{$1}++;
            }
            my @chains = sort { $a cmp $b } keys %lchains;
            $uid{$uid}{chain} = join ",", @chains;
        }
        else {
            $uid{$uid}{chain} = $chain;
        }
    }

    return ( \%uid, \%uid_lookup );
}

sub count_domain_nodes {
    $_[0]->findnodes('.//domain')->size();
}

sub count_rep_domain_nodes {
    $_[0]->findnodes(qq{.//domain/cluster[\@level="$_[1]"][\@method="$_[2]"][\@domain_rep="true"]})->size();
}

sub assemble_pf_hhm_db {

    my ( $ecod_xml_doc, $tag, $pf_fasta_dir, $method ) = @_;

    foreach my $pf_group ( find_nodes( $ecod_xml_doc, '//pf_group' ) ) {
        my $pf_id = get_pf_id($pf_group);
        if ( count_domain_nodes($pf_group) > 0 ) {

            if ( $tag && count_rep_domain_nodes( $pf_group, $tag, $method ) < 2 ) {
                next;
            }

            my $pf_id_hhm_db;
            $tag
              ? $pf_id_hhm_db = "$pf_fasta_dir/$pf_id/$pf_id.hhm.db"
              : $pf_id_hhm_db = "$pf_fasta_dir/$pf_id/$pf_id.$tag.hhm_db";

            print "pf: $pf_id_hhm_db\n";

            open my $out_fh, ">", $pf_id_hhm_db
              or die "ERROR! Could not open $pf_id_hhm_db for writing:$!\n";

          DOMAIN:
            foreach my $domain ( find_domain_nodes($pf_group) ) {
                if ( is_domain_rep( $domain, $tag, $method ) ) {
                    my ( $uid, $ecod_domain_id ) = get_ids($domain);
                    my $short_uid = substr( $uid, 2, 5 );

                    my $hhm_file = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.hhm";

                    if ( -f $hhm_file ) {
                        open my $fh, "<", $hhm_file
                          or die "ERROR! Could not open $hhm_file for reading:$!\n";
                        while ( my $ln = <$fh> ) {
                            print $out_fh $ln;
                        }
                    }
                    else {
                        print "WARNING! No hhm file found for $pf_id:$uid\n";
                    }
                }
            }
            close $out_fh;
        }
    }
}

sub cluster_strip { 
	my ($ecod_xml_doc, $tag, $method) = @_;
	foreach my $domain_node (find_domain_nodes($ecod_xml_doc)) { 
		if ($domain_node->exists("cluster[\@level='$tag'][\@method='$method']")) { 
			foreach my $cluster_node ($domain_node->findnodes("cluster[\@level='$tag'][\@method='$method']")) { 
				$cluster_node->unbindNode;
			}
		}
	}
}



sub generate_hhsearch_jobs_for_domain_reps {

    my ( $ecod_xml_doc, $tag, $pf_fasta_dir, $global_opt, $method ) = @_;

    my %pf_rep_count;
    foreach my $pf_group ( find_nodes( $ecod_xml_doc, "//pf_group" ) ) {
        #my $pf_id = get_pf_id($pf_group);
		my $pf_id = $pf_group->findvalue('@pf_id');
        my $rep_count = count_rep_domain_nodes( $pf_group, $tag, $method );
        $pf_rep_count{$pf_id} = $rep_count;
    }
    my @cmds;
    foreach my $domain_rep ( find_domain_reps( $ecod_xml_doc, $tag, $method ) ) {
        my $pf_id = get_pf_id($domain_rep);
        if ($pf_rep_count{$pf_id} < 2) { 
			#Not warning-worthy. All singleton pf_groups are skipped	
			#print "WARNING! $pf_id rep count $pf_rep_count{$pf_id}\n";
			next;
		}

        my ( $uid, $ecod_domain_id ) = get_ids($domain_rep);
        my $short_uid = substr( $uid, 2, 5 );
        my $hhm_file = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.hhm";

        if ( !-f $hhm_file ) {
            print "WARNING! $hhm_file not found, skipping...\n";
            next;
        }

        my $result_fn = "$uid.$tag.hhr";
        if ( !-d "$pf_fasta_dir/$pf_id" ) {
            mkdir("$pf_fasta_dir/$pf_id");
        }
        my $ref_db;
        my $output_fn;
        if ( $mode eq 'F_ANALYZE' ) {
            $ref_db    = "$pf_fasta_dir/$pf_id/$pf_id.$tag.hhm_db";
            $output_fn = "$pf_fasta_dir/$pf_id/$result_fn";
            if ( !-f $ref_db ) {
                print "WARNING! $ref_db not found, skipping...\n";
                next;
            }
        }
        elsif ( $mode eq 'ALL_REPS' ) {
            use ECOD::Reference;
            load_references();
            my $current_version = $LATEST_REFERENCE;
            $output_fn = "$pf_fasta_dir/all/$result_fn";
            $ref_db    = "$pf_fasta_dir/all/ecod.$current_version.$tag.db";
        }
        else {
            die "ERROR! Unk mode $mode\n";
        }
        if ( !-f $output_fn || $$global_opt{force_overwrite} ) {
            my $cmd = "$HHSEARCH_EXE -i $hhm_file -o $output_fn -d $ref_db -cpu 8 -B 2500 -Z 2500";
            push( @cmds, $cmd );
        }
        if ( $mode eq 'F_ANALYZE' ) {
            if ( -f "$pf_fasta_dir/$pf_id/$result_fn" ) {

            }
        }
    }
    return \@cmds;
}

sub verify {
    my ( $pf_node, $pf_fasta_fn ) = @_;

    open my $fh, "<", $pf_fasta_fn
      or die "ERROR! Could not open $pf_fasta_fn for reading:$!";
    my %fasta_uids;
    while ( my $ln = <$fh> ) {
        if ( $ln =~ /^>(\d{9})/ ) {
            $fasta_uids{$1}++;
        }

    }
    my @f1 = sort keys %fasta_uids;
    my %pf_uids;
    foreach my $domain_node ( find_domain_nodes($pf_node) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
        $pf_uids{$uid}++;
    }
    my @f2 = sort keys %pf_uids;
    if ( array_eq( \@f1, \@f2 ) ) {
        return 1;
    }
    else {
        return 0;
    }

}

sub array_eq {

    use List::AllUtils qw( each_arrayref );

    my ( $aref1, $aref2 ) = @_;

    return 0 if scalar @$aref1 != scalar @$aref2;

    my $it = each_arrayref( $aref1, $aref2 );

    while ( my ( $x, $y ) = $it->() ) {
        return 0 unless $x eq $y;
    }
    return 1;

}

sub gen_pf_fasta {
    my ( $pf_node, $pf_fasta_fn ) = @_;
	return if $pf_node->findnodes('.//domain')->size() == 0;
    open my $out_fh, ">", $pf_fasta_fn
          or die "ERROR! Could not open $pf_fasta_fn for appending:$!\n";
    foreach my $domain_node ( find_domain_nodes($pf_node) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);

        my $domain_fasta_fn = get_domain_fasta_fn($uid);

        open my $fh, "<", $domain_fasta_fn
          or die "ERROR! Could not open $domain_fasta_fn for reading:$!\n";

        while ( my $ln = <$fh> ) {
            if ( $ln =~ /^>.*(\d{9})/ ) {
                print $out_fh ">$1\n";
            }
            else {
                print $out_fh $ln;
            }
        }
        close $fh;
    }
	close $out_fh;
}

sub get_domain_fasta_fn {
    my ($uid) = @_;

    my $s = substr( $uid, 2, 5 );

    my $path = "$DOMAIN_DATA_DIR/$s/$uid/$uid.fa";
    if ( -f $path && -s $path > 0 ) {
        return $path;
    }
    else {
        warn "WARNING! $path not found\n";
        return 0;
    }
}

sub is_domain_rep {
    $_[0]->findvalue(qq{cluster[\@level='$_[1]'][\@method='$_[2]']/\@domain_rep}) eq 'true';
}



1;
