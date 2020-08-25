#!/usr/bin/perl

#use warnings;
use warnings FATAL => 'all';
use File::Slurp;

unless(scalar(@ARGV) == 3){
    print "Give files with iPAO1 rxns, mets and genes\n";
    print "Ex: $0 iPAO1_rxns.csv iPAO1_mets.csv iPAO1_genes.csv\n";
    exit(0);
}
$ipao1_rxns_file=$ARGV[0];
$ipao1_mets_file=$ARGV[1];
$ipao1_genes_file=$ARGV[2];

$gene_prot_file="ProteinTable.txt";
read_PAO1_geneprot_table($gene_prot_file);
    
$ccbh_gbk_file="CCBH4851.gb1";
read_CCBH_gb($ccbh_gbk_file);
$tcdb_ortho_file="orthologues_blast_ccbh_tcdb.txt";
find_orthologues_tcdb($tcdb_ortho_file,0.000001);
#exit(0);

create_ipao1_sbml($ipao1_rxns_file,$ipao1_mets_file,$ipao1_genes_file);
print "Finished parsing iPAO1\n";
exit(0);

sub create_ipao1_sbml{
    my $xls_rxns=$_[0];
    my $xls_mets=$_[1];
    my $xls_genes=$_[2];

    $gapfill_rxn{"rxn00695"}=1;
    
    $sbml_header="<\?xml version=\"1.0\" encoding=\"UTF-8\"\?>\n";
    $sbml_header.="<sbml xmlns=\"http\://www.sbml.org/sbml/level3/version1/core\" level=\"3\" version=\"1\" xmlns\:fbc=\"http\:\/\/www\.sbml\.org\/sbml\/level3\/version1\/fbc\/version2\" fbc\:required=\"false\">\n";
    $sbml_ipao1_header=$sbml_header."<model id=\"iPAO1\" timeUnits=\"time\" fbc\:strict=\"true\" name=\"Reconstruction iPAO1 of multidrug\-resistant Pseudomonas aeruginosa\">\n";
    $sbml_ccbh_header=$sbml_header."<model id=\"iCCBH4851\" timeUnits=\"time\" fbc\:strict=\"true\" name=\"Reconstruction iCCBH4851 of multidrug\-resistant Pseudomonas aeruginosa CCBH4851\">\n";
    $sbml_header.="<model id=\"iPAO1\" timeUnits=\"time\" fbc\:strict=\"true\" name=\"Reconstruction iPAO1 of multidrug\-resistant Pseudomonas aeruginosa\">\n";
    $sbml_footer="\<\/model\>\n\<\/sbml\>\n";

    $obj_rxn="bio00991";
    $sbml_objective= " <fbc\:listOfObjectives fbc\:activeObjective=\"objective1\">\n";
    $sbml_objective.="  <fbc\:objective fbc\:type=\"maximize\" fbc\:id=\"objective1\">\n";
    $sbml_objective.="  <fbc\:listOfFluxObjectives>\n";
    $sbml_objective.="   <fbc:fluxObjective fbc:reaction=\"".$obj_rxn."\" fbc\:coefficient=\"1\.0\"\/>\n";
    $sbml_objective.="  <\/fbc\:listOfFluxObjectives>\n";
    $sbml_objective.="  <\/fbc\:objective>\n";
    $sbml_objective.=" <\/fbc\:listOfObjectives>\n";

    $sbml_units= " <listOfUnitDefinitions>\n";
    $sbml_units.="  <unitDefinition id=\"mmol_per_gDW_per_hr\">\n";
    $sbml_units.="  <listOfUnits>\n";
    $sbml_units.="   <unit scale=\"-3\" exponent=\"1\" multiplier=\"1\" kind=\"mole\"\/>\n";
    $sbml_units.="   <unit scale=\"0\" exponent=\"-1\" multiplier=\"1\" kind=\"gram\"\/>\n";
    $sbml_units.="   <unit scale=\"0\" exponent=\"-1\" multiplier=\"3600\" kind=\"second\"\/>\n";
    $sbml_units.="  <\/listOfUnits>\n";
    $sbml_units.="  <\/unitDefinition>\n";
    $sbml_units.=" <\/listOfUnitDefinitions>\n";

    $sbml_compartments= " <listOfCompartments>\n";
    $sbml_compartments.="  <compartment id=\"c\" name=\"Cytoplasm\" spatialDimensions=\"3\" constant=\"false\"\/>\n";
    $sbml_compartments.="  <compartment id=\"p\" name=\"Periplasm\" spatialDimensions=\"3\" constant=\"false\"\/>\n";
    $sbml_compartments.="  <compartment id=\"e\" name=\"Extracellular\" spatialDimensions=\"3\" constant=\"false\"\/>\n";
    $sbml_compartments.=" <\/listOfCompartments>\n";

    $sbml_param_header=" <listOfParameters>\n";
    $sbml_param_footer=" <\/listOfParameters>\n";
    $sbml_parameters="";
    $nbounds=0;
    
    $sbml_mets_header =" <listOfSpecies>\n";
    $sbml_mets_footer =" <\/listOfSpecies>\n";

    $sbml_genes_header=" <fbc:listOfGeneProducts>\n";
    $sbml_genes_footer=" <\/fbc:listOfGeneProducts>\n";

    $sbml_rxns_header =" <listOfReactions>\n";
    $sbml_rxns_footer =" <\/listOfReactions>\n";
    $sbml_rxn_footer="\<\/reaction\>\n";
    $sbml_from_header="<listOfReactants>\n";
    $sbml_from_footer="<\/listOfReactants>\n";
    $sbml_to_header="<listOfProducts>\n";
    $sbml_to_footer="<\/listOfProducts>\n";
    
    @component=read_file($xls_mets);
    $nmetabolites_ipao1=0;
    $sbml_metabolites="";
    for(my $im=0;$im<scalar(@component);$im++){
	chomp($component[$im]);
	unless(($component[$im] =~ /^\#/)||($component[$im] !~ /\w+/)){
	    @info=split(/\@/,$component[$im]);
	    $met_id=$info[1];
	    chomp($met_id);
	    $met_code=$info[2];
	    chomp($met_code);	   
	    $met_code =~ s/^\s+//;
	    $met_code =~ s/\s+$//;
	    $met_cmp=$info[3];
	    chomp($met_cmp);
	    $met_name=$info[5];
	    chomp($met_name);
	    $met_name =~ s/^\s+//;
	    $met_name =~ s/\s+$//;
	    $met_formula=$info[7];
	    chomp($met_formula);
	    $met_charge=$info[8];
	    chomp($met_charge);
	    if(defined $ipao1_met{$met_code}){
		print "$met_code already defined\n";
		exit(0);
	    }else{
		$nmetabolites_ipao1++;
		$ipao1_met{$met_code}=$met_name;
		print "($nmetabolites_ipao1) $met_code = $met_name (cmp $met_cmp)\n";
		$sbml_metabolites.="  <species id=\"$met_code\" name=\"$met_name\" compartment=\"$met_cmp\" hasOnlySubstanceUnits=\"false\" boundaryCondition=\"false\" constant=\"false\"\>\n";
		$sbml_metabolites.="   <notes>\n";
		$sbml_metabolites.="    <body xmlns=\"http\:\/\/www\.w3\.org\/1999\/xhtml\">\n";
		$sbml_metabolites.="     <p>FORMULA\: $met_formula<\/p>\n";
		$sbml_metabolites.="     <p>CHARGE\: $met_charge<\/p>\n";
		$sbml_metabolites.="    <\/body>\n";
		$sbml_metabolites.="   <\/notes>\n";
		$sbml_metabolites.="  <\/species>\n";
	    }
	}
    }   


    @component=read_file($xls_genes);
    $ngenes_ipao1=0;
    $sbml_genes="";
    for(my $ig=0;$ig<scalar(@component);$ig++){
	chomp($component[$ig]);
	unless($component[$ig] =~ /^\#/){
	    @info=split(/\@/,$component[$ig]);
	    chomp($info[0]);
	    $gid=$info[0];
	    $gid=~s/^\s+//;
	    $gid=~s/\s+$//;
	    if(defined $ipao1_gene{$gid}){
		print "$gid already defined\n";
		exit(0);
	    }else{
		$ngenes_ipao1++;
		$ipao1_gene{$gid}=$ngenes_ipao1;
		$sbml_genes.="  <fbc\:geneProduct fbc:id=\"$gid\" label=\"g\_$gid\"\/>\n";
		if(defined $ccbh_gid{$gid}){
		    $gene_in{$gid}=1;
		}else{
		    $gene_in{$gid}=0;
		}

	    }
	}
    }

# warning: there are rxns where spontaneity is indicated not on GENES col (ex:line 4316, rxn 9984)
#      <reaction id="R_rxn03253" name="acyl-CoA dehydrogenase (decanoyl-CoA)" fbc:lowerFluxBound="flux_b2" fbc:upperFluxBound="flux_b1" reversible="false">
#        <listOfReactants>
#          <speciesReference species="M_cpd00015_c" stoichiometry="1"/>
#          <speciesReference species="M_cpd03128_c" stoichiometry="1"/>
#        </listOfReactants>
#        <listOfProducts>
#          <speciesReference species="M_cpd00982_c" stoichiometry="1"/>
#          <speciesReference species="M_cpd03129_c" stoichiometry="1"/>
#        </listOfProducts>
#        <notes>
#          <body xmlns="http://www.w3.org/1999/xhtml">
#            <p>GENE_ASSOCIATION: ((PA0879) or (PA0506) or (PA2550))</p>
#            <p>SUBSYSTEM: Fatty acid metabolism</p>
#            <p>EC Number: </p>
#            <p>Confidence Level: </p>
#            <p>AUTHORS: </p>
#            <p/>
#          </body>
#        </notes>
#        <fbc:geneProductAssociation>
#          <fbc:or>
#          <fbc:geneProductRef fbc:geneProduct="PA0879"/>
#          <fbc:geneProductRef fbc:geneProduct="PA0506"/>
#          <fbc:geneProductRef fbc:geneProduct="PA2550"/>
#</fbc:or>
#        </fbc:geneProductAssociation>
#      </reaction>

    
    @component=read_file($xls_rxns);
    $nrxns_ipao1=0;
    $sbml_reactions="";
    $sbml_ipao1_reactions="";
    $sbml_ccbh_reactions="";
    $sbml_ccbh_reactions_mm="";
    $sbml_ccbh_reactions_uptake="";
    $sbml_ccbh_gapfill_reactions="";
    for(my $ir=0;$ir<scalar(@component);$ir++){ ## Parse each rxn from iPAO1 xls file
	chomp($component[$ir]);
#      <reaction id="R_rxn03253" name="acyl-CoA dehydrogenase (decanoyl-CoA)" fbc:lowerFluxBound="flux_b2" fbc:upperFluxBound="flux_b1" reversible="false">
	unless($component[$ir] =~ /^\#/){
	    @info=split(/\@/,$component[$ir]);
	    $info[1]=~s/\s+$//; # rxn_id
	    $info[1]=~s/^\s+//; 
	    $info[2]=~s/\s+$//; # rxn_name
	    $info[2]=~s/^\s+//; 
	    $info[3]=~s/\s+$//; # rxn_formula
	    $info[3]=~s/^\s+//; 
	    $info[4]=~s/\s+$//; # rxn_GPR
	    $info[4]=~s/^\s+//; 
	    $info[7]=~s/\s+$//; # rxn_EC_TC_numbers
	    $info[7]=~s/^\s+//; 
	    $info[8]=~s/\s+$//; # rxn_subsystem
	    $info[8]=~s/^\s+//; 
	    $info[9]=~s/\s+$//; # rxn_reversibility
	    $info[9]=~s/^\s+//; 
	    $info[10]=~s/\s+$//; # rxn_flux_lb
	    $info[10]=~s/^\s+//; 
	    $info[11]=~s/\s+$//; # rxn_flux_ub
	    $info[11]=~s/^\s+//; 
	    $rxn_id=$info[1];
	    $rxn_name=$info[2];
	    $rxn_formula=$info[3];
	    $rxn_gpr=$info[4];
	    $rxn_EC_TC=$info[7];
	    $rxn_system=$info[8];
	    $rxn_rev=$info[9];
	    $rxn_lb=$info[10];
	    $rxn_ub=$info[11];

	    if($rxn_rev == 1){
		$rxn_rev="true";
		$rxn_id_rev=$rxn_id."_REV";
		$rxn_name_rev=$rxn_name."_REV";
	    }else{	    
		$rxn_rev="false";
	    }

	    if($rxn_id=~/EX\_/i){
		print "XXXSP ($rxn_id) $component[$ir]\n";
		$exchange_rxn{$rxn_id}=1;
		if($rxn_rev eq "true"){
		    $exchange_rxn{$rxn_id_rev}=1;
		}
	    }
	    if($component[$ir]=~/spontaneous/i){
		print "XXXSP ($rxn_id) $component[$ir]\n";
		$spontaneous_rxn{$rxn_id}=1;
		if($rxn_rev eq "true"){
		    $spontaneous_rxn{$rxn_id_rev}=1;
		}
	    }elsif(($rxn_gpr =~ /null/i)&&($rxn_EC_TC =~ /null/i)){
		print "XXXNULL $rxn_id ($rxn_name)\n";
		$null_rxn{$rxn_id}=1;
		if($rxn_rev eq "true"){
		    $null_rxn{$rxn_id_rev}=1;
		}
	    }

	    unless(($rxn_EC_TC =~ /spontaneous/i)||($rxn_EC_TC =~ /null/i)){
		if($rxn_EC_TC =~ /\d\.[a-z]/i){ # TC numbers
		    if($rxn_EC_TC =~ /\|/){
			@rxn_TCs= split('\|', $rxn_EC_TC);
			for(my $rt=0;$rt<scalar(@rxn_TCs);$rt++){
			    $rxn_TCs[$rt]=~s/\s//g;
			    if(defined $ccbh_tc{$rxn_TCs[$rt]}){
				print "XXXTC iPAO1 ($rxn_id) $rxn_TCs[$rt]\n";
				$ccbh_tc_in{$rxn_id}=1;
				if($rxn_rev eq "true"){
				    $ccbh_tc_in{$rxn_id_rev}=1;
				}	
			    }	
			}
		    }else{
			$rxn_TC=$rxn_EC_TC;
			$rxn_TC=~s/\s//g;
			if(defined $ccbh_tc{$rxn_TC}){
			    print "XXXTC iPAO1 ($rxn_id) $rxn_TC\n";
			    $ccbh_tc_in{$rxn_id}=1;
			    if($rxn_rev eq "true"){
				$ccbh_tc_in{$rxn_id_rev}=1;
			    }	
			}	
		    }		   
		}elsif($rxn_EC_TC =~ /\d+\.\d+\./i){ # EC numbers
		    if($rxn_EC_TC =~ /\|/){
			@rxn_ECs= split('\|', $rxn_EC_TC);
			for(my $re=0;$re<scalar(@rxn_ECs);$re++){
			    $rxn_ECs[$re]=~s/\s//g;
			    if(defined $ccbh_ec_gbk{$rxn_ECs[$re]}){
				print "XXXEC iPAO1 ($rxn_id) $rxn_ECs[$re]\n";
				$ccbh_ec_in{$rxn_id}=1;
				if($rxn_rev eq "true"){
				    $ccbh_ec_in{$rxn_id_rev}=1;
				}	
			    }
			}
		    }else{
			$rxn_EC=$rxn_EC_TC;
			$rxn_EC=~s/\s//g;
			if(defined $ccbh_ec_gbk{$rxn_EC}){
			    print "XXXEC iPAO1 ($rxn_id) $rxn_EC\n";
			    $ccbh_ec_in{$rxn_id}=1;
			    if($rxn_rev eq "true"){
				$ccbh_ec_in{$rxn_id_rev}=1;
			    }	
			}
		    }
		}else{ 
		    print "Waaaat? ($rxn_id) EC_TC ($rxn_EC_TC)\n";
		    exit(0);
		}
	
	    }	    

	    $sbml_gpr="";
	    if($rxn_gpr=~/PA\d{4}/g){
		    print "XXXB $rxn_id check GPR[$rxn_gpr]\n";		    
		    $ccbh_gpr=check_geneset($rxn_gpr,$rxn_id); # Check if GPR is true with ccbh genes
		    if($ccbh_gpr =~ m/yes/){
			$ccbh_gpr_in{$rxn_id}=1;
			if($rxn_rev eq "true"){
			    $ccbh_gpr_in{$rxn_id_rev}=1;
			}
		    }
		    $sbml_gpr=make_fbc_GPR($rxn_gpr);
		    print "-----------++++++++++++++++++++++-------------------++++++++++++++++++\n";
		    print "CCBH GPR RXN:: $ccbh_gpr\n";
		    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		    print "$sbml_gpr";
		    print "-----------++++++++++++++++++++++-------------------++++++++++++++++++\n";
	    }
	    
	    if($rxn_formula =~ /(.*)(\<\=\>|\-\>)(.*)/){
		$from=$1;
		$dir=$2;
		$to=$3;
#		print "RXN $rxn_id :: [[$from]] $dir [[$to]]\n RXN $rxn_id :: ";
	    }else{
		print "XXXF Why? $rxn_formula\n";
		exit(0);
	    }
	    undef @met_from;
	    $nfrom=0;
	    if(!defined $from){
		$sbml_from_str="";
	    }elsif($from =~ /\+/){
		$sbml_from_str="";
		my @met_from= split('\+', $from);
		my $nm=scalar(@met_from);
		$nfrom=$nm;
		for(my $m=0;$m<$nm;$m++){
		    chomp($met_from[$m]);
		    $met_from[$m] =~ s/^\s+//;
		    $met_from[$m] =~ s/\s+$//;
		    if($met_from[$m] =~ /(.+?)\s(cpd.*)/){
			$s=$1;
			$met=$2;
			$s =~ s/\s+$//;
		    }elsif($met_from[$m] =~ /(cpd.*)/){
			$met=$1;
			$s=1;
		    }else{
			print "What? $met_from[$m]\n";
			exit(0);
		    }
		    $met =~ s/^\s+//;
		    $met =~ s/\s+$//;
		    unless(defined $ipao1_met{$met}){
			print "MM What? $met\n";
			exit(0);
		    }
		    $sbml_from_str.="\<speciesReference species=\"$met\" stoichiometry=\"$s\"\/\>\n";
#		    print " [[($s) $met]]";
		}
	    }elsif($from =~ /cpd/){
		$nfrom=1;
		$sbml_from_str="";
		chomp($from);
		$from =~ s/^\s+//;
		$from =~ s/\s+$//;
		if($from =~ /(.+?)\s(cpd.*)/){
		    $s=$1;
		    $met=$2;
		    $s =~ s/\s+$//;
		}elsif($from =~ /(cpd.*)/){
		    $met=$1;
		    $s=1;
		}else{
		    print "What? $from\n";
		    exit(0);
		}
		$met =~ s/^\s+//;
		$met =~ s/\s+$//;
		unless(defined $ipao1_met{$met}){
		    print "MM What? $met\n";
		    exit(0);
		}
		$sbml_from_str.="\<speciesReference species=\"$met\" stoichiometry=\"$s\"\/\>\n";
#		print " [[($s) $met]]";
	    }else{
		$sbml_from_str="";
		print "from What? $from\n";
		exit(0);
	    }
	    $dir =~ s/^\s+//;
	    $dir =~ s/\s+$//;
#	    print "\<\<$dir\>\>";

	    undef @met_to;
	    $nto=0;
	    if(!defined $to){
		$sbml_to_str="";
	    }elsif($to =~ /\+/){
		$sbml_to_str="";
		my @met_to= split('\+', $to);
		my $nm=scalar(@met_to);
		$nto=$nm;
		for(my $m=0;$m<$nm;$m++){
		    chomp($met_to[$m]);
		    $met_to[$m] =~ s/^\s+//;
		    $met_to[$m] =~ s/\s+$//;
		    if($met_to[$m] =~ /(.+?)\s(cpd.*)/){
			$s=$1;
			$met=$2;
			$s =~ s/\s+$//;
		    }elsif($met_to[$m] =~ /(cpd.*)/){
			$met=$1;
			$s=1;
		    }else{
			print "What? $met_to[$m]\n";
			exit(0);
		    }
		    $met =~ s/^\s+//;
		    $met =~ s/\s+$//;
		    unless(defined $ipao1_met{$met}){
			print "MM What? $met\n";
			exit(0);
		    }
		    $sbml_to_str.="\<speciesReference species=\"$met\" stoichiometry=\"$s\"\/\>\n";
#		    print " [[($s) $met]]";
		}
	    }elsif($to =~ /cpd/){
		$sbml_to_str="";
		$nto=1;
		chomp($to);
		$to =~ s/^\s+//;
		$to =~ s/\s+$//;
		if($to =~ /(.+?)(cpd.*)/){
		    $s=$1;
		    $met=$2;
		    $s =~ s/\s+$//;
		}elsif($to =~ /(cpd.*)/){
		    $met=$1;
		    $s=1;
		}else{
		    print "To What? ($rxn_formula) $to\n";
		    exit(0);
		}
		$met =~ s/^\s+//;
		$met =~ s/\s+$//;
		unless(defined $ipao1_met{$met}){
		    print "MM What? $met\n";
		    exit(0);
		}
		$sbml_to_str.="\<speciesReference species=\"$met\" stoichiometry=\"$s\"\/\>\n";
#		print " [[($s) $met]]";
	    }else{
		$sbml_to_str="";
#		print "To nw? $to\n";
#		exit(0);
	    }
#	    print "\n";	    
	    
	    if($rxn_rev eq "true"){

		print "XXXBLR ($rxn_id) $rxn_lb\n";
		print "XXXBUR ($rxn_id) $rxn_ub\n";
		if($rxn_lb > 0.0){
		    $rxn_ub_rev=0.0;
		    $rxn_lb_rev=0.0;
		}else{
		    if(defined $rxn_lb){
			$rxn_ub_rev=(-1.0)*$rxn_lb;
			$rxn_lb=0.0;
		    }else{
			print "Undefined rxn_lb on $rxn_id\n";
			exit(0);
		    }			
		    if($rxn_ub < 0.0){
			if(defined $rxn_ub){
			    $rxn_lb_rev=(-1.0)*$rxn_ub;
			    $rxn_ub=0.0;
			}else{
			    print "Undefined rxn_ub on $rxn_id\n";
			    exit(0);
			}
		    }else{
			$rxn_lb_rev=0.0;
		    }
		}
		print "XXXBLR1 ($rxn_id_rev) $rxn_lb_rev\n";
		print "XXXBUR1 ($rxn_id_rev) $rxn_ub_rev\n";

		unless(defined $fbc_bound{$rxn_ub_rev}){
		    $nbounds++;
		    $fbc_bound{$rxn_ub_rev}="flux_b"."$nbounds";		
		    $sbml_parameters.="\<parameter id=\"$fbc_bound{$rxn_ub_rev}\" constant=\"true\"  value=\"$rxn_ub_rev\" units=\"mmol\_per\_gDW\_per\_hr\"\/\>\n";
		}
		$ub_rev_str=$fbc_bound{$rxn_ub_rev};
		$reac_ub_rev_str="fbc\:upperFluxBound=\"".$ub_rev_str."\"";
		
		unless(defined $fbc_bound{$rxn_lb_rev}){
		    $nbounds++;
		    $fbc_bound{$rxn_lb_rev}="flux_b"."$nbounds";		
		    $sbml_parameters.="\<parameter id=\"$fbc_bound{$rxn_lb_rev}\" constant=\"true\"  value=\"$rxn_lb_rev\" units=\"mmol\_per\_gDW\_per\_hr\"\/\>\n";
		}
		$lb_rev_str=$fbc_bound{$rxn_lb_rev};
		$reac_lb_rev_str="fbc\:lowerFluxBound=\"".$lb_rev_str."\"";

		$sbml_rxn_rev_header="\<reaction id=\"$rxn_id_rev\" name=\"$rxn_name_rev\" $reac_lb_rev_str $reac_ub_rev_str reversible=\"false\"\>\n";
		if($nto>0){
		    $sbml_from_rev=$sbml_from_header.$sbml_to_str.$sbml_from_footer;
		}else{
		    $sbml_from_rev="";
		}
		if($nfrom>0){
		    $sbml_to_rev=$sbml_to_header.$sbml_from_str.$sbml_to_footer;
		}else{
		    $sbml_to_rev="";
		}
		$sbml_rxn_rev=$sbml_rxn_rev_header.$sbml_from_rev.$sbml_to_rev.$sbml_gpr.$sbml_rxn_footer;
		$sbml_reactions.=$sbml_rxn_rev;
		print "XXXBr  $sbml_rxn_rev";
		print "---------------------------------------------------------- \n";
# iPAO1 rxn is also a CCBH rxn if EC or TC match or GPR is true with CCBH genes
		if((defined $ccbh_tc_in{$rxn_id_rev})||(defined $ccbh_ec_in{$rxn_id_rev})||(defined $ccbh_gpr_in{$rxn_id_rev})||(defined $spontaneous_rxn{$rxn_id_rev})||(defined $null_rxn{$rxn_id_rev})||(defined $exchange_rxn{$rxn_id_rev})||(defined $gapfill_rxn{$rxn_id_rev})){
		    $sbml_ccbh_reactions.=$sbml_rxn_rev;
		    if(defined $exchange_rxn{$rxn_id_rev}){
			$sbml_ccbh_reactions_uptake.=$sbml_rxn_rev;
		    }else{
			$sbml_ccbh_reactions_mm.=$sbml_rxn_rev;		    
		    }
		}else{
		    $sbml_ccbh_gapfill_reactions.=$sbml_rxn_rev;
		}		    
	    }
	    print "XXXBL ($rxn_id) $rxn_lb\n";
	    print "XXXBU ($rxn_id) $rxn_ub\n";
	    	    
	    unless(defined $fbc_bound{$rxn_ub}){
		$nbounds++;
		$fbc_bound{$rxn_ub}="flux_b"."$nbounds";		
		$sbml_parameters.="\<parameter id=\"$fbc_bound{$rxn_ub}\" constant=\"true\"  value=\"$rxn_ub\" units=\"mmol\_per\_gDW\_per\_hr\"\/\>\n";
	    }
	    $ub_str=$fbc_bound{$rxn_ub};
	    $reac_ub_str="fbc\:upperFluxBound=\"".$ub_str."\"";
	    unless(defined $fbc_bound{$rxn_lb}){
		$nbounds++;
		$fbc_bound{$rxn_lb}="flux_b"."$nbounds";		
		$sbml_parameters.="\<parameter id=\"$fbc_bound{$rxn_lb}\" constant=\"true\"  value=\"$rxn_lb\" units=\"mmol\_per\_gDW\_per\_hr\"\/\>\n";
	    }
	    $lb_str=$fbc_bound{$rxn_lb};
	    $reac_lb_str="fbc\:lowerFluxBound=\"".$lb_str."\"";

	    $sbml_rxn_header="\<reaction id=\"$rxn_id\" name=\"$rxn_name\" $reac_lb_str $reac_ub_str reversible=\"false\"\>\n";
	    if($nfrom>0){
		$sbml_from=$sbml_from_header.$sbml_from_str.$sbml_from_footer;
	    }else{
		$sbml_from="";
	    }
	    if($nto>0){
		$sbml_to=$sbml_to_header.$sbml_to_str.$sbml_to_footer;
	    }else{
		$sbml_to="";
	    }
	    $sbml_rxn=$sbml_rxn_header.$sbml_from.$sbml_to.$sbml_gpr.$sbml_rxn_footer;
	    $sbml_reactions.=$sbml_rxn;
	    print "XXXB ($rxn_id) $sbml_rxn\n";
#	    print "$sbml_rxn";
	    print "---------------------------------------------------------- \n";
	    if($rxn_id eq $obj_rxn){
		$sbml_ccbh_reactions.=$sbml_rxn;
		$sbml_ccbh_reactions_mm.=$sbml_rxn;		    
#		$sbml_ccbh_gapfill_reactions.=$sbml_rxn;
	    }elsif((defined $ccbh_tc_in{$rxn_id})||(defined $ccbh_ec_in{$rxn_id})||(defined $ccbh_gpr_in{$rxn_id})||(defined $spontaneous_rxn{$rxn_id})||(defined $null_rxn{$rxn_id})||(defined $exchange_rxn{$rxn_id})||(defined $gapfill_rxn{$rxn_id})){
		$sbml_ccbh_reactions.=$sbml_rxn;
		$sbml_ccbh_reactions_mm.=$sbml_rxn;		    
	    }else{
		$sbml_ccbh_gapfill_reactions.=$sbml_rxn;
	    }		    
		    
	}
    }

    $sbml_genes=$sbml_genes_header.$sbml_genes.$sbml_genes_footer;
    $sbml_metabolites=$sbml_mets_header.$sbml_metabolites.$sbml_mets_footer;
    $sbml_reactions=$sbml_rxns_header.$sbml_reactions.$sbml_rxns_footer;
    $sbml_ccbh_reactions=$sbml_rxns_header.$sbml_ccbh_reactions.$sbml_rxns_footer;
    $sbml_ccbh_reactions_mm=$sbml_rxns_header.$sbml_ccbh_reactions_mm.$sbml_rxns_footer;
    $sbml_ccbh_reactions_uptake=$sbml_rxns_header.$sbml_ccbh_reactions_uptake.$sbml_rxns_footer;
    $sbml_ccbh_gapfill_reactions=$sbml_rxns_header.$sbml_ccbh_gapfill_reactions.$sbml_rxns_footer;
    $sbml_parameters=$sbml_param_header.$sbml_parameters.$sbml_param_footer;

    $sbml_model=$sbml_header.$sbml_objective.$sbml_units.$sbml_compartments.$sbml_parameters.$sbml_metabolites.$sbml_genes.$sbml_reactions.$sbml_footer;
    $sbml_ccbh_model=$sbml_ccbh_header.$sbml_objective.$sbml_units.$sbml_compartments.$sbml_parameters.$sbml_metabolites.$sbml_genes.$sbml_ccbh_reactions.$sbml_footer;
    $sbml_ccbh_model_mm=$sbml_ccbh_header.$sbml_objective.$sbml_units.$sbml_compartments.$sbml_parameters.$sbml_metabolites.$sbml_genes.$sbml_ccbh_reactions_mm.$sbml_footer;
    $sbml_ccbh_model_uptake=$sbml_ccbh_header.$sbml_units.$sbml_compartments.$sbml_parameters.$sbml_metabolites.$sbml_genes.$sbml_ccbh_reactions_uptake.$sbml_footer;
    $sbml_ccbh_gapfill_model=$sbml_ccbh_header.$sbml_units.$sbml_compartments.$sbml_parameters.$sbml_metabolites.$sbml_genes.$sbml_ccbh_gapfill_reactions.$sbml_footer;
    my $ipao1_sbml_file="ipao1.sbml";
    write_file($ipao1_sbml_file,\$sbml_model);

    my $ccbh_sbml_file="ccbh4851.sbml";
    write_file($ccbh_sbml_file,\$sbml_ccbh_model);
    my $ccbh_mm_sbml_file="ccbh4851_mm.sbml";
    write_file($ccbh_mm_sbml_file,\$sbml_ccbh_model_mm);
    my $ccbh_uptake_sbml_file="ccbh4851_uptake.sbml";
    write_file($ccbh_uptake_sbml_file,\$sbml_ccbh_model_uptake);
    

    my $ccbh_sbml_gapfill_file="ccbh4851_gapfill.sbml";
    write_file($ccbh_sbml_gapfill_file,\$sbml_ccbh_gapfill_model);

}
##############################################################

sub make_fbc_GPR{
    my $orig_gpr=$_[0];    
    undef %boolean_str;
    undef %fbc_istr;
    undef %sstr;
    undef %fbc_rxn; # Identify fbc_str for each rxn_id

    my $reac_gpr=$orig_gpr;
    do{
	undef @queue;
	$queue[0] = $reac_gpr;
	$regex = qr/
	    (                   # start of bracket 1
		\(                   # match an opening angle bracket
			(?:               
			 [^\(\)]++     # one or more non angle brackets, non backtracking
			 |                  
			 (?1)        # recurse to bracket 1
			)*                 
			\)                   # match a closing angle bracket
	    )                   # end of bracket 1
	    /x;

	$" = "\n\t";

	while( @queue ){
	    my $string = shift @queue;
	    my @groups = $string =~ m/$regex/g;
	    my $ng=scalar(@groups);
	    if($ng > 0){
		for(my $k=0;$k<$ng;$k++){
		    unless(($groups[$k] =~ / or /)&&($groups[$k] =~ / and /)){
			unless(defined $sstr{$groups[$k]}){
			    $nsstr++;
			    $istr="str".$nsstr;
			    $sstr{$groups[$k]}=$istr;
			    $boolean_str{$istr}=$groups[$k];
			}
			$reac_gpr =~ s/$groups[$k]/$sstr{$groups[$k]}/g;
		    }
		}
	    }
	    unshift @queue, map { s/^\(//; s/\)$//; $_ } @groups;
	}
	if($reac_gpr=~/\(str\d+\)/){
	    my @sp=($reac_gpr=~/\((str\d+)\)/g);
	    $nsp=scalar(@sp);
	    for($j=0;$j<$nsp;$j++){
		$reac_gpr =~ s/\($sp[$j]\)/$sp[$j]/g;
	    }
	}
    }while(($reac_gpr =~ / or /)&&($reac_gpr =~ / and /));

    print "XXXOIII $reac_gpr\n";
    $reaction_gpr=$reac_gpr;
    $fbc_str="";
    if(($reaction_gpr =~ / or /)&&($reaction_gpr !~ / and /)){
	my @tmp_or= split(' or ', $reaction_gpr);
	$fbc_str.="        <fbc:or>\n";
	for(my $ior=0;$ior<scalar(@tmp_or);$ior++){
	    $tmp_or[$ior]=~s/\(+//g;
	    $tmp_or[$ior]=~s/\)+//g;
	    $fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$tmp_or[$ior]."\"\/>\n";
	}
	$fbc_str.="        <\/fbc:or>\n";
    }elsif(($reaction_gpr =~ / and /)&&($reaction_gpr !~ / or /)){
	my @tmp_and= split(' and ', $reaction_gpr);
	$fbc_str.="        <fbc:and>\n";
	for(my $iand=0;$iand<scalar(@tmp_and);$iand++){
	    $tmp_and[$iand]=~s/\(+//g;
	    $tmp_and[$iand]=~s/\)+//g;
	    $fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$tmp_and[$iand]."\"\/>\n";
	}
	$fbc_str.="        <\/fbc:and>\n";
    }elsif(($reaction_gpr !~ / and /)&&($reaction_gpr !~ / or /)){
	if($reaction_gpr =~ m/(PA\d{4}\.d{1})/){
	    $rg=$1;
	}elsif($reaction_gpr =~ m/(PA\d{4})/){
	    $rg=$1;
	}else{
	    print "Waaat? GPR $reaction_gpr\n";
	    exit(0);
	}
	$fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$rg."\"\/>\n";
    }
    $fbc_rxn=$fbc_str;


    foreach my $istr(keys %boolean_str){
	my $fbc_str="";
	$reaction_gpr=$boolean_str{$istr};
	if(($reaction_gpr =~ / or /)&&($reaction_gpr !~ / and /)){
	    print "XXXOR $reaction_gpr ";
	    my @tmp_or= split(' or ', $reaction_gpr);
#	    my @tmp_or=($reaction_gpr =~ /(\w.*\d+)/g);
	    $fbc_str.="        <fbc:or>\n";
	    for(my $ior=0;$ior<scalar(@tmp_or);$ior++){
#		unless($tmp_or[$ior] =~ /or/){
		    $tmp_or[$ior]=~s/\(+//g;
		    $tmp_or[$ior]=~s/\)+//g;
		    print "($tmp_or[$ior])";
		    $fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$tmp_or[$ior]."\"\/>\n";
#		}
	    }
	    print "\n";
	    $fbc_str.="        <\/fbc:or>\n";
	}elsif(($reaction_gpr =~ / and /)&&($reaction_gpr !~ / or /)){
	    print "XXXAND $reaction_gpr ";
	    $fbc_str.="        <fbc:and>\n";
	    my @tmp_and= split(' and ', $reaction_gpr);
#	    my @tmp_and=($reaction_gpr =~ /(\w+)/g);
	    for(my $iand=0;$iand<scalar(@tmp_and);$iand++){
#		unless($tmp_and[$iand] =~ /and/){
		    $tmp_and[$iand]=~s/\(+//g;
		    $tmp_and[$iand]=~s/\)+//g;
		    print "($tmp_and[$iand])";
		    $fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$tmp_and[$iand]."\"\/>\n";
#		}
	    }
	    print "\n";
	    $fbc_str.="        <\/fbc:and>\n";
	}elsif(($reaction_gpr !~ / and /)&&($reaction_gpr !~ / or /)){
	    if($reaction_gpr =~ m/(PA\d{4}\.d{1})/){
		$rg=$1;
	    }elsif($reaction_gpr =~ m/(PA\d{4})/){
		$rg=$1;
	    }else{
		print "Waaat? GPR $reaction_gpr\n";
		exit(0);
	    }
	    $fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$rg."\"\/>\n";

#	    my @tmp_s = ($reaction_gpr =~ /(\w+)/);
#	    unless(scalar(@tmp_s)>0){
#		print "WAAAT $reaction_gpr\n";
#		exit(0);
#	    }
#	    $fbc_str.="          <fbc:geneProductRef fbc:geneProduct=\"".$tmp_s[0]."\"\/>\n";
	}
	$fbc_istr{$istr}=$fbc_str;
    }

    $fbc_str="        <fbc:geneProductAssociation>\n".$fbc_rxn;
    undef @istr_tmp;
    @istr_tmp=($fbc_str=~/<fbc:geneProductRef fbc:geneProduct=\"(str\d+)\"\/>\n/g);
    while(scalar(@istr_tmp)>0){
	for(my $l=0;$l<scalar(@istr_tmp);$l++){
	    $lstr=$istr_tmp[$l];
	    $fbc_str=~s/          <fbc:geneProductRef fbc:geneProduct=\"$lstr\"\/>\n/$fbc_istr{$lstr}/g;
	}
	@istr_tmp=($fbc_str=~/<fbc:geneProductRef fbc:geneProduct=\"(str\d+)\"\/>\n/g);
    }
    $fbc_str.="        <\/fbc:geneProductAssociation>\n";
    $fbc_rxn=$fbc_str;

#    print "GPR($orig_gpr)==>($fbc_rxn)\n";
   return($fbc_rxn);
}
##################a#################################################

sub check_geneset{

    my $entry=$_[0];
    my $rid=$_[1];    
    my $gpr_ok="no";
    chomp($entry);
    print "XXXZ ($rid) Checking GPR |**|$entry|**|\n";
    
    $missing_gene="";
    my $perl_expr = $entry;

    @geneset = ($perl_expr =~ m/(PA\d{4}(?!\.))/g);
    $ng=scalar(@geneset);
    if($ng>0){
	print "XXXZ1b ($rid) ($ng): "; 
	for(my $i=0;$i<$ng;$i++){
	    $geneset[$i] =~ s/\s+//g;
	    if(defined $ccbh_gid{$geneset[$i]}){
		$gene_in{$geneset[$i]}=1;
	    }else{
		$gene_in{$geneset[$i]}=0;
	    }		
	    print "<<$geneset[$i]>><<($gene_in{$geneset[$i]})>>";
	    if(!defined $ipao1_gene{$geneset[$i]}){
		$ngenes_ipao1++;
		$ipao1_gene{$geneset[$i]}=$ngenes_ipao1;
		$sbml_genes.="  <fbc\:geneProduct fbc:id=\"$geneset[$i]\" label=\"g\_$geneset[$i]\"\/>\n";
		print "\n XXXNG <<$geneset[$i] not defined ($gene_in{$geneset[$i]})>>\n";
#		exit(0);
	    }
	    if(defined $gene_in{$geneset[$i]}){
		print "($gene_in{$geneset[$i]}) ";
	    }else{
		print "($geneset[$i] not checked) ";
		exit(0);
		$gene_in{$geneset[$i]}=0;
	    }	    
 	    if($gene_in{$geneset[$i]} == 0){
		$missing_gene.=" $geneset[$i]";
	    }
#	    $perl_expr =~ s/$geneset[$i]/$gene_in{$geneset[$i]}/eg;
	}
	$perl_expr =~ s/(PA\d{4}(?!\.))/$gene_in{$1}/eg;
	print "\n";
	print "XXXZ1b ($rid) perl_expr = $perl_expr\n";
    }

    @geneset = ($perl_expr =~ m/(PA\d{4}\.\d{1})/g);
    $ng=scalar(@geneset);
    if($ng>0){
	print "XXXZ1a ($rid) ($ng): ";
	for(my $i=0;$i<$ng;$i++){
#	    print "$geneset[$i]";
	    if(defined $ccbh_gid{$geneset[$i]}){
		$gene_in{$geneset[$i]}=1;
	    }else{
		$gene_in{$geneset[$i]}=0;
	    }		
	    print "<<$geneset[$i]>><<($gene_in{$geneset[$i]})>>";
	    if(!defined $ipao1_gene{$geneset[$i]}){
		$ngenes_ipao1++;
		$ipao1_gene{$geneset[$i]}=$ngenes_ipao1;
		$sbml_genes.="  <fbc\:geneProduct fbc:id=\"$geneset[$i]\" label=\"g\_$geneset[$i]\"\/>\n";
		print "\n XXXNG <<$geneset[$i] not defined ($gene_in{$geneset[$i]})>>\n";
#		exit(0);
	    }
	    if(defined $gene_in{$geneset[$i]}){
		print "($gene_in{$geneset[$i]}) ";
	    }else{
		print "($geneset[$i] not checked) ";
		exit(0);
		$gene_in{$geneset[$i]}=0;
	    }	    
# 	    if(!defined $ccbh_gene_fasta{$geneset[$i]}){
 	    if($gene_in{$geneset[$i]} == 0){
		$missing_gene.=" $geneset[$i]";
	    }
	}
	$perl_expr =~ s/(PA\d{4}\.\d{1})/$gene_in{$1}/eg;
	print "\n";
	print "XXXZ1a ($rid) perl_expr = $perl_expr\n";
    }
    
    $perl_expr =~ s/ and /\*/g;
    $perl_expr =~ s/ or /\|/g;
    print "XXXZ2 ($rid) perl_expr = $perl_expr\n";
    my $result = eval $perl_expr;
    print "XXXZ3 ($rid) result = $result\n";
    if($result == 1){
	$gpr_ok="yes";
    }
    if($missing_gene =~ /\w+/){
	$gpr_ok.=" \(missing $missing_gene\)";
    }
    print "XXXZ4 Done! ($rid) GPR=$gpr_ok\n";
    return($gpr_ok);
}

# ===================================================================================
sub read_CCBH_gb{
## Corrected file CCBH4851.gb
## New name CCBH4851.gb1
# before:  
#           sequence:RefSeq:protein TadG
# after:
#           sequence:RefSeq:WP_015648121.1
    
    my $gfile=$_[0];

    my $ec_file_gbk="ec_ccbh4851_gbk.data";
    my $ec_file_gbk_transfer="ec_ccbh4851_gbk_transfer.data";
    if(-e $ec_file_gbk){
	$old_ec_file_gbk=$ec_file_gbk."\.old";
	$str="mv $ec_file_gbk $old_ec_file_gbk";
	print "       Moving existing $ec_file_gbk to $old_ec_file_gbk\n";
	system($str);
    }	
    $nccbh_ec_brenda=0;
    $nccbh_ec_transfer=0;
    my $ec_file_brenda="ec_ccbh4851_BRENDA.data";
    my $ec_file_brenda_transfer="ec_ccbh4851_BRENDA_transfer.data";
    if(-e $ec_file_brenda){
	$old_ec_file_brenda=$ec_file_brenda."\.old";
	$str="mv $ec_file_brenda $old_ec_file_brenda";
	print "       Moving existing $ec_file_brenda to $old_ec_file_brenda\n";
	system($str);
    }	
    write_file($ec_file_gbk,"");
    write_file($ec_file_gbk_transfer,"");
    write_file($ec_file_brenda,"");

    $ffile=$gfile.".fasta";
    print "-------------------------------------------------------------------------------------\n";
    print "CCBH4851 GBK:\n";
    print "    Parsing gbk file $gfile and creating fasta file $ffile\n";
    $aastr="";
    write_file($ffile,\$aastr);
    
    $CDS_re_from= "     CDS             ";
    $CDS_re_to=   "\@\@END\@\@\@\n";
    $CDS_re_to.=  "&&&&&CDS#############";
    $gene_re_from="     gene            ";
    $gene_re_to=  "\@\@END\@\@\@\n";
    $gene_re_to.= "&&&&&gene############";
    $aa_end_re_from="BASE COUNT";
    $aa_end_re_to="\@\@END\@\@\@\n";
    $aa_end_re_to.="BASE COUNT";
    $tRNA_re_from="     tRNA            ";
    $tRNA_re_to=  "\@\@END\@\@\@\n";
    $tRNA_re_to.= "&&&&&tRNA############";
    
    $gstr=read_file($gfile);
    $gstr =~ s/$CDS_re_from/$CDS_re_to/g;
    $gstr =~ s/$gene_re_from/$gene_re_to/g;
    $gstr =~ s/$aa_end_re_from/$aa_end_re_to/g;
    $gstr =~ s/$tRNA_re_from/$tRNA_re_to/g;
    $gstr =~ s/\t+/ /g;
    $gstr =~ s/\n+/ /g;
    $gstr =~ s/\r+/ /g;
    $gstr =~ s/\s+/ /g;
    $gstr =~ s/BASE COUNT.*//g;
    
    @cdsinfo=($gstr=~/&&&&&CDS#############(.+?)\@\@END\@\@\@/g);
    $ncds=0;
    $ntr=0;
    $nccbh_genes_with_gid=0;
    $nrefseq=0;
    for(my $i=0;$i<scalar(@cdsinfo);$i++){
	if(($cdsinfo[$i] =~ /(^complement\(.*?\))/)||($cdsinfo[$i] =~ /(^join\(.*?\))/)||($cdsinfo[$i] =~ /(^\d+\.\.\d+)/)){
	    $cds_pos=$1;
	    $cds_pos=~s/\.\./\_\_/g;
	    $cds_pos=~s/\,/\_/g;
	    $cds_pos=~s/\s+//g;
	    if(defined $cds_id{$cds_pos}){
		print "CDS $cds_id{$cds_pos} already defined\n";
		exit(0);
	    }else{
		$ncds++;
		$cds_id{$cds_pos}=$ncds;
#		print "XXXCDS $i=$cds_pos ($cdsinfo[$i])\n";
	    }
	}else{
	    print "No info in $cdsinfo[$i] ?\n";
	    exit(0);
	}
	$fasta_header=">gnl\|ccbh4851";
	if($cdsinfo[$i] =~ /\/translation=\"(.*?)\"/){
	    $aaseq=$1;
	    $aaseq =~ s/\s+//g;
	}else{
	    print "No aastr in $cdsinfo[$i] ?\n";
	    exit(0);
	}
	$fasta_CDS=$cds_pos;

	$id_code=">gnl\|ccbh4851\|CDS\:$cds_pos";
	if($cdsinfo[$i] =~ /\/EC_number=\"(.*?)\"/){
	    $ec_str_gbk=$1;
	    @ecl=($ec_str_gbk =~ /(\d+\.\d+\.\d+\.\d+)/g);
	    if(scalar(@ecl)>0){
		$ec_gbk_str="--------------- GBK EC\:";
		$id_code.="\|EC\:";
		$fasta_ec_gbk="EC GBK\:";
		for(my $e=0;$e<scalar(@ecl);$e++){
		    print "XXXEC $ecl[$e]\n";
		    $ccbh_ec_gbk{$ecl[$e]}=1;
		    if($ecl[$e] =~ /(\d+\.\d+\.)(\d+\.)\d+/){
			$ecs1=$1."\-\.\-";
			$ecs2=$1.$2."-";
			$ccbh_ec_gbk{$ecs1}=1;
			$ccbh_ec_gbk{$ecs2}=1;
			$ccbh_ec_gbk_s{$ecs1}=1;
			$ccbh_ec_gbk_s{$ecs2}=1;
		    }
#		    my $ec_tr=check_ec_transfer($ecl[$e]);
		    my $ec_tr=$ecl[$e];
		    if($ec_tr ne $ecl[$e]){
			$ccbh_ec_gbk{$ec_tr}=1;
			if($ec_tr =~ /(\d+\.\d+\.)(\d+\.)\d+/){
			    $ecs1=$1."\-\.\-";
			    $ecs2=$1.$2."-";
			    $ccbh_ec_gbk{$ecs1}=1;
			    $ccbh_ec_gbk{$ecs2}=1;
			    $ccbh_ec_gbk_s{$ecs1}=1;
			    $ccbh_ec_gbk_s{$ecs2}=1;
			}
			$ec_gbk_str.="\[$ecl[$e](transferred to $ec_tr)\]";
			$fasta_ec_gbk.="\[$ecl[$e](transferred to $ec_tr)\]";
			$transf_str="$ecl[$e] transferred to $ec_tr\n";
			append_file($ec_file_gbk_transfer,\$transf_str);			
		    }else{
			$ec_str="\[$ecl[$e]\]\n";
			$ec_gbk_str.="\[$ecl[$e]\]";
			$fasta_ec_gbk.="\[$ecl[$e]\]";
			$ec_str="\[$ecl[$e]\]";
			append_file($ec_file_gbk,\$ec_str);		       
		    }
		    $id_code.="$ecl[$e]";
		}
	    }else{
		$ec_gbk_str="--------------- NO GBK EC";
		$fasta_ec_gbk="[NO GBK EC]";
	    }
	}else{
	    $ec_gbk_str="--------------- NO GBK EC";
	    $fasta_ec_gbk="[NO GBK EC]";
	}
#	$fasta_ec_brenda=get_ec_BRENDA($aaseq,$cds_pos,$ec_file_brenda,$ec_file_brenda_transfer);
#	$fasta_ec_brenda="No BRENDA EC";
#	$fasta_ec="$fasta_ec_gbk\|$fasta_ec_brenda";
	$fasta_ec="$fasta_ec_gbk";

	$fasta_prod="no product";
	$gbk_prod="";
	if($cdsinfo[$i]=~/\/product=\"(.*?)\"/){
	    $prod=$1;
	    $gbk_prod="($prod)";
	    $prod =~ s/\s+/\_/g;
	    $fasta_prod=$prod;
	    $id_code.="\|$prod";
	}

	$fasta_gene="\[no gene\]";
	$gbk_gene="";
	if($cdsinfo[$i]=~/\/gene=\"(.*?)\"/){
	    $gene=$1;
	    $gene=~s/\s+//g;
	    $gbk_gene="($gene)";
	    $fasta_gene=$gene;
	    $id_code.="\|$gene";
	}
	
	$fasta_refseq="(no RefSeq)";
	if($cdsinfo[$i] =~ /RefSeq\:(.*?)\"/){
	    $rstr=$1;
	    $rstr=~s/\s+//g;
#	    $rstr=~s/\s+$//;
	    if(defined $ref{$rstr}){
		print "XXXREF Why? $rstr\n";
#		exit(0);
	    }else{
		$nrefseq++;
		$ref{$rstr}=$nrefseq;
	    }
	    $id_code.=" RefSeq\:$rstr";
	    my $refseq="$rstr";
	    if(defined $gid_refseq{$rstr}){
		$ccbh_gid{$gid_refseq{$rstr}}=1;
		$refseq.="\/$gid_refseq{$rstr}";
		$id_code.="\($gid_refseq{$rstr}\)";
		unless(defined $ccbh_gene_fasta{$gid_refseq{$rstr}}){
		    $aastr=$fasta_header."\|".$rstr."\|".$gid_refseq{$rstr}."\n";
		    $aastr.="$aaseq\n";
		    $ccbh_gene_fasta{$gid_refseq{$rstr}}=$aastr;
		    $ccbh_gene_str{$gid_refseq{$rstr}}="\[$rstr\]\[$gid_refseq{$rstr}\] $gbk_prod $gbk_gene\n";
		    $nccbh_genes_with_gid++;
		}
	    }else{
		$ccbh_gene_no_gid{$rstr}="\[$rstr\] $gbk_prod $gbk_gene\n";
	    }
	    $fasta_refseq="RefSeq\:$refseq";
	    $ntr++;
	}
	    
	$fasta_prod="no product";
	if($cdsinfo[$i]=~/\/product=\"(.*?)\"/){
	    $prod=$1;
	    $prod =~ s/\s+/\_/g;
	    $fasta_prod=$prod;
	    $id_code.="\|$prod";
	}

	$fasta_gene="\[no gene\]";
	if($cdsinfo[$i]=~/\/gene=\"(\w+)\"/){
	    my $gene=$1;
	    $gene=~s/\s+//g;
	    $fasta_gene=$gene;
	    $id_code.="\|$gene";

	    if(defined $gid_gene{uc($gene)}){
		if(defined $multi_gid_gene{uc($gene)}){
		    print "XXXMGID mgid $gene ";
		    foreach my $mgid(keys %{$gid_genes{uc($gene)}}){
			$ccbh_gid{$mgid}=1;
			print "($mgid)";
		    }
		    print "\n";
		}else{
		    $ccbh_gid{$gid_gene{uc($gene)}}=1;
		}
	    }
	}

	$aastr=$fasta_header."\|".$fasta_CDS."\|".$fasta_refseq."\|".$fasta_prod."\n";
	$aastr.="$aaseq\n";

	$fasta_info=$fasta_header."\|".$fasta_CDS."\|".$fasta_refseq."\|".$fasta_ec."\|".$fasta_prod."\|".$fasta_gene."\n";
	$fasta_info.="$aaseq\n";
	append_file($ffile,\$fasta_info);
#	append_file($ffile,\$aastr);
#	print "=================================================================================\n";
#	print "$fasta_info";
    }   
    print "Found $ncds coding sequences (CDS) and respective AA patterns.\n";
}

	
# gene 
# A region of biological interest identified as a gene and for
# which a name has been assigned. The base span for the gene feature
# is dependent on the furthest 5' and 3' features. 
# Entrez Search Field: Feature Key [FKEY] 
# Example:  gene    687..3158 
# The feature extends from base 687 through base 3158 in the sequence shown
# <     indicates partial on the 5' end
# >     indicates partial on the 3' end
# (complement)  indicates that the feature is on the complementary strand

# CDS 
# Coding sequence; region of nucleotides that corresponds with the
# sequence of amino acids in a protein (location includes start and
# stop codons). The CDS feature includes an amino acid translation.

# translation 
# The amino acid translation corresponding to the nucleotide
# coding sequence (CDS). In many cases, the translations are conceptual.

# ================================================================================= #
sub read_PAO1_geneprot_table{
# Gene set from PAO1 with ncbi accession number IDs
# http://www.ncbi.nlm.nih.gov/genome/proteins/187?genome_assembly_id=165496

# Info on file:
#Locus Tag@RefSeq Accession@GI Number@UniProtKB Accession@UniProtKB ID@Uniparc@UniRef100 ID@UniRef90 ID@UniRef50 ID@NCBI Locus Tag@NCBI Old Locus Tag@Affymetrix@Feature Type@Coordinates@Gene Length@Gene Name@Product Description@Product Name Confidence@Estimated MW (kDa)@Estimated Isoelectric Point (pI)@Estimated Charge (ph7)

#PA0034@NP_248724.1@15595232@P24908@TRPO_PSEAE@UPI000013763B@UniRef100_P24908@UniRef90_P24908@UniRef50_P24908@PA0034@@PA0034_at@CDS@36278..36901(-)@624@@probable two-component response regulator@Class 3@22.9@10.15@6.92    

    my $table_file=$_[0];
    print "----------------------------------------------------------------------------------------\n";
    print "PROTEINS:\n";
    print "    Reading protein table $table_file.\n";
    my @component=read_file($table_file);
    $nc=scalar(@component);
#
# Maybe should add a mapping between also genes and locus tags here to see further on
# if genes present in the CCBH annotation maps to locus tags on iPAO1 GPR relations
#    
    for($i=0;$i<$nc;$i++){
	unless($component[$i] =~ /^\#/){
	    chomp($component[$i]);
	    undef @info;
	    @info=split/@/,$component[$i];
	    if($info[0]=~/\w+/){
		chomp($info[0]);
		print "GID=$info[0]: ";
		$locus_tag=uc($info[0]);
		if(defined $info[1]){
		if($info[1]=~/\w+/){
		    chomp($info[1]);
		    print " RefSeq=$info[1]";
		    my $refseq=uc($info[1]);
#		    if(defined $gid_refseq{$refseq}){
#			if($gid_refseq{$refseq} ne $locus_tag){
#			    print "  Nooooo ($refseq) $gid_refseq{$refseq} ($locus_tag)\n";
#			    exit(0);
#			}
#		    }
		    $gid_refseq{$refseq}=$locus_tag;
		}
		}
		if(defined $info[15]){
		if($info[15]=~/\w+/){
		    chomp($info[15]);
		    print " gene=$info[15]";
		    my $gname=uc($info[15]);
		    if((defined $gid_gene{$gname})&&($gid_gene{$gname} ne $locus_tag)){ # there are repeated lines on file
			unless(defined $multi_gid_gene{$gname}){
			    $multi_gid_gene{$gname}=1;
			    $gid_genes{$gname}{$gid_gene{$gname}}=1;
			}
			$gid_genes{$gname}{$locus_tag}=1;
#			if($gid_gene{$gname} ne $locus_tag){
#			    print "  Nooooo (gene $gname, locus_gid $gid_gene{$gname}) also at locus_gid $locus_tag\n";
#			    exit(0);
#			}
		    }else{
			$gid_gene{$gname}=$locus_tag;
		    }
		}
		}
		print "\n";
	    }
	}
    }
}

###################################################################
sub get_ec_BRENDA{
    my $aaseq=$_[0];
    my $idprot=$_[1];
    my $ec_file=$_[2];
    my $ec_file_transfer=$_[3];
    my $aa_str="";
    my $aa_file=$database_dir."id-info_".$idprot.".html";

    undef %ec_list;
    
    if(-e $aa_file){
	$aa_str = read_file($aa_file);
    }
    unless((-e $aa_file)||($aa_str =~ /NOURL/)){
	my $url = $brenda_front.$aaseq.$brenda_rear; # Get list of ECs for gene idprot at BRENDA
	my $response = $browser->get("$url");
	if($response->{success}){
	    $aa_str=$response->{content};
	    print "Dump result in file $aa_file\n";	
	    write_file($aa_file,\$aa_str);
	}else{
	    $aa_str="NOURL";
#	    die "Can't download $url -- ", $response->{status_line};
	}
    }
    my @aa_ec = ($aa_str =~ /enzyme\.php\?ecno\=(.*?)\"/g);
    my $ec_str_brenda="EC BRENDA\:";
    if(scalar(@aa_ec)>0){
	for(my $i=0;$i<scalar(@aa_ec);$i++){
	    chomp($aa_ec[$i]);
	    my $ec_number=$aa_ec[$i];
	    unless((defined $ccbh_ec_brenda{$ec_number})||(defined $ccbh_ec_transfer{$ec_number})){
		my $ec_transf=check_ec_transfer($ec_number);
		if($ec_transf ne $ec_number){
		    $nccbh_ec_transfer++;
		    $ccbh_ec_transfer{$ec_number}=$ec_transf;
		    $ccbh_ec_brenda{$ec_number}=$ec_transf;
		    append_file($ec_file_transfer,"$ec_transf ((((replaces $ec_number))))\n");
		    unless(defined $ccbh_ec_brenda{$ec_transf}){
			$nccbh_ec_brenda++;
			$ccbh_ec_brenda{$ec_transf}=$ec_number;   
			if($ec_transf =~ /(\d+\.\d+\.)(\d+\.)\d+/){
			    $ecs1=$1."\-\.\-";
			    $ecs2=$1.$2."-";
			    $ccbh_ec_brenda{$ecs1}=1;
			    $ccbh_ec_brenda{$ecs2}=1;
			    $ccbh_ec_brenda_s{$ecs1}=1;
			    $ccbh_ec_brenda_s{$ecs2}=1;
			}
		    }
		}else{
		    $nccbh_ec_brenda++;
		    $ccbh_ec_brenda{$ec_number}=$ec_number;   
		    append_file($ec_file,"$ec_number\n");
		    if($ec_number =~ /(\d+\.\d+\.)(\d+\.)\d+/){
			$ecs1=$1."\-\.\-";
			$ecs2=$1.$2."-";
			$ccbh_ec_brenda{$ecs1}=1;
			$ccbh_ec_brenda{$ecs2}=1;
			$ccbh_ec_brenda_s{$ecs1}=1;
			$ccbh_ec_brenda_s{$ecs2}=1;
		    }
		}
	    }
	    unless(defined $ec_list{$ccbh_ec_brenda{$ec_number}}){
		$ec_list{$ccbh_ec_brenda{$ec_number}}=1;
		if(defined $ccbh_ec_transfer{$ec_number}){
		    $ec_str_brenda.="\[$ec_number (transferred to $ccbh_ec_brenda{$ec_number})\]";
		}else{
		    $ec_str_brenda.="\[$ccbh_ec_brenda{$ec_number}\]";
		}
	    }
	}
	$out_str=$ec_str_brenda;
    }else{
	$out_str="\[NO BRENDA EC\]";
    }
    return($out_str);
}
###################################################################
sub check_ec_transfer{
    my $ec_number=$_[0];
    my $ec_transfer=$ec_number;		
    my $ec_str_tr="";
    my $ec_file = $database_dir."ec-info_".$ec_number.".html";
    if(-e $ec_file){
	$ec_str_tr = read_file($ec_file);
    }else{
	my $url = $EC_to_rxn_front.$ec_number; # Fetch info at KEGG db for ec_number
	my $response = $browser->get("$url"); 
	if($response->{success}){
	    $ec_str_tr=$response->{content};
	    print "Dump result in file $ec_file\n";	
	    write_file($ec_file,\$ec_str_tr);
	}
    }
    if($ec_str_tr =~ m/Transferred to \<a href\=\"\/dbget\-bin\/www\_bget\?ec\:.*?\"/){
	my @aa_tr=($ec_str_tr=~ m/Transferred to \<a href\=\"\/dbget\-bin\/www\_bget\?ec\:(.*?)\"/g);
	my $ntr=scalar(@aa_tr);
	if($ntr!=1){
	    print "oooow $ec_number\n";
	    exit(0);
	}else{
	    $ec_transfer=$aa_tr[0];
	}
    }
    return($ec_transfer);
}

#############################################################################################
sub find_orthologues_tcdb{
    my $ortho_file=$_[0];
    my $elim=$_[1];
    my $tc_file="tc_ccbh4851.data";
    write_file($tc_file,"");
    $ezero=2.0e-308;
    $ntc_ccbh=0;
# evalue bitscore qseqid sseqid qlen slen length qstart qend sstart send
# 0.051	35.0	gnl|ccbh4851|483__2027|RefSeq:NP_064721.1/PA0001|[NO	gnl|TC-DB|F7X0P7|9.B.110.2.1	514	968	58	213	270	768	825    
    @blast_data=read_file($ortho_file);
    for($i=0;$i<scalar(@blast_data);$i++){
	chomp($blast_data[$i]);
	unless($blast_data[$i]=~/^\#/){
	    @info = split/\s+/,$blast_data[$i];
	    if($info[0]>0.0){
		$evalue=(1.0)*$info[0];
	    }else{
		$evalue=$ezero;
	    }
	    chomp($info[2]);
	    chomp($info[3]);
	    if($evalue < $elim){	    
		unless((defined $tcdb_ortho{$info[2]}{$info[3]})&&($tcdb_ortho{$info[2]}{$info[3]}<$evalue)){
		    $tcdb_ortho{$info[2]}{$info[3]}=$evalue;
		    if(defined $tcdb_ortho{$info[3]}{$info[2]}){
			if($info[2] =~ /ccbh4851/){
			    my @tci=split/\|/,$info[3];
			    $tcstr=$tci[3];
			}elsif($info[3] =~ /ccbh4851/){
			    my @tci=split/\|/,$info[2];
			    $tcstr=$tci[3];
			}else{
			    print "Error, no gene: $component[$i]";
			    exit;
			}
			if($tcstr=~/(\d+\.\w+\.\d+)(\.\d+)(\.\d+)/){
			    $s1=$1;
			    $s2=$2;
			    $s3=$3;
			    $tc=$s1.$s2.$s3;
			    $stc1=$s1;
			    $stc2=$s1.$s2;
			    unless(defined $ccbh_tc{$tc}){
				$ntc_ccbh++;
				$ccbh_tc_specific{$stc1}{$tc}=1;
				$ccbh_tc_specific{$stc2}{$tc}=1;
				$ccbh_tc_specific{$tc}{$tc}=1;
				$ccbh_tc{$stc1}=1;
				$ccbh_tc{$stc2}=1;
				$ccbh_tc{$tc}=1;
				print "XXXTC ($tc)\n";
			    }
			}else{
			    print "Not specific $component[$i] ?\n";
			    exit;
			}
		    }
		}
	    }
	}
    }
    print "CCBH4851 annotation have $ntc_ccbh TC numbers (e-value < $elim with tcdb genes)\n";
}


###################################################################
    

