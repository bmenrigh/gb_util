#!/usr/bin/perl

use strict;
use warnings;

use CGI;


# Global vars
my $comp_edge_same = 0;
my $comp_face2_same = 0;
my $disas_allow_xn = 0;
my $disas_allow_concat = 0;
my $max_concat_trials = 100000;
my $cur_concat_trials = 0;

my $cgi_var = new CGI;

print $cgi_var->header('text/html');
print $cgi_var->start_html('Sequence utility for Gelatinbrain\'s Virtual Magic Polyhedra');

my %gb_to_kd = (
    'A' => 'F',
    'B' => 'U',
    'C' => 'R',
    'D' => 'DR',
    'E' => 'DL',
    'F' => 'L',
    'G' => 'BR',
    'H' => 'UR',
    'I' => 'UL',
    'J' => 'BL',
    'K' => 'D',
    'L' => 'B'
    );

my %gb_to_kvtd = (
    'BHI' => 'U',
    'BCH' => 'UR',
    'ABC' => 'UFR',
    'ABF' => 'UFL',
    'BFI' => 'UL',
    'FIJ' => 'L',
    'EFJ' => 'FDL',
    'AEF' => 'FL',
    'ADE' => 'F',
    'ACD' => 'FR',
    'CDG' => 'FDR',
    'CGH' => 'R',
    'EJK' => 'DL',
    'DEK' => 'D',
    'DGK' => 'DR',
    'HIL' => 'B',
    'IJL' => 'BL',
    'JKL' => 'DBL',
    'GKL' => 'DBR',
    'GHL' => 'BR'
    );

sub html_warn {
    my $warnstr = shift;

    #print '<p><b><font color="red">WARN: ', $warnstr, '</font></b></p>', "\n";
}

sub parse_sequence {
    my $move_input = shift;

    my @parsed_moves;

    while ($move_input =~ m/(?:\/\*\d+\*\/)?\s*([A-Z]+\'?\d*(&\d+)?)\s*,?/ig) {
        push @parsed_moves, $1;
    }

    return @parsed_moves;
}


sub seq_disas {
    my $mref = shift;

    my @moves = @{$mref};
    my $mlen = scalar @moves;

    html_warn('head of seq_disas with [' . join(', ', @moves) . ']');

    # There is a need to minimize the number of concatentation operations
    # otherwise sequences that barely use concatenation might be decoded
    # in a way that suggests that they use it a lot, missing a lot of their
    # structure.
    # To solve that, for each depth of the recursion, only return a result
    # if either that results uses no concatenation or if it uses less concatenation
    # than any other result.
    my $candidate_solution;


    # Case 1: we only have a single move
    if (scalar $mlen == 1) {
	html_warn('Got to 1-move: ' . $moves[0]);

	return [(1, 1, 0, $mlen, $moves[0])];
    }

    # Case 2: check for commutator
    if (($mlen % 2 == 0) && ($mlen >= 4)) {

	for (my $xlen = 1; ($mlen - (2 * $xlen)) >= 2; $xlen++) {

	    my $ylen = ($mlen - (2 * $xlen)) / 2;

	    my $xref = [@moves[0 .. ($xlen - 1)]];
	    my $xiref = [@moves[($xlen + $ylen) .. (((2 * $xlen) + $ylen) - 1)]];

	    my $yref = [@moves[$xlen .. (($xlen + $ylen) - 1)]];
	    my $yiref = [@moves[((2 * $xlen) + $ylen) .. ((2 * ( $xlen + $ylen)) - 1)]];

	    html_warn('Trying commutator with ' .
		      '[' . join(', ', @{$xref}) . ']' . ' ' .
		      '[' . join(', ', @{$yref}) . ']' . ' ' .
		      '[' . join(', ', @{$xiref}) . ']' . ' ' .
		      '[' . join(', ', @{$yiref}) . ']');

	    # Did we get x and x' ?
	    next unless (are_inverses($xref, $xiref));

	    # Did we get y and y' ?
	    next unless (are_inverses($yref, $yiref));

	    # Okay now we need to be able to disassemble x and y and get something meaningful
	    my ($xres, $xstd, $xconcat, $xfrm, $xstr) = @{seq_disas($xref)};
	    my ($yres, $ystd, $yconcat, $yfrm, $ystr) = @{seq_disas($yref)};

	    my $this_std = 1;

	    $this_std = 0 unless (($xstd == 1) && ($ystd == 1));

	    if (($xres == 1) && ($yres == 1)) {

		# Minimize concat
		my $total_concat = $xconcat + $yconcat;
		if ($total_concat == 0) {

		    # Great, we got a commutator, return it
		    return [(1, $this_std, $total_concat, '[' . $xfrm . ',' . $yfrm . ']', '[' . $xstr . ',' . $ystr . ']')];
		}
		else {
		    if ((not defined $candidate_solution) ||
			($total_concat < $candidate_solution->[2])) {
			$candidate_solution = [(1, $this_std, $total_concat, '[' . $xfrm . ',' . $yfrm . ']', '[' . $xstr . ',' . $ystr . ']')];
		    }

		    # Check if we should abandon backtracking
		    $cur_concat_trials++;
		    if ($cur_concat_trials >= $max_concat_trials) {
			return $candidate_solution;
		    }
		}
	    }
	    else {
		# Either x or y can't be disassembled
		next;
	    }
	}
    }


    # Case 3: check for conjugate
    if ($mlen >= 3) {

	for (my $xlen = 1; ($mlen - (2 * $xlen)) >= 1; $xlen++) {

	    my $ylen = ($mlen - (2 * $xlen));

	    my $xref = [@moves[0 .. ($xlen - 1)]];
	    my $xiref = [@moves[($xlen + $ylen) .. (((2 * $xlen) + $ylen) - 1)]];

	    my $yref = [@moves[$xlen .. (($xlen + $ylen) - 1)]];

	    html_warn('Trying conjugate with ' .
		      '[' . join(', ', @{$xref}) . ']' . ' ' .
		      '[' . join(', ', @{$yref}) . ']' . ' ' .
		      '[' . join(', ', @{$xiref}) . ']');

	    # Did we get x and x' ?
	    next unless (are_inverses($xref, $xiref));

	    # Okay now we need to be able to disassemble x and y and get something meaningful
	    my ($xres, $xstd, $xconcat, $xfrm, $xstr) = @{seq_disas($xref)};
	    my ($yres, $ystd, $yconcat, $yfrm, $ystr) = @{seq_disas($yref)};

	    my $this_std = 1;

	    $this_std = 0 unless (($xstd == 1) && ($ystd == 1));

	    if (($xres == 1) && ($yres == 1)) {

		# Minimize concat
		my $total_concat = $xconcat + $yconcat;
		if ($total_concat == 0) {

		    # Great, we got a conjugate, return it
		    return [(1, $this_std, $total_concat, '[' . $xfrm . ':' . $yfrm . ']', '[' . $xstr . ':' . $ystr . ']')];
		}
		else {
		    if ((not defined $candidate_solution) ||
			($total_concat < $candidate_solution->[2])) {
			$candidate_solution = [(1, $this_std, $total_concat, '[' . $xfrm . ':' . $yfrm . ']', '[' . $xstr . ':' . $ystr . ']')];
		    }

		    # Check if we should abandon backtracking
		    $cur_concat_trials++;
		    if ($cur_concat_trials >= $max_concat_trials) {
			return $candidate_solution;
		    }
		}
	    }
	    else {
		# Either x or y can't be disassembled
		next;
	    }
	}
    }


    # Case 4: check for []xN (if allowed)
    if (($disas_allow_xn == 1) && ($mlen >= 2)) {

	for (my $i = 1; $i < $mlen; $i++) {

	    if ($mlen % $i == 0) {
		my $n = $mlen / $i;

		my $all_same = 1;
		for (my $j = 0; $j < ($n - 1); $j++) {

		    my $x1ref = [@moves[($j * $i) .. ((($j + 1) * $i) - 1)]];
		    my $x2ref = [@moves[(($j + 1) * $i) .. ((($j + 2) * $i) - 1)]];

		    $all_same = 0 unless (are_same($x1ref, $x2ref));

		}

		if ($all_same == 1) {
		    my ($xres, $xstd, $xconcat, $xfrm, $xstr) = @{seq_disas([@moves[0 .. ($i - 1)]])};

		    if ($xres == 1) {

			# Minimize concat
			#my $total_concat = $xconcat * $n;  # I can't decide if * $n is better or not
			my $total_concat = $xconcat;
			if ($total_concat == 0) {

			    # Great, we []xN, return it
			    return [(1, 0, $total_concat, $xfrm . 'x' . $n, $xstr . 'x' . $n)];
			}
			else {
			    if ((not defined $candidate_solution) ||
				($total_concat < $candidate_solution->[2])) {
				$candidate_solution = [(1, 0, $total_concat, $xfrm . 'x' . $n, $xstr . 'x' . $n)];
			    }

			    # Check if we should abandon backtracking
			    $cur_concat_trials++;
			    if ($cur_concat_trials >= $max_concat_trials) {
				return $candidate_solution;
			    }

			}
		    }
		}
	    }
	}
    }

    # Case 5: check for concatenation (if allowed)
    if (($disas_allow_concat == 1) && ($mlen >= 2)) {

	for (my $i = 1; ($mlen - $i) > 0; $i++) {

	    my $xref = [@moves[0 .. ($i - 1)]];
	    my $yref = [@moves[$i .. ($mlen - 1)]];

	    my ($xres, $xstd, $xconcat, $xfrm, $xstr) = @{seq_disas($xref)};
	    my ($yres, $ystd, $yconcat, $yfrm, $ystr) = @{seq_disas($yref)};

	    if (($xres == 1) && ($yres == 1)) {

		# Minimize concat
		my $total_concat = $xconcat + $yconcat + 1;

		if ((not defined $candidate_solution) ||
		    ($total_concat < $candidate_solution->[2])) {
		    $candidate_solution = [(1, 0, ($xconcat + $yconcat + 1), '[' . $xfrm . ' ' . $yfrm . ']', '[' . $xstr . ' ' . $ystr . ']')];
		}

		# Check if we should abandon backtracking
		$cur_concat_trials++;
		if ($cur_concat_trials >= $max_concat_trials) {
		    return $candidate_solution;
		}
	    }
	}
    }

    if (defined $candidate_solution) {
	return $candidate_solution;
    }
    else {
	# Case default: nothing matched
	return [(0, 0, 0, '', '')];
    }

}


sub are_inverses {
    my ($lref1, $lref2) = @_;

    return 0 if (scalar @{$lref1} != scalar @{$lref2});

    my $len = scalar @{$lref1};
    for (my $i = 0; $i < $len; $i++) {
	return 0 unless (liberal_inv_comp($lref1->[$i], $lref2->[($len - $i) - 1]));
    }

    return 1;
}


sub are_same {
    my ($lref1, $lref2) = @_;

    return 0 if (scalar @{$lref1} != scalar @{$lref2});

    my $len = scalar @{$lref1};
    for (my $i = 0; $i < $len; $i++) {
	return 0 unless (liberal_comp($lref1->[$i], $lref2->[$i]));
    }

    return 1;
}


# This function compares moves to check if they are inverses of each other
# allowing for fuzzy things like edge-turns to be self-inverses
sub liberal_inv_comp {
    my ($m1, $m2) = @_;

    my $norm1 = norm_name($m1);
    my $norm2 = norm_name($m2);

    my ($name1, $dir1, $amt1, $slice1, $nlen1);
    my ($name2, $dir2, $amt2, $slice2, $nlen2);

    if ($norm1 =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name1, $dir1, $amt1, $slice1, $nlen1) = ($1, $2, $3, $4, length $1);
    }
    else {
	return 0;
    }

    if ($norm2 =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name2, $dir2, $amt2, $slice2, $nlen2) = ($1, $2, $3, $4, length $1);
    }
    else {
	return 0;
    }

    # Fuzzy-check faces
    if (($comp_face2_same == 1) && ($nlen1 == 1) && ($amt1 eq '2')) {
	return (($name1 eq $name2) && ($amt1 eq $amt2) && ($slice1 eq $slice2));
    }

    # Fuzzy-check edges
    if (($comp_edge_same == 1) && ($nlen1 == 2)) {
	return (($name1 eq $name2) && ($amt1 eq $amt2) && ($slice1 eq $slice2));
    }

    # Well neither of those panned out so just do the regular check
    return (($name1 eq $name2) && ($dir1 ne $dir2) && ($amt1 eq $amt2) && ($slice1 eq $slice2));

}


sub liberal_comp {
    my ($m1, $m2) = @_;

    my $norm1 = norm_name($m1);
    my $norm2 = norm_name($m2);

    my ($name1, $dir1, $amt1, $slice1, $nlen1);
    my ($name2, $dir2, $amt2, $slice2, $nlen2);

    if ($norm1 =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name1, $dir1, $amt1, $slice1, $nlen1) = ($1, $2, $3, $4, length $1);
    }
    else {
	return 0;
    }

    if ($norm2 =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name2, $dir2, $amt2, $slice2, $nlen2) = ($1, $2, $3, $4, length $1);
    }
    else {
	return 0;
    }

    # Fuzzy-check faces
    if (($comp_face2_same == 1) && ($nlen1 == 1) && ($amt1 eq '2')) {
	return (($name1 eq $name2) && ($amt1 eq $amt2) && ($slice1 eq $slice2));
    }

    # Fuzzy-check edges
    if (($comp_edge_same == 1) && ($nlen1 == 2)) {
	return (($name1 eq $name2) && ($amt1 eq $amt2) && ($slice1 eq $slice2));
    }

    # Well neither of those panned out so just do the regular check
    return (($name1 eq $name2) && ($dir1 eq $dir2) && ($amt1 eq $amt2) && ($slice1 eq $slice2));

}


sub norm_name {
    my $move = shift;

    my ($name, $dir, $amt, $slice);

    if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name, $dir, $amt, $slice) = ($1, $2, $3, $4);

	my $clean_name = join('', (sort split(//, lc($name))));

	return $clean_name . $dir . $amt . $slice;
    }
    else {
	return $move;
    }
}


sub m_inv {
    my $move = shift;

    my ($name, $dir, $amt, $slice);

    if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name, $dir, $amt, $slice) = ($1, $2, $3, $4);
    }
    else {
	print 'Unable to parse move: ', $move, "\n";
	return;
    }

    return $name . ($dir eq '\''? '' : '\'') . $amt . $slice;
}


sub m_mir_cube1 {
    my $move = shift;

    my ($name, $dir, $amt, $slice);

    if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name, $dir, $amt, $slice) = ($1, $2, $3, $4);
    }
    else {
	#print 'Unable to parse move: ', $move, "\n";
	return;
    }

    $name =~ tr/UFLBRDuflbrd/UFRBLDufrbld/;

    return $name . ($dir eq '\''? '' : '\'') . $amt . $slice;
}


sub m_mir_cube2 {
    my $move = shift;

    my ($name, $dir, $amt, $slice);

    if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name, $dir, $amt, $slice) = ($1, $2, $3, $4);
    }
    else {
	#print 'Unable to parse move: ', $move, "\n";
	return;
    }

    $name =~ tr/UFLBRDuflbrd/URBLFDurblfd/;

    return $name . ($dir eq '\''? '' : '\'') . $amt . $slice;
}


sub m_mir_dodeca {
    my $move = shift;

    my ($name, $dir, $amt, $slice);

    if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	($name, $dir, $amt, $slice) = ($1, $2, $3, $4);
    }
    else {
	#print 'Unable to parse move: ', $move, "\n";
	return;
    }

    $name =~ tr/ABCDEFGHIJKLabcdefghijkl/ABFEDCJIHGKLabfedcjihgkl/;

    return $name . ($dir eq '\''? '' : '\'') . $amt . $slice;
}


sub all_cube_rotations {
    my $seq = shift;

    my %all_cube_rotations_list = ();

    all_cube_r(\%all_cube_rotations_list, 'UFRBLD', $seq);

    my $all_str = '';
    foreach my $or (sort keys %all_cube_rotations_list) {
	$all_str .= sprintf('%s -> %s : %s',
			    'UFRBLD',
			    $or,
			    $all_cube_rotations_list{$or}) . "\n";

    }

    return $all_str;
}


sub all_cube_r {
    my $crhref = shift;
    my $or = shift;
    my $seq = shift;

    my $new_or = cube_rot1($or);
    my $new_seq = cube_rot1($seq);

    unless (exists $crhref->{$new_or}) {
	$crhref->{$new_or} = $new_seq;

	all_cube_r($crhref, $new_or, $new_seq);
    }

    $new_or = cube_rot2($or);
    $new_seq = cube_rot2($seq);

    unless (exists $crhref->{$new_or}) {
	$crhref->{$new_or} = $new_seq;

	all_cube_r($crhref, $new_or, $new_seq);
    }

}


sub cube_rot1 {
    my $seq = shift;

    $seq =~ tr/UFLBRDuflbrd/URFLBDurfldb/;

    return $seq;
}

sub cube_rot2 {
    my $seq = shift;

    $seq =~ tr/UFLBRDuflbrd/RFUBDLrfudbl/;

    return $seq;
}


sub format_output {
    my ($to, $move) = @_;

    if ($to eq 'kd') {

	my ($name, $dir, $amt, $slice);
	if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	    ($name, $dir, $amt, $slice) = ($1, $2, $3, $4);

	    #print 'Name ', $name, ' converted to: ', $gb_to_kd{$name};

	    if (exists $gb_to_kd{uc($name)}) {
		$name = $gb_to_kd{uc($name)};
		#print 'new name: ', $name, "\n";
	    }

	    if ($slice eq '&2') {
		$name = lc($name);
	    }

	    return $name . $dir . $amt;
	}
	else {
	    return $move;
	}

    }
    elsif ($to eq 'kvtd') {

	my ($name, $dir, $amt, $slice);
	if ($move =~ m/([A-Z]+)(\'?)(\d*)((?:&\d+)?)/i) {
	    ($name, $dir, $amt, $slice) = ($1, $2, $3, $4);

	    my $norm_name = norm_name($name);
	    #print 'Name ', $name, ' converted to: ', $gb_to_kd{$name};

	    if (exists $gb_to_kvtd{uc($norm_name)}) {
		$name = $gb_to_kvtd{uc($norm_name)};
		#print 'new name: ', $name, "\n";
	    }

	    if ($slice eq '&2') {
		$name = lc($name);
	    }

	    return $name . $dir . $amt;
	}
	else {
	    return $move;
	}

    }
    else {
	return $move;
    }

}


########
if (defined $cgi_var->param('action_button')) {
    if ($cgi_var->param('action_button') eq 'Reset') {
	$cgi_var->delete_all();
    }
}

print $cgi_var->h4('Enter move sequence<sup>*</sup>:');
print $cgi_var->start_form();
print '<p>', $cgi_var->textarea('input', '', 12, 64), '</p>';
print '<p>Puzzle is: ', $cgi_var->scrolling_list('in_puzzle',
				      ['in_cube',
				       'in_dodeca',
				       'in_other'],
				      ['Cube / Octahedron'],
				      1, 0,
				      {'in_cube'=>'Cube / Octahedron',
				       'in_dodeca'=>'Dodecahedron / Icosahedron',
				       'in_other'=>'Other'}, 0), '</p>', "\n";
print '<hr>', "\n";
print '<p><b>Output options:</b></p>', "\n";
print '<p>',
    $cgi_var->checkbox('use_commas', 'checked', 'on', 'Use commas'), ' | ',
    $cgi_var->checkbox('use_spaces', 'checked', 'on', 'Use spaces'), ' | ',
    $cgi_var->checkbox('use_brackets', 'checked', 'on', 'Use brackets'), ' | ',
    $cgi_var->checkbox('trailing_comma', '', 'on', 'Use trailing comma'),
    '</p>', "\n";
print '<p>Output notation: ', $cgi_var->scrolling_list('output_format',
				      ['gb',
				       'kd',
				       'kvtd'],
				      ['Gelatinbrain\'s'],
				      1, 0,
				      {'gb'=>'Gelatinbrain\'s',
				       'kd'=>'Konrad\'s (dodecahedron)',
				       'kvtd'=>'Konrad\'s (Vertex-Turning-Dodecahedron)'}, 0), '</p>', "\n";
print '<hr>', "\n";
print '<p><b>Sequence disassembly options:</b></p>', "\n";
print '<p>', $cgi_var->checkbox('assume_edge_inv', 'checked', 'on', 'Assume edge-looking moves are self-inverses (e.g. FU == FU\')'), '</p>', "\n";
print '<p>', $cgi_var->checkbox('assume_face2_inv', '', 'on', 'Assume half-face-turn-looking moves are self-inverses (e.g. U2 == U\'2)'), '</p>', "\n";
print '<p>', $cgi_var->checkbox('disas_allow_xn', 'checked', 'on', 'Allow the []xN operator (e.g. [[A, B]:C]x3)'), '</p>', "\n";
print '<p>', $cgi_var->checkbox('disas_allow_concat', '', 'on', 'Allow sequence concatenation (e.g. [A [[B:C] D]])'), '&nbsp;<sup>&dagger;</sup></p>', "\n";
print '<hr>', "\n";
print $cgi_var->submit('action_button', 'Process Sequence');
print $cgi_var->submit('action_button', 'Reset');
print $cgi_var->end_form();


if (defined $cgi_var->param('input')) {
    my $input = $cgi_var->param('input');

    my $input_format = 'in_other';
    if ((defined  $cgi_var->param('in_puzzle')) &&
	($cgi_var->param('in_puzzle') eq 'in_cube')) {
	$input_format = 'in_cube';
    }
    elsif ((defined  $cgi_var->param('in_puzzle')) &&
	($cgi_var->param('in_puzzle') eq 'in_dodeca')) {
	$input_format = 'in_dodeca';
    }

    # output options
    my $join_char = '';
    if (((defined  $cgi_var->param('use_commas'))) &&
	($cgi_var->param('use_commas') eq 'on')) {
	$join_char .= ',';
    }

    if (((defined  $cgi_var->param('use_spaces'))) &&
	($cgi_var->param('use_spaces') eq 'on')) {
	$join_char .= ' ';
    }

    my ($start_bracket, $end_bracket) = ('[', ']');
    if (((not defined  $cgi_var->param('use_brackets'))) ||
	($cgi_var->param('use_brackets') ne 'on')) {
	($start_bracket, $end_bracket) = ('', '');
    }

    my $trailing_char = '';
    if ((defined  $cgi_var->param('trailing_comma')) &&
	($cgi_var->param('trailing_comma') eq 'on')) {
	$trailing_char = ',';
    }

    my $output_format = 'gb';
    if ((defined  $cgi_var->param('output_format')) &&
	($cgi_var->param('output_format') eq 'kd')) {
	$output_format = 'kd';
    }
    elsif ((defined  $cgi_var->param('output_format')) &&
	($cgi_var->param('output_format') eq 'kvtd')) {
	$output_format = 'kvtd';
    }


    # disassembly options
    if ((defined  $cgi_var->param('assume_edge_inv')) &&
	($cgi_var->param('assume_edge_inv') eq 'on')) {
	$comp_edge_same = 1;
    }

    if ((defined  $cgi_var->param('assume_face2_inv')) &&
	($cgi_var->param('assume_face2_inv') eq 'on')) {
	$comp_face2_same = 1;
    }

    if ((defined  $cgi_var->param('disas_allow_xn')) &&
	($cgi_var->param('disas_allow_xn') eq 'on')) {
	$disas_allow_xn = 1;
    }

    if ((defined  $cgi_var->param('disas_allow_concat')) &&
	($cgi_var->param('disas_allow_concat') eq 'on')) {
	$disas_allow_concat = 1;
    }

    my @moves = parse_sequence($input);

    print '<hr>', "\n";

    my ($mret, $mstd, $mconcat, $mfrm, $mstr) = @{seq_disas(\@moves)};
    if ($mret == 1) {
	print '<p><b>Sequence disassembled:</b><br>';
	print $mstr, '</p>', "\n";
	print '<p><b>Sequence form (shorthand):</b><br>';
	print $mfrm, '</p>', "\n";
	print '<hr>', "\n";
    }
    else {
	print '<p><b>Unable to disassemble sequence.  You may need to adjust the disassembly options or',
	' the sequence may be of a non-standard form.', '</p>', "\n";
	print '<hr>', "\n";
    }

    print '<p><b>Original:</b><br>';
    print $start_bracket, join($join_char, map {format_output($output_format, $_)} @moves), $end_bracket, $trailing_char, '</p>', "\n";

    print '<p><b>Inverted:</b><br>';
    print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_inv($_)} reverse @moves), $end_bracket, $trailing_char, '</p>', "\n";

    if ($input_format eq 'in_cube') {
	print '<p><b>Mirror (cube; plane bisecting U, F, D, B faces):</b><br>';
	print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_mir_cube1($_)} @moves), $end_bracket, $trailing_char, '</p>', "\n";

	print '<p><b>Inverted Mirror (cube; plane bisecting U, F, D, B faces):</b><br>';
	print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_inv(m_mir_cube1($_))} reverse @moves), $end_bracket, $trailing_char, '</p>', "\n";

	print '<p><b>Mirror (cube; plane containing FR, BL edges):</b><br>';
	print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_mir_cube2($_)} @moves), $end_bracket, $trailing_char, '</p>', "\n";

	print '<p><b>Inverted Mirror (cube; plane containing FR, BL edges):</b><br>';
	print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_inv(m_mir_cube2($_))} reverse @moves), $end_bracket, $trailing_char, '</p>', "\n";
    }

    if ($input_format eq 'in_dodeca') {
	print '<p><b>Mirror (dodecahedron):</b><br>';
	print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_mir_dodeca($_)} @moves), $end_bracket, $trailing_char, '</p>', "\n";

	print '<p><b>Inverted Mirror (dodecahedron):</b><br>';
	print $start_bracket, join($join_char, map {format_output($output_format, $_)} map {m_inv(m_mir_dodeca($_))} reverse @moves), $end_bracket, $trailing_char, '</p>', "\n";
    }

    if ($input_format eq 'in_dodeca') {
	print '<hr>', "\n";

	if ($output_format eq 'gb') {
	    print '<p><b>Gelatinbrain\'s dodecahedron notation:</b></p>', "\n";
	    print '<p><img src="http://www.brandonenright.net/twistypuzzles/other/gb_notation.png"></p>', "\n";
	}
	elsif ($output_format eq 'kd') {
	    print '<p><b>Konrad\'s dodecahedron notation:</b></p>', "\n";
	    print '<p><img src="http://www.brandonenright.net/twistypuzzles/other/kd_notation.png"></p>', "\n";
	}
	elsif ($output_format eq 'kvtd') {
	    print '<p><b>Konrad\'s Vertex-Turning-Dodecahedron notation:</b></p>', "\n";
	    print '<p><img src="http://www.brandonenright.net/twistypuzzles/other/kvtd_notation.png"></p>', "\n";
	}
    }
    elsif ($input_format eq 'in_cube') {
	print '<hr>', "\n";

	my $orig = $start_bracket .
	    join($join_char, map {format_output($output_format, $_)} @moves) .
	    $end_bracket .
	    $trailing_char;

	print '<p><b>All cube rotations:</b></p>', "\n";
	print '<p>', $cgi_var->textarea('all_cube_r',
					all_cube_rotations($orig),
					6, 64), '</p>';

    }

}

print '<hr>', "\n";
print '<p><i><sup>*</sup> The sequence format is somewhat flexible and the /*0000*/ comments and such are acceptable.  The input notation should be Gelatinbrain\'s.</i></p>';
print '<p><i><sup>&dagger;</sup> If you allow concatenation the disassembled sequence may have other, more compact (or more standard) representations that won\'t be found.  Backtracking is used to try to minimize the use of the concatenation operator but due to the CPU load, the number of backtracking trials is currently limited to ', $max_concat_trials, '.</i></p>';

print $cgi_var->end_html();
