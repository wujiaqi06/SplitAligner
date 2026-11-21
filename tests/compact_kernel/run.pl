#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use utf8;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use Test::More;
use JSON::PP;
use Scalar::Util qw(weaken);
use File::Spec;
use SplitAligner::Bitset qw(zero_bits ids_to_bits bits_to_ids bits_to_hex hex_to_bits
    validate_bits and_bits or_bits xor_bits min_bit_id canonical_bits);
use SplitAligner::Newick qw(parse_newick assign_branch_ids retained_primitive_axis
    tip_labels split_taxon_sets canonical_split edge_records gene_split_records serialize_plain);
use SplitAligner::CompactKernel qw(prepare_reference compute_gene compute_gene_materialized materialize_gene construct_cell);

sub rejected {
    my ($code, $pattern, $name) = @_;
    my $ok = eval { $code->(); 1 };
    ok(!$ok && $@ =~ $pattern, $name) or diag($@ || 'unexpected success');
}

for my $n (2,3,7,8,9,63,64,65,127,128,129,302) {
    my $full = ids_to_bits($n, [0 .. $n - 1]);
    my $a = ids_to_bits($n, [0]);
    my $b = xor_bits($a, $full);
    is(length($full), int(($n + 7)/8), "fixed width $n");
    is_deeply(bits_to_ids(hex_to_bits(bits_to_hex($a,$n),$n),$n),[0],"hex roundtrip $n");
    is(canonical_bits($a,$full,0),canonical_bits($b,$full,0),"complement canonical identity $n");
    is(and_bits($a,$b),zero_bits($n),"disjoint sides $n");
    is(min_bit_id(ids_to_bits($n,[$n-1]),$n),$n-1,"last bit anchor $n");
    rejected(sub { ids_to_bits($n,[$n]) },qr/Invalid or duplicate/,'reject out-of-range bit');
}
my ($numeric_a,$numeric_b)=(pack('C*',49,50),pack('C*',51,52));
my $coercion=0+$numeric_a;
is(unpack('H*',xor_bits($numeric_a,$numeric_b)),'0206','numeric scalar slot cannot change byte XOR');
my $upgraded=zero_bits(8); utf8::upgrade($upgraded);
rejected(sub { validate_bits($upgraded,8) },qr/byte string/,'reject UTF8 flagged bit vector');
rejected(sub { hex_to_bits('ff',7) },qr/tail bits/,'reject nonzero tail');
rejected(sub { hex_to_bits('0g',7) },qr/hexadecimal/,'reject invalid hex');
rejected(sub { hex_to_bits('00ff',8) },qr/hexadecimal/,'reject wrong width');
rejected(sub { and_bits("a","bb") },qr/equal-length/,'reject unequal operands');
rejected(sub { canonical_bits("\x04","\x03",0) },qr/not a subset/,'reject foreign side bits');

# Optional frozen-code oracle: read its module into a separate package in this
# process. The original file is untouched. This is a bounded fixture test, not
# an installation, production run, or full-panel benchmark.
my $oracle='SplitAligner::Newick';
if (my $baseline=$ENV{SPLITALIGNER_BASELINE_DIR}) {
    my $path=File::Spec->catfile($baseline,'lib','SplitAligner','Newick.pm');
    open my $fh,'<:raw',$path or die "Cannot read frozen oracle $path: $!";
    local $/; my $code=<$fh>; close $fh;
    $code =~ s/\Apackage SplitAligner::Newick;/package FrozenNewick;/
        or die 'Unexpected frozen module package';
    eval $code; die $@ if $@;
    $oracle='FrozenNewick';
    pass('loaded frozen Newick module as independent oracle');
} else {
    note('SPLITALIGNER_BASELINE_DIR unset: using unchanged legacy algorithm helpers in candidate');
}
sub oracle_call {
    my ($name,@args)=@_;
    no strict 'refs';
    return &{"${oracle}::$name"}(@args);
}

sub compare_case {
    my ($reference,$newicks,$label)=@_;
    my $root=parse_newick($reference);
    assign_branch_ids($root);
    my $axis=retained_primitive_axis($root);
    my $oldroot=oracle_call('parse_newick',$reference);
    oracle_call('assign_branch_ids',$oldroot);
    is_deeply($axis,oracle_call('retained_primitive_axis',$oldroot),"$label reference identities unchanged");
    is(serialize_plain($root),oracle_call('serialize_plain',$oldroot),"$label serialization order unchanged");
    my $ctx=prepare_reference($axis,[tip_labels($root)]);
    for my $tree (@{$newicks}) {
        my $gene=parse_newick($tree);
        my $oldgene=oracle_call('parse_newick',$tree);
        my $record=compute_gene($gene,$ctx);
        my ($paired_record, $paired_data)=compute_gene_materialized($gene,$ctx);
        is_deeply($paired_record,$record,"$label paired API preserves source record");
        is_deeply($paired_data,materialize_gene($record,$ctx),"$label paired API reuses equivalent materialization");
        my @single_record=compute_gene($gene,$ctx);
        is(scalar(@single_record),1,"$label original API still returns one record in list context");
        my $raw=JSON::PP->new->canonical(1)->encode($record);
        my $decoded=JSON::PP->new->decode($raw);
        my $m=materialize_gene($decoded,$ctx);
        is_deeply($record,$decoded,"$label ledger JSON roundtrip");
        my $old=oracle_call('gene_split_records',$oldgene);
        my %old=map { $_->{split} => $_ } @{$old};
        my $tg=hex_to_bits($record->{tg},$ctx->{n});
        my %present=map { $_=>1 } oracle_call('tip_labels',$oldgene);
        my %groups;
        for my $coordinate (@{$axis}) {
            my ($a,$b)=split_taxon_sets($coordinate->{split});
            my @a=grep { $present{$_} } @{$a};
            my @b=grep { $present{$_} } @{$b};
            next unless @a && @b;
            push @{$groups{canonical_split(\@a,\@b)}},$coordinate->{branch_id};
        }
        my (@truth,@state,@primitive,@expected_groups);
        my %expected;
        for my $split (sort keys %groups) {
            my @members=sort { substr($a,1)<=>substr($b,1) } @{$groups{$split}};
            my $obs=$old{$split};
            $expected{$_}=[$#members>0?'NA_fuse':$obs?'Mapped':'not_displayed',
                $obs?($#members>0?'F':'D'):'U',@members==1&&$obs?$obs->{value}:'NA'] for @members;
            push @expected_groups,{split=>$split,members=>join('|',@members),
                status=>$obs?'recovered':'not_displayed',length=>$obs?$obs->{value}:'NA'} if @members>1;
        }
        for my $coordinate (@{$axis}) {
            my $row=$expected{$coordinate->{branch_id}} || ['NA_struct','S','NA'];
            push @truth,$row->[0]; push @state,$row->[1]; push @primitive,$row->[2];
        }
        is_deeply($m->{truth},\@truth,"$label reference truth");
        is_deeply($m->{state},\@state,"$label compatibility states");
        is_deeply($m->{primitive},\@primitive,"$label primitive numeric tokens");
        is_deeply($m->{groups},\@expected_groups,"$label reference groups incl absent observations");
        is(scalar(@{$record->{observations}}),scalar(@{$old}),"$label all input observations retained");
        for my $obs (@{$record->{observations}}) {
            my $side=hex_to_bits($obs->{key},$ctx->{n});
            my $split=canonical_split([map { $ctx->{taxa}[$_] } @{bits_to_ids($side,$ctx->{n})}],
                [map { $ctx->{taxa}[$_] } @{bits_to_ids(xor_bits($side,$tg),$ctx->{n})}]);
            ok(exists($old{$split}),"$label observation identity");
            is($obs->{length},$old{$split}{value},"$label observation raw or pooled length");
            is_deeply($obs->{source_nodes},$old{$split}{source_nodes},"$label preorder sources");
        }
        my $bad=JSON::PP->new->decode($raw); $bad->{truth}[0]='broken';
        rejected(sub { materialize_gene($bad,$ctx) },qr/Stored truth differs/,'tampered truth rejected');
        $bad=JSON::PP->new->decode($raw); $bad->{components}[0]{members}='B999999';
        rejected(sub { materialize_gene($bad,$ctx) },qr/Stored reference components differ/,'tampered component rejected');
        $bad=JSON::PP->new->decode($raw); $bad->{observations}[0]{source_nodes}[0]=0;
        rejected(sub { materialize_gene($bad,$ctx) },qr/source node/,'invalid source reference rejected');
        $bad=JSON::PP->new->decode($raw); $bad->{reference_id}='wrong';
        rejected(sub { materialize_gene($bad,$ctx) },qr/reference identity/,'foreign dictionary/axis rejected');
        $bad=JSON::PP->new->decode($raw); $bad->{observations}[0]{length}='1e-999';
        rejected(sub { materialize_gene($bad,$ctx) },qr/binary64 range/,'underflowing saved value rejected');
    }
    my @incomplete=grep { $_->{branch_id} ne 'B1' } @{$axis};
    my $axis_error=@incomplete ? qr/missing terminal/ : qr/nonempty axis/;
    rejected(sub { prepare_reference(\@incomplete,[tip_labels($root)]) },$axis_error,'incomplete reference axis rejected');
    return $ctx;
}

my $ctx=compare_case('(((A,B),(C,D)),((E,F),(G,H)));',[
    '(((A:1.00,B:0),(C:2,D:3):4):5,((E:6,F:7):8,(G:9,H:10):11):12);',
    '(A:1,B:2,(G:3,H:4):5);',
    '((A:1,G:2),(B:3,H:4):0);',
    '((A:0,B:0):0,(E:0,H:0):0);',
    '(A,B);',
    '((A:1e308,B:3):1.5e308,(E:4,H:5):0);',
],'eight-tip');
compare_case('((Alpha..x,B|eta),(猫,Δ));',[
    '(Alpha..x:1,猫:2,Δ:3);',
    '((Alpha..x:1.00,B|eta:0):2,(猫:3,Δ:4):5);',
],'structured-Unicode');
compare_case('(A,B);',['(A:2,B:5);'],'two-tip');
rejected(sub { compute_gene(parse_newick('((A:1,B:2):1.2e308,(E:4,H:5):1.2e308);'),$ctx) },
    qr/Pooled unrooted branch length is not finite/,'input observation overflow keeps legacy error');
rejected(sub { compute_gene(parse_newick('(A:1,Z:2);'),$ctx) },qr/Unknown taxon 'Z'/,'unknown taxon named');
my $duplicate_tree=parse_newick('(A:1,B:2);');
$duplicate_tree->{children}[1]{label}='A';
rejected(sub { compute_gene($duplicate_tree,$ctx) },qr/Duplicate terminal taxon 'A'/,'duplicate taxon named');

my $b1={id=>'B1',size=>1,recovered=>1,length=>'1'};
my $b23={id=>'B2|B3',size=>2,recovered=>1,length=>'2'};
my $owner={B1=>$b1,B2=>$b23,B3=>$b23};
is_deeply([construct_cell([qw(B1 B2 B3)],$owner,'whole')],
    ['decomposable','recovered','available','selected_component_edge_sum','B1,B2|B3','','','3'],
    'whole uses complete smaller fusion once');
is((construct_cell([qw(B1 B2 B3)],$owner,'direct'))[2],'direct_only','direct suppresses decomposable sum');
is((construct_cell([qw(B1 B2)],$owner,'whole'))[0],'unavailable','cannot cut a larger component');
is((construct_cell([qw(B1 B4)],$owner,'whole'))[0],'unavailable','cannot ignore deleted member');
$b1->{length}='NA';
is((construct_cell([qw(B1 B2 B3)],$owner,'whole'))[2],'missing_length','no partial sum with missing length');
is((construct_cell([qw(B1 B2 B3)],$owner,'direct'))[2],'direct_only','direct takes priority over missing length');
$b1->{length}='1.5e308';$b23->{length}='1.5e308';
is((construct_cell([qw(B1 B2 B3)],$owner,'whole'))[2],'nonfinite_sum','whole overflow classified');
is((construct_cell([qw(B1 B2 B3)],$owner,'direct'))[2],'direct_only','direct takes priority over overflow');
$b23->{recovered}=0;
is((construct_cell([qw(B1 B2 B3)],$owner,'direct'))[1],'not_recovered','recovery precedes direct suppression');
$b23->{recovered}=1;$b1->{length}='0';$b23->{length}='0';
is((construct_cell([qw(B1 B2 B3)],$owner,'whole'))[-1],'0','valid zero retained');

my ($weak_root,$weak_leaf);
{
    my $tree=parse_newick('((A:1,B:2):3,C:4);');
    $weak_root=$tree;weaken($weak_root);
    $weak_leaf=$tree->{children}[0]{children}[0];weaken($weak_leaf);
    assign_branch_ids($tree); my $edges=edge_records($tree); my $plain=serialize_plain($tree);
}
ok(!defined($weak_root)&&!defined($weak_leaf),'legacy tree traversals release root and descendants');
done_testing();
