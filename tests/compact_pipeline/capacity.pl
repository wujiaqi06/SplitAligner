#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use utf8;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use Test::More;
use File::Temp qw(tempdir);
use File::Spec;
use Cwd qw(getcwd);
use Digest::SHA qw(sha256_hex);
use JSON::PP;
use SplitAligner::TextIO qw(encode_utf8_strict);
use SplitAligner::Provenance qw(write_axis_ledger file_sha256);
use SplitAligner::Newick qw(parse_newick assign_branch_ids retained_primitive_axis);
use SplitAligner::CompactPipeline qw(run_compact_matrix visit_input_lines read_ledger_block compact_limits);

my $tmp=tempdir('compact-capacity-XXXXXX',TMPDIR=>1,CLEANUP=>1);
my $start=getcwd();
my $json=JSON::PP->new->canonical(1)->utf8(1);
sub raw_file {
    my ($name,$raw)=@_;
    my $path=File::Spec->catfile($tmp,$name);
    open my $fh,'>:raw',$path or die $!;
    print {$fh} $raw or die $!;
    close $fh or die $!;
    return $path;
}
sub rejected {
    my ($code,$pattern,$message)=@_;
    my $ok=eval { $code->(); 1 };
    ok(!$ok && $@ =~ $pattern,$message) or diag($@ || 'unexpected success');
}

my $reference='((A,B),((C,D),(E,F)));';
my $genes="g1(A:1,C:2,E:3);\ng2((A:1,B:2):3,(C:4,D:5,E:6,F:7):8);\n";
my $species_path=raw_file('reference.nwk',$reference);
my $gene_path=raw_file('genes.nwk',$genes);
my $root=parse_newick($reference); assign_branch_ids($root);
my $axis=retained_primitive_axis($root);
my $axis_path=File::Spec->catfile($tmp,'axis.tsv'); write_axis_ledger($axis_path,$axis);
my %input_hash=map { $_=>file_sha256($_) } ($species_path,$gene_path,$axis_path);
my $public=File::Spec->catdir($tmp,'public'); mkdir $public or die $!;
my $sentinel=File::Spec->catfile($public,'existing-result.txt');
open my $old,'>:raw',$sentinel or die $!; print {$old} "keep old result\n";close $old or die $!;
my $old_hash=file_sha256($sentinel);

my %common=(label=>'panel',axis_file=>$axis_path,species_root=>$root,gene_path=>$gene_path,
    view=>'whole',parser_sha256=>file_sha256("$RealBin/../../lib/SplitAligner/Newick.pm"),
    input_sha256=>{species_tree=>file_sha256($species_path),gene_trees=>file_sha256($gene_path)},
    input_bytes=>{species_tree=>length($reference),gene_trees=>length($genes)});

# A capacity exception must prevent the matrix stage from returning a result
# eligible for the controller's publication step. Partial files are permitted
# only inside this private staging directory; transaction publication itself is
# covered by the CLI/ownership suites, not simulated by this internal API test.
for my $kind (sort keys %{compact_limits()}) {
    my $work=File::Spec->catdir($tmp,"stage-$kind");mkdir $work or die $!;
    chdir $work or die $!;
    my $result;
    my $ok=eval { $result=run_compact_matrix(%common,limits=>{$kind=>1});1 };
    my $error=$@;
    chdir $start or die $!;
    ok(!$ok && $error =~ /\[Capacity\].*\Q$kind\E/,"$kind low limit rejects") or diag($error||'unexpected success');
    ok(!defined($result),"$kind returns no publishable result");
    is(file_sha256($sentinel),$old_hash,"$kind leaves outside result unchanged");
}
for my $path (sort keys %input_hash) {
    is(file_sha256($path),$input_hash{$path},'all capacity failures preserve source input');
}
opendir my $pub,$public or die $!; my @public_files=sort grep { $_ ne '.'&&$_ ne '..' } readdir $pub;closedir $pub;
is_deeply(\@public_files,['existing-result.txt'],'private capacity failures create no public artifacts');

my $bounded=raw_file('long-no-newline', 'x' x 33);
my $visits=0;
rejected(sub { visit_input_lines($bounded,32,sub {++$visits}) },qr/Capacity.*input_line_bytes/,
    'unterminated line is rejected before callback');
is($visits,0,'no callback for oversized unterminated line');
my $boundary=raw_file('exact-no-newline','abcde');
my @got;
is(visit_input_lines($boundary,5,sub {push @got,[@_]}),1,'exact byte limit accepted without trailing newline');
is_deeply(\@got,[['abcde',1,0,5]],'raw offsets and lengths include original bytes');
$boundary=raw_file('newline-over-limit',"abcde\n");
rejected(sub {visit_input_lines($boundary,5,sub {})},qr/Capacity.*input_line_bytes/,
    'newline counts towards original byte limit');
my $unicode=encode_utf8_strict("甲\r\n\n乙",'test');
my $unicode_path=raw_file('unicode-lines',$unicode);
@got=();
is(visit_input_lines($unicode_path,32,sub {push @got,[@_]}),3,'CRLF blank line and EOF are visited');
is_deeply(\@got,[["甲",1,0,5],['',2,5,1],["乙",3,6,3]],'Unicode offsets are bytes rather than characters');
my $chunked=raw_file('utf8-crossing-read',('a' x 65535).encode_utf8_strict("猫\n",'test'));
@got=();
visit_input_lines($chunked,70000,sub {push @got,[$_[1],$_[2],$_[3],substr($_[0],-1)]});
is_deeply(\@got,[[1,0,65539,"猫"]],'UTF8 character crossing raw-read boundary remains intact');
my $bom=raw_file('bom',"\xEF\xBB\xBFg1(A,B);\n");
rejected(sub {visit_input_lines($bom,100,sub {})},qr/BOM is not permitted/,'BOM rejected');
my $invalid=raw_file('bad-utf8',"\xff\n");
rejected(sub {visit_input_lines($invalid,100,sub {})},qr/Invalid UTF-8/,'invalid UTF8 rejected');

my $record={gene_id=>'g1',input_record_id=>1,input_line=>2,observations=>[],components=>[]};
my $raw=$json->encode($record)."\n";
sub block_entry {
    my ($bytes,$offset)=@_;
    return {gene_id=>'g1',input_record_id=>1,input_line=>2,observations=>0,components=>0,
        ledger_offset=>$offset//0,ledger_bytes=>length($bytes),sha256=>sha256_hex($bytes)};
}
sub read_case {
    my ($name,$content,$entry,$maximum)=@_;
    my $path=raw_file($name,$content);
    open my $fh,'<:raw',$path or die $!;
    my $result=read_ledger_block($fh,$entry,$maximum//4096);
    close $fh or die $!;
    return $result;
}
is_deeply(read_case('valid-ledger','PREFIX'.$raw.'SUFFIX',block_entry($raw,6)),$record,
    'bounded ledger block read uses exact offset and length');
my $entry=block_entry($raw);$entry->{sha256}='0'x64;
rejected(sub {read_case('wrong-hash',$raw,$entry)},qr/fingerprint differs/,'ledger block hash mismatch rejected');
$entry=block_entry($raw);$entry->{ledger_bytes}++;
rejected(sub {read_case('short-ledger',$raw,$entry)},qr/Short or failed read/,'short ledger block rejected');
for my $bad_offset (-1,'oops',undef) {
    my $bad=block_entry($raw);$bad->{ledger_offset}=$bad_offset;
    rejected(sub {read_case('bad-offset',$raw,$bad)},qr/block boundary/,'invalid ledger offset rejected');
}
for my $bad_length (0,-1,'1.5',undef,5000) {
    my $bad=block_entry($raw);$bad->{ledger_bytes}=$bad_length;
    rejected(sub {read_case('bad-length',$raw,$bad)},qr/block boundary/,'invalid or over-budget ledger length rejected');
}
my $no_newline=$json->encode($record);
rejected(sub {read_case('no-terminator',$no_newline,block_entry($no_newline))},qr/newline-terminated/,
    'ledger requires final newline');
my $multi=$raw.$raw;
rejected(sub {read_case('multiple-records',$multi,block_entry($multi))},qr/newline-terminated/,
    'ledger block cannot contain two records');
my $bad_json="not-json\n";
rejected(sub {read_case('invalid-json',$bad_json,block_entry($bad_json))},qr/Invalid compact ledger JSON/,
    'invalid ledger JSON rejected after fingerprint check');
for my $field (qw(gene_id input_record_id input_line observations components)) {
    my $bad=block_entry($raw);$bad->{$field}= $field eq 'gene_id' ? 'other' : 9;
    rejected(sub {read_case('wrong-identity',$raw,$bad)},qr/identity or record counts differ/,
        "$field index-record mismatch rejected");
}
done_testing();
