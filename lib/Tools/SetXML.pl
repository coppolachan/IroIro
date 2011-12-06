#!/usr/bin/perl -w
use warnings;
use strict;
use XML::Twig;

# usage
# SetXML <input_XML> [<xml path>] [<new value>] [<-o out_XML>]
# modify XML on browser?

my $filename  = $ARGV[0];
my $xml_path  = $ARGV[1];
my $new_value = $ARGV[2];

my $twig = XML::Twig->new(
    twig_handlers => 
    {$xml_path => sub {+$_->set_text($new_value)} },
    pretty_print => 'indented',
    comments => 'keep');
$twig->parsefile($filename);
$twig->print_to_file('out.xml',pretty_print => 'indented');

