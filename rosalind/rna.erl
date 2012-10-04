-module(rna).

-export([translate_codon_to_protein/1]).

translate_codon_to_protein(S) when is_list(S) ->
    translate_codon_to_protein(list_to_binary(S));
translate_codon_to_protein(B) when is_binary(B) ->
    Trans = translate(B, <<>>),
    binary_to_list(Trans).

translate(<<>>, Trans) ->
    Trans;
translate(<<"UUU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "F">>);
translate(<<"CUU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "L">>);
translate(<<"AUU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "I">>);
translate(<<"GUU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "V">>);
translate(<<"UUC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "F">>);
translate(<<"CUC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "L">>);
translate(<<"AUC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "I">>);
translate(<<"GUC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "V">>);
translate(<<"UUA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "L">>);
translate(<<"CUA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "L">>);
translate(<<"AUA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "I">>);
translate(<<"GUA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "V">>);
translate(<<"UUG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "L">>);
translate(<<"CUG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "L">>);
translate(<<"AUG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "M">>);
translate(<<"GUG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "V">>);
translate(<<"UCU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "S">>);
translate(<<"CCU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "P">>);
translate(<<"ACU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "T">>);
translate(<<"GCU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "A">>);
translate(<<"UCC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "S">>);
translate(<<"CCC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "P">>);
translate(<<"ACC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "T">>);
translate(<<"GCC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "A">>);
translate(<<"UCA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "S">>);
translate(<<"CCA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "P">>);
translate(<<"ACA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "T">>);
translate(<<"GCA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "A">>);
translate(<<"UCG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "S">>);
translate(<<"CCG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "P">>);
translate(<<"ACG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "T">>);
translate(<<"GCG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "A">>);
translate(<<"UAU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "Y">>);
translate(<<"CAU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "H">>);
translate(<<"AAU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "N">>);
translate(<<"GAU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "D">>);
translate(<<"UAC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "Y">>);
translate(<<"CAC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "H">>);
translate(<<"AAC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "N">>);
translate(<<"GAC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "D">>);
translate(<<"UAA", _Rest/binary>>, <<Trans/binary>>) ->
    Trans; % stop
translate(<<"CAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "Q">>);
translate(<<"AAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "K">>);
translate(<<"GAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "E">>);
translate(<<"UAG", _Rest/binary>>, <<Trans/binary>>) ->
    Trans; % stop
translate(<<"CAG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "Q">>);
translate(<<"AAG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "K">>);
translate(<<"GAG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "E">>);
translate(<<"UGU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "C">>);
translate(<<"CGU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "R">>);
translate(<<"AGU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "S">>);
translate(<<"GGU", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "G">>);
translate(<<"UGC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "C">>);
translate(<<"CGC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "R">>);
translate(<<"AGC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "S">>);
translate(<<"GGC", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "G">>);
translate(<<"UGA", _Rest/binary>>, <<Trans/binary>>) ->
    Trans; % stop
translate(<<"CGA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "R">>);
translate(<<"AGA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "R">>);
translate(<<"GGA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "G">>);
translate(<<"UGG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "W">>);
translate(<<"CGG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "R">>);
translate(<<"AGG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "R">>);
translate(<<"GGG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "G">>).
 
