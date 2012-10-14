-module(rna).

-compile(export_all).

-include("rna.hrl").

% Inferring mRNA from Protein
% http://rosalind.info/problems/mrna/
% Given a string of amino acid symbols, return the number of
% possible RNA strings that could have generated it.
count_codon_combinations(S) when is_list(S) ->
    Bin = list_to_binary(S),
    count_codon_combinations(Bin);
count_codon_combinations(B) when is_binary(B) ->
    count_codon_combinations(B, 1);
count_codon_combinations(_) ->
    {error, invalid_argument}.

count_codon_combinations(<<>>, C) ->
    StopC = symbol_lookup(<<"Stop">>),
    C * length(StopC);
count_codon_combinations(<<S:8, T/binary>>, C) ->
    Codons = symbol_lookup(S),
    C2 = C * length(Codons),
    count_codon_combinations(T, C2).

symbol_lookup(<<"Stop">>) ->
    [<<"UAA">>,<<"UAG">>,<<"UGA">>];
symbol_lookup($F) ->
    [<<"UUU">>,<<"UUC">>];
symbol_lookup($L) ->
    [<<"CUU">>,<<"CUC">>,<<"UUA">>,<<"CUA">>,<<"UUG">>,<<"CUG">>];
symbol_lookup($I) ->
    [<<"AUU">>,<<"AUC">>,<<"AUA">>];
symbol_lookup($V) ->
    [<<"GUU">>,<<"GUC">>,<<"GUA">>,<<"GUG">>];
symbol_lookup($M) ->
    [<<"AUG">>];
symbol_lookup($S) ->
    [<<"UCU">>,<<"UCC">>,<<"UCA">>,<<"UCG">>,<<"AGU">>,<<"AGC">>];
symbol_lookup($P) ->
    [<<"CCU">>,<<"CCC">>,<<"CCA">>,<<"CCG">>];
symbol_lookup($T) ->
    [<<"ACU">>,<<"ACC">>,<<"ACA">>,<<"ACG">>];
symbol_lookup($A) ->
    [<<"GCU">>,<<"GCC">>,<<"GCA">>,<<"GCG">>];
symbol_lookup($Y) ->
    [<<"UAU">>,<<"UAC">>];
symbol_lookup($H) ->
    [<<"CAU">>,<<"CAC">>];
symbol_lookup($N) ->
    [<<"AAU">>,<<"AAC">>];
symbol_lookup($D) ->
    [<<"GAU">>,<<"GAC">>];
symbol_lookup($Q) ->
    [<<"CAA">>,<<"CAG">>];
symbol_lookup($K) ->
    [<<"AAA">>,<<"AAG">>];
symbol_lookup($E) ->
    [<<"GAA">>,<<"GAG">>];
symbol_lookup($C) ->
    [<<"UGU">>,<<"UGC">>];
symbol_lookup($R) ->
    [<<"CGU">>,<<"CGC">>,<<"CGA">>,<<"AGA">>,<<"CGG">>,<<"AGG">>];
symbol_lookup($G) ->
    [<<"GGU">>,<<"GGC">>,<<"GGA">>,<<"GGG">>];
symbol_lookup($W) ->
    [<<"UGG">>].

translate_codon_to_protein(S) when is_list(S) ->
    translate_codon_to_protein(list_to_binary(S));
translate_codon_to_protein(B) when is_binary(B) ->
    translate(B, <<>>).

% Translate stop to a dash which helps with some problems.
translate(<<>>, Trans) ->
    Trans;
translate(<<"UAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, ?StopAA>>); % stop
translate(<<"UAG", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, ?StopAA>>); % stop
translate(<<"UGA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, ?StopAA>>); % stop
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
translate(<<"CAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "Q">>);
translate(<<"AAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "K">>);
translate(<<"GAA", Rest/binary>>, <<Trans/binary>>) ->
    translate(Rest, <<Trans/binary, "E">>);
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
    translate(Rest, <<Trans/binary, "G">>);
translate(_Codon, <<Trans/binary>>) ->
    Trans. % unknown codon - we're done

% RNA Transcription
% http://rosalind.info/problems/rna/
% Transcribe T to U.
transcribe_rna(S) when is_list(S) ->
    transcribe_rna(list_to_binary(S));
transcribe_rna(B) when is_binary(B) ->
    Trans = transcribe_rna(B, <<>>),
    binary_to_list(Trans).

transcribe_rna(<<>>, Trans) ->
    Trans;
transcribe_rna(<<"T", Rest/binary>>, <<Trans/binary>>) ->
    transcribe_rna(Rest, <<Trans/binary, "U">>);
transcribe_rna(<<N:8, Rest/binary>>, <<Trans/binary>>) ->
    transcribe_rna(Rest, <<Trans/binary, N>>).

% RNA reverse complement.
rev_comp(S) when is_list(S) ->
    rev_comp(list_to_binary(S));
rev_comp(B) when is_binary(B) ->
    RC = rev_comp(B, <<>>),
    binary_to_list(RC).

rev_comp(<<>>, RC) ->
    RC;
rev_comp(<<"A", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"U", RC/binary>>);
rev_comp(<<"U", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"A", RC/binary>>);
rev_comp(<<"C", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"G", RC/binary>>);
rev_comp(<<"G", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"C", RC/binary>>).

