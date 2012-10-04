-module(user_default).

-export([
    count_nuc/1, rna_trans/1, rev_comp/1, count_point_mut/2, pro_trans/1
]).

% Protein Translation
% http://rosalind.info/problems/prot/
pro_trans(S) ->
    rna:translate_codon_to_protein(S).

% Counting Point Mutations
% http://rosalind.info/problems/hamm/
% Given two strings of equal length, count the number of
% positions with a difference.  This is the Hamming distance.
count_point_mut(S1, S2) when is_list(S1), is_list(S2) ->
    count_point_mut(list_to_binary(S1), list_to_binary(S2));
count_point_mut(B1, B2) when is_binary(B1), is_binary(B2) ->
    count_point_mut(B1, B2, 0).

count_point_mut(<<>>, <<>>, Count) ->
    Count;
count_point_mut(<<N1:8, Rest1/binary>>, <<N2:8, Rest2/binary>>, Count) when N1 =:= N2 ->
    count_point_mut(Rest1, Rest2, Count);
count_point_mut(<<_N1:8, Rest1/binary>>, <<_N2:8, Rest2/binary>>, Count) ->
    count_point_mut(Rest1, Rest2, Count+1).

% Reverse Complement
% http://rosalind.info/problems/revc/
% Reverse the string and complement each nucleotide:
% A -> T, C -> G
rev_comp(S) when is_list(S) ->
    rev_comp(list_to_binary(S));
rev_comp(B) when is_binary(B) ->
    RC = rev_comp(B, <<>>),
    binary_to_list(RC).

rev_comp(<<>>, RC) ->
    RC;
rev_comp(<<"A", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"T", RC/binary>>);
rev_comp(<<"T", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"A", RC/binary>>);
rev_comp(<<"C", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"G", RC/binary>>);
rev_comp(<<"G", Rest/binary>>, <<RC/binary>>) ->
    rev_comp(Rest, <<"C", RC/binary>>).

%complement(N) when is_list(N) ->
%    complement(list_to_binary(N));
%complement(<<"A">>) ->
%    <<"T">>;
%complement(<<"T">>) ->
%    <<"A">>;
%complement(<<"C">>) ->
%    <<"G">>;
%complement(<<"G">>) ->
%    <<"C">>;
%complement(_) ->
%    {error, invalid}.

% RNA Transcription
% http://rosalind.info/problems/rna/
% Transcribe T to U.
rna_trans(S) when is_list(S) ->
    rna_trans(list_to_binary(S));
rna_trans(B) when is_binary(B) ->
    Trans = rna_trans(B, <<>>),
    binary_to_list(Trans).

rna_trans(<<>>, Trans) ->
    Trans;
rna_trans(<<"T", Rest/binary>>, <<Trans/binary>>) ->
    rna_trans(Rest, <<Trans/binary, "U">>);
rna_trans(<<N:8, Rest/binary>>, <<Trans/binary>>) ->
    rna_trans(Rest, <<Trans/binary, N>>).

% Counting Nucleotides
% http://rosalind.info/problems/dna/
% Count occurrences of each nucleotide.
count_nuc(S) when is_list(S) ->
    count_nuc(list_to_binary(S));
count_nuc(B) when is_binary(B) ->
    {A, C, G, T} = count_nuc(B, {0, 0, 0, 0}),
    io:format("~p ~p ~p ~p~n", [A, C, G, T]).

count_nuc(<<>>, Acc) ->
    Acc;
count_nuc(<<"A", Rest/binary>>, {A, C, G, T}) ->
    count_nuc(Rest, {A+1, C, G, T});
count_nuc(<<"C", Rest/binary>>, {A, C, G, T}) ->
    count_nuc(Rest, {A, C+1, G, T});
count_nuc(<<"G", Rest/binary>>, {A, C, G, T}) ->
    count_nuc(Rest, {A, C, G+1, T});
count_nuc(<<"T", Rest/binary>>, {A, C, G, T}) ->
    count_nuc(Rest, {A, C, G, T+1}).

