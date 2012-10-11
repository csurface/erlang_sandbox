-module(user_default).

-compile(export_all).

-define(FORMAT_LIST(L), [io:format("~p ", [P]) || P <- L], io:format("~n")).

% Finding a Shared Motif
% http://rosalind.info/problems/lcs/
% Find one longest common substring of all given strings.
% MinLength gives a performance improvement.
find_shared_motif(Filename, MinLength) ->
    {ok, Handle} = file:open(Filename, read),
    Strings = read_strings(Handle),
    file:close(Handle),
    find_shared_motif_strings(Strings, MinLength).

read_strings(Handle) ->
    lists:reverse(read_strings(Handle, [])).

read_strings(Handle, Acc) ->
    case file:read_line(Handle) of
        {ok, String} ->
            case string:strip(String) of
                S when length(S) > 0 ->
                    read_strings(Handle, [string:strip(S, right, $\n)|Acc]);
                _ ->
                    Acc
            end;
        eof ->
            Acc;
        Error ->
            Error
    end.

find_shared_motif_strings(Strings, MinLength) ->
    % Sort ascending so that we select motifs from the shortest string.
    SortAsc = fun(A, B) -> length(A) =< length(B) end,
    Sorted = lists:sort(SortAsc, Strings),
    % Reverse the motifs so that we are searching for longest ones first.
    Motifs = lists:reverse(all_substrings(lists:nth(1, Sorted), MinLength)),
    io:format("Found ~p motifs~n", [length(Motifs)]),
    Remainder = lists:nthtail(1, Sorted),
    find_in_all(Motifs, Remainder).

% Return a subset of strings in Motifs that exist in all Strings.
find_in_all([], _) ->
    not_found;
find_in_all([Motif|T], Strings) ->
    case is_motif_in_strings(Motif, Strings) of
        true ->
            Motif;
        _ ->
            find_in_all(T, Strings)
    end.

is_motif_in_strings([], []) ->
    false;
is_motif_in_strings(_Motif, []) ->
    true;
is_motif_in_strings(Motif, [S|Tail]) ->
    case find_motif(Motif, S) of
        [] ->
            false;
        L when is_list(L), length(L) > 0 ->
            is_motif_in_strings(Motif, Tail)
    end.

all_substrings(S, MinLength) ->
    all_substrings(S, length(S), MinLength, []).

all_substrings(_S, 0, _MinLength, Acc) ->
    Acc;
all_substrings(_S, Length, MinLength, Acc) when Length < MinLength ->
    Acc;
all_substrings(S, Length, MinLength, Acc) when length(S) =:= Length ->
    all_substrings(S, Length-1, MinLength, [S|Acc]);
all_substrings(S, Length, MinLength, Acc) ->
    Acc2 = find_all_substrings(S, Length, lists:seq(1, length(S) - Length + 1), Acc),
    all_substrings(S, Length-1, MinLength, Acc2).

find_all_substrings(_S, _Length, [], Acc) ->
    Acc;
find_all_substrings(S, Length, [Position|Tail], Acc) ->
    Sub = lists:sublist(S, Position, Length),
    find_all_substrings(S, Length, Tail, [Sub|Acc]).

% Overlap Graphs
% http://rosalind.info/problems/grph/
overlap_graphs(Filename, K) ->
    Fasta = read_fasta_file(Filename),
    AL = adjacency_list(Fasta, K, length(Fasta), []),
    [io:format("~s ~s~n", [N1, N2]) || {N1, N2} <- AL].

adjacency_list(_, _, 0, Acc) ->
    Acc;
adjacency_list(Fasta, K, N, Acc) ->
    {Name, String} = lists:nth(N, Fasta),
    Remainder = remainder(Fasta, N),
    Acc2 = is_adjacent({Name, String}, K, Remainder, []),
    adjacency_list(Fasta, K, N-1, lists:flatten([Acc2|Acc])).

is_adjacent(_, _, [], Acc) ->
    lists:reverse(Acc);
is_adjacent({Name1, String1}, K, [{Name2, String2}|Tail], Acc) ->
    N = length(String1),
    Suffix = lists:nthtail(N-K, String1),
    Acc2 = case lists:prefix(Suffix, String2) of
        true ->
            [{Name1, Name2}|Acc];
        _ ->
            Acc
    end,
    is_adjacent({Name1, String1}, K, Tail, Acc2).

remainder([_|T], 1) ->
    T;
remainder(L, N) ->
    A = lists:sublist(L, N-1),
    B = lists:nthtail(N, L),
    lists:flatten([A|B]).

read_fasta_file(Filename) ->
    {ok, Handle} = file:open(Filename, read),
    Fasta = read_fasta(Handle),
    file:close(Handle),
    Fasta.

% GC Content
% http://rosalind.info/problems/gc/
% Read a file in FASTA format and return the GC content of the highest
% string:
%
% >Name1
% ACGT...
% ACGT...
% >Name2
% TCGA...
%
% Return:
% NameM
% XX.YY%
gc_content(Filename) ->
    {ok, Handle} = file:open(Filename, read),
    Fasta = read_fasta(Handle),
    file:close(Handle),
    {Name, GC} = compute_max_gc(Fasta, {none, 0}),
    io:format("~s~n~p%~n", [Name, GC]).

compute_max_gc([], {Name, Value}) ->
    {Name, erlang:round(Value * 10000) / 100};
compute_max_gc([{Name, DnaString}|Tail], {MaxName, MaxValue}) ->
    GcValue = compute_gc(DnaString),
    NextMax = case (GcValue > MaxValue) of
        true ->
            {Name, GcValue};
        _ ->
            {MaxName, MaxValue}
    end,
    compute_max_gc(Tail, NextMax).

compute_gc(S) ->
    {_A, C, G, _T} = count_nuc(list_to_binary(S), {0, 0, 0, 0}),
    L = length(S),
    (G + C) / L.

% Return tuples of {Name, DNA_String}.
read_fasta(Handle) ->
    lists:reverse(read_fasta(Handle, {}, [])).

read_fasta(Handle, {}, Acc) ->
    case file:read_line(Handle) of
        {ok, [S|Tail]} ->
            case S of
                $> ->
                    % This line contains a name.
                    read_fasta(Handle, {string:strip(Tail, right, $\n), []}, Acc);
                _ ->
                    {error, unexpected_value}
            end;
        _ ->
            {error, invalid_format1}
    end;
read_fasta(Handle, {Name, DNA}, Acc) ->
    case file:read_line(Handle) of
        {ok, Data} ->
            [S|Tail] = Data,
            case S of
                $> ->
                    Tuple = {Name, lists:flatten(lists:reverse(DNA))},
                    read_fasta(Handle, {string:strip(Tail, right, $\n), []}, [Tuple|Acc]);
                _ ->
                    read_fasta(Handle, {Name, [string:strip(Data, right, $\n)|DNA]}, Acc)
            end;
        eof ->
            Tuple = {Name, lists:flatten(lists:reverse(DNA))},
            [Tuple|Acc];
        _ ->
            {error, invalid_format2}
    end.

% Consensus and Profile
% http://rosalind.info/problems/cons/
%
% Build a profile matrix and consensus string from a set of strings.
% Each string has length N.  The profile matrix has dimensions 4 x N.
% Each element of the matrix contains the number of times the given
% letter (A, T, C, G) is found at that position in the set of strings.
% The output should look like:
%
% <consensus string>
% A: . . . 
% C: . . . 
% G: . . . 
% T: . . . 
%
% The consensus string has length N and contains the most common symbol
% at each position.  There can be multiple consensus strings for each 
% set of input strings.  Return at least one.
cons_and_prof(Strings) ->
    N = length(lists:nth(1, Strings)),
    {A, C, G, T} = profile_matrix(Strings, N, {[], [], [], []}),
    Cons = consensus({A, C, G, T}, []),
    io:format("~s~n", [Cons]),
    io:format("A: "),
    ?FORMAT_LIST(A),
    io:format("C: "),
    ?FORMAT_LIST(C),
    io:format("G: "),
    ?FORMAT_LIST(G),
    io:format("T: "),
    ?FORMAT_LIST(T).

profile_matrix(_, 0, Mx) ->
    Mx;
profile_matrix(Strings, N, {A, C, G, T}) ->
    CA = occurs_at_pos($A, Strings, N),
    CC = occurs_at_pos($C, Strings, N),
    CG = occurs_at_pos($G, Strings, N),
    CT = occurs_at_pos($T, Strings, N),
    profile_matrix(Strings, N-1, {[CA|A], [CC|C], [CG|G], [CT|T]}).

consensus({[], [], [], []}, Acc) ->
    lists:reverse(Acc);
consensus({[A|ARest], [C|CRest], [G|GRest], [T|TRest]}, Acc) ->
    M1 = max_term({A, "A"}, {C, "C"}),
    M2 = max_term(M1, {G, "G"}),
    {_, S} = max_term(M2, {T, "T"}),
    consensus({ARest, CRest, GRest, TRest}, [S|Acc]).
    
max_term({A, AD}, {B, BD}) ->
    case erlang:max(A, B) of
        A ->
            {A, AD};
        B ->
            {B, BD}
    end.

% Count occurrences of Symbol at position N across each string in Strings.
occurs_at_pos(Symbol, Strings, N) ->
    F = 
    fun(String, Count) ->
        case lists:nth(N, String) of
            Symbol ->
                Count+1;
            _ ->
                Count
        end
    end,
    lists:foldl(F, 0, Strings).

% Finding a Motif in DNA
% http://rosalind.info/problems/subs/
% Find positions of T in S (1-indexed).
find_motif(T, S) when is_list(T), is_list(S) ->
    find_motif(list_to_binary(T), list_to_binary(S));
find_motif(T, S) when is_binary(T), is_binary(S), size(T) =< size(S) ->
    find_motif_at(T, S, 1, []).

find_motif_formatted(T, S) ->
    L = find_motif(T, S),
    ?FORMAT_LIST(L).

% Call do_find_motif for every value in S.
find_motif_at(_, <<>>, _, Acc) ->
    lists:reverse(Acc);
find_motif_at(T, S, Pos, Acc) ->
    Acc2 = case do_find_motif(T, S) of
        true ->
            [Pos|Acc];
        _ ->
            Acc
    end,
    <<_:8, Rest/binary>> = S,
    find_motif_at(T, Rest, Pos+1, Acc2).

do_find_motif(<<>>, _) ->
    true;
do_find_motif(<<T:8, TRest/binary>>, <<S:8, SRest/binary>>) when T =:= S ->
    do_find_motif(TRest, SRest);
do_find_motif(_, _) ->
    false.

%    find_motif(T, S, T, 1, S, 1, 0, S, []).
%
%% Match attempt on first position.
%find_motif(T, S, <<TN:8>>, 1, <<SN:8, SRest/binary>>, SI, _SStart, _SRestAtStart, Found) when TN =:= SN ->
%    find_motif(T, S, T, 1, SRest, SI+1, 0, SRest, [SI|Found]);
%
%% First position match of this attempt: T[1] =:= S[SI]
%% This call sets up the SRestAtStart value which is needed when attempting after a failed match.
%find_motif(T, S, <<TN:8, TRest/binary>>, 1, <<SN:8, SRest/binary>>, SI, _SStart, _SRestAtStart, Found) when TN =:= SN ->
%    io:format("1: ~p ~p ~p ~p~n", [TN, TRest, SN, SRest]),
%    find_motif(T, S, TRest, 2, SRest, SI+1, SI, SRest, Found);
%
%% This attempt is a match.  Start another attempt at SSTart+1
%find_motif(T, S, <<TN:8>>, _TI, <<SN:8, SRest/binary>>, _SI, SStart, SRestAtStart, Found) when TN =:= SN ->
%    io:format("2: ~p <<>> ~p ~p~n", [TN, SN, SRest]),
%    find_motif(T, S, T, 1, SRestAtStart, SStart+1, 0, SRest, [SStart|Found]);
%
%% Another position match of this attempt: T[TI] =:= S[SI]
%find_motif(T, S, <<TN:8, TRest/binary>>, TI, <<SN:8, SRest/binary>>, SI, SStart, SRestAtStart, Found) when TN =:= SN ->
%    io:format("3: ~p ~p ~p ~p~n", [TN, TRest, SN, SRest]),
%    find_motif(T, S, TRest, TI+1, SRest, SI+1, SStart, SRestAtStart, Found);
%
%% No match.  Start another attempt at SStart+1.
%find_motif(T, S, TRest, TI, <<_:8, SRest/binary>>, SI, SStart, _, Found) ->
%    io:format("4: ~p ~p ~p ~p~n", [TI, TRest, SI, SRest]),
%    find_motif(T, S, T, 1, SRest, SStart+1, SStart+1, SRest, Found);
%
%% Done searching.
%find_motif(_, _, TRest, TI, <<>>, SI, _, _, Found) ->
%    Found;
%
%find_motif(_, _, TRest, TI, SRest, SI, _, _, Found) ->
%    {error, TRest, TI, SRest, SI, Found}.

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

