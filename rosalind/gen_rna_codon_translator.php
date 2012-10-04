<?php

if ($f = fopen("rna_codon_table.csv", "r")) {

    $s = "";
    $count = 0;

    while ($row = fgetcsv($f)) {
        $count++;
        $codon = $row[0];
        $symbol = $row[1];
        $s = append($s, $codon, $symbol);
    }

    fclose($f);

    echo $s . "\n";
}

function append($s, $codon, $symbol) {
    $s .= "translate(<<\"$codon\", Rest/binary>>, <<Trans/binary>>) ->\n";
    $s .= "    translate(Rest, <<Trans/binary, \"$symbol\">>);\n";
    return $s;
}

