<?php

if ($f = fopen("rna_codon_table.csv", "r")) {

    $s = "";
    $count = 0;

    $aa = array();

    while ($row = fgetcsv($f)) {

        $count++;
        $codon = $row[0];
        $symbol = $row[1];

        // Create a lookup of symbols to codons.
        if (isset($aa[$symbol])) {
            $aa[$symbol][] = $codon;
        }
        else {
            $aa[$symbol] = array();
            $aa[$symbol][] = $codon;
        }


        //$s = gen_translate($s, $codon, $symbol);
    }

    fclose($f);

    //print_r($aa);

    $s = gen_lookup($s, $aa);

    echo $s . "\n";
}

function gen_lookup($s, $aa) {
    foreach ($aa as $symbol => $codons) {
        $s .= "symbol_lookup($$symbol) ->\n";
        $s .= "    [";
        foreach ($codons as $codon) {
            $s .= "<<\"$codon\">>,";
        }
        $s = trim($s, ",");
        $s .= "];\n";
    }
    return $s;
}

function gen_translate($s, $codon, $symbol) {
    $s .= "translate(<<\"$codon\", Rest/binary>>, <<Trans/binary>>) ->\n";
    $s .= "    translate(Rest, <<Trans/binary, \"$symbol\">>);\n";
    return $s;
}

