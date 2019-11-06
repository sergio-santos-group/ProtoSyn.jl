for f in $(ls -1 *pdb | sort); do
    RES=$(echo $f | sed 's/.pdb//' | tr '[:lower:]' '[:upper:]')
    echo "reslib[\"$RES\"] = Residue("
    echo "    \"$RES\"",
    echo "    ["
    grep ATOM $f | awk '{printf("        Atom(%2s, \"%4s\", \"%4s\", %8s, %8s, %8s),\n",$2,$11,$3,$6,$7,$8)}'
    echo "    ],"
    echo "    ConnectGraphByName("
    grep CONECT $f | awk '{print "       ", $2, "=>", "[", $0, "],"}' | sed 's/CONECT//'
    echo "    )"
    echo ")"
done
