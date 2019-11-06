let
using ProtoSyn.Calculators.Forcefield

if !isdefined(@__MODULE__, :ffamber03)
    ffamber03 = ProtoSyn.Calculators.Forcefield.ForcefieldSpec(
        name = "amber03",
        fudgeLJ = 0.5,
        fudgeQQ = 0.83333,
        genpairs = true,
        exclusion_depth = 3)
end


angles = get!(ffamber03.components, :angles, Dict())
improper = get!(ffamber03.components, :improper, Dict())

@ffdef ffamber03.components :atoms AtomType begin
    # type    Z        mass    sigma         epsilon
    ("H0",    1,      1.008,  2.47135e-01,  6.56888e-02),
    ("Br",   35,     79.90 ,  0.00000e+00,  0.00000e+00),
    ( "C",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CA",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CB",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CC",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CK",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CM",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CN",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CQ",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CR",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CT",    6,     12.01 ,  3.39967e-01,  4.57730e-01),
    ("CV",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("CW",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("C*",    6,     12.01 ,  3.39967e-01,  3.59824e-01),
    ("C0",   20,     40.08 ,  3.05240e-01,  1.92376e+00),
    ( "F",    9,     19.00 ,  3.11815e-01,  2.55224e-01),
    ( "H",    1,      1.008,  1.06908e-01,  6.56888e-02),
    ("HC",    1,      1.008,  2.64953e-01,  6.56888e-02),
    ("H1",    1,      1.008,  2.47135e-01,  6.56888e-02),
    ("H2",    1,      1.008,  2.29317e-01,  6.56888e-02),
    ("H3",    1,      1.008,  2.11499e-01,  6.56888e-02),
    ("HA",    1,      1.008,  2.59964e-01,  6.27600e-02),
    ("H4",    1,      1.008,  2.51055e-01,  6.27600e-02),
    ("H5",    1,      1.008,  2.42146e-01,  6.27600e-02),
    ("HO",    1,      1.008,  0.00000e+00,  0.00000e+00),
    ("HS",    1,      1.008,  1.06908e-01,  6.56888e-02),
    ("HW",    1,      1.008,  0.00000e+00,  0.00000e+00),
    ("HP",    1,      1.008,  1.95998e-01,  6.56888e-02),
    ( "I",   53,    126.9  ,  4.18722e-01,  1.67360e+00),
    ("Cl",   17,     35.45 ,  4.40104e-01,  4.18400e-01),
    ("Na",   11,     22.99 ,  3.32840e-01,  1.15897e-02),
    ("IB",    0,    131.0  ,  8.90899e-01,  4.18400e-01),
    ("MG",   12,     24.305,  1.41225e-01,  3.74342e+00),
    ( "N",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ("NA",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ("NB",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ("NC",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ("N2",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ("N3",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ("N*",    7,     14.01 ,  3.25000e-01,  7.11280e-01),
    ( "O",    8,     16.00 ,  2.95992e-01,  8.78640e-01),
    ("OW",    8,     16.00 ,  3.15061e-01,  6.36386e-01),
    ("OH",    8,     16.00 ,  3.06647e-01,  8.80314e-01),
    ("OS",    8,     16.00 ,  3.00001e-01,  7.11280e-01),
    ("O2",    8,     16.00 ,  2.95992e-01,  8.78640e-01),
    ( "P",   15,     30.97 ,  3.74177e-01,  8.36800e-01),
    ( "S",   16,     32.06 ,  3.56359e-01,  1.04600e+00),
    ("SH",   16,     32.06 ,  3.56359e-01,  1.04600e+00),
    ("CU",   29,     63.55 ,  3.39967e-01,  3.59824e-01),
    ("FE",   26,     55.00 ,  0.00000e+00,  0.00000e+00),
    ( "K",   19,     39.10 ,  4.73602e-01,  1.37235e-03),
    ("Rb",   37,     85.47 ,  5.26699e-01,  7.11280e-04),
    ("Cs",   55,    132.91 ,  6.04920e-01,  3.37230e-04)
end
#@ffdef ffamber03 :atoms AtomType ("XXX",   99,    999.9 ,  99.9,  9.9)

@ffdef ffamber03.components :bonds HarmonicBondType begin
    # at1   at2     r0       k
    ("CT", "H0",    0.1090,  284512.0),  #  03GLY changed from 331 bsd on NMA nmodes; AA, SUGARS
    ( "C",  "C",    0.1525,  259408.0),  #  new99
    ( "C", "OS",    0.1323,  376560.0),  #  new99
    ( "C", "H4",    0.1080,  307105.6),  #  new99
    ( "C", "H5",    0.1080,  307105.6),  #  new99
    ("CA", "OH",    0.1364,  376560.0),  #  new99
    ("CM", "OS",    0.1240,  401664.0),  #  new99
    ("Cl", "CT",    0.1766,  194137.6),  #  new99
    ("Br", "CT",    0.1944,  133051.2),  #  new99
    ( "I", "CT",    0.2166,  123846.4),  #  new99
    ( "F", "CA",    0.1359,  323004.8),  #  new99
    ("Cl", "CA",    0.1727,  161502.4),  #  new99
    ( "I", "CA",    0.2075,  143092.8),  #  new99
    ("Br", "CA",    0.1890,  143929.6),  #  new99
    ("OW", "HW",   0.09572,  462750.4),  #  P water
    ("HW", "HW",   0.15136,  462750.4),  #  P water
    ( "C", "CA",   0.14090,  392459.2),  #  7,(1986),230; TYR
    ( "C", "CB",   0.14190,  374049.6),  #  7,(1986),230; GUA
    ( "C", "CM",   0.14440,  343088.0),  #  7,(1986),230; THY,URA
    ( "C", "CT",   0.15220,  265265.6),  #  7,(1986),230; AA
    ( "C", "NA",   0.13880,  349782.4),  #  7,(1986),230; GUA.URA
    ( "C", "NC",   0.13580,  382417.6),  #  7,(1986),230; CYT
    ( "C",  "O",   0.12290,  476976.0),  #  7,(1986),230; AA,CYT,GUA,THY,URA
    ( "C", "O2",   0.12500,  548940.8),  #  7,(1986),230; GLU,ASP
    ( "C", "OH",   0.13640,  376560.0),  #  7,(1986),230; TYR
    ("CA", "CA",   0.14000,  392459.2),  #  7,(1986),230; BENZENE,PHE,TRP,TYR
    ("CA", "CB",   0.14040,  392459.2),  #  7,(1986),230; ADE,TRP
    ("CA", "CM",   0.14330,  357313.6),  #  7,(1986),230; CYT
    ("CA", "CT",   0.15100,  265265.6),  #  7,(1986),230; PHE,TYR
    ("CA", "HA",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; PHE,TRP,TYR
    ("CA", "H4",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; no assigned
    ("CA", "N2",   0.13400,  402500.8),  #  7,(1986),230; ARG,CYT,GUA
    ("CA", "NA",   0.13810,  357313.6),  #  7,(1986),230; GUA
    ("CA", "NC",   0.13390,  404174.4),  #  7,(1986),230; ADE,CYT,GUA
    ("CB", "CB",   0.13700,  435136.0),  #  7,(1986),230; ADE,GUA
    ("CB", "NB",   0.13910,  346435.2),  #  7,(1986),230; ADE,GUA
    ("CB", "NC",   0.13540,  385764.8),  #  7,(1986),230; ADE,GUA
    ("CK", "H5",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; ADE,GUA
    ("CK", "NB",   0.13040,  442667.2),  #  7,(1986),230; ADE,GUA
    ("CM", "CM",   0.13500,  459403.2),  #  7,(1986),230; CYT,THY,URA
    ("CM", "CT",   0.15100,  265265.6),  #  7,(1986),230; THY
    ("CM", "HA",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; CYT,URA
    ("CM", "H4",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; CYT,URA
    ("CM", "H5",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; not assigned
    ("CQ", "H5",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; ADE
    ("CQ", "NC",   0.13240,  420073.6),  #  7,(1986),230; ADE
    ("CT", "CT",   0.15260,  259408.0),  #  7,(1986),230; AA, SUGARS
    ("CT", "HC",   0.10900,  284512.0),  #  ged from 331 bsd on NMA nmodes; AA, SUGARS
    ("CT", "H1",   0.10900,  284512.0),  #  ged from 331 bsd on NMA nmodes; AA, RIBOSE
    ("CT", "H2",   0.10900,  284512.0),  #  ged from 331 bsd on NMA nmodes; SUGARS
    ("CT", "H3",   0.10900,  284512.0),  #  ged from 331 bsd on NMA nmodes; not assigned
    ("CT", "HP",   0.10900,  284512.0),  #  nged from 331; AA-lysine, methyl ammonium cation
    ("CT", "N2",   0.14630,  282001.6),  #  7,(1986),230; ARG
    ("CT", "OH",   0.14100,  267776.0),  #  7,(1986),230; SUGARS
    ("CT", "OS",   0.14100,  267776.0),  #  7,(1986),230; NUCLEIC ACIDS
    ( "H", "N2",   0.10100,  363171.2),  #  7,(1986),230; ADE,CYT,GUA,ARG
    ( "H", "NA",   0.10100,  363171.2),  #  7,(1986),230; GUA,URA,HIS
    ("HO", "OH",   0.09600,  462750.4),  #  7,(1986),230; SUGARS,SER,TYR
    ("HO", "OS",   0.09600,  462750.4),  #  7,(1986),230; NUCLEOTIDE ENDS
    ("O2",  "P",   0.14800,  439320.0),  #  7,(1986),230; NA PHOSPHATES
    ("OH",  "P",   0.16100,  192464.0),  #  7,(1986),230; NA PHOSPHATES
    ("OS",  "P",   0.16100,  192464.0),  #  7,(1986),230; NA PHOSPHATES
    ( "C",  "N",   0.13350,  410032.0),  #  7,(1986),230; AA
    ("CA", "CN",   0.14000,  392459.2),  #  7,(1986),230; TRP
    ("CB", "CN",   0.14190,  374049.6),  #  7,(1986),230; TRP
    ("CC", "CT",   0.15040,  265265.6),  #  7,(1986),230; HIS
    ("CC", "CV",   0.13750,  428441.6),  #  7,(1986),230; HIS(delta)
    ("CC", "CW",   0.13710,  433462.4),  #  7,(1986),230; HIS(epsilon)
    ("CC", "NA",   0.13850,  353129.6),  #  7,(1986),230; HIS
    ("CC", "NB",   0.13940,  343088.0),  #  7,(1986),230; HIS
    ("CN", "NA",   0.13800,  358150.4),  #  7,(1986),230; TRP
    ("CR", "H5",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes;HIS
    ("CR", "NA",   0.13430,  399153.6),  #  7,(1986),230; HIS
    ("CR", "NB",   0.13350,  408358.4),  #  7,(1986),230; HIS
    ("CT",  "N",   0.14490,  282001.6),  #  7,(1986),230; AA
    ("CT", "N3",   0.14710,  307105.6),  #  7,(1986),230; LYS
    ("CT",  "S",   0.18100,  189953.6),  #  ged from 222.0 based on dimethylS nmodes
    ("CT", "SH",   0.18100,  198321.6),  #  ged from 222.0 based on methanethiol nmodes
    ("CV", "H4",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes; HIS
    ("CV", "NB",   0.13940,  343088.0),  #  7,(1986),230; HIS
    ("CW", "H4",   0.10800,  307105.6),  #  ged from 340. bsd on C6H6 nmodes;HIS(epsilon,+)
    ("CW", "NA",   0.13810,  357313.6),  #  7,(1986),230; HIS,TRP
    ( "H",  "N",   0.10100,  363171.2),  #  7,(1986),230; AA
    ( "H", "N3",   0.10100,  363171.2),  #  7,(1986),230; LYS
    ("HS", "SH",   0.13360,  229283.2),  #  7,(1986),230; CYS
    ( "S",  "S",   0.20380,  138908.8),  #  7,(1986),230; CYX   (SCHERAGA)
    ("CT",  "F",   0.13800,  307105.6)   #  13,(1992),963;CF4; R0=1.332 FOR CHF3
end

@ffdef ffamber03.components :angles HarmonicAngleType begin
    ("H0" , "CT" , "H0"  , deg2rad(109.500),  292.880),  #  03GLY
    ("H0" , "CT" , "N"   , deg2rad(109.500),  418.400),  #  03GLY AA general changed based on NMA nmodes
    ("C"  , "CT" , "H0"  , deg2rad(109.500),  418.400),  #  03GLY AA general changed based on NMA nmodes
    ("HW" , "OW" , "HW"  , deg2rad(104.520),  836.800),  #  TIP3P water
    ("HW" , "HW" , "OW"  , deg2rad(127.740),    0.000),  #  (found in crystallographic water with 3 bonds)
    ("C"  , "C"  , "O"   , deg2rad(120.000),  669.440),  #  new99
    ("C"  , "C"  , "OH"  , deg2rad(120.000),  669.440),  #  new99
    ("CT" , "C"  , "CT"  , deg2rad(117.000),  527.184),  #  new99
    ("CT" , "C"  , "OS"  , deg2rad(115.000),  669.440),  #  new99
    ("O"  , "C"  , "OS"  , deg2rad(125.000),  669.440),  #  new99
    ("H4" , "C"  , "C"   , deg2rad(120.000),  418.400),  #  new99
    ("H4" , "C"  , "CM"  , deg2rad(115.000),  418.400),  #  new99
    ("H4" , "C"  , "CT"  , deg2rad(115.000),  418.400),  #  new99
    ("H4" , "C"  , "O"   , deg2rad(120.000),  418.400),  #  new99
    ("H4" , "C"  , "OH"  , deg2rad(120.000),  418.400),  #  new99
    ("H5" , "C"  , "N"   , deg2rad(120.000),  418.400),  #  new99
    ("H5" , "C"  , "O"   , deg2rad(119.000),  418.400),  #  new99
    ("H5" , "C"  , "OH"  , deg2rad(107.000),  418.400),  #  new99
    ("H5" , "C"  , "OS"  , deg2rad(107.000),  418.400),  #  new99
    ("CA" , "CA" , "OH"  , deg2rad(120.000),  585.760),  #  new99
    ("CA" , "OH" , "HO"  , deg2rad(113.000),  418.400),  #  new99
    ("F"  , "CA" , "CA"  , deg2rad(121.000),  585.760),  #  new99
    ("Cl" , "CA" , "CA"  , deg2rad(118.800),  585.760),  #  new99
    ("Br" , "CA" , "CA"  , deg2rad(118.800),  585.760),  #  new99
    ("I"  , "CA" , "CA"  , deg2rad(118.800),  585.760),  #  new99
    ("CM" , "CM" , "OS"  , deg2rad(125.000),  669.440),  #  new99
    ("H4" , "CM" , "OS"  , deg2rad(113.000),  418.400),  #  new99
    ("HA" , "CM" , "HA"  , deg2rad(120.000),  292.880),  #  new99
    ("HA" , "CM" , "CT"  , deg2rad(120.000),  418.400),  #  new99
    ("H1" , "CT" , "CM"  , deg2rad(109.500),  418.400),  #  new99
    ("HC" , "CT" , "CM"  , deg2rad(109.500),  418.400),  #  new99
    ("C"  , "CT" , "OS"  , deg2rad(109.500),  502.080),  #  new99
    ("CM" , "CT" , "CT"  , deg2rad(111.000),  527.184),  #  new99
    ("CM" , "CT" , "OS"  , deg2rad(109.500),  418.400),  #  new99
    ("CT" , "CT" , "CA"  , deg2rad(114.000),  527.184),  #  new99
    ("OS" , "CT" , "OS"  , deg2rad(101.000),  502.080),  #  new99
    ("F"  , "CT" , "CT"  , deg2rad(109.000),  418.400),  #  new99
    ("F"  , "CT" , "H2"  , deg2rad(109.500),  418.400),  #  new99
    ("Cl" , "CT" , "CT"  , deg2rad(108.500),  418.400),  #  new99
    ("Cl" , "CT" , "H1"  , deg2rad(108.500),  418.400),  #  new99
    ("Br" , "CT" , "CT"  , deg2rad(108.000),  418.400),  #  new99
    ("Br" , "CT" , "H1"  , deg2rad(106.500),  418.400),  #  new99
    ("I"  , "CT" , "CT"  , deg2rad(106.000),  418.400),  #  new99
    ("CB" , "C"  , "NA"  , deg2rad(111.300),  585.760),  #  NA
    ("CB" , "C"  , "O"   , deg2rad(128.800),  669.440),  #
    ("CM" , "C"  , "NA"  , deg2rad(114.100),  585.760),  #
    ("CM" , "C"  , "O"   , deg2rad(125.300),  669.440),  #
    ("CT" , "C"  , "O"   , deg2rad(120.400),  669.440),  #
    ("CT" , "C"  , "O2"  , deg2rad(117.000),  585.760),  #
    ("CT" , "C"  , "OH"  , deg2rad(110.000),  669.440),  #  new99
    ("N*" , "C"  , "NA"  , deg2rad(115.400),  585.760),  #
    ("N*" , "C"  , "NC"  , deg2rad(118.600),  585.760),  #
    ("N*" , "C"  , "O"   , deg2rad(120.900),  669.440),  #
    ("NA" , "C"  , "O"   , deg2rad(120.600),  669.440),  #
    ("NC" , "C"  , "O"   , deg2rad(122.500),  669.440),  #
    ("CT" , "C"  , "N"   , deg2rad(116.600),  585.760),  #  AA                   general
    ("N"  , "C"  , "O"   , deg2rad(122.900),  669.440),  #  AA                   general
    ("O"  , "C"  , "O"   , deg2rad(126.000),  669.440),  #  AA                   COO-        terminal          residues
    ("O2" , "C"  , "O2"  , deg2rad(126.000),  669.440),  #  AA                   GLU         (SCH              JPC          79,2379)
    ("O"  , "C"  , "OH"  , deg2rad(120.000),  669.440),  #
    ("CA" , "C"  , "CA"  , deg2rad(120.000),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      tyr
    ("CA" , "C"  , "OH"  , deg2rad(120.000),  585.760),  #  AA                   tyr
    ("C"  , "CA" , "CA"  , deg2rad(120.000),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes
    ("CA" , "CA" , "CA"  , deg2rad(120.000),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes
    ("CA" , "CA" , "CB"  , deg2rad(120.000),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes
    ("CA" , "CA" , "CT"  , deg2rad(120.000),  585.760),  #
    ("CA" , "CA" , "HA"  , deg2rad(120.000),  418.400),  #  new99
    ("CA" , "CA" , "H4"  , deg2rad(120.000),  418.400),  #  new99
    ("CB" , "CA" , "HA"  , deg2rad(120.000),  418.400),  #  new99
    ("CB" , "CA" , "H4"  , deg2rad(120.000),  418.400),  #  new99
    ("CB" , "CA" , "N2"  , deg2rad(123.500),  585.760),  #
    ("CB" , "CA" , "NC"  , deg2rad(117.300),  585.760),  #
    ("CM" , "CA" , "N2"  , deg2rad(120.100),  585.760),  #
    ("CM" , "CA" , "NC"  , deg2rad(121.500),  585.760),  #
    ("N2" , "CA" , "NA"  , deg2rad(116.000),  585.760),  #
    ("N2" , "CA" , "NC"  , deg2rad(119.300),  585.760),  #
    ("NA" , "CA" , "NC"  , deg2rad(123.300),  585.760),  #
    ("C"  , "CA" , "HA"  , deg2rad(120.000),  418.400),  #  new99                tyr
    ("N2" , "CA" , "N2"  , deg2rad(120.000),  585.760),  #  AA                   arg
    ("CN" , "CA" , "HA"  , deg2rad(120.000),  418.400),  #  new99                trp
    ("CA" , "CA" , "CN"  , deg2rad(120.000),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      trp
    ("C"  , "CB" , "CB"  , deg2rad(119.200),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  NA      gua
    ("C"  , "CB" , "NB"  , deg2rad(130.000),  585.760),  #
    ("CA" , "CB" , "CB"  , deg2rad(117.300),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  NA      ade
    ("CA" , "CB" , "NB"  , deg2rad(132.400),  585.760),  #
    ("CB" , "CB" , "N*"  , deg2rad(106.200),  585.760),  #
    ("CB" , "CB" , "NB"  , deg2rad(110.400),  585.760),  #
    ("CB" , "CB" , "NC"  , deg2rad(127.700),  585.760),  #
    ("N*" , "CB" , "NC"  , deg2rad(126.200),  585.760),  #
    ("C*" , "CB" , "CA"  , deg2rad(134.900),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      trp
    ("C*" , "CB" , "CN"  , deg2rad(108.800),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      trp
    ("CA" , "CB" , "CN"  , deg2rad(116.200),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      trp
    ("H5" , "CK" , "N*"  , deg2rad(123.050),  418.400),  #  new99
    ("H5" , "CK" , "NB"  , deg2rad(123.050),  418.400),  #  new99
    ("N*" , "CK" , "NB"  , deg2rad(113.900),  585.760),  #
    ("C"  , "CM" , "CM"  , deg2rad(120.700),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  NA      thy
    ("C"  , "CM" , "CT"  , deg2rad(119.700),  585.760),  #
    ("C"  , "CM" , "HA"  , deg2rad(119.700),  418.400),  #  new99
    ("C"  , "CM" , "H4"  , deg2rad(119.700),  418.400),  #  new99
    ("CA" , "CM" , "CM"  , deg2rad(117.000),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  NA      cyt
    ("CA" , "CM" , "HA"  , deg2rad(123.300),  418.400),  #  new99
    ("CA" , "CM" , "H4"  , deg2rad(123.300),  418.400),  #  new99
    ("CM" , "CM" , "CT"  , deg2rad(119.700),  585.760),  #
    ("CM" , "CM" , "HA"  , deg2rad(119.700),  418.400),  #  new99
    ("CM" , "CM" , "H4"  , deg2rad(119.700),  418.400),  #  new99
    ("CM" , "CM" , "N*"  , deg2rad(121.200),  585.760),  #
    ("H4" , "CM" , "N*"  , deg2rad(119.100),  418.400),  #  new99
    ("H5" , "CQ" , "NC"  , deg2rad(115.450),  418.400),  #  new99
    ("NC" , "CQ" , "NC"  , deg2rad(129.100),  585.760),  #
    ("CM" , "CT" , "HC"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("CT" , "CT" , "CT"  , deg2rad(109.500),  334.720),  #
    ("CT" , "CT" , "HC"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("CT" , "CT" , "H1"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("CT" , "CT" , "H2"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("CT" , "CT" , "HP"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("CT" , "CT" , "N*"  , deg2rad(109.500),  418.400),  #
    ("CT" , "CT" , "OH"  , deg2rad(109.500),  418.400),  #
    ("CT" , "CT" , "OS"  , deg2rad(109.500),  418.400),  #
    ("HC" , "CT" , "HC"  , deg2rad(109.500),  292.880),  #
    ("H1" , "CT" , "H1"  , deg2rad(109.500),  292.880),  #
    ("HP" , "CT" , "HP"  , deg2rad(109.500),  292.880),  #  AA                   lys,        ch3nh4+
    ("H2" , "CT" , "N*"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("H1" , "CT" , "N*"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("H1" , "CT" , "OH"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("H1" , "CT" , "OS"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("H2" , "CT" , "OS"  , deg2rad(109.500),  418.400),  #  changed              based       on                NMA          nmodes
    ("N*" , "CT" , "OS"  , deg2rad(109.500),  418.400),  #
    ("H1" , "CT" , "N"   , deg2rad(109.500),  418.400),  #  AA                   general     changed           based        on           NMA           nmodes
    ("C"  , "CT" , "H1"  , deg2rad(109.500),  418.400),  #  AA                   general     changed           based        on           NMA           nmodes
    ("C"  , "CT" , "HP"  , deg2rad(109.500),  418.400),  #  AA                   zwitterion  changed           based        on           NMA           nmodes
    ("H1" , "CT" , "S"   , deg2rad(109.500),  418.400),  #  AA                   cys         changed           based        on           NMA           nmodes
    ("H1" , "CT" , "SH"  , deg2rad(109.500),  418.400),  #  AA                   cyx         changed           based        on           NMA           nmodes
    ("CT" , "CT" , "S"   , deg2rad(114.700),  418.400),  #  AA                   cyx         (SCHERAGA         JPC          79,1428)
    ("CT" , "CT" , "SH"  , deg2rad(108.600),  418.400),  #  AA                   cys
    ("H2" , "CT" , "H2"  , deg2rad(109.500),  292.880),  #  AA                   lys
    ("H1" , "CT" , "N2"  , deg2rad(109.500),  418.400),  #  AA                   arg         changed           based        on           NMA           nmodes
    ("HP" , "CT" , "N3"  , deg2rad(109.500),  418.400),  #  AA                   lys,        ch3nh3+,          changed      based        on            NMA      nmodes
    ("CA" , "CT" , "CT"  , deg2rad(114.000),  527.184),  #  AA                   phe         tyr               (SCH         JPC          79,2379)
    ("C"  , "CT" , "HC"  , deg2rad(109.500),  418.400),  #  AA                   gln         changed           based        on           NMA           nmodes
    ("C"  , "CT" , "N"   , deg2rad(110.100),  527.184),  #  AA                   general
    ("CT" , "CT" , "N2"  , deg2rad(111.200),  669.440),  #  AA                   arg         (JCP              76,          1439)
    ("CT" , "CT" , "N"   , deg2rad(109.700),  669.440),  #  AA                   ala,        general           (JACS        94,          2657)
    ("C"  , "CT" , "CT"  , deg2rad(111.100),  527.184),  #  AA                   general
    ("CA" , "CT" , "HC"  , deg2rad(109.500),  418.400),  #  AA                   tyr         changed           based        on           NMA           nmodes
    ("CT" , "CT" , "N3"  , deg2rad(111.200),  669.440),  #  AA                   lys         (JCP              76,          1439)
    ("CC" , "CT" , "CT"  , deg2rad(113.100),  527.184),  #  AA                   his
    ("CC" , "CT" , "HC"  , deg2rad(109.500),  418.400),  #  AA                   his         changed           based        on           NMA           nmodes
    ("C"  , "CT" , "N3"  , deg2rad(111.200),  669.440),  #  AA                   amino       terminal          residues
    ("C*" , "CT" , "CT"  , deg2rad(115.600),  527.184),  #  AA                   trp
    ("C*" , "CT" , "HC"  , deg2rad(109.500),  418.400),  #  AA                   trp         changed           based        on           NMA           nmodes
    ("CT" , "CC" , "NA"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CT" , "CC" , "CV"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CT" , "CC" , "NB"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CV" , "CC" , "NA"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CW" , "CC" , "NA"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CW" , "CC" , "NB"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CT" , "CC" , "CW"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("H5" , "CR" , "NA"  , deg2rad(120.000),  418.400),  #  new99                his
    ("H5" , "CR" , "NB"  , deg2rad(120.000),  418.400),  #  new99                his
    ("NA" , "CR" , "NA"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("NA" , "CR" , "NB"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CC" , "CV" , "H4"  , deg2rad(120.000),  418.400),  #  new99                his
    ("CC" , "CV" , "NB"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("H4" , "CV" , "NB"  , deg2rad(120.000),  418.400),  #  new99                his
    ("CC" , "CW" , "H4"  , deg2rad(120.000),  418.400),  #  new99                his
    ("CC" , "CW" , "NA"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("H4" , "CW" , "NA"  , deg2rad(120.000),  418.400),  #  new99                his
    ("C*" , "CW" , "H4"  , deg2rad(120.000),  418.400),  #  new99                trp
    ("C*" , "CW" , "NA"  , deg2rad(108.700),  585.760),  #  AA                   trp
    ("CT" , "C*" , "CW"  , deg2rad(125.000),  585.760),  #  AA                   trp
    ("CB" , "C*" , "CT"  , deg2rad(128.600),  585.760),  #  AA                   trp
    ("CB" , "C*" , "CW"  , deg2rad(106.400),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      trp
    ("CA" , "CN" , "NA"  , deg2rad(132.800),  585.760),  #  AA                   trp
    ("CB" , "CN" , "NA"  , deg2rad(104.400),  585.760),  #  AA                   trp
    ("CA" , "CN" , "CB"  , deg2rad(122.700),  527.184),  #  changed              from        85.0              bsd          on           C6H6          nmodes;  AA      trp
    ("C"  , "N"  , "CT"  , deg2rad(121.900),  418.400),  #  AA                   general
    ("C"  , "N"  , "H"   , deg2rad(120.000),  418.400),  #  new99                general,    gln,              asn,changed  based        on            NMA      nmodes
    ("CT" , "N"  , "H"   , deg2rad(118.040),  418.400),  #  new99                general,    changed           based        on           NMA           nmodes
    ("CT" , "N"  , "CT"  , deg2rad(118.000),  418.400),  #  AA                   pro         (DETAR            JACS         99,1232)
    ("H"  , "N"  , "H"   , deg2rad(120.000),  292.880),  #  ade,cyt,gua,gln,asn  **
    ("C"  , "N*" , "CM"  , deg2rad(121.600),  585.760),  #
    ("C"  , "N*" , "CT"  , deg2rad(117.600),  585.760),  #
    ("C"  , "N*" , "H"   , deg2rad(119.200),  418.400),  #  new99
    ("CB" , "N*" , "CK"  , deg2rad(105.400),  585.760),  #
    ("CB" , "N*" , "CT"  , deg2rad(125.800),  585.760),  #
    ("CB" , "N*" , "H"   , deg2rad(125.800),  418.400),  #  new99
    ("CK" , "N*" , "CT"  , deg2rad(128.800),  585.760),  #
    ("CK" , "N*" , "H"   , deg2rad(128.800),  418.400),  #  new99                for         unmethylated      n.a.         bases,chngd  bsd           NMA      nmodes
    ("CM" , "N*" , "CT"  , deg2rad(121.200),  585.760),  #
    ("CM" , "N*" , "H"   , deg2rad(121.200),  418.400),  #  new99                for         unmethylated      n.a.         bases,chngd  bsd           NMA      nmodes
    ("CA" , "N2" , "H"   , deg2rad(120.000),  418.400),  #  new99
    ("H"  , "N2" , "H"   , deg2rad(120.000),  292.880),  #
    ("CT" , "N2" , "H"   , deg2rad(118.400),  418.400),  #  new99                arg
    ("CA" , "N2" , "CT"  , deg2rad(123.200),  418.400),  #  AA                   arg
    ("CT" , "N3" , "H"   , deg2rad(109.500),  418.400),  #  AA                   lys,        changed           based        on           NMA           nmodes
    ("CT" , "N3" , "CT"  , deg2rad(109.500),  418.400),  #  AA                   pro/nt
    ("H"  , "N3" , "H"   , deg2rad(109.500),  292.880),  #  AA                   lys,        AA(end)
    ("C"  , "NA" , "C"   , deg2rad(126.400),  585.760),  #
    ("C"  , "NA" , "CA"  , deg2rad(125.200),  585.760),  #
    ("C"  , "NA" , "H"   , deg2rad(116.800),  418.400),  #  new99
    ("CA" , "NA" , "H"   , deg2rad(118.000),  418.400),  #  new99
    ("CC" , "NA" , "CR"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CC" , "NA" , "H"   , deg2rad(120.000),  418.400),  #  new99                his
    ("CR" , "NA" , "CW"  , deg2rad(120.000),  585.760),  #  AA                   his
    ("CR" , "NA" , "H"   , deg2rad(120.000),  418.400),  #  new99                his
    ("CW" , "NA" , "H"   , deg2rad(120.000),  418.400),  #  new99                his
    ("CN" , "NA" , "CW"  , deg2rad(111.600),  585.760),  #  AA                   trp
    ("CN" , "NA" , "H"   , deg2rad(123.100),  418.400),  #  new99                trp
    ("CB" , "NB" , "CK"  , deg2rad(103.800),  585.760),  #
    ("CC" , "NB" , "CR"  , deg2rad(117.000),  585.760),  #  AA                   his
    ("CR" , "NB" , "CV"  , deg2rad(117.000),  585.760),  #  AA                   his
    ("C"  , "NC" , "CA"  , deg2rad(120.500),  585.760),  #
    ("CA" , "NC" , "CB"  , deg2rad(112.200),  585.760),  #
    ("CA" , "NC" , "CQ"  , deg2rad(118.600),  585.760),  #
    ("CB" , "NC" , "CQ"  , deg2rad(111.000),  585.760),  #
    ("C"  , "OH" , "HO"  , deg2rad(113.000),  418.400),  #  new99
    ("CT" , "OH" , "HO"  , deg2rad(108.500),  460.240),  #
    ("HO" , "OH" , "P"   , deg2rad(108.500),  376.560),  #
    ("CT" , "OS" , "CT"  , deg2rad(109.500),  502.080),  #
    ("CT" , "OS" , "P"   , deg2rad(120.500),  836.800),  #
    ("P"  , "OS" , "P"   , deg2rad(120.500),  836.800),  #
    ("O2" , "P"  , "OH"  , deg2rad(108.230),  376.560),  #
    ("O2" , "P"  , "O2"  , deg2rad(119.900), 1171.520),  #
    ("O2" , "P"  , "OS"  , deg2rad(108.230),  836.800),  #
    ("OH" , "P"  , "OS"  , deg2rad(102.600),  376.560),  #
    ("OS" , "P"  , "OS"  , deg2rad(102.600),  376.560),  #
    ("CT" , "S"  , "CT"  , deg2rad( 98.900),  518.816),  #  AA                   met
    ("CT" , "S"  , "S"   , deg2rad(103.700),  569.024),  #  AA                   cyx         (SCHERAGA         JPC          79,1428)
    ("CT" , "SH" , "HS"  , deg2rad( 96.000),  359.824),  #  changed              from        44.0              based        on           methanethiol  nmodes
    ("HS" , "SH" , "HS"  , deg2rad( 92.070),  292.880),  #  AA                   cys
    ("F"  , "CT" , "F"   , deg2rad(109.100),  644.336),  #  JCC,13,(1992),963;
    ("F"  , "CT" , "H1"  , deg2rad(109.500),  418.400),  #  new99
    ("N"  , "C"  , "N"   , deg2rad(120.000),  585.760)   #  Added                for         Urea              (same        as           N2-CA-N2)     -        EJS
end


improper["CA:CA:CA:OH"] = CosineDihedralType("CA:CA:CA:OH",  180.00,   4.60240,  2)  #  new99
improper["H5:O:C:OH"] = CosineDihedralType("H5:O:C:OH",  180.00,   4.60240,  2)  #  new99
improper["H5:O:C:OS"] = CosineDihedralType("H5:O:C:OS",  180.00,   4.60240,  2)  #  new99
improper["CM:CT:CM:HA"] = CosineDihedralType("CM:CT:CM:HA",  180.00,   4.60240,  2)  #  new99
improper["CA:CA:CA:Br"] = CosineDihedralType("CA:CA:CA:Br",  180.00,   4.60240,  2)  #  new99
improper["CM:H4:C:O"] = CosineDihedralType("CM:H4:C:O",  180.00,   4.60240,  2)  #  new99
improper[ "C:CT:N:H"] = CosineDihedralType( "C:CT:N:H",  180.00,   4.60240,  2)  #  new99
improper[ "C:CT:N:O"] = CosineDihedralType( "C:CT:N:O",  180.00,   4.60240,  2)  #  new99
improper["CB:CK:N*:CT"] = CosineDihedralType("CB:CK:N*:CT",  180.00,   4.18400,  2)  # 
improper["CK:CB:N*:CT"] = CosineDihedralType("CK:CB:N*:CT",  180.00,   4.18400,  2)  # 
improper[ "C:CM:N*:CT"] = CosineDihedralType( "C:CM:N*:CT",  180.00,   4.18400,  2)  #  dac guess, 9/94
improper["CM:C:CM:CT"] = CosineDihedralType("CM:C:CM:CT",  180.00,   4.60240,  2)  # 
improper["CT:O:C:OH"] = CosineDihedralType("CT:O:C:OH",  180.00,  43.93200,  2)  # 
improper["NA:CV:CC:CT"] = CosineDihedralType("NA:CV:CC:CT",  180.00,   4.60240,  2)  # 
improper["NB:CW:CC:CT"] = CosineDihedralType("NB:CW:CC:CT",  180.00,   4.60240,  2)  # 
improper["NA:CW:CC:CT"] = CosineDihedralType("NA:CW:CC:CT",  180.00,   4.60240,  2)  # 
improper["CW:CB:C*:CT"] = CosineDihedralType("CW:CB:C*:CT",  180.00,   4.60240,  2)  # 
improper["CA:CA:CA:CT"] = CosineDihedralType("CA:CA:CA:CT",  180.00,   4.60240,  2)  # 
improper[ "C:CM:CM:CT"] = CosineDihedralType( "C:CM:CM:CT",  180.00,   4.60240,  2)  #  dac guess, 9/94
improper["NC:CM:CA:N2"] = CosineDihedralType("NC:CM:CA:N2",  180.00,   4.60240,  2)  #  dac guess, 9/94
improper["CB:NC:CA:N2"] = CosineDihedralType("CB:NC:CA:N2",  180.00,   4.60240,  2)  #  dac, 10/94
improper["NA:NC:CA:N2"] = CosineDihedralType("NA:NC:CA:N2",  180.00,   4.60240,  2)  #  dac, 10/94
improper["CA:CA:C:OH"] = CosineDihedralType("CA:CA:C:OH",  180.00,   4.60240,  2)  # 
improper["CT:CV:CC:NA"] = CosineDihedralType("CT:CV:CC:NA",  180.00,   4.60240,  2)  # 
improper["CT:CW:CC:NB"] = CosineDihedralType("CT:CW:CC:NB",  180.00,   4.60240,  2)  # 
improper["CT:CW:CC:NA"] = CosineDihedralType("CT:CW:CC:NA",  180.00,   4.60240,  2)  # 
improper["CB:CT:C*:CW"] = CosineDihedralType("CB:CT:C*:CW",  180.00,   4.60240,  2)  # 
improper["CM:N2:CA:NC"] = CosineDihedralType("CM:N2:CA:NC",  180.00,   4.60240,  2)  # 
improper["CB:N2:CA:NC"] = CosineDihedralType("CB:N2:CA:NC",  180.00,   4.60240,  2)  # 
improper["N2:NA:CA:NC"] = CosineDihedralType("N2:NA:CA:NC",  180.00,   4.60240,  2)  # 
improper[ "N:N:C:O"] = CosineDihedralType( "N:N:C:O",  180.00,  43.93200,  2)  #  urea
improper[ "X:O2:C:O2"] = CosineDihedralType( "X:O2:C:O2",  180.00,  43.93200,  2)  #  JCC,7,(1986),230
improper[ "X:N2:CA:N2"] = CosineDihedralType( "X:N2:CA:N2",  180.00,  43.93200,  2)  #  JCC,7,(1986),230
improper[ "X:CT:N:CT"] = CosineDihedralType( "X:CT:N:CT",  180.00,   4.18400,  2)  #  JCC,7,(1986),230
improper[ "X:X:C:O"] = CosineDihedralType( "X:X:C:O",  180.00,  43.93200,  2)  #  JCC,7,(1986),230
improper[ "X:X:N:H"] = CosineDihedralType( "X:X:N:H",  180.00,   4.18400,  2)  #  JCC,7,(1986),230
improper[ "X:X:N2:H"] = CosineDihedralType( "X:X:N2:H",  180.00,   4.18400,  2)  #  JCC,7,(1986),230
improper[ "X:X:NA:H"] = CosineDihedralType( "X:X:NA:H",  180.00,   4.18400,  2)  #  JCC,7,(1986),230
improper[ "X:X:CA:HA"] = CosineDihedralType( "X:X:CA:HA",  180.00,   4.60240,  2)  #  bsd.on C6H6 nmodes
improper[ "X:X:CW:H4"] = CosineDihedralType( "X:X:CW:H4",  180.00,   4.60240,  2)  # 
improper[ "X:X:CR:H5"] = CosineDihedralType( "X:X:CR:H5",  180.00,   4.60240,  2)  # 
improper[ "X:X:CV:H4"] = CosineDihedralType( "X:X:CV:H4",  180.00,   4.60240,  2)  # 
improper[ "X:X:CQ:H5"] = CosineDihedralType( "X:X:CQ:H5",  180.00,   4.60240,  2)  # 
improper[ "X:X:CK:H5"] = CosineDihedralType( "X:X:CK:H5",  180.00,   4.60240,  2)  # 
improper[ "X:X:CM:H4"] = CosineDihedralType( "X:X:CM:H4",  180.00,   4.60240,  2)  # 
improper[ "X:X:CM:HA"] = CosineDihedralType( "X:X:CM:HA",  180.00,   4.60240,  2)  # 
improper[ "X:X:CA:H4"] = CosineDihedralType( "X:X:CA:H4",  180.00,   4.60240,  2)  #  bsd.on C6H6 nmodes
improper[ "X:X:CA:H5"] = CosineDihedralType( "X:X:CA:H5",  180.00,   4.60240,  2)  #  bsd.on C6H6 nmodes

# @ffdef ffamber03 :proper NCosineDihedralType begin
#     ("CT", "CT", "OS", "CT",  [0.0, 180.0],  [1.60247, 0.41840],  [3, 2]),
#     ("CT", "CT", "OS", "CT",  [0.0, 180.0],  [1.60247, 0.41840],  [3, 2]),
#     ("CT", "CT", "OS", "CT",  [0.0, 180.0],  [1.60247, 0.41840],  [3, 2]),
#     ("CT", "CT", "OS", "CT",  v1,  v2,  [3, 2]),
#     ("CT", "CT", "OS", "CT",  v1,  v2,  [3, 2]),
#     ("CT", "CT", "OS", "CT",  v1,  v2,  [3, 2]),
# end


# proper["CT:CT:OS:CT"] = NCosineDihedralType("CT:CT:OS:CT",    0.0,   1.60247,  3)  # 
# proper["CT:CT:OS:CT"] = NCosineDihedralType("CT:CT:OS:CT",  180.0,   0.41840,  2)  # 
# proper[ "N:CT:C:N"] = NCosineDihedralType( "N:CT:C:N",  180.0,   2.86144,  1)  #  Amber03
# proper[ "N:CT:C:N"] = NCosineDihedralType( "N:CT:C:N",  180.0,   6.08228,  2)  #  Amber03
# proper[ "N:CT:C:N"] = NCosineDihedralType( "N:CT:C:N",  180.0,   1.93092,  3)  #  Amber03
# proper[ "C:N:CT:C"] = NCosineDihedralType( "C:N:CT:C",    0.0,   4.25053,  1)  #  Amber03
# proper[ "C:N:CT:C"] = NCosineDihedralType( "C:N:CT:C",  180.0,   1.44390,  2)  #  Amber03
# proper[ "C:N:CT:C"] = NCosineDihedralType( "C:N:CT:C",    0.0,   0.94517,  3)  #  Amber03
# proper["CT:CT:C:N"] = NCosineDihedralType("CT:CT:C:N",  180.0,   3.25683,  1)  #  Amber03
# proper["CT:CT:C:N"] = NCosineDihedralType("CT:CT:C:N",  180.0,   0.27489,  2)  #  Amber03
# proper["CT:CT:C:N"] = NCosineDihedralType("CT:CT:C:N",    0.0,   0.23430,  3)  #  Amber03
# proper[ "C:N:CT:CT"] = NCosineDihedralType( "C:N:CT:CT",  180.0,   1.47988,  1)  #  Amber03
# proper[ "C:N:CT:CT"] = NCosineDihedralType( "C:N:CT:CT",  180.0,   3.69698,  2)  #  Amber03
# proper[ "C:N:CT:CT"] = NCosineDihedralType( "C:N:CT:CT",  180.0,   0.94977,  3)  #  Amber03
# proper[ "C:N:CT:H0"] = NCosineDihedralType( "C:N:CT:H0",    0.0,   1.91418,  1)  #  Amber03 GLY
# proper[ "C:N:CT:H0"] = NCosineDihedralType( "C:N:CT:H0",  180.0,   5.25427,  2)  #  Amber03 GLY
# proper["H0:CT:C:N"] = NCosineDihedralType("H0:CT:C:N",  180.0,   4.43797,  1)  #  Amber03 GLY
# proper["H0:CT:C:N"] = NCosineDihedralType("H0:CT:C:N",    0.0,   0.04602,  2)  #  Amber03 GLY
# proper[ "H:N:C:O"] = NCosineDihedralType( "H:N:C:O",  180.0,  10.46000,  2)  #  JCC,7,(1986),230
# proper[ "H:N:C:O"] = NCosineDihedralType( "H:N:C:O",    0.0,   8.36800,  1)  #  J.C.cistrans-NMA DE
# proper["CT:S:S:CT"] = NCosineDihedralType("CT:S:S:CT",    0.0,  14.64400,  2)  #  JCC,7,(1986),230
# proper["CT:S:S:CT"] = NCosineDihedralType("CT:S:S:CT",    0.0,   2.51040,  3)  #  JCC,7,(1986),230
# proper["OS:CT:CT:OS"] = NCosineDihedralType("OS:CT:CT:OS",    0.0,   0.60250,  3)  #  parm98, TC,PC,PAK
# proper["OS:CT:CT:OS"] = NCosineDihedralType("OS:CT:CT:OS",    0.0,   4.91620,  2)  #  Piotr et al.
# proper["OS:CT:CT:OH"] = NCosineDihedralType("OS:CT:CT:OH",    0.0,   0.60250,  3)  #  parm98, TC,PC,PAK
# proper["OS:CT:CT:OH"] = NCosineDihedralType("OS:CT:CT:OH",    0.0,   4.91620,  2)  #  parm98, TC,PC,PAK
# proper["OH:CT:CT:OH"] = NCosineDihedralType("OH:CT:CT:OH",    0.0,   0.60250,  3)  #  parm98, TC,PC,PAK
# proper["OH:CT:CT:OH"] = NCosineDihedralType("OH:CT:CT:OH",    0.0,   4.91620,  2)  #  parm98, TC,PC,PAK
# proper["OH:P:OS:CT"] = NCosineDihedralType("OH:P:OS:CT",    0.0,   1.04600,  3)  #  JCC,7,(1986),230
# proper["OH:P:OS:CT"] = NCosineDihedralType("OH:P:OS:CT",    0.0,   5.02080,  2)  #  gg&gt ene.631g*/mp2
# proper["OS:P:OS:CT"] = NCosineDihedralType("OS:P:OS:CT",    0.0,   1.04600,  3)  #  JCC,7,(1986),230
# proper["OS:P:OS:CT"] = NCosineDihedralType("OS:P:OS:CT",    0.0,   5.02080,  2)  #  gg&gt ene.631g*/mp2
# proper["OS:CT:N*:CK"] = NCosineDihedralType("OS:CT:N*:CK",    0.0,  10.46000,  1)  #  parm98, TC,PC,PAK
# proper["OS:CT:N*:CM"] = NCosineDihedralType("OS:CT:N*:CM",    0.0,  10.46000,  1)  #  parm98, TC,PC,PAK
# proper["H1:CT:C:O"] = NCosineDihedralType("H1:CT:C:O",    0.0,   3.34720,  1)  #  Junmei et al, 1999
# proper["H1:CT:C:O"] = NCosineDihedralType("H1:CT:C:O",  180.0,   0.33472,  3)  #  Junmei et al, 1999
# proper["HC:CT:C:O"] = NCosineDihedralType("HC:CT:C:O",    0.0,   3.34720,  1)  #  Junmei et al, 1999
# proper["HC:CT:C:O"] = NCosineDihedralType("HC:CT:C:O",  180.0,   0.33472,  3)  #  Junmei et al, 1999
# proper["HC:CT:CT:HC"] = NCosineDihedralType("HC:CT:CT:HC",    0.0,   0.62760,  3)  #  Junmei et al, 1999
# proper["HC:CT:CT:CT"] = NCosineDihedralType("HC:CT:CT:CT",    0.0,   0.66944,  3)  #  Junmei et al, 1999
# proper["HC:CT:CM:CM"] = NCosineDihedralType("HC:CT:CM:CM",  180.0,   1.58992,  3)  #  Junmei et al, 1999
# proper["HC:CT:CM:CM"] = NCosineDihedralType("HC:CT:CM:CM",    0.0,   4.81160,  1)  #  Junmei et al, 1999
# proper["HO:OH:CT:CT"] = NCosineDihedralType("HO:OH:CT:CT",    0.0,   0.66944,  3)  #  Junmei et al, 1999
# proper["HO:OH:CT:CT"] = NCosineDihedralType("HO:OH:CT:CT",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["HO:OH:C:O"] = NCosineDihedralType("HO:OH:C:O",  180.0,   9.62320,  2)  #  Junmei et al, 1999
# proper["HO:OH:C:O"] = NCosineDihedralType("HO:OH:C:O",    0.0,   7.94960,  1)  #  Junmei et al, 1999
# proper["CM:CM:C:O"] = NCosineDihedralType("CM:CM:C:O",  180.0,   9.10020,  2)  #  Junmei et al, 1999
# proper["CM:CM:C:O"] = NCosineDihedralType("CM:CM:C:O",    0.0,   1.25520,  3)  #  Junmei et al, 1999
# proper["CT:CM:CM:CT"] = NCosineDihedralType("CT:CM:CM:CT",  180.0,  27.82360,  2)  #  Junmei et al, 1999
# proper["CT:CM:CM:CT"] = NCosineDihedralType("CT:CM:CM:CT",  180.0,   7.94960,  1)  #  Junmei et al, 1999
# proper["CT:CT:CT:CT"] = NCosineDihedralType("CT:CT:CT:CT",    0.0,   0.75312,  3)  #  Junmei et al, 1999
# proper["CT:CT:CT:CT"] = NCosineDihedralType("CT:CT:CT:CT",  180.0,   1.04600,  2)  #  Junmei et al, 1999
# proper["CT:CT:CT:CT"] = NCosineDihedralType("CT:CT:CT:CT",  180.0,   0.83680,  1)  #  Junmei et al, 1999
# proper["CT:CT:OS:C"] = NCosineDihedralType("CT:CT:OS:C",    0.0,   1.60247,  3)  #  Junmei et al, 1999
# proper["CT:CT:OS:C"] = NCosineDihedralType("CT:CT:OS:C",  180.0,   3.34720,  1)  #  Junmei et al, 1999
# proper["CT:OS:CT:OS"] = NCosineDihedralType("CT:OS:CT:OS",    0.0,   0.41840,  3)  #  Junmei et al, 1999
# proper["CT:OS:CT:OS"] = NCosineDihedralType("CT:OS:CT:OS",  180.0,   3.55640,  2)  #  Junmei et al, 1999
# proper["CT:OS:CT:OS"] = NCosineDihedralType("CT:OS:CT:OS",  180.0,   5.64840,  1)  #  Junmei et al, 1999
# proper[ "O:C:OS:CT"] = NCosineDihedralType( "O:C:OS:CT",  180.0,  11.29680,  2)  #  Junmei et al, 1999
# proper[ "O:C:OS:CT"] = NCosineDihedralType( "O:C:OS:CT",  180.0,   5.85760,  1)  #  Junmei et al, 1999
# proper[ "F:CT:CT:F"] = NCosineDihedralType( "F:CT:CT:F",  180.0,   5.02080,  1)  #  Junmei et al, 1999
# proper["Cl:CT:CT:Cl"] = NCosineDihedralType("Cl:CT:CT:Cl",  180.0,   1.88280,  1)  #  Junmei et al, 1999
# proper["Br:CT:CT:Br"] = NCosineDihedralType("Br:CT:CT:Br",    0.0,   0.00000,  0)  #  Junmei et al, 1999
# proper["H1:CT:CT:OS"] = NCosineDihedralType("H1:CT:CT:OS",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["H1:CT:CT:OH"] = NCosineDihedralType("H1:CT:CT:OH",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["H1:CT:CT:F"] = NCosineDihedralType("H1:CT:CT:F",    0.0,   0.79496,  1)  #  Junmei et al, 1999
# proper["H1:CT:CT:Cl"] = NCosineDihedralType("H1:CT:CT:Cl",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["H1:CT:CT:Br"] = NCosineDihedralType("H1:CT:CT:Br",    0.0,   2.30120,  1)  #  Junmei et al, 1999
# proper["HC:CT:CT:OS"] = NCosineDihedralType("HC:CT:CT:OS",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["HC:CT:CT:OH"] = NCosineDihedralType("HC:CT:CT:OH",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["HC:CT:CT:F"] = NCosineDihedralType("HC:CT:CT:F",    0.0,   0.79496,  1)  #  Junmei et al, 1999
# proper["HC:CT:CT:Cl"] = NCosineDihedralType("HC:CT:CT:Cl",    0.0,   1.04600,  1)  #  Junmei et al, 1999
# proper["HC:CT:CT:Br"] = NCosineDihedralType("HC:CT:CT:Br",    0.0,   2.30120,  1)  #  Junmei et al, 1999
# proper["CT:OS:CT:N*"] = NCosineDihedralType("CT:OS:CT:N*",    0.0,   1.60247,  3)  #  parm98.dat, TC,PC,PAK
# proper["CT:OS:CT:N*"] = NCosineDihedralType("CT:OS:CT:N*",    0.0,   2.71960,  2)  #  Piotr et al.
# proper[ "X:C:C:X"] = NCosineDihedralType( "X:C:C:X",  180.0,  15.16700,  2)  #  Junmei et al, 1999
# proper[ "X:C:O:X"] = NCosineDihedralType( "X:C:O:X",  180.0,  11.71520,  2)  #  Junmei et al, 1999
# proper[ "X:C:OS:X"] = NCosineDihedralType( "X:C:OS:X",  180.0,  11.29680,  2)  #  Junmei et al, 1999
# proper[ "X:CA:OH:X"] = NCosineDihedralType( "X:CA:OH:X",  180.0,   3.76560,  2)  #  Junmei et al, 99
# proper[ "X:CM:OS:X"] = NCosineDihedralType( "X:CM:OS:X",  180.0,   4.39320,  2)  #  Junmei et al, 1999
# proper[ "X:C:CA:X"] = NCosineDihedralType( "X:C:CA:X",  180.0,  15.16700,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:C:CB:X"] = NCosineDihedralType( "X:C:CB:X",  180.0,  12.55200,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:C:CM:X"] = NCosineDihedralType( "X:C:CM:X",  180.0,   9.10020,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:C:N*:X"] = NCosineDihedralType( "X:C:N*:X",  180.0,   6.06680,  2)  #  JCC,7,(1986),230
# proper[ "X:C:NA:X"] = NCosineDihedralType( "X:C:NA:X",  180.0,   5.64840,  2)  #  JCC,7,(1986),230
# proper[ "X:C:NC:X"] = NCosineDihedralType( "X:C:NC:X",  180.0,  16.73600,  2)  #  JCC,7,(1986),230
# proper[ "X:C:OH:X"] = NCosineDihedralType( "X:C:OH:X",  180.0,   9.62320,  2)  #  Junmei et al, 1999
# proper[ "X:C:CT:X"] = NCosineDihedralType( "X:C:CT:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CA:CA:X"] = NCosineDihedralType( "X:CA:CA:X",  180.0,  15.16700,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CA:CB:X"] = NCosineDihedralType( "X:CA:CB:X",  180.0,  14.64400,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CA:CM:X"] = NCosineDihedralType( "X:CA:CM:X",  180.0,  10.66920,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CA:CT:X"] = NCosineDihedralType( "X:CA:CT:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CA:N2:X"] = NCosineDihedralType( "X:CA:N2:X",  180.0,  10.04160,  2)  #  reinterpolated 93'
# proper[ "X:CA:NA:X"] = NCosineDihedralType( "X:CA:NA:X",  180.0,   6.27600,  2)  #  JCC,7,(1986),230
# proper[ "X:CA:NC:X"] = NCosineDihedralType( "X:CA:NC:X",  180.0,  20.08320,  2)  #  JCC,7,(1986),230
# proper[ "X:CB:CB:X"] = NCosineDihedralType( "X:CB:CB:X",  180.0,  22.80280,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CB:N*:X"] = NCosineDihedralType( "X:CB:N*:X",  180.0,   6.90360,  2)  #  JCC,7,(1986),230
# proper[ "X:CB:NB:X"] = NCosineDihedralType( "X:CB:NB:X",  180.0,  10.66920,  2)  #  JCC,7,(1986),230
# proper[ "X:CB:NC:X"] = NCosineDihedralType( "X:CB:NC:X",  180.0,  17.36360,  2)  #  JCC,7,(1986),230
# proper[ "X:CK:N*:X"] = NCosineDihedralType( "X:CK:N*:X",  180.0,   7.11280,  2)  #  JCC,7,(1986),230
# proper[ "X:CK:NB:X"] = NCosineDihedralType( "X:CK:NB:X",  180.0,  41.84000,  2)  #  JCC,7,(1986),230
# proper[ "X:CM:CM:X"] = NCosineDihedralType( "X:CM:CM:X",  180.0,  27.82360,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CM:CT:X"] = NCosineDihedralType( "X:CM:CT:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CM:N*:X"] = NCosineDihedralType( "X:CM:N*:X",  180.0,   7.74040,  2)  #  JCC,7,(1986),230
# proper[ "X:CQ:NC:X"] = NCosineDihedralType( "X:CQ:NC:X",  180.0,  28.45120,  2)  #  JCC,7,(1986),230
# proper[ "X:CT:CT:X"] = NCosineDihedralType( "X:CT:CT:X",    0.0,   0.65084,  3)  #  JCC,7,(1986),230
# proper[ "X:CT:N:X"] = NCosineDihedralType( "X:CT:N:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CT:N*:X"] = NCosineDihedralType( "X:CT:N*:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CT:N2:X"] = NCosineDihedralType( "X:CT:N2:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CT:OH:X"] = NCosineDihedralType( "X:CT:OH:X",    0.0,   0.69733,  3)  #  JCC,7,(1986),230
# proper[ "X:CT:OS:X"] = NCosineDihedralType( "X:CT:OS:X",    0.0,   1.60387,  3)  #  JCC,7,(1986),230
# proper[ "X:OH:P:X"] = NCosineDihedralType( "X:OH:P:X",    0.0,   1.04600,  3)  #  JCC,7,(1986),230
# proper[ "X:OS:P:X"] = NCosineDihedralType( "X:OS:P:X",    0.0,   1.04600,  3)  #  JCC,7,(1986),230
# proper[ "X:C:N:X"] = NCosineDihedralType( "X:C:N:X",  180.0,  10.46000,  2)  #  AA,NMA
# proper[ "X:CT:N3:X"] = NCosineDihedralType( "X:CT:N3:X",    0.0,   0.65084,  3)  #  JCC,7,(1986),230
# proper[ "X:CT:S:X"] = NCosineDihedralType( "X:CT:S:X",    0.0,   1.39467,  3)  #  JCC,7,(1986),230
# proper[ "X:CT:SH:X"] = NCosineDihedralType( "X:CT:SH:X",    0.0,   1.04600,  3)  #  JCC,7,(1986),230
# proper[ "X:CA:CN:X"] = NCosineDihedralType( "X:CA:CN:X",  180.0,  15.16700,  2)  #  reinterpolated 93'
# proper[ "X:CB:CN:X"] = NCosineDihedralType( "X:CB:CN:X",  180.0,  12.55200,  2)  #  reinterpolated 93'
# proper[ "X:CC:CT:X"] = NCosineDihedralType( "X:CC:CT:X",    0.0,   0.00000,  0)  #  JCC,7,(1986),230
# proper[ "X:CC:CV:X"] = NCosineDihedralType( "X:CC:CV:X",  180.0,  21.54760,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CC:CW:X"] = NCosineDihedralType( "X:CC:CW:X",  180.0,  22.48900,  2)  #  intrpol.bsd.on C6H6
# proper[ "X:CC:NA:X"] = NCosineDihedralType( "X:CC:NA:X",  180.0,   5.85760,  2)  #  JCC,7,(1986),230
# proper[ "X:CC:NB:X"] = NCosineDihedralType( "X:CC:NB:X",  180.0,  10.04160,  2)  #  JCC,7,(1986),230
# proper[ "X:CN:NA:X"] = NCosineDihedralType( "X:CN:NA:X",  180.0,   6.38060,  2)  #  reinterpolated 93'
# proper[ "X:CR:NA:X"] = NCosineDihedralType( "X:CR:NA:X",  180.0,   9.72780,  2)  #  JCC,7,(1986),230
# proper[ "X:CR:NB:X"] = NCosineDihedralType( "X:CR:NB:X",  180.0,  20.92000,  2)  #  JCC,7,(1986),230
# proper[ "X:CV:NB:X"] = NCosineDihedralType( "X:CV:NB:X",  180.0,  10.04160,  2)  #  JCC,7,(1986),230
# proper[ "X:CW:NA:X"] = NCosineDihedralType( "X:CW:NA:X",  180.0,   6.27600,  2)  #  JCC,7,(1986),230

ffamber03

end
