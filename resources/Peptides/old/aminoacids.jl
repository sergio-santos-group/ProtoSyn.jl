let

using .ProtoSyn: Atom, Residue, ResidueLib, ConnectGraph, ConnectGraphByName

reslib = ResidueLib()

# ==============================================================================
reslib["BKB"] = Residue(
    "BKB",
    [
        Atom(1,  "H", "H", 0.3909, 0.0724, 0.0000),
        Atom(2,  "N", "N", 0.3326, 0.1548, 0.0000),
        Atom(3, "CA", "C", 0.3970, 0.2846, 0.0000),
        Atom(4,  "C", "C", 0.5486, 0.2705, 0.0000),
        Atom(5,  "O", "O", 0.6009, 0.1593, 0.0000),
        
    ],
    ConnectGraph(
        1 => [2],
        2 => [1, 3],
        3 => [2, 4],
        4 => [3, 5],
        5 => [4]
    )
)

# ==============================================================================
reslib["BKX"] = Residue(
    "BKX",
    [
        Atom(1,  "H", "H", 0.3909, 0.0724, 0.0000),
        Atom(1,  "N", "N", 0.3326, 0.1548, 0.0000),
        Atom(1, "CA", "C", 0.3970, 0.2846, 0.0000),
        Atom(1,  "C", "C", 0.5486, 0.2705, 0.0000),
        Atom(1,  "O", "O", 0.6009, 0.1593, 0.0000)
    ],
    ConnectGraphByName(
        "H"  => ["N"],
        "N"  => ["H", "CA"],
        "CA" => ["N", "C"],
        "C"  => ["CA", "O"],
        "O"  => ["C"]
    )
)


# # ==============================================================================
# reslib["ALA"] = Residue(
#     "ALA",
#     [
#         Atom(id= 1, name=  "N", symbol="N",  x=0.3326, y=0.1548, z= 0.0000),
#         Atom(id= 2, name=  "H", symbol="H",  x=0.3909, y=0.0724, z= 0.0000),
#         Atom(id= 3, name= "CA", symbol="C",  x=0.3970, y=0.2846, z= 0.0000),
#         Atom(id= 4, name= "HA", symbol="H",  x=0.3672, y=0.3400, z=-0.0890),
#         Atom(id= 5, name= "CB", symbol="C",  x=0.3577, y=0.3654, z= 0.1232),
#         Atom(id= 6, name="HB1", symbol="H",  x=0.3877, y=0.3116, z= 0.2131),
#         Atom(id= 7, name="HB2", symbol="H",  x=0.4075, y=0.4623, z= 0.1206),
#         Atom(id= 8, name="HB3", symbol="H",  x=0.2497, y=0.3801, z= 0.1241),
#         Atom(id= 9, name=  "C", symbol="C",  x=0.5486, y=0.2705, z= 0.0000),
#         Atom(id=10, name=  "O", symbol="O",  x=0.6009, y=0.1593, z= 0.0000),
#     ],
#     ConnectGraph(
#         1 => [2, 3],
#         2 => [1],
#         3 => [1, 4, 5, 9],
#         4 => [3],
#         5 => [3, 6, 7, 8],
#         6 => [5],
#         7 => [5],
#         8 => [5],
#         9 => [3, 10],
#        10 => [9]
#     )
# )

# # residues["BALA"] = Residue("BKB", deepcopy(residues["BKB"].atoms), residues["BKB"].bonds)

reslib
end