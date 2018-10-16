# inspired by from https://github.com/JuliaAudio/SampledSignals.jl/blob/master/src/WAVDisplay.jl

using .Common
using Random
using .Print

function embed_javascript()
    # js_path = joinpath(dirname(dirname(@__FILE__)), "deps", "ChemDoodleWeb.js")
    # js_text = open(js_path) do io
    #     read(io, String)
    # end
    # # the javascript file contains the code to add itself to the require module
    # # cache under the name 'wavesurfer'
    # display("text/html", """
    # <script charset="utf-8" type='text/javascript'>
    # $js_text
    # console.log("Loaded ChemDoodleWeb");
    # let CHEMDOOD_WIDGET_WIDTH=600;
    # let CHEMDOOD_WIDGET_HEIGHT=400;
    # function set_chemdoodle_size(w,h) {
    #     CHEMDOOD_WIDGET_WIDTH = w;
    #     CHEMDOOD_WIDGET_HEIGHT = h;
    # }
    # </script>
    # """)
    display("text/html", """
    <script src="https://hub.chemdoodle.com/cwc/latest/ChemDoodleWeb.js"></script>
    <script charset="utf-8" type='text/javascript'>
    let CHEMDOOD_WIDGET_WIDTH=600;
    let CHEMDOOD_WIDGET_HEIGHT=400;
    function set_chemdoodle_size(w,h) {
        CHEMDOOD_WIDGET_WIDTH = w;
        CHEMDOOD_WIDGET_HEIGHT = h;
    }
    </script>
    """)
end

function Base.show(io::IO, ::MIME"text/html", state::Common.State)
    rndstr = randstring()
    divid = string("div_", rndstr)
    canvasid = string("canvas_", rndstr)

    # println(io, """
    # <div id="$divid">
    #     <h4>ProtoSyn display requires javascript</h4>
    #     <p>To enable for the whole notebook select "Trust Notebook" from the
    #     "File" menu. You can also trust this cell by re-running it. You may
    #     also need to re-run `using SampledSignals` if the module is not yet
    #     loaded in the Julia kernel, or `ProtoSyn.embed_javascript()`
    #     if the Julia module is loaded but the javascript isn't initialized.</p>
    #     <canvas id="$canvasid">
    #         Your browser doesn't support HTML5 canvas.
    #     </canvas>
    # </div>""")
    
    iobuffer = IOBuffer()
    Print.as_xyz(state, ostream=iobuffer)
    xyz = replace(String(take!(iobuffer)), "\n" => "\\n")

    println(io, """
    <div id="$divid">
    <canvas id="$canvasid">Your browser doesn't support HTML5 canvas.</canvas>
    <script>
        var transformBallAndStick = new ChemDoodle.TransformCanvas3D('$canvasid', CHEMDOOD_WIDGET_WIDTH, CHEMDOOD_WIDGET_HEIGHT);
        transformBallAndStick.specs.set3DRepresentation('Ball and Stick');
        //transformBallAndStick.specs.backgroundColor = 'black';
        //var molFile = '3036\\n  CHEMDOOD12280913053D\\n\\n 28 29  0     0  0  0  0  0  0999 V2000\\n    0.0456    1.0544   -1.9374 Cl  0  0  0  0  0  0  0  0  0  0  0  0\\n   -0.7952   -1.7026   -1.7706 Cl  0  0  0  0  0  0  0  0  0  0  0  0\\n    0.6447   -0.8006   -4.1065 Cl  0  0  0  0  0  0  0  0  0  0  0  0\\n    1.8316   -0.9435    4.4004 Cl  0  0  0  0  0  0  0  0  0  0  0  0\\n    6.9949    1.1239   -3.9007 Cl  0  0  0  0  0  0  0  0  0  0  0  0\\n    1.9032   -1.0692   -1.6001 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.8846   -1.0376   -0.1090 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    3.2176   -0.5035   -2.1949 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    0.5585   -0.6223   -2.3126 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    2.2670    0.1198    0.5688 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    4.3480   -1.2638   -2.0859 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.4856   -2.1660    0.6075 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    3.1719    0.7242   -2.7939 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    2.2506    0.1490    1.9633 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    5.5313   -0.7541   -2.6203 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.4691   -2.1369    2.0020 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    4.3552    1.2340   -3.3284 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.8515   -0.9793    2.6800 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    5.5350    0.4948   -3.2417 C   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.9777   -2.1366   -1.8749 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    2.5727    1.0177    0.0401 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    4.3513   -2.2356   -1.6034 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.1951   -3.0814    0.0991 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    2.3077    1.3562   -2.8879 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    2.5491    1.0585    2.4783 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    6.4431   -1.3411   -2.5451 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    1.1584   -3.0244    2.5473 H   0  0  0  0  0  0  0  0  0  0  0  0\\n    4.3449    2.2098   -3.8075 H   0  0  0  0  0  0  0  0  0  0  0  0\\n  1  9  1  0  0  0  0\\n  2  9  1  0  0  0  0\\n  3  9  1  0  0  0  0\\n  4 18  1  0  0  0  0\\n  5 19  1  0  0  0  0\\n  6  7  1  0  0  0  0\\n  6  8  1  0  0  0  0\\n  6  9  1  0  0  0  0\\n  6 20  1  0  0  0  0\\n  7 10  2  0  0  0  0\\n  7 12  1  0  0  0  0\\n  8 11  2  0  0  0  0\\n  8 13  1  0  0  0  0\\n 10 14  1  0  0  0  0\\n 10 21  1  0  0  0  0\\n 11 15  1  0  0  0  0\\n 11 22  1  0  0  0  0\\n 12 16  2  0  0  0  0\\n 12 23  1  0  0  0  0\\n 13 17  2  0  0  0  0\\n 13 24  1  0  0  0  0\\n 14 18  2  0  0  0  0\\n 14 25  1  0  0  0  0\\n 15 19  2  0  0  0  0\\n 15 26  1  0  0  0  0\\n 16 18  1  0  0  0  0\\n 16 27  1  0  0  0  0\\n 17 19  1  0  0  0  0\\n 17 28  1  0  0  0  0\\nM  END\\n';
        //var molecule = ChemDoodle.readMOL(molFile, 1);
        var xyzFile = '$xyz';
        var molecule = ChemDoodle.readXYZ(xyzFile, 1);
        transformBallAndStick.loadMolecule(molecule);
    </script>
    </div>""")

end