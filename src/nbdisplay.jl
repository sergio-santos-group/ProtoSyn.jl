# inspired by https://github.com/JuliaAudio/SampledSignals.jl/blob/master/src/WAVDisplay.jl

using Random

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
    # let CWC_WIDTH=600;
    # let CWC_HEIGHT=400;
    # function set_chemdoodle_size(w,h) {
    #     CWC_WIDTH = w;
    #     CWC_HEIGHT = h;
    # }
    # </script>
    # """)
    display("text/html", """
    <script src="https://hub.chemdoodle.com/cwc/latest/ChemDoodleWeb.js"></script>
    <script charset="utf-8" type='text/javascript'>
    window.CWC_WIDTH=600;
    window.CWC_HEIGHT=400;
    function set_chemdoodle_size(w,h) {
        window.CWC_WIDTH = w;
        window.CWC_HEIGHT = h;
    }
    </script>
    """)
end

function Base.show(io::IO, ::MIME"text/html", molstate::Tuple{Molecule, State})
    rndstr = randstring()
    divid = string("div_", rndstr)
    canvasid = string("canvas_", rndstr)
    cwcid = string("cwc_", rndstr)  # ChemDoodleWeb component ID

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

    mol, state = molstate
    iobuffer = IOBuffer()
    write(iobuffer, mol, state, XYZ)
    mol = replace(String(take!(iobuffer)), "\n" => "\\n")
    # println(mol)
    close(iobuffer)

    # var w = document.getElementById('$divid').offsetWidth
    # var $cwcid = new ChemDoodle.TransformCanvas3D('$canvasid', w, w/1.5);
    # $(cwcid).specs.backgroundColor = 'black';
        
    println(io, """
    <div id="$divid">
    <canvas id="$canvasid">Your browser doesn't support HTML5 canvas.</canvas>
    <script>
        var $cwcid = new ChemDoodle.TransformCanvas3D('$canvasid', CWC_WIDTH, CWC_HEIGHT);
        $(cwcid).specs.set3DRepresentation('Ball and Stick');
        var molfile = '$mol';
        var molecule = ChemDoodle.readXYZ(molfile, 1);
        $(cwcid).loadMolecule(molecule);
    </script>
    </div>""")

end