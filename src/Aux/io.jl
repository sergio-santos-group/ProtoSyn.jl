#TO DO:
# 1) DOCUMENT THE FUNCTION
function read_JSON(i_file::String)::Dict

    open(i_file, "r") do f
        json_txt = read(f, String)
        return JSON.parse(json_txt)
    end
end