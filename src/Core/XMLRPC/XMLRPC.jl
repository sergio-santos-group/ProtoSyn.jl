module XMLRPC

using HTTP

export ServerProxy

"""
    ServerProxy(host::String, port::Int)

Instantiate a server proxy listening on `host`and  `port`. This
object is meant to act as a proxy to a remote service implementing
the XML/RPC protocol, and was designed to work with PyMOL.
"""
struct ServerProxy
    host::String
    port::Int
    url::String
    ServerProxy(host::String, port::Int) = begin
        p = new(host, port, "$host:$port")
    end
end

Base.getproperty(p::ServerProxy, k::Symbol) = begin
    if k in (:host, :port, :url)
        return getfield(p, k)
    end
    url = getfield(p, :url)
    (args...) -> rpccall(url, string(k), args...)
end


function rpccall(url::String, methodname::String, args...)
    headers = Dict("Content-type" => "text/xml")
    body = xml(methodname, args...)
    r = HTTP.post(url, headers, body)
end


function xml(methodname::String, args...)
    io = IOBuffer()
    write(io, "<?xml version=\"1.0\"?>\n")
    write(io, "<methodCall>\n")
    write(io, "<methodName>", methodname, "</methodName>\n")
    write(io, "<params>\n")
    for arg in args
        write(io, "<param>\n")
        xml(io, arg)
        write(io, "</param>\n")
    end
    write(io, "</params>\n")
    write(io, "</methodCall>\n")
    body = String(take!(io))
    close(io)
    return body
end

xml(io::IO, x::Int) = write(io, "<value><int>$(x)</int></value>\n")
xml(io::IO, x::Bool) = write(io, "<value><boolean>$(x ? 1 : 0)</boolean></value>\n")
xml(io::IO, x::String) = write(io, "<value><string>$(x)</string></value>\n")
xml(io::IO, x::T) where {T<:AbstractFloat} = write(io, "<value><double>$(Float64(x))</double></value>\n")

xml(io::IO, d::Dict) = begin
    write(io, "<value><struct>\n")
    foreach(x->xml(io, x), d)
    write(io, "</struct></value>\n")
end

xml(io::IO, p::Pair) = begin
    write(io, "<member>")
    write(io, "<name>", string(p.first), "</name>\n")
    xml(io, p.second)
    write(io, "</member>\n")
end

xml(io::IO, v::Vector) = begin
    write(io, "<value><array><data>\n")
    foreach(x->xml(io, x), v)
    write(io, "</data></array></value>\n")
end

xml(io::IO, ar::Array{T,2}) where {T} = begin
    write(io, "<value><array><data>\n")
    foreach(row->xml(io, Vector(row)), eachrow(ar))
    write(io, "</data></array></value>\n")
end

xml(x) = xml(stdout, x)

# proxy = ServerProxy("http://localhost", 9123)

# # proxy.ola("eu")
# # proxy.adeus("tu",1,[2,4,5])
# # proxy.adeus("tu",1,Dict(10=>20))
# println("HOST: ", proxy.host)
# println("PORT: ", proxy.port)
# println("URL: ", proxy.url)
# # proxy.delete("all")
# # proxy.fetch("1ctf")
# # proxy.show("lines", "1ctf")

# xml(stdout, Array(reshape(1:12,3,4)'))

# https://blog.papercut.com/write-xml-rpc-clients/
# https://wiki.python.org/moin/XmlRpc


end