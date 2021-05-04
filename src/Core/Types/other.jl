"""
    ResidueName(s::String) <: AbstractString

A string overload to accomodate different known and expected denominations of
certain residues ("HIS" == "HIE" for example).

# Example
```jldoctest
julia> residue.name = ResidueName("HIS")
```
"""
mutable struct ResidueName <: AbstractString
    content::String
end

Base.show(io::IO, rn::ResidueName) = begin
    print(rn.content)
end

Base.:(==)(s1::String, s2::ResidueName)      = s2 == s1
Base.:(==)(s1::ResidueName, s2::ResidueName) = s1 == s2.content
Base.:(==)(s1::ResidueName, s2::String)      = begin
    if s1.content in ["HIE", "HIS"] && s2 in ["HIE", "HIS"]
        return true
    else
        return s1.content == s2
    end
end

Base.iterate(rn::ResidueName) = Base.iterate(rn.content)
Base.iterate(rn::ResidueName, n::Int64) = Base.iterate(rn.content, n)

# * ---- XML-RPC ---------------------------------------------------------------

module XMLRPC

    using HTTP

    export ClientProxy

    """
        ClientProxy(host::String, port::Int)

    Instantiate a client proxy listening on `host`and  `port`. This object is
    meant to act as a proxy to a remote service implementing the XML/RPC
    protocol, and was designed to work with Python and PyMOL. For more detailed
    information, please read:

    https://blog.papercut.com/write-xml-rpc-clients/

    https://wiki.python.org/moin/XmlRpc

    # Examples
    ```jldoctest
    julia> ProtoSyn.XMLRPC.ClientProxy("http://localhost", 50000)
    ProtoSyn.XMLRPC.ClientProxy("http://localhost", 50000, "http://localhost:50000")
    ```
    """
    struct ClientProxy
        host::String
        port::Int
        url::String
        ClientProxy(host::String, port::Int) = begin
            p = new(host, port, "$host:$port")
        end
    end # struct

    Base.getproperty(p::ClientProxy, k::Symbol) = begin
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

end # module