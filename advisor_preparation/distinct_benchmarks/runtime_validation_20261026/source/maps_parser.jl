module StrictProcMaps
export parse_proc_maps, read_proc_maps_complete

"""Read the proc stream line by line until EOF; pseudo-file stat size is unused."""
function read_proc_maps_complete(path::AbstractString)
    lines=String[]
    open(path,"r") do io
        while !eof(io)
            push!(lines,readline(io;keep=true))
        end
    end
    return join(lines)
end

"""Strict syntax only. Identity, provider and completion gates are separate."""
function parse_proc_maps(raw::AbstractString)
    !isempty(raw) && endswith(raw,"\n") || error("Empty or truncated process maps")
    rows=Dict{String,String}[]
    pattern=r"^([0-9A-Fa-f]+)-([0-9A-Fa-f]+)\s+([r-][w-][x-][ps])\s+([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+:[0-9A-Fa-f]+)\s+([0-9]+)(?:\s+(.*))?$"
    for (number,line) in enumerate(split(chop(raw;tail=1),'\n';keepempty=true))
        m=match(pattern,line)
        m!==nothing || error("Malformed process map line "*string(number))
        first,last,perms,offset,device,inode,name=m.captures
        parse(UInt64,first;base=16)<parse(UInt64,last;base=16) || error("Nonpositive mapped address span")
        name=name===nothing ? "" : strip(name)
        endswith(name," (deleted)") && error("Deleted mapped provider")
        isempty(name) || startswith(name,"/") || (startswith(name,"[") && endswith(name,"]")) || error("Unrecognized non-file mapping name")
        push!(rows,Dict("address"=>first*"-"*last,"permissions"=>perms,
            "offset"=>offset,"device"=>device,"inode"=>inode,"path"=>name))
    end
    !isempty(rows) && any(startswith(v["path"],"/") for v in rows) || error("Zero file-backed mappings")
    return rows
end
end
