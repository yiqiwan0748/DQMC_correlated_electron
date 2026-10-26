module FullFileIdentity
using TOML, SHA
export verify_stat_utility, stat_batch, mapped_identity_matches, parse_stat_output

file_sha(path)=bytes2hex(sha256(read(path)))

function verify_stat_utility(utility)
    utility["schema"]=="P17_PINNED_GNU_STAT_IDENTITY_UTILITY_V3" || error("Stat utility schema")
    utility["path"]=="/usr/bin/stat" && utility["realpath"]=="/usr/bin/stat" || error("Stat utility path")
    realpath(utility["path"])==utility["realpath"] || error("Stat utility realpath")
    filesize(utility["path"])==utility["bytes"] && file_sha(utility["path"])==utility["sha256"] || error("Stat utility bytes")
    nothing
end

"""Pure parser for exact GNU output and request order; no file/library work."""
function parse_stat_output(raw::AbstractString,paths::Vector{String})
    !isempty(paths) && length(unique(paths))==length(paths) || error("Nonempty unique parser paths required")
    !isempty(raw) && endswith(raw,"\n") || error("Empty/truncated GNU stat output")
    lines=split(chop(raw;tail=1),'\n';keepempty=true)
    length(lines)==length(paths) || error("GNU stat output record count/order mismatch")
    result=Dict{String,Dict{String,Any}}()
    rows=Dict{String,Any}[]
    for (index,(path,line)) in enumerate(zip(paths,lines))
        occursin(r"^[0-9]+\|[0-9]+\|[0-9]+$",line) || error("Malformed/noninteger GNU stat output")
        parts=split(line,'|')
        device,inode,size=parse.(UInt64,parts)
        row=Dict{String,Any}("request_index"=>index,"path"=>path,
            "device_decimal"=>String(parts[1]),"inode_decimal"=>String(parts[2]),
            "size_decimal"=>String(parts[3]),"device_hex"=>string(device;base=16),
            "full_width_decimal_parse_verified"=>true)
        push!(rows,row)
        result[path]=row
    end
    return result,rows
end

"""One pinned GNU-stat batch, direct argv, exact order/count, decimal UInt64.

File identities are not obtained from Julia's narrowed jl_stat_dev/ino API.
All raw utility output and the command/exit/utility identity receipt are saved
before rejection. This module calls no numerical or scientific API.
"""
function stat_batch(paths::Vector{String},utility,output_prefix::String)
    !isempty(paths) && length(unique(paths))==length(paths) || error("Nonempty unique stat paths required")
    all(isabspath,paths) && all(realpath(path)==path for path in paths) || error("Exact absolute realpaths required")
    outpath=output_prefix*".stdout_PRIVATE.txt"
    errpath=output_prefix*".stderr_PRIVATE.txt"
    receiptpath=output_prefix*".receipt_PRIVATE.toml"
    any(ispath,(outpath,errpath,receiptpath)) && error("Fresh stat receipt paths required")
    argv=vcat([utility["path"],"--dereference","--printf=%d|%i|%s\n","--"],paths)
    receipt=Dict{String,Any}("schema"=>"P17_FULL_WIDTH_STAT_BATCH_V3","status"=>"UNRESOLVED",
        "argv"=>argv,"utility_sha256"=>utility["sha256"],"paths_in_request_order"=>paths,
        "expected_record_count"=>length(paths),"integer_encoding"=>"unsigned decimal strings parsed as UInt64; no floating conversion",
        "rows"=>Any[])
    result=Dict{String,Dict{String,Any}}()
    try
        verify_stat_utility(utility)
        receipt["utility_before_verified"]=true
        process=open(outpath,"w") do stdout
            open(errpath,"w") do stderr
                run(pipeline(ignorestatus(Cmd(argv));stdout=stdout,stderr=stderr))
            end
        end
        receipt["exit_code"]=process.exitcode
        verify_stat_utility(utility)
        receipt["utility_after_verified"]=true
        receipt["stdout_sha256"]=file_sha(outpath)
        receipt["stderr_sha256"]=file_sha(errpath)
        process.exitcode==0 || error("GNU stat batch failed")
        raw=read(outpath,String)
        result,rows=parse_stat_output(raw,paths)
        receipt["rows"]=rows
        receipt["status"]="PASS_METADATA_ONLY"
    catch e
        receipt["error"]=sprint(showerror,e)
    finally
        open(receiptpath,"w") do io;TOML.print(io,receipt);end
    end
    receipt["status"]=="PASS_METADATA_ONLY" || error(get(receipt,"error","GNU stat unresolved"))
    return result,receipt
end

function mapped_identity_matches(segment,identity)
    device=parse(UInt64,identity["device_decimal"])
    inode=parse(UInt64,identity["inode_decimal"])
    map_major,map_minor=parse.(UInt64,split(segment["device"],':');base=16)
    # These glibc functions return the decoded major/minor. Their INPUT is the
    # full UInt64 GNU/kernel device value, never Julia's narrowed stat.device.
    major=UInt64(ccall(:gnu_dev_major,Cuint,(Culonglong,),device))
    minor=UInt64(ccall(:gnu_dev_minor,Cuint,(Culonglong,),device))
    return major==map_major && minor==map_minor && inode==parse(UInt64,segment["inode"])
end
end
