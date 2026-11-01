# Actual helper regressions: fabricated filesystem metadata only. Run solely
# inside a separately admitted compute preflight, never during preparation.
function run_full_identity_fixtures()
    checks=Dict{String,Any}[]
    function record(name,passed)
        push!(checks,Dict("name"=>name,"passed"=>passed))
        passed || error("Full-width identity metadata regression failed: "*name)
    end
    major=UInt64(0x1234);minor=UInt64(0x56789)
    device=UInt64(ccall(:gnu_dev_makedev,Culonglong,(Cuint,Cuint),major,minor))
    high=(UInt64(1)<<57)+UInt64(17)
    segment=Dict("device"=>string(major;base=16)*":"*string(minor;base=16),"inode"=>string(high))
    identity=Dict("device_decimal"=>string(device),"inode_decimal"=>string(high),"size_decimal"=>"42")
    record("Full2pow57plus17inode accepted by actual helper",mapped_identity_matches(segment,identity))
    low=copy(identity);low["inode_decimal"]="17"
    record("Low32alias17 REJECTED against same full mapped inode",!mapped_identity_matches(segment,low))
    maximal=copy(segment);maximal["inode"]=string(typemax(UInt64))
    maxid=copy(identity);maxid["inode_decimal"]=string(typemax(UInt64))
    wrong=copy(maxid);wrong["inode_decimal"]=string(typemax(UInt64)-UInt64(1))
    record("UInt64maximum retained; distinct near-maximum inode rejected",
        mapped_identity_matches(maximal,maxid) && !mapped_identity_matches(maximal,wrong))
    wrongdev=copy(identity);wrongdev["device_decimal"]=string(device ⊻ UInt64(1))
    record("Wrong full device rejected",!mapped_identity_matches(segment,wrongdev))

    # Exercise the actual pure parser, not a duplicate implementation. Its
    # output attaches records to the original DIRECT-ARGV request order.
    paths=["/metadata/first","/metadata/second"]
    raw=string(device)*"|"*string(high)*"|42\n"*string(device)*"|"*string(typemax(UInt64))*"|0\n"
    parsed,rows=parse_stat_output(raw,paths)
    record("Actual parser preserves full UInt64 and exact request order",
        length(rows)==2 && rows[1]["path"]==paths[1] && rows[2]["path"]==paths[2] &&
        parsed[paths[1]]["inode_decimal"]==string(high) &&
        parse(UInt64,parsed[paths[2]]["inode_decimal"])==typemax(UInt64))
    function rejects(text,request)
        try
            parse_stat_output(text,request)
            return false
        catch
            return true
        end
    end
    record("Actual parser rejects missing and extra records",
        rejects("1|2|3\n",paths) && rejects("1|2|3\n4|5|6\n7|8|9\n",paths))
    invalid=["1|-2|3\n","1|2.0|3\n","1|18446744073709551616|3\n","1|2|NaN\n"]
    record("Actual parser rejects negative/floating/overflow/noninteger metadata",
        all(text->rejects(text,[paths[1]]),invalid))
    record("Actual parser rejects empty and truncated output",
        rejects("",[paths[1]]) && rejects("1|2|3",[paths[1]]))
    return checks
end
