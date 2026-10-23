#!/usr/bin/env julia
"""Independent returned-CSV arithmetic audit, without loading ITensors or MPS objects."""
module ALMDMRGReturnAudit20260911
using LinearAlgebra, Serialization, SHA, TOML, Dates
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const DEFAULT_ATTEMPT=joinpath(ROOT,"results/literature_transfer/alm_dmrg_benchmark_20260911/attempt01")
const DEFAULT_REFERENCE=joinpath(ROOT,"results/literature_transfer/alm_exact_locking_20260911/dmrg_reference_target.jls")
const DEFAULT_OUTPUT=joinpath(ROOT,"results/literature_transfer/alm_method_benchmark_20260911/dmrg_return_audit.toml")
const SEEDS=(314159,271828)
const CAPS=(32,64,128,256,256)
const COMPARE_ATOL=5e-12 # Roundtrip Float64 CSV and 8x8 reductions; not an accuracy gate.
const ENERGY_ATOL=1e-8
const VARIANCE_ATOL=1e-8
const SPIN_ATOL=1e-6
const NUMBER_ATOL=1e-10
check(condition,message)=condition || error(message)
sha(path)=bytes2hex(sha256(read(path)))
function csv_fields(line)
    chars=collect(chomp(line));fields=String[];buffer=IOBuffer();quoted=false;i=1
    while i<=length(chars)
        c=chars[i]
        if c=='"'
            if quoted && i<length(chars) && chars[i+1]=='"'
                print(buffer,'"');i+=1
            else
                quoted=!quoted
            end
        elseif c==',' && !quoted
            push!(fields,String(take!(buffer)))
        elseif c!='\r'
            print(buffer,c)
        end
        i+=1
    end
    check(!quoted,"Unterminated CSV quoted field")
    push!(fields,String(take!(buffer)));fields
end
function read_csv(path)
    lines=readlines(path);check(!isempty(lines),"Empty CSV: $path")
    header=csv_fields(first(lines));check(length(unique(header))==length(header),"Duplicate CSV header")
    rows=Dict{String,String}[]
    for (i,line) in enumerate(lines[2:end])
        fields=csv_fields(line);check(length(fields)==length(header),"Malformed row $(i+1) in $path")
        push!(rows,Dict(zip(header,fields)))
    end
    rows
end
number(row,key)=parse(Float64,row[key])
integer(row,key)=parse(Int,row[key])
function boolean(row,key)
    row[key]=="true" && return true
    row[key]=="false" && return false
    error("Expected Boolean $key, got $(row[key])")
end
stagekey(row)=(integer(row,"Sz"),integer(row,"seed"),integer(row,"stage"))
function compare(actual,expected,label,residuals)
    check(isfinite(actual) && isfinite(expected),"Nonfinite compared value: $label")
    difference=abs(actual-expected);residuals[label]=max(get(residuals,label,0.),difference)
    check(difference<=COMPARE_ATOL,"$label mismatch: recomputed=$actual exported=$expected")
end
function grouped(rows)
    groups=Dict{Tuple{Int,Int,Int},Vector{Dict{String,String}}}()
    for row in rows;push!(get!(groups,stagekey(row),Dict{String,String}[]),row);end
    groups
end
function write_report(path,payload)
    check(!ispath(path),"Refusing audit overwrite: $path");mkpath(dirname(path))
    temporary=path*".partial";check(!ispath(temporary),"Existing partial audit: $temporary")
    open(temporary,"w") do io;TOML.print(io,payload);end
    mv(temporary,path)
end

function compute_audit(attempt_dir,reference_path)
    files=["stages.csv","spin_matrices.csv","local_observables.csv","initialization_agreement.csv",
        "summary.toml","metadata.toml","status.txt","run_benchmark.jl","Project.toml","Manifest.toml"]
    for file in files;check(isfile(joinpath(attempt_dir,file)),"Missing returned artifact $file");end
    hashes=Dict(file=>sha(joinpath(attempt_dir,file)) for file in files)
    metadata=TOML.parsefile(joinpath(attempt_dir,"metadata.toml"));summary=TOML.parsefile(joinpath(attempt_dir,"summary.toml"))
    check(sha(reference_path)==metadata["reference_sha256"],"Original ED reference hash mismatch")
    check(hashes["run_benchmark.jl"]==metadata["runner_sha256"],"Returned runner hash mismatch")
    for file in ("Project.toml","Manifest.toml")
        check(hashes[file]==metadata[file*"_sha256"],"Returned environment hash mismatch")
    end
    check(metadata["seeds"]==collect(SEEDS) && metadata["caps"]==collect(CAPS[1:4]),"Registered seeds/caps mismatch")
    check(metadata["physical_Sz"]==[0,1] && metadata["mode"]=="run","Registered sector/mode mismatch")
    for (key,value) in (("energy_tolerance",ENERGY_ATOL),("variance_tolerance",VARIANCE_ATOL),
        ("spin_entry_tolerance",SPIN_ATOL),("particle_Sz_tolerance",NUMBER_ATOL),
        ("sweeps_per_cap",6),("extra_zero_noise_sweeps",6),("cutoff",1e-12),
        ("eigsolve_tol",1e-11),("eigsolve_krylovdim",10),("eigsolve_maxiter",8),
        ("julia_threads",1),("blas_threads",1))
        check(metadata[key]==value,"Registered metadata mismatch: $key")
    end
    check(metadata["ITensors_version"]=="0.9.30" && metadata["ITensorMPS_version"]=="0.4.1","Package pins mismatch")
    ref=deserialize(reference_path)
    check(ref.N==6 && ref.U==20. && ref.t4==-.65 && ref.cells==4,"Reference model mismatch")
    check(Set(keys(ref.sectors))==Set((0,1)),"Reference sector coverage mismatch")
    stages=read_csv(joinpath(attempt_dir,"stages.csv"))
    spins=read_csv(joinpath(attempt_dir,"spin_matrices.csv"))
    localrows=read_csv(joinpath(attempt_dir,"local_observables.csv"))
    check(length(stages)==20 && length(spins)==1280 && length(localrows)==160,"Incomplete or extra raw rows")
    expected=Set((M,seed,stage) for M in 0:1 for seed in SEEDS for stage in 1:5)
    check(Set(stagekey(r) for r in stages)==expected && length(unique(stagekey.(stages)))==20,"Repeated or missing stages")
    spin_groups=grouped(spins);local_groups=grouped(localrows)
    check(Set(keys(spin_groups))==expected && Set(keys(local_groups))==expected,"Raw table stage coverage mismatch")
    residuals=Dict{String,Float64}();audit_rows=Dict{String,Any}[];recomputed=Dict{Tuple{Int,Int,Int},Any}()
    odd=[isodd(p) ? 1. : -1. for p in 1:8]
    for row in sort(stages;by=stagekey)
        key=stagekey(row);M,seed,stage=key;cap=CAPS[stage];sec=ref.sectors[M]
        check(integer(row,"cap")==cap,"Cap mismatch at $key")
        matrixrows=spin_groups[key];occupationrows=local_groups[key]
        check(length(matrixrows)==64 && length(occupationrows)==8,"Raw matrix/vector length mismatch at $key")
        C=zeros(8,8);Z=zeros(8,8);seen=Set{Tuple{Int,Int}}()
        for r in matrixrows
            p,q=integer(r,"p"),integer(r,"q")
            check(1<=p<=8 && 1<=q<=8 && !((p,q) in seen),"Duplicate/invalid spin matrix index")
            check(integer(r,"cap")==cap,"Spin-matrix cap mismatch")
            push!(seen,(p,q));C[p,q]=number(r,"spin_dot");Z[p,q]=number(r,"SzSz")
        end
        up=zeros(8);down=zeros(8);density=zeros(8);doubles=zeros(8);sz=zeros(8);modes=Set{Int}()
        for r in occupationrows
            p=integer(r,"mode");check(1<=p<=8 && !(p in modes),"Duplicate/invalid local mode")
            check(integer(r,"cap")==cap,"Local-observable cap mismatch")
            push!(modes,p);up[p]=number(r,"Nup");down[p]=number(r,"Ndown")
            density[p]=number(r,"density");doubles[p]=number(r,"double_occupancy");sz[p]=number(r,"Sz_local")
        end
        check(all(isfinite,C) && all(isfinite,Z) && all(isfinite,vcat(up,down,density,doubles,sz)),"Nonfinite raw data")
        compare(maximum(abs,density-up-down),0.,"raw_density_identity",residuals)
        compare(maximum(abs,sz-(up-down)/2),0.,"raw_local_Sz_identity",residuals)
        compare(maximum(abs,diag(Z)-(density-2doubles)/4),0.,"raw_longitudinal_contact_identity",residuals)
        compare(maximum(abs,C-C'),0.,"raw_scalar_Hermiticity",residuals)
        # These are recomputed from raw CSV values and the original plain ED reference.
        quantities=Dict(
            "energy_error"=>abs(number(row,"energy")-sec.energy),
            "spin_matrix_error"=>maximum(abs,C-sec.spin_dot),
            "total_spin_squared"=>sum(C),
            "spin_squared_error"=>abs(sum(C)-sec.total_spin_squared),
            "Sodd_z"=>dot(odd,Z*odd)/ref.cells,
            "Sodd_scalar"=>dot(odd,C*odd)/(3ref.cells),
            "Delta"=>dot(odd,sz)/ref.cells,
            "particle_error"=>max(abs(sum(up)-sec.Nup),abs(sum(down)-sec.Ndown)),
            "sz_error"=>abs(sum(sz)-sec.Sz),
            "contact_error"=>maximum(abs,diag(C)-.75*(density-2doubles)),
            "density_error"=>maximum(abs,density-sec.density),
            "double_occupancy_error"=>maximum(abs,doubles-sec.double_occupancy))
        for (name,value) in quantities;compare(value,number(row,name),name,residuals);end
        # Only the reported variance/internal diagnostics are available in CSV;
        # they are gate-checked, not represented as independent MPS measurements.
        variance=number(row,"variance");fockvariance=number(row,"fock_variance")
        actionerror=number(row,"energy_action_error");leakage=number(row,"sector_leakage")
        check(all(isfinite,(variance,fockvariance,actionerror,leakage)),"Nonfinite reported diagnostic")
        check(fockvariance>=0 && actionerror>=0 && leakage>=0,"Negative squared/absolute diagnostic")
        gates=Dict("energy_pass"=>quantities["energy_error"]<=ENERGY_ATOL,
            "variance_pass"=>(-VARIANCE_ATOL<=variance<=VARIANCE_ATOL && fockvariance<=VARIANCE_ATOL),
            "spin_pass"=>(quantities["spin_matrix_error"]<=SPIN_ATOL && quantities["spin_squared_error"]<=SPIN_ATOL),
            "number_pass"=>(quantities["particle_error"]<=NUMBER_ATOL && quantities["sz_error"]<=NUMBER_ATOL),
            "internal_pass"=>(actionerror<=ENERGY_ATOL && quantities["contact_error"]<=SPIN_ATOL && leakage<=NUMBER_ATOL))
        gates["passed"]=all(values(gates))
        for (name,pass) in gates;check(boolean(row,name)==pass,"Exported stage gate $name disagrees at $key");end
        check(boolean(row,"ensemble_sensitive")== (sec.ground_degeneracy>1),"Degeneracy flag mismatch")
        data=(;C,Z,E=number(row,"energy"),variance,fockvariance,fidelity=number(row,"fidelity"),quantities,gates)
        recomputed[key]=data
        push!(audit_rows,merge(Dict{String,Any}("Sz"=>M,"seed"=>seed,"stage"=>stage,"cap"=>cap),
            Dict{String,Any}(quantities),Dict{String,Any}(gates),
            Dict("reported_variance"=>variance,"reported_fock_variance"=>fockvariance,
                 "reported_energy_action_error"=>actionerror,"reported_sector_leakage"=>leakage)))
    end
    agreements=read_csv(joinpath(attempt_dir,"initialization_agreement.csv"))
    check(length(agreements)==2 && Set(integer(r,"Sz") for r in agreements)==Set((0,1)),"Two-start agreement coverage mismatch")
    agreement_pass=true;agreement_rows=Dict{String,Any}[]
    for row in agreements
        M=integer(row,"Sz");a=recomputed[M,SEEDS[1],5];b=recomputed[M,SEEDS[2],5]
        dE=abs(a.E-b.E);dC=maximum(abs,a.C-b.C);dS=abs(sum(a.C)-sum(b.C))
        compare(dE,number(row,"energy_difference"),"two_start_energy_difference",residuals)
        compare(dC,number(row,"spin_matrix_difference"),"two_start_spin_matrix_difference",residuals)
        compare(dS,number(row,"S2_difference"),"two_start_S2_difference",residuals)
        pass=dE<=2ENERGY_ATOL && dC<=2SPIN_ATOL && dS<=2SPIN_ATOL
        check(boolean(row,"within_twice_reference_tolerance")==pass,"Two-start gate mismatch")
        agreement_pass &= pass
        push!(agreement_rows,Dict("Sz"=>M,"energy_difference"=>dE,"spin_matrix_difference"=>dC,"S2_difference"=>dS,"passed"=>pass))
    end
    final=[recomputed[M,seed,5] for M in 0:1 for seed in SEEDS]
    final_passes=count(r->r.gates["passed"],final)
    ensemble_sensitive=any(ref.sectors[M].ground_degeneracy>1 for M in 0:1)
    passed=final_passes==4 && agreement_pass
    expected_status=passed ? "PASS" : ensemble_sensitive ? "INCONCLUSIVE_ENSEMBLE_SENSITIVE" : "FAIL_REGISTERED_GATES"
    check(strip(read(joinpath(attempt_dir,"status.txt"),String))==expected_status && summary["status"]==expected_status,"Overall status mismatch")
    check(summary["completed_stages"]==20 && summary["completed_arms"]==4,"Summary stage/arm count mismatch")
    check(summary["initialization_agreement_pass"]==agreement_pass && summary["ensemble_sensitive"]==ensemble_sensitive,"Summary flags mismatch")
    maxima=Dict("max_final_energy_error"=>maximum(r.quantities["energy_error"] for r in final),
        "max_final_spin_matrix_error"=>maximum(r.quantities["spin_matrix_error"] for r in final),
        "max_final_S2_error"=>maximum(r.quantities["spin_squared_error"] for r in final),
        "max_final_variance"=>maximum(r.variance for r in final),
        "max_final_fock_variance"=>maximum(r.fockvariance for r in final),
        "minimum_final_ground_subspace_fidelity"=>minimum(r.fidelity for r in final))
    for (name,value) in maxima;compare(value,summary[name],"summary_"*name,residuals);end
    maxima["max_absolute_final_reported_variance"]=maximum(abs(r.variance) for r in final)
    maxima["max_final_particle_error"]=maximum(r.quantities["particle_error"] for r in final)
    maxima["max_final_Sz_error"]=maximum(r.quantities["sz_error"] for r in final)
    check(all(sha(joinpath(attempt_dir,file))==digest for(file,digest) in hashes),"Returned input changed during audit")
    Dict{String,Any}("audit_status"=>"PASS","benchmark_status"=>expected_status,
        "created_utc"=>string(now(UTC))*"Z","audit_julia_version"=>string(VERSION),
        "attempt_directory"=>abspath(attempt_dir),"original_ED_reference"=>abspath(reference_path),
        "original_ED_reference_sha256"=>sha(reference_path),"audit_source_sha256"=>sha(@__FILE__),
        "stage_count"=>20,"final_count"=>4,"final_pass_count"=>final_passes,
        "raw_spin_rows"=>1280,"raw_local_rows"=>160,"comparison_absolute_tolerance"=>COMPARE_ATOL,
        "maximum_exported_value_discrepancy"=>maximum(values(residuals)),
        "input_sha256"=>hashes,"comparison_residuals"=>residuals,"final_extrema"=>maxima,
        "stage_audits"=>audit_rows,"two_start_agreement"=>agreement_rows,
        "scalar_convention"=>"Sodd_scalar = a' Cdot a/(3*Ncells); ED scalar/Ncells is three times this field.",
        "independently_recomputed_from_raw_CSV"=>["energy error against original ED","full spin-dot matrix error",
            "total S2","orbital odd scalar and z structure factors","Delta","particle and Sz errors",
            "contact terms","density and double-occupancy errors","two-start differences","stage and final gate arithmetic"],
        "reported_only_gate_checked"=>["MPO moment-subtraction variance","dense-H residual variance",
            "energy-action disagreement","sector leakage"],
        "limitation"=>"This audit independently recomputes arithmetic from returned CSV and the original ED reference. It does not reconstruct MPS states, recompute energies/variances from MPS coefficients, rerun MPO action, or prove a bulk magnetic phase.")
end
function audit_return(;attempt_dir=DEFAULT_ATTEMPT,reference_path=DEFAULT_REFERENCE,output_path=DEFAULT_OUTPUT)
    check(!ispath(output_path),"Refusing existing audit output")
    try
        report=compute_audit(attempt_dir,reference_path);write_report(output_path,report)
        println("CSV audit PASS; benchmark=",report["benchmark_status"],"; stages20; final gates ",report["final_pass_count"],"/4")
        println("Max recomputed/exported discrepancy=",report["maximum_exported_value_discrepancy"])
        println("Saved: ",output_path);report
    catch err
        write_report(output_path,Dict("audit_status"=>"FAILED","error"=>sprint(showerror,err),
            "attempt_directory"=>abspath(attempt_dir),"original_ED_reference"=>abspath(reference_path),
            "audit_source_sha256"=>sha(@__FILE__),"created_utc"=>string(now(UTC))*"Z"))
        rethrow()
    end
end
function selftest()
    check(csv_fields("a,b,c")==["a","b","c"],"CSV simple parser fixture")
    check(csv_fields("a,\"b,c\",\"d\"\"e\"")==["a","b,c","d\"e"],"CSV quoted parser fixture")
    check(csv_fields("a,,c")==["a","","c"],"CSV empty field fixture")
    check(csv_fields("a,b\r\n")==["a","b"],"CSV CRLF fixture")
    check(csv_fields("\"π\",\"-0.65\"")==["π","-0.65"],"CSV Unicode fixture")
    println("Five CSV parser fixtures passed")
end
function main(args=ARGS)
    if args==["--selftest"];selftest();return;end
    isempty(args) || error("Usage: alm_dmrg_return_audit_20260911.jl [--selftest]; use API for explicit paths")
    audit_return()
end
end
if abspath(PROGRAM_FILE)==@__FILE__;ALMDMRGReturnAudit20260911.main();end
