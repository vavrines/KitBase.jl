"""
$(SIGNATURES)

Generate variables from dictionary
"""
function generate_vars(dict::AbstractDict)
    for key in keys(dict)
        s = key
        @eval $s = $(dict[key])
    end
end


"""
$(SIGNATURES)

Run command and collect output from shell
"""
function collect_run(cmd::Cmd)
    out = Pipe()
    err = Pipe()

    process = run(pipeline(ignorestatus(cmd), stdout = out, stderr = err))
    close(out.in)
    close(err.in)
    stdout = @async String(read(out))
    stderr = @async String(read(err))

    return (stdout = String(read(out)), stderr = String(read(err)), code = process.exitcode)
end
