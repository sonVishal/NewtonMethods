# NLEQ-Options for Nonlinear Equations solvers

import Base: show

"""macro for importing OptionsNLEQ and option handling.
macro import_options()
    :(
    using NLEQ: OptionsNLEQ, getOption, setOption!, setOptions!,
    copyOptions!
    )
end

super-type for all types storing options for NLEQ solvers.

NLEQ solvers often have serveral parameters for fine-tuning them.
In this NLEQInterface this parameters are called 'options' and
they are stored in key/value paris. For the key a
`ASCIIString` is used. The value can be `Any`-thing.
The key is often called the option-name.

All types for this purpose have this abstract type as super-type.

Required fields are: `name`, `lastchanged`, `options`
"""
abstract AbstractOptionsNLEQ <: Any

"""
Stores options for NLEQ Solver(s) together with a name.
Additionally the time of the last change is saved.

Options can be set at construction time, e.g.

opt=OptionsNLEQ("test",
"loglevel" => NLEQ.LOG_ALL,
"logio"    => STDERR)

or later. For changing single options

oldValue = setOption!(opt,"myopt","new value")
oldValue = setOption!(opt,"myopt" => "new value")

and for changing many options at once:

oldValues = setOption!(opt, "myopt" => "new value",
    "oldopt" => 56)

see also: `setOption!`, `setOptions!`
"""
type OptionsNLEQ <: AbstractOptionsNLEQ
    name        :: AbstractString
    lastchanged :: DateTime
    options     :: Dict{ASCIIString,Any}

    function OptionsNLEQ(name::AbstractString="")
        obj = new(name,now(),Dict{ASCIIString,Any}())
        return obj
    end
end

function OptionsNLEQ(name::AbstractString,copyOptionsFrom::AbstractOptionsNLEQ)
    opt = OptionsNLEQ(name)
    copyOptions!(opt,copyOptionsFrom)
    return opt
end

function OptionsNLEQ(copyOptionsFrom::AbstractOptionsNLEQ)
    return OptionsNLEQ("",copyOptionsFrom)
end

function OptionsNLEQ(name::AbstractString,pairs::Pair...)
    opt = OptionsNLEQ(name)
    setOptions!(opt,pairs...)
    return opt
end

function OptionsNLEQ(pairs::Pair...)
    return OptionsNLEQ("",pairs...)
end

"""
function getOption(opt::AbstractOptionsNLEQ,name::ASCIIString,
default::Any=nothing)

get saved value of option with given `name` or `default`
if option is unknown.
"""
function getOption(opt::AbstractOptionsNLEQ,name::ASCIIString,
    default::Any=nothing)
    return haskey(opt.options,name)?opt.options[name]:default
end

"""
function setOption!(opt::AbstractOptionsNLEQ,name::ASCIIString,value::Any)

set NLEQ-Option with given `name` and return old value
(or `nothing` if there was no old value).
"""
function setOption!(opt::AbstractOptionsNLEQ,name::ASCIIString,value::Any)
    oldValue = getOption(opt,name)
    opt.options[name]=value
    opt.lastchanged=now()
    return oldValue
end

"""
function setOption!(opt::AbstractOptionsNLEQ,pair::Pair)

set NLEQ-Option with given (`name`,`value`) pair and return old value
(or `nothing` if there was no old value).
"""
function setOption!(opt::AbstractOptionsNLEQ,pair::Pair)
    return setOption!(opt,pair.first,pair.second)
end

"""
function setOptions!(opt::AbstractOptionsNLEQ,pairs::Pair...)

set many NLEQ-Options and return an array with the old option values.
"""
function setOptions!(opt::AbstractOptionsNLEQ,pairs::Pair...)
    oldValues=Any[]
    for (name,value) in pairs
        push!(oldValues, setOption!(opt,name,value))
    end
    return oldValues
end

"""
function initOption!(opt::AbstractOptionsNLEQ,name::ASCIIString,value::Any)

initialize NLEQ-Option with given `name` only if it does not exist already and return the old value if changed
(or `nothing` if there was an old value).
"""
function initOption!(opt::AbstractOptionsNLEQ,name::ASCIIString,value::Any)
    return haskey(opt.options,name)?nothing:setOption!(opt,name,value);
end

"""
function initOption!(opt::AbstractOptionsNLEQ,pair::Pair)

initialize NLEQ-Option with given (`name`,`value`) pair only if it does not exist already and return the old value if changed
(or `nothing` if there was an old value).
"""
function initOption!(opt::AbstractOptionsNLEQ,pair::Pair)
    return initOption!(opt,pair.first,pair.second);
end

"""
function initOptions!(opt::AbstractOptionsNLEQ,pairs::Pair...)

initialize many NLEQ-Options and return an array with the old option values if the options did not exist before
"""
function initOptions!(opt::AbstractOptionsNLEQ,pairs::Pair...)
    oldValues=Any[]
    for (name,value) in pairs
        push!(oldValues, initOption!(opt,name,value))
    end
    return oldValues
end

"""
function copyOptions!(dest::AbstractOptionsNLEQ,source::AbstractOptionsNLEQ)

copy all options from other NLEQ-Option object.
"""
function copyOptions!(dest::AbstractOptionsNLEQ,source::AbstractOptionsNLEQ)
    merge!(dest.options,source.options)
    dest.lastchanged=now()
    return dest
end

function show(io::IO, opt::AbstractOptionsNLEQ)
    print(io,typeof(opt)," ")
    isempty(opt.name) || print(io,"'",opt.name,"' ")
    len=length(opt.options)
    print(io,"with ",len," option",len!=1?"s":"",len>0?":":"."); println(io)
    if len>0
        maxLen=2+max( 0,map(length,keys(opt.options))... )
        for key in sort(collect(keys(opt.options)))
            print(io,lpad(key,maxLen),": ")
            show(io,opt.options[key]); println(io)
        end
    end
    print(io,"lasted changed "); show(io,opt.lastchanged);
end
