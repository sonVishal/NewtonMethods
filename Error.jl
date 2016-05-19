import Base: showerror

# TODO: Do not know why they are throwing a load error alongwith the normal error

# Error when the option has invalid value(s)
type InvalidOption <: Exception
    optName::AbstractString
    value::Any
end

Base.showerror(io::IO, e::InvalidOption) = print(io, "The value \"", e.value, "\" is not valid for the option ", e.optName,
".\nPlease assign the correct value and try again!\n")

# Error when the function evaluation is not successful
type EvaluationError <: Exception
    funcName::Function
end

# TODO: Display why the function failed

Base.showerror(io::IO, e::EvaluationError) = print(io, "The function \"",
string(e.funcName),"\" could not be evaluated.\n",
"Please check that the function is correct!\n")

# Error when the option 1 is incompatible with option 2
type IncompatibleOption <: Exception
    optName1::AbstractString
    optName2::AbstractString
end

Base.showerror(io::IO, e::InvalidOption) = print(io, "The option \"",
string(e.optName1), "\" is not compatible with option ", string(e.optName2),
".\n Please assign compatible options and try again!\n")
