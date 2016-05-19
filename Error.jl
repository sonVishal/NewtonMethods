import Base: showerror

type InvalidOption <: Exception
    optName::AbstractString
    value::Any
end

# TODO: Do not know why this is throwing a load error alongwith the normal error

Base.showerror(io::IO, e::InvalidOption) = print(io, "The value \"", e.value, "\" is not valid for the option ", e.optName,
".\n\tPlease assign the correct value and try again!\n")
