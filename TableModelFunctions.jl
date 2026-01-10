using SpectralFitting, CFITSIO, ProgressMeter #not sure if you want to use something other than CFITSIO for writing files, also ProgressMeter could be removed.
module TableModels #used module definition to be able to change my struct definitions 
using SpectralFitting

mutable struct TableParameter{T} #parameter used in the table model
    _param::FitParam{T} 
    name::String 
    num_vals::Int #number of values evaluated for in the table model
    logged::Bool #if true the paraters evalated are over a logged range
end
TableParameter(_param, name; num_vals=10, logged=false) = TableParameter(_param, name, num_vals, logged) 

mutable struct TableModel #table model struct.
    free_params::Vector
    frozen_params::Vector
    Spectra_Units::String
    Out_Path::String #ouput file path/name
    Redshift::Bool
    Escale::Bool
    E_units::String #energy units
    E_min::Float64  #minmum energy of spectrum
    E_max::Float64 #maximum energy of spectrum
    E_bins::Int #number of energy bins
end
TableModel(free_params, frozen_params; Spectra_Units="photons/cm^2/s", Out_Path="Table_Model.fits", Redshift=false, Escale=false, E_units="keV", E_min=0.5, E_max=70.0, E_bins = 1000) = TableModel(free_params, frozen_params, Spectra_Units, Out_Path, Redshift, Escale, E_units, E_min, E_max, E_bins)
end

function ReplaceGreek(s, mapping) #this function repleaces any greek letters in paramter names to regular letters to comply with the fits standard 
    return join([get(mapping, c, string(c)) for c in s])
end

function addparams(A,B) #for creating the paramters table PARAMVAL column. 'Adds' two vectors of paramters together to create a single vector iterating through all their combinations
    out = []
    for a in A
        for b in B
            push!(out,[a; b])
        end
    end
    return out
end

function multiplyparams(V) #for creating the paramters table PARAMVAL column. Takes the vector of vectors of parameter values and 'multiplies' out to ctrate a single vector of all their possible combinations
    out = V[1]
    for i in 2:length(V)
        out = addparams(out,V[i])
    end
    return out
end

function MakeTable(model) #takes any spectral fitting model and creates the initial table model structure.
    full_model_vals , model_names = SpectralFitting._all_parameters_with_names(model)

if typeof(model).name.name == :CompositeModel #should handle composite or non composite models. 
    full_model_names = model_names
else
    full_model_names = String.(stack(split.(model_names,"."))[2,:])
end

deleteat!(full_model_vals, findall(x-> occursin(r"(?<!\.)\bK\b|\.K",x),full_model_names)) #removes normalisation paramters as they aren't needed
deleteat!(full_model_names, findall(x-> occursin(r"(?<!\.)\bK\b|\.K",x),full_model_names))

NamestoValues = Dict(full_model_names.=>full_model_vals) 

free_param_names = filter(x -> SpectralFitting.isfree(NamestoValues[x]), full_model_names)
free_param_values = filter(x -> SpectralFitting.isfree(x), full_model_vals)

frozen_param_names = filter(x -> !SpectralFitting.isfree(NamestoValues[x]), full_model_names)
frozen_param_values = filter(x -> !SpectralFitting.isfree(x), full_model_vals)

free_table_params = TableModels.TableParameter.(free_param_values,free_param_names)
frozen_table_params = TableModels.TableParameter.(frozen_param_values,frozen_param_names) 

# this is where I originaly had a print statment similar to when a model is called in SpectralFitting.
    
return TableModels.TableModel(free_table_params,frozen_table_params)

end


function OutputTable(model,TableModel) #this writes the table model to the fits file. It currently takes as an argument the model and the TableModel struct. I couldn'f find a way to take just the TableModel
    Model_Name = splitpath(TableModel.Out_Path)[end]
    file_name = splitpath(TableModel.Out_Path)[end] #this section hopfully handles whatever is used as Out_Path for writing the fits file to 
    if length(file_name) <5 
        file_name *= ".fits"
        TableModel.Out_Path *= ".fits"
    elseif file_name[end-4:end] !== ".fits"
        file_name = file_name*".fits"
        TableModel.Out_Path *= ".fits"
    elseif file_name[end-4:end] == ".fits"
        Model_Name = Model_Name[1:end-5]
    end
    
    if modelkind(model) == Additive()
        AddModel = "T"
    else
        AddModel = "F"
    end

    if TableModel.Redshift
        REDSHIFT = "T"
    else
        REDSHIFT = "F"
    end

    if TableModel.Escale
        ESCALE = "T"
    else
        ESCALE = "F"
    end

    #dictionary to convert from greek letters to english text. If you can find a better way to do this please do it. 
    greek_map = Dict(
        'α' => "alpha",  'β' => "beta",  'γ' => "gamma",  'δ' => "delta",
        'ϵ' => "epsilon",  'ζ' => "zeta",  'η' => "eta",  'θ' => "theta",
        'ι' => "iota",  'κ' => "kappa",  'λ' => "lambda",  'μ' => "mu",
        'ν' => "nu",  'ξ' => "xi",  'ο' => "omicron",  'π' => "pi",
        'ρ' => "rho",  'σ' => "sigma",  'τ' => "tau",  'υ' => "upsilon",
        'ϕ' => "phi", 'χ' => "chi", 'ψ' => "psi", 'ω' => "omega",
        'Α' => "Alpha",  'Β' => "Beta",  'Γ' => "Gamma",  'Δ' => "Delta",
        'Ε' => "Epsilon",  'Ζ' => "Zeta",  'Η' => "Eta",  'Θ' => "Theta",
        'Ι' => "Iota",  'Κ' => "Kappa",  'Λ' => "Lambda",  'Μ' => "Mu",
        'Ν' => "Nu",  'Ξ' => "Xi",  'Ο' => "Omega",  'Π' => "Pi",
        'Ρ' => "Rho",  'Σ' => "Sigma",  'Τ' => "Tau",  'Υ' => "Upsilon",
        'Φ' => "Phi", 'Χ' => "Chi", 'Ψ' => "Psi", 'Ω' => "Omega")

# Writing the Primary Header 

f = fits_clobber_file(TableModel.Out_Path)
fits_create_empty_img(f)
fits_write_key(f,"MODLNAME", Model_Name, "the name of the model")
fits_write_key(f,"MODLUNIT", TableModel.Spectra_Units, "the units for the model")
fits_write_key(f,"REDSHIFT", REDSHIFT, "whether the model contains redshift as a parameter")
fits_write_key(f,"ESCALE", ESCALE, " whether escale is to be a parameter")
fits_write_key(f,"ADDMODEL", AddModel, "whether the model is additive or not")
fits_write_key(f,"LOELIMIT", 0, "the model value for energies below those tabulated")
fits_write_key(f,"HIELIMIT", 0, "the model value for energies above those tabulated")

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUVERS", "1.1.0", "format version")

# Parameters 

free_param_names = [ReplaceGreek(x.name,greek_map) for x in TableModel.free_params]
free_param_values = [x._param for x in TableModel.free_params]
logged = [Int(x.logged) for x in TableModel.free_params]
NumbVals = [x.num_vals for x in TableModel.free_params]

params=[]
for P in TableModel.free_params
    if P.logged
        push!(params,collect(logrange(SpectralFitting.get_lowerlimit(P._param),SpectralFitting.get_upperlimit(P._param),P.num_vals)))
    else
        push!(params,collect(range(SpectralFitting.get_lowerlimit(P._param),SpectralFitting.get_upperlimit(P._param),P.num_vals)))
    end
end
# Writing the Parameters Table
PARAMETERS_colsdef = [("NAME", "12A", ""),("METHOD", "J", ""),("INITIAL", "E", ""),("DELTA", "E", ""),("MINIMUM", "E", ""),("BOTTOM", "E", ""),("TOP", "E", ""),("MAXIMUM", "E", ""),("NUMBVALS", "J", ""),("VALUE", string(maximum(NumbVals))*"E", "")]
fits_create_binary_tbl(f, 0, PARAMETERS_colsdef, "PARAMETERS")
fits_write_col(f, 1, 1, 1, free_param_names)
fits_write_col(f, 2, 1, 1, logged) 
fits_write_col(f, 3, 1, 1, SpectralFitting.get_value.(free_param_values))
fits_write_col(f, 4, 1, 1, SpectralFitting.get_error.(free_param_values))
fits_write_col(f, 5, 1, 1, SpectralFitting.get_lowerlimit.(free_param_values))
fits_write_col(f, 6, 1, 1, SpectralFitting.get_lowerlimit.(free_param_values))
fits_write_col(f, 7, 1, 1, SpectralFitting.get_upperlimit.(free_param_values))
fits_write_col(f, 8, 1, 1, SpectralFitting.get_upperlimit.(free_param_values))
fits_write_col(f, 9, 1, 1, NumbVals)
for i in eachindex(params)
    fits_write_col(f, 10, i, 1,params[i])
end

for P in TableModel.frozen_params
    fits_write_key(f,ReplaceGreek(P.name,greek_map),SpectralFitting.get_value(P._param),"physical parameter held constant") 
end
fits_write_key(f,"NINTPARM", length(free_param_names), "the number of interpolated parameters")
fits_write_key(f,"NADDPARM", 0, "the number of additional parameters")
fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "PARAMETERS", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")

# Writing the Parameters Table

ENERGIES_colsdef = [("ENERG_LO", "E", TableModel.E_units),("ENERG_HI", "E", TableModel.E_units)]
fits_create_binary_tbl(f, TableModel.E_bins, ENERGIES_colsdef, "ENERGIES")
Energies = collect(logrange(TableModel.E_min,TableModel.E_max,TableModel.E_bins+1))
E_low = Energies[1:end-1]
E_high = Energies[2:end]
fits_write_col(f, 1, 1, 1, E_low)
fits_write_col(f, 2, 1, 1, E_high)

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "ENERGIES", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")

# Writing the Spectrums Table

iter_params = multiplyparams(params)

SPECTRA_colsdef = [("PARAMVAL", string(length(TableModel.free_params))*"E", ""),("INTPSPEC", string(TableModel.E_bins)*"E", TableModel.Spectra_Units)]
fits_create_binary_tbl(f, prod(length.(params)), SPECTRA_colsdef, "SPECTRA")
fits_write_col(f, 1, 1, 1, vec(stack(iter_params)))

@showprogress for j in eachindex(iter_params)
    ps = iter_params[j]
    for i in eachindex(ps)
        symbol_name = split(TableModel.free_params[i].name,".")
        if length(symbol_name) == 2
            setproperty!(getproperty(model,Symbol(symbol_name[1])), Symbol(symbol_name[2]), ps[i])
        else
            setproperty!(model,Symbol(symbol_name[1]),ps[i])
        end
        global spec = invokemodel(Energies,model)
        fits_write_col(f, 2, j, 1, spec.parent[:,1]) #CHECK THIS <<<<------------<<<<
    end
end

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "MODEL SPECTRA", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")
close(f)

end

##

# this is the rough way to make the table model. First call MakeTable() on the model. The model should have desired paramters unfrozen and limits on those apropriatly set

TableModel = MakeTable(model)
#this is where a print out allowing you to easly change parameters would be usefull

# then call the TableModel and the model to OutputTable() to write the fits file
OutputTable(model,TableModel)
