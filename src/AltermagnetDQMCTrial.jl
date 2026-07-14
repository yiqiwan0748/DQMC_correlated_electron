module AltermagnetDQMCTrial

export ModelParameters,
       DEFAULT_PARAMETERS,
       bloch_hamiltonian,
       band_energies,
       BandGridSamples,
       band_path,
       sample_band_energies,
       sorted_band_energies,
       density_of_states,
       zero_temperature_filling,
       estimate_van_hove_peaks,
       infer_doping_convention

include("noninteracting.jl")

end
