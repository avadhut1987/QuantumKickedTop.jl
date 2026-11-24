module QuantumKickedTop

    # Classical system
    include("classical_dynamics/ClassicalMap.jl")
    include("classical_dynamics/ClassicalUtils.jl")
    using .ClassicalMap, .ClassicalUtils

    # Quantum system
    include("quantum_dynamics/FloquetSystem.jl")
    include("quantum_dynamics/PhiStates.jl")
    using .FloquetSystem, .PhiStates

    # Quantum measures
    include("quantum_information_measures/QuantumUtils.jl")
    include("quantum_information_measures/QuantumUtilsDiscord.jl")
    include("quantum_information_measures/EntanglementDynamics.jl")
    include("eigenstate_entanglement_studies/EigenstateEntanglement.jl")
    using .QuantumUtils, .QuantumUtilsDiscord, .EntanglementDynamics, .EigenstateEntanglement

    # Parameter scans
    include("parameter_scans/LinearEntropyVsKR.jl")
    using .LinearEntropyVsKR

    # Analysis tools
    include("visualization/HeatmapEntropy.jl")
    include("visualization/HeatmapTimeAveragedEntropy.jl")
    include("visualization/HeatmapDiscord.jl")
    include("visualization/HeatmapLinearEntropy.jl")
    include("visualization/HeatmapLLE.jl")
    using .HeatmapEntropy, .HeatmapTimeAveragedEntropy, .HeatmapDiscord, .HeatmapLinearEntropy, .HeatmapLLE

end
