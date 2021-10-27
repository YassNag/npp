########### Importations ###########

# Système
using Dates
using Random

# Affichage graphique
using PyCall
using PyPlot
pygui(:qt5)
Path = PyPlot.matplotlib.path.Path

# Modélisation
using JuMP

# Solveurs libres
using GLPK
using Cbc
using SCIP

# Solveurs commerciaux
#using Gurobi
using CPLEX
using MosekTools

#Graph library
using LightGraphs

# Divers
include("ascii_art.jl")

# Essentiels
include("structures.jl")
include("parser.jl")
include("utilitaires.jl")
include("visualisation.jl")
include("controles_saisies.jl")
include("parametres_solveurs.jl")

# Modèle mono-objectif
#include("modele1.jl")
#include("modele2.jl")
#include("modele3.jl")
#include("modele4.jl")
#include("modele5.jl")
# include("modele_symFlowModified.jl")
include("modelcapacity.jl")


########### Programme principal ###########

function main()
    # Nettoyage de la console
    clear()
    ascii_art_pmcn()

    # Choix de l'instance
    instance = lecture_instance()

    # Affichage du contexte
    choix = choix_binaire("\n --> Souhaitez-vous afficher le contexte (o/n) ? ")
    if choix == "o"
        affichage_contexte(instance)
    end

    # Visualisation du scenario
    choix = choix_binaire("\n --> Souhaitez-vous visualiser le scénario (o/n) ? ")
    if choix == "o"
        visualisation_scenario(instance)
    end

    # Poursuivre vers la résolution
    choix = choix_binaire("\n --> Souhaitez-vous poursuivre vers la résolution (o/n) ? ")
    if choix == "o"
        # Résolution de l'instance
        solution = resolution_modele(instance)

        # Visualisation de la solution
        if solution.statut == 1
            choix = choix_binaire("\n --> Souhaitez-vous visualiser la solution (o/n) ? ")
            if choix == "o"
                visualisation(instance, solution)
            end
        end
    end
    # Fin
    ascii_art_fin()
end

main()
